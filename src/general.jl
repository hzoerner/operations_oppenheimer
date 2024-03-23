### this file is supposed to contain the general model without any specifications to a plan case

using JuMP
using GLPK
using DataFrames
using CSV

using HiGHS

include("utils.jl")

path = "../operations_oppenheimer/ExtendedNuclearData.xlsx"

years = 2030:2099

transport_details = get_transport_details(path)

transport_costs = transport_details[2]
TRANSPORT_POSSIBLE = transport_details[1]

snf = get_snf_at_reactors(path)
reactors = keys(snf)

production_capacities = get_hot_cell_capacities(path)
hot_cells = keys(production_capacities)

storage_cap = get_storage_capacities(path)
nodes = keys(storage_cap)

esf_costs = get_end_storage_transport_costs(path)

cisf_costs = get_cisf_costs(path)
interim_storages = keys(cisf_costs)

general = get_general_data(path)

CISF_BUILDING_COSTS = general["CISF_building_costs"]
CISF_OPERATING_COSTS = general["CISF_operation_costs"]
REACTOR_OPERATING_COSTS = general["Reaktor_operation_costs"]

BIG_trans = 10000000

CASK_DECAY = 0.05
BIG_cap = 1000000
BIG_cisf = 1000000

#create model and attach solver
model = Model()
set_optimizer(model, HiGHS.Optimizer)

#create all variables we need
@variable(model, SNF_t[nodes, nodes, years] >= 0)
@variable(model, SNF_s[nodes, (minimum(years) - 1):maximum(years)] >= 0)
@variable(model, NC_t[nodes, nodes, years] >= 0)
@variable(model, NC_s[nodes, years] >= 0)
@variable(model, NC_t_end[nodes] >= 0)
@variable(model, B[interim_storages], Bin)

#TODO add cost parameters to obj function!!!
#TODO are costs for SNF and NC the same? what measurments do we use?
#TODO two transport modes? really nessecary?
@objective(model, Min, 
    sum(
        sum(transport_costs[string(n, "-", m)] * (SNF_t[n, m, y] + NC_t[n, m, y]) for n in nodes, m in nodes) + 
        REACTOR_OPERATING_COSTS * sum(SNF_s[r, y] + NC_s[r, y] for r in reactors) + 
        CISF_OPERATING_COSTS * sum( SNF_s[i, y] + NC_s[i, y] for i in interim_storages) for y in years) + 
    CISF_BUILDING_COSTS * sum(B[i] for i in interim_storages) + 
    sum(NC_t_end[n] * esf_costs[n] for n in nodes)
)

# mass balances
@constraint(
    model, 
    mass_balance_SNF[n = nodes, y = years; n ∉ hot_cells], 
    sum(SNF_t[m, n, y] for m in nodes) - 
    sum(SNF_t[n, m, y] for m in nodes) + SNF_s[n, y - 1] == SNF_s[n, y]
)

@constraint(
    model, 
    mass_balance_NC[n = nodes, y = years; n ∉ hot_cells], 
    sum(NC_t[m, n, y] for m in nodes) - 
    sum(NC_t[n, m, y] for m in nodes) + 
    get_conditional_variable(model, n, y, y_ -> y_ > minimum(years), NC_s) ==  NC_s[n, y]
)
# TODO can we move casks from a cisf (before moving to end storage)?
# TODO can we store SNF in a cisf?

# hot cell shit here (Kosten für Umladen notwendig?)
@constraint(
    model, 
    mass_balance_hc[hc = hot_cells, y = years], 
    sum(SNF_t[n, hc, y] for n in nodes) == sum(NC_t[hc, n, y] for n in nodes)
)

# hot cell processing capacities
@constraint(
    model, 
    hc_production_cap[hc = hot_cells, y = years], 
    sum(SNF_t[n, hc, y] for n in nodes) <= production_capacities[hc]
)

# No new casks to hot cell
@constraint(
    model, 
    no_casks_to_hc[hc = hot_cells, y = years], 
    sum(NC_t[n, hc, y] for n in nodes) == 0)

# no "old" cask from hot cell
@constraint(
    model, 
    no_snf_from_hc[hc = hot_cells, y = years], 
    sum(SNF_t[hc, n, y] for n in nodes) == 0
)

# storage capacities
@constraint(
    model, 
    storage[n = nodes, y = years], 
    SNF_s[n, y] + NC_s[n, y] <= storage_cap[n]
) # are capacities fix?

# transport is possible
@constraint(
    model, 
    transport_possibility[n = nodes, m = nodes, y = years], 
    SNF_t[n, m, y] + NC_t[n, m, y] <= TRANSPORT_POSSIBLE[string(n, "-", m)] * BIG_trans
)

# transport to end storage facility
@constraint(model, transport_to_end[n = nodes], NC_t_end[n] == NC_s[n, maximum(years)])


# SNF clearing condition
@constraint(
    model, 
    SNF_clear[n = nodes], 
    SNF_s[n, maximum(years)] == 0
)

# dont store before you build
# TODO double check!
@constraint(
    model, 
    cisf_build_before_store[i = interim_storages, y = years], 
    SNF_s[i, y] + NC_s[i, y] + 
    sum(NC_t[n, i, y] + 
    SNF_t[n, i, y] for n in nodes) <= B[i] * BIG_cisf
)

# constraints to prevent model from using cisfs as transport node
@constraint(
    model, 
    dont_move_nc_from_cisf[i = interim_storages, y = years], 
    sum(NC_t[i, n, y] for n in nodes) == 0
)

@constraint(
    model, 
    dont_move_snf_to_cisf[i = interim_storages, y = years], 
    sum(SNF_t[n, i, y] for n in nodes) == 0
)

# cask decay
@constraint(
    model, 
    cask_decay[y = (minimum(years)+1):maximum(years), 
    n = nodes; n ∉ hot_cells], 
    sum(SNF_t[n, hc, y] for hc in hot_cells) >= CASK_DECAY * SNF_s[n, y - 1]
)

# initialize snf at reactors
@constraint(
    model, 
    initial_snf[n = nodes], 
    SNF_s[n, minimum(years) - 1] == get_snf_at_node(snf, n)
)


optimize!(model)
obj_value = objective_value(model)
println("Total costs are ", obj_value)

df_SNF_s = DataFrame(node = String[], year = Int[], SNF = Float32[])
df_NC_s = DataFrame(node = String[], year = Int[], NC = Float32[])

for n in nodes, y in years
    append!(df_SNF_s, DataFrame(node = n, year = y, SNF = round(value.(SNF_s[n,y]))))
    append!(df_NC_s, DataFrame(node = n, year = y, NC = round(value.(NC_s[n,y]))))
end

CSV.write("snf_stored.csv", df_SNF_s)
CSV.write("nc_stored.csv", df_NC_s)

df_SNF_t = DataFrame(from = String[], to = String[], year = Int[], SNF = Float32[])
df_NC_t = DataFrame(from = String[], to = String[], year = Int[], NC = Float32[])

for n in nodes, m in nodes, y in years
    append!(df_SNF_t, DataFrame(from = n, to = m, year = y, SNF = round(value.(SNF_t[n,m,y]))))
    append!(df_NC_t, DataFrame(from = n, to = m, year = y, NC = round(value.(NC_t[n,m,y]))))
end

CSV.write("snf_shipped.csv", df_SNF_t)
CSV.write("nc_shipped.csv", df_NC_t)

df_Bi = DataFrame(cisf = String[], build = Int[])
for i in interim_storages
    append!(df_Bi, DataFrame(cisf = i, build = round(value.(B[i]))))
end

CSV.write("cisf_build.csv", df_Bi)

df_to_end = DataFrame(node = String[], NC = Int[])
for i in nodes
    append!(df_to_end, DataFrame(node = i,NC = round(value.(NC_t_end[i]))))
end

CSV.write("end_transport.csv", df_to_end)

println("Done!")