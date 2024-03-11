### this file is supposed to contain the general model without any specifications to a plan case

using JuMP
using GLPK

include("utils.jl")

path = "src/NuclearData.xlsx"

years = 2030:2100

transport_details = get_transport_details(path)

transport_costs = transport_details[2]
TRANSPORT_POSSIBLE = transport_details[1]

snf = get_snf_at_reactors(path)
reactors = keys(snf)

production_capacities = get_hot_cell_capacities(path)
hot_cells = keys(production_capacities)

storage_cap = get_storage_capacities(path)
nodes = keys(storage_cap)

for n in nodes
    if !haskey(snf, n)
        snf[n] = 0
        end
end

cisf_costs = get_cisf_costs(path)
interim_storages = keys(cisf_costs)

BIG_trans = 10000000

CASK_DECAY = 0.05
BIG_cap = 1000000
BIG_cisf = BIG_cap ^2

#create model and attach solver
model = Model()
set_optimizer(model, GLPK.Optimizer)

#create all variables we need
@variable(model, SNF_t[nodes, nodes, years] >= 0)
@variable(model, SNF_s[nodes, (minimum(years) - 1):maximum(years)] >= 0)
@variable(model, NC_t[nodes, nodes, years] >= 0)
@variable(model, NC_s[nodes, years] >= 0)
@variable(model, NC_t_end[nodes] >= 0)
@variable(model, B[interim_storages, years], Bin)

#TODO add cost parameters to obj function!!!
#TODO are costs for SNF and NC the same? what measurments do we use?
#TODO two transport modes? really nessecary?
@objective(model, Min, 
            sum(sum(transport_costs[string(n, "-", m)] * (SNF_t[n, m, y] + NC_t[n, m, y]) for n in nodes, m in nodes) + 
            sum(SNF_s[n, y] + NC_s[n, y] for n in nodes) + sum(B[i, y] for i in interim_storages) + 
            sum(SNF_t[n, hc, y] for hc in hot_cells, n in nodes) for y in years) + 
            sum(NC_t_end[i] for i in interim_storages)
)

# mass balances
@constraint(model, 
            mass_balance_SNF[n = nodes, y = years; n ∉ hot_cells], 
            sum(SNF_t[m, n, y] for m in nodes) - # incoming
            sum(SNF_t[n, m, y] for m in nodes) + # outgoing
            SNF_s[n, y - 1] == SNF_s[n, y])
@constraint(model, 
            mass_balance_NC[n = nodes, y = years; n ∉ hot_cells], 
            sum(NC_t[m, n, y] for m in nodes) - sum(NC_t[n, m, y] for m in nodes) + 
            get_conditional_NC_s(model, n, y, y_ -> y_ > minimum(years)) ==  NC_s[n, y])
# TODO can we move casks from a cisf (before moving to end storage)?
# TODO can we store SNF in a cisf?

# hot cell shit here
@constraint(model, mass_balance_hc[hc = hot_cells, y = years], sum(SNF_t[n, hc, y] for n in nodes) == sum(NC_t[hc, n, y] for n in nodes))

@constraint(model, hc_production_cap[hc = hot_cells, y = years], sum(SNF_t[n, hc, y] for n in nodes) <= production_capacities[hc])

@constraint(model, no_casks_to_hc[hc = hot_cells, y = years], sum(NC_t[n, hc, y] for n in nodes) == 0)

@constraint(model, no_snf_from_hc[hc = hot_cells, y = years], sum(SNF_t[hc, n, y] for n in nodes) == 0)

# storage capacities
@constraint(model, storage[n = nodes, y = years], SNF_s[n, y] + NC_s[n, y] <= storage_cap[n]) # are capacities fix?

# transport is possible
@constraint(model, transport_possibility[n = nodes, m = nodes, y = years], SNF_t[n, m, y] + NC_t[n, m, y] <= TRANSPORT_POSSIBLE[string(n, "-", m)] * BIG_trans)

# transport to end storage facility
@constraint(model, transport_to_end[n = nodes], NC_t_end[n] == NC_s[n, maximum(years)])

# SNF clearing condition
@constraint(model, SNF_clear[n = nodes], SNF_s[n, maximum(years)] == 0)

# contraint for building cisf only once
@constraint(model, cisf_build[i = interim_storages], sum(B[i, y] for y in years) <= 1)

# dont store before you build
@constraint(model, cisf_build_before_store[i = interim_storages, y = years], sum(SNF_s[i, cur] + NC_s[i, cur] + sum(SNF_t[n, i, cur] + SNF_t[i, n, cur] + NC_t[i, n, cur] + NC_t[n, i, cur] for n in nodes) for cur in minimum(years): y) <= sum(B[i, cur] for cur in minimum(years): y) * BIG_cisf)

# cask decay
@constraint(model, cask_decay[y = (minimum(years)+1):maximum(years), n = nodes; n ∉ hot_cells], sum(SNF_t[n, hc, y] for hc in hot_cells) >= sum(CASK_DECAY * SNF_s[n, y - 1]))

# initialize snf at reactors
@constraint(model, initial_snf[n = nodes], SNF_s[n, minimum(years) - 1] == snf[n])

optimize!(model)
obj_value = objective_value(model)
println("Total costs are ", obj_value)

using DataFrames
using CSV

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