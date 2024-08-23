### this file is supposed to contain the general model without any specifications to a plan case
using Pkg
Pkg.add("JuMP")
Pkg.add("GLPK")
Pkg.add("HiGHS")
Pkg.add("DataFrames")
Pkg.add("CSV")
Pkg.add("XLSX")
using JuMP
using GLPK
using DataFrames
using CSV

using HiGHS

include("utils.jl")

path = ""

if isempty(ARGS)
    path = "dummy/path/ExtendedNuclearData.xlsx" #individual absolute filepath
else
    path = ARGS[1]
end


years = 2030:2060

transport_details = get_transport_details(path)

transport_costs = transport_details[2]
TRANSPORT_POSSIBLE = transport_details[1]

production_capacities = get_hot_cell_capacities(path)
hot_cells = keys(production_capacities)

storage_cap = get_storage_capacities(path)
nodes = keys(storage_cap)

esf_costs = get_end_storage_transport_costs(path)

cisf_costs = get_cisf_costs(path)
interim_storages = keys(cisf_costs)

general = get_general_data(path)

versions = [Version("small", 150, 96444026.38), Version("medium", 375, 142844431.51), Version("large", 500, 189244836.64)]

CISF_BUILDING_COSTS = general["CISF_building_costs"]
CISF_OPERATING_COSTS = general["CISF_operation_costs"]
REACTOR_OPERATING_COSTS = general["Reaktor_operation_costs"]
HC_BUILDING_COSTS = 300000000
FIX_COST_RATE = 100000
HOT_CELL_COSTS = 100 # just for model logic reasons
build_time = 20
snf = get_snf_at_reactors(path)
reactors = keys(snf)

reactor_costs = Dict()

for r in reactors
    reactor_costs[r] = REACTOR_OPERATING_COSTS
end

CASK_DECAY = 0.05

DISCOUNT = 1.03

INITIAL_VOLUME = sum(values(snf))

#create model and attach solver
model = Model()
set_optimizer(model, HiGHS.Optimizer)

#create all variables we need
@variable(model, SNF_t[nodes, nodes, years] >= 0)
@variable(model, SNF_s[nodes, (minimum(years) - 1):maximum(years)] >= 0)
@variable(model, NC_t[nodes, nodes, years] >= 0)
@variable(model, NC_s[nodes, years] >= 0)
@variable(model, NC_t_end[nodes] >= 0)
@variable(model, B[interim_storages, map(v -> v.size, versions), minimum(years)-1:maximum(years)], Bin)
@variable(model, A[reactors, years], Bin)
@variable(model, HC[hot_cells, minimum(years)-1:maximum(years)], Bin)

cost_factor = 1/1000

@objective(model, Min, 
    sum(
        sum(cost_factor * transport_costs[string(n, "-", m)] * (SNF_t[n, m, y] + NC_t[n, m, y]) for n in nodes, m in nodes) + 
        sum(cost_factor * reactor_costs[r] * (SNF_s[r, y] + NC_s[r, y]) for r in reactors) +
        cost_factor * HOT_CELL_COSTS * sum(SNF_t[r, hc, y] for r in reactors, hc in hot_cells) +
        cost_factor * CISF_OPERATING_COSTS * sum( SNF_s[i, y] + NC_s[i, y] for i in interim_storages) +
        cost_factor * FIX_COST_RATE * sum(A[r, y] for r in reactors) + 
        cost_factor * FIX_COST_RATE * sum(B[d, v.size, y] for v in versions, d in interim_storages) for y in years) + 

    #sum(LICENSE_EXTENSION_COSTS * A[r, minimum(years) + p * MAXIMUM_RUNTIME] for p in 1:EXTENSION_PERIOD) +
    sum(
        cost_factor * sum(v.costs * B[d, v.size, y] - B[d, v.size, y - 1] for v in versions, d in interim_storages) +
        cost_factor * HC_BUILDING_COSTS * sum(HC[d, y] - HC[d, y - 1] for d in hot_cells) for y in years )
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
    get_conditional_variable(n, y, y_ -> y_ > minimum(years), NC_s) ==  NC_s[n, y]
)

@constraint(
    model,
    no_waste_on_the_road[y = years],
    INITIAL_VOLUME == sum(NC_s[n, y] + SNF_s[n, y] for n in nodes)
)

 ### Hot Cell constraints
 @constraint(
     model, 
     hc_production_cap[hc = hot_cells, y = years], 
     sum(SNF_t[n, hc, y] for n in nodes) <= HC[hc, y] * 100
 )

 @constraint(
    model,
    hc_operating[hc = hot_cells, y = (minimum(years) - 1):(maximum(years) - 1)],
    HC[hc, y] <= HC[hc, y + 1]
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

# hot cell mass balance
@constraint(
    model, 
    mass_balance_hc[hc = hot_cells, y = years], 
    sum(SNF_t[n, hc, y] for n in nodes) == sum(NC_t[hc, n, y] for n in nodes)
)

# storage capacities
@constraint(
    model, 
    storage[n = nodes, y = years], 
    SNF_s[n, y] + NC_s[n, y] <= storage_cap[n]
) 

### transport capacity restrictions
@constraint(
    model, 
    transport_possibility[n = nodes, m = nodes, y = years], 
    SNF_t[n, m, y] + NC_t[n, m, y] <= TRANSPORT_POSSIBLE[string(n, "-", m)] * 50
)

@constraint(
    model,
    countrywide_transport_cap[y=years],
    sum(SNF_t[n, m, y] + NC_t[n, m, y] for n in nodes, m in nodes) <= 5*50
)

# permanent cisf storing constraint
@constraint(model, no_waste_at_reactors[r = reactors], NC_s[r, maximum(years)] == 0)

# SNF clearing condition
@constraint(
    model, 
    SNF_clear[n = nodes], 
    SNF_s[n, maximum(years)] == 0
)

### Reactor operation logic
@constraint(
    model,
    reactor_operating[r = reactors, y = years],
    SNF_s[r, y] + NC_s[r, y] <= 500 * A[r, y]
)

@constraint(
    model,
    operation_logic[r = reactors, y = minimum(years):(maximum(years) - 1)],
    A[r, y] >= A[r, y + 1]
)

# CISF related constraints
@constraint(
    model, 
    cisf_build_before_store[i = interim_storages, y = years], 
    sum(SNF_t[n, i, y] + NC_t[n, i, y] for n in nodes) <= sum(B[i, v.size, y] * v.capacity for v in versions) #TODO set cisf capacity instead of BIG
)

# CISF capacity constraint
@constraint(
    model,
    cisf_storage_capacity[i=interim_storages, y=years],
    SNF_s[i, y] + NC_s[i, y] <= sum(B[i, v.size, y] * v.capacity for v in versions)
)

@constraint(
    model,
    cisf_logic[d = interim_storages, y = minimum(years):(maximum(years) - 1), v=versions],
    B[d, v.size, y] <= B[d, v.size, y + 1]
)

@constraint(
    model,
    cisf_init[d = interim_storages, v = versions, y = minimum(years):(minimum(years) + build_time)],
    B[d, v.size, y] == 0
)

@constraint(
    model,
    number_cisf_built,
    sum(B[d, v.size, maximum(years)] for v in versions, d in interim_storages) == 5
)

@constraint(
    model,
    cisf_build_only_once[d = interim_storages, y = years],
    sum(B[d, v.size, y] for v in versions) <= 1
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

# Variablen und deren Werte ausgeben
for v in all_variables(model)
    if value(v) != 0
            println("Variable: ", v, " = ", value(v))
    end
end


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

df_Bi = DataFrame(cisf = String[], size=String[], year=Int[], build = Int[])
# for i in interim_storages
#     append!(df_Bi, DataFrame(cisf = i, build = round(value.(B[i]))))
# end

for i in interim_storages, y in years, v in versions
    # println(i, y)
    # println(round(value.(B[i,y])))
    append!(df_Bi, DataFrame(cisf=i, size=v.size, year=y, build=round(value.(B[i,v.size,y]))))
end

CSV.write("cisf_build.csv", df_Bi)

df_to_end = DataFrame(node = String[], NC = Int[])
for i in nodes
    append!(df_to_end, DataFrame(node = i,NC = round(value.(NC_t_end[i]))))
end

CSV.write("end_transport.csv", df_to_end)

println("Done!")

for y in years
    println(DISCOUNT ^ (y - minimum(years)) * CISF_BUILDING_COSTS * sum(value.(B[i, y]) - value.(B[i, y - 1]) for i in interim_storages))
end

transport_eval = Dict()
reactor_eval = Dict()
hc_eval = Dict()
cisf_eval = Dict()
reactor_fix = Dict()
cisf_fix = Dict()
cisf_investment = Dict()

for y in years
    transport_eval[y] = (DISCOUNT ^ (y - minimum(years)) * sum(transport_costs[string(n, "-", m)] * (value.(SNF_t[n, m, y]) + value.(NC_t[n, m, y])) for n in nodes, m in nodes))
    reactor_eval[y] = DISCOUNT ^ (y - minimum(years)) * sum(reactor_costs[r] * value.((SNF_s[r, y]) + value.(NC_s[r, y])) for r in reactors)
    hc_eval[y] = DISCOUNT ^ (y - minimum(years)) * HOT_CELL_COSTS * sum(value.(SNF_t[r, hc, y]) for r in reactors, hc in hot_cells)
    cisf_eval[y] = DISCOUNT ^ (y - minimum(years)) * CISF_OPERATING_COSTS * sum( value.(SNF_s[i, y]) + value.(NC_s[i, y]) for i in interim_storages)
    reactor_fix[y] = DISCOUNT ^ (y - minimum(years)) * FIX_COST_RATE * sum(value.(A[r, y]) for r in reactors)
    cisf_fix[y] = DISCOUNT ^ (y - minimum(years)) * FIX_COST_RATE * sum(value.(B[i, y]) for i in interim_storages)
    cisf_investment[y] = DISCOUNT ^ (y - minimum(years)) * CISF_BUILDING_COSTS * sum(value.(B[i, y]) - value.(B[i, y - 1]) for i in interim_storages)
end

yearly_costs = DataFrame(year=Int[], transport=Int[], reactor=Int[], hc=Int[], cisf=Int[], reactor_fix=Int[], cisf_fix=Int[], cisf_investment=Int[])
for y in years
    append!(yearly_costs, DataFrame(year=y, transport=round(transport_eval[y]), reactor=round(reactor_eval[y]),
            hc=round(hc_eval[y]), cisf=round(cisf_eval[y]), reactor_fix=round(reactor_fix[y]), cisf_fix=round(cisf_fix[y]), cisf_investment=round(cisf_investment[y])))
end

CSV.write("yearly_costs.csv", yearly_costs)