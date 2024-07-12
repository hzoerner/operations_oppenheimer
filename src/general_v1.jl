### this file is supposed to contain the general model without any specifications to a plan case

using JuMP
using DataFrames
using CSV

using HiGHS

include("utils.jl")

path = "C:/Users/ACER/Desktop/Uni/9_WiSe_23_24/OR-INF/operations_oppenheimer/operations_oppenheimer/data/NuclearDummyData.xlsx"

years = 2030:2050

transport_details = get_transport_details(path)
transport_costs = transport_details[2]

hc, cisfs, districts = get_hot_cells_and_cisf(path)

production_capacities = get_hot_cell_capacities(path)
hot_cells = keys(production_capacities)

snf = get_snf_at_reactors(path)
reactors = keys(snf)

storage_cap = get_storage_capacities(path)
nodes = vcat(districts, collect(reactors))

cisf_costs = get_cisf_costs(path)
interim_storages = keys(cisf_costs)

general = get_general_data(path)

versions = [Version("small", 150, 50000000), Version("medium", 375, 100000000), Version("large", 500, 150000000)]

HC_BUILDING_COSTS = 300000000
CISF_BUILDING_COSTS = general["CISF_building_costs"]
CISF_OPERATING_COSTS = general["CISF_operation_costs"]
REACTOR_OPERATING_COSTS = general["Reaktor_operation_costs"]
FIX_COST_RATE = 100000
HOT_CELL_COSTS = 100 # just for model logic reasons
build_time = 2

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
@variable(model, SNF_last_mile[districts, ["hc", "cisf"], years] >= 0)
@variable(model, SNF_first_mile[ ["hc", "cisf"], districts, years] >= 0)
@variable(model, SNF_s[nodes, (minimum(years) - 1):maximum(years)] >= 0)
@variable(model, NC_t[nodes, nodes, years] >= 0)
@variable(model, NC_last_mile[districts, ["hc", "cisf"], years] >= 0)
@variable(model, NC_first_mile[ ["hc", "cisf"], districts, years] >= 0)
@variable(model, NC_s[nodes, years] >= 0)
@variable(model, B[districts, map(v -> v.size, versions), minimum(years)-1:maximum(years)], Bin)
@variable(model, A[reactors, years], Bin)
@variable(model, HC[districts, minimum(years)-1:maximum(years)], Bin)

cost_factor = 1/1000

@objective(model, Min, 
    sum(
        DISCOUNT ^ (y - minimum(years)) * sum(cost_factor * transport_costs[string(n, "-", m)] * (SNF_t[n, m, y] + NC_t[n, m, y]) for n in nodes, m in nodes) + 
        DISCOUNT ^ (y - minimum(years)) * sum(cost_factor * reactor_costs[r] * (SNF_s[r, y] + NC_s[r, y]) for r in reactors) +
        DISCOUNT ^ (y - minimum(years)) * cost_factor * HOT_CELL_COSTS * sum(SNF_last_mile[d, "hc", y] for d in districts) +
        DISCOUNT ^ (y - minimum(years)) * cost_factor * CISF_OPERATING_COSTS * sum( SNF_s[d, y] + NC_s[d, y] for d in districts) +
        DISCOUNT ^ (y - minimum(years)) * cost_factor * FIX_COST_RATE * sum(A[r, y] for r in reactors) +
        DISCOUNT ^ (y - minimum(years)) * cost_factor * FIX_COST_RATE * sum(B[d, v.size, y] for v in versions, d in districts) for y in years ) + 
    sum(
        DISCOUNT ^ (y - minimum(years)) * cost_factor * sum(v.costs * B[d, v.size, y] - B[d, v.size, y - 1] for v in versions, d in districts) +
        DISCOUNT ^ (y - minimum(years)) * cost_factor * HC_BUILDING_COSTS * sum(HC[d, y] - HC[d, y - 1] for d in districts) for y in years )
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
    distribution_SNF[d = districts, y = years],
    sum(SNF_t[n, d, y] for n in nodes) == SNF_last_mile[d, "hc", y] + SNF_last_mile[d, "cisf", y]
)

@constraint(
    model,
    collection_SNF[d = districts, y = years],
    sum(SNF_t[d, n, y] for n in nodes) == SNF_first_mile["hc", d, y] + SNF_first_mile["cisf", d, y]
)

@constraint(
    model,
    distribution_NC[d = districts, y = years],
    sum(NC_t[n, d, y] for n in nodes) == NC_last_mile[d, "hc", y] + NC_last_mile[d, "cisf", y]
)

@constraint(
    model,
    collection_NC[d = districts, y = years],
    sum(NC_t[d, n, y] for n in nodes) == NC_first_mile["hc", d, y] + NC_first_mile["cisf", d, y]
)

@constraint(
    model,
    no_waste_on_the_road[y = years],
    INITIAL_VOLUME == sum(NC_s[n, y] + SNF_s[n, y] for n in nodes)
)

# hot cell shit here (Kosten für Umladen notwendig?)
@constraint(
    model, 
    mass_balance_hc[hc = interim_storages, y = years], 
    sum(SNF_last_mile[d, "hc", y] for d in districts) == sum(NC_first_mile["hc", d, y] for d in districts)
)

 # hot cell processing capacities
 #@constraint(
 #    model, 
 #    hc_production_cap[hc = hot_cells, y = years], 
 #    sum(SNF_t[n, hc, y] for n in nodes) <= production_capacities[hc]
 #)

 # No new casks to hot cell
 @constraint(
     model, 
     no_casks_to_hc[d = districts, y = years], 
     NC_last_mile[d, "hc", y] == 0
     )

# no "old" cask from hot cell
 @constraint(
     model, 
     no_snf_from_hc[d = districts, y = years], 
     SNF_first_mile["hc", d, y] == 0
)

# dont refill before you build
@constraint(
    model, 
    hc_build_before_refill[d = districts, y = years], 
    SNF_last_mile[d, "hc", y] <= HC[d, y] * 500
)

# storage capacities
#@constraint(
#    model, 
#    storage[n = nodes, y = years], 
#    SNF_s[n, y] + NC_s[n, y] <= storage_cap[n]
#) # are capacities fix?

# transport is possible
#TODO: (Gorleben,Hot Cell 1/2,), (Ahaus, HotCell 1/2) rausnehmen?
@constraint(
    model, 
    transport_possibility[n = nodes, m = nodes, y = years], 
    SNF_t[n, m, y] + NC_t[n, m, y] <= is_transport_possible(n, m)  * 50 #TODO set BIG as small as possible
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

# dont store before you build
@constraint(
    model, 
    cisf_build_before_store[d = districts, y = years], 
    sum(SNF_last_mile[d, "cisf", y] + NC_last_mile[d, "cisf", y] for n in nodes) <= sum(B[d, v.size, y] * v.capacity for v in versions)
)

@constraint(
    model,
    cisf_build_only_once[d = districts, y = years],
    sum(B[d, v.size, y] for v in versions) <= 1
)

@constraint(
    model,
    cisf_logic[d = districts, y = minimum(years):(maximum(years) - 1)],
    sum(B[d, v.size, y] for v in versions) <= sum(B[d, v.size, y + 1] for v in versions)
)

@constraint(
    model,
    number_cisf_built,
    sum(B[d, v.size, maximum(years)] for v in versions, d in districts) == 3
)

@constraint(
    model,
    cisf_init[d = districts, v = versions],
    B[d, v.size, (minimum(years) + build_time)] == 0
)

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

# cask decay
#@constraint(
#    model, 
#    cask_decay[y = (minimum(years)+1):maximum(years), 
#    n = nodes; n ∉ hot_cells], 
#    sum(SNF_t[n, hc, y] for hc in hot_cells) >= CASK_DECAY * SNF_s[n, y - 1]
#)

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

df_Bi = DataFrame(cisf = String[], year=Int[], build = Int[])
# for i in interim_storages
#     append!(df_Bi, DataFrame(cisf = i, build = round(value.(B[i]))))
# end

for i in interim_storages, y in years
    # println(i, y)
    # println(round(value.(B[i,y])))
    append!(df_Bi, DataFrame(cisf=i, year=y, build=round(value.(B[i,y]))))
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