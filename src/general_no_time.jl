### this file is supposed to contain the general model without any specifications to a plan case
using Pkg
# Pkg.add("JuMP")
# Pkg.add("GLPK")
# Pkg.add("HiGHS")
# Pkg.add("DataFrames")
# Pkg.add("CSV")
# Pkg.add("XLSX")
using JuMP
using DataFrames
using CSV

using HiGHS

include("utils.jl")

path = ""

if isempty(ARGS)
    path = "C:/Users/ACER/Desktop/Uni/9_WiSe_23_24/OR-INF/operations_oppenheimer/operations_oppenheimer/data/ExtendedNuclearData.xlsx" #individual absolute filepath
else
    path = ARGS[1]
end


transport_details = get_transport_details(path)

transport_costs = transport_details[2]
TRANSPORT_POSSIBLE = transport_details[1]

production_capacities = get_hot_cell_capacities(path)
hot_cells = keys(production_capacities)

storage_cap = get_storage_capacities(path)
nodes = keys(storage_cap)

cisf_costs = get_cisf_costs(path)
interim_storages = keys(cisf_costs)

general = get_general_data(path)

versions = [Version("small", 150, 96444026.38), Version("medium", 375, 142844431.51), Version("large", 500, 189244836.64)]

HC_BUILDING_COSTS = 300000000
FIX_COST_RATE = 100000
HOT_CELL_COSTS = 2000000 # just for model logic reasons
snf = get_snf_at_reactors(path)
reactors = keys(snf)

reactor_costs = Dict()

for r in reactors
    reactor_costs[r] = REACTOR_OPERATING_COSTS
end

INITIAL_VOLUME = sum(values(snf))

#create model and attach solver
model = Model()
set_optimizer(model, HiGHS.Optimizer)

#create all variables we need
@variable(model, SNF_t[nodes, nodes] >= 0)
@variable(model, NC_t[nodes, nodes] >= 0)
@variable(model, B[interim_storages, map(v -> v.size, versions)], Bin)
@variable(model, HC[hot_cells], Bin)

cost_factor = 1/1000

@objective(model, Min,
        # transportation costs
        sum(cost_factor * transport_costs[string(n, "-", m)] * (SNF_t[n, m] + NC_t[n, m]) for n in nodes, m in nodes) +
        # hot cell repacking costs
        cost_factor * HOT_CELL_COSTS * sum(SNF_t[n, hc] for n in nodes, hc in hot_cells) + 
        # CISF construction costs
        sum(v.costs * B[d, v.size] for v in versions, d in interim_storages) +
        # hot cell construction costs
        HC_BUILDING_COSTS * sum(HC[hc] for hc in hot_cells)
)

 ### Hot Cell constraints
 @constraint(
     model, 
     hc_production_cap[hc = hot_cells], 
     sum(SNF_t[n, hc] for n in nodes) <= HC[hc] * 2000
 )

 # No new casks to hot cell
 @constraint(
     model, 
     no_casks_to_hc[hc = hot_cells], 
     sum(NC_t[n, hc] for n in nodes) == 0)

# no "old" cask from hot cell
 @constraint(
     model, 
     no_snf_from_hc[hc = hot_cells], 
     sum(SNF_t[hc, n] for n in nodes) == 0
)

# hot cell mass balance
@constraint(
    model, 
    mass_balance_hc[hc = hot_cells], 
    sum(SNF_t[n, hc] for n in nodes) == sum(NC_t[hc, n] for n in nodes)
)

@constraint(
    model,
    hot_cell_count,
    sum(HC[hc] for hc in hot_cells) <= 5
)

### transport capacity restrictions
@constraint(
    model, 
    transport_possibility[n = nodes, m = nodes], 
    SNF_t[n, m] + NC_t[n, m] <= TRANSPORT_POSSIBLE[string(n, "-", m)] * 5000
)

#NC to CISF
@constraint(
    model, 
    target_condition, 
    sum(NC_t[n, i] for n in nodes, i in interim_storages) == INITIAL_VOLUME
)

# CISF related constraints
@constraint(
    model, 
    cisf_build_before_store[i = interim_storages], 
    sum(SNF_t[n, i] + NC_t[n, i] for n in nodes) <= sum(B[i, v.size] * v.capacity for v in versions)
)

@constraint(
    model,
    number_cisf_built,
    sum(B[d, v.size] for v in versions, d in interim_storages) == 5
)

@constraint(
    model,
    cisf_build_only_once[d = interim_storages],
    sum(B[d, v.size] for v in versions) <= 1
)

@constraint(
    model,
    no_snf_to_cisf,
    sum(SNF_t[n, i] for n in nodes, i in interim_storages) == 0
)

# initialize snf at reactors
@constraint(
    model, 
    initial_snf[n = nodes], 
    sum(SNF_t[n, m] for m in nodes) == get_snf_at_node(snf, n)
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

df_SNF_t = DataFrame(from = String[], to = String[], SNF = Float32[])
df_NC_t = DataFrame(from = String[], to = String[], NC = Float32[])

for n in nodes, m in nodes, y in years
    append!(df_SNF_t, DataFrame(from = n, to = m, SNF = round(value.(SNF_t[n,m]))))
    append!(df_NC_t, DataFrame(from = n, to = m, NC = round(value.(NC_t[n,m]))))
end

CSV.write("snf_shipped_static.csv", df_SNF_t)
CSV.write("nc_shipped_static.csv", df_NC_t)

df_Bi = DataFrame(cisf = String[], size=String[], build = Int[])

for i in interim_storages, v in versions
    # println(i, y)
    # println(round(value.(B[i,y])))
    append!(df_Bi, DataFrame(cisf=i, size=v.size, build=round(value.(B[i,v.size]))))
end

CSV.write("cisf_build.csv", df_Bi)

df_Ai = DataFrame(reactor = String[], year=Int[], build = Int[])
for r in reactors, y in years
    append!(df_Ai, DataFrame(reactor=r, year=y, build=round(value.(A[r,y]))))
end
CSV.write("reactor_operating.csv", df_Ai)

df_to_end = DataFrame(node = String[], NC = Int[])
for i in nodes
    append!(df_to_end, DataFrame(node = i,NC = round(value.(NC_t_end[i]))))
end

CSV.write("end_transport.csv", df_to_end)

println("Done!")