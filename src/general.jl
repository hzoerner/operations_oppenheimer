### this file is supposed to contain the general model without any specifications to a plan case

using JuMP
using GLPK

path = "C:/Users/ACER/Documents/Downloads/NuclearData.xlsx"

years = 2031:2101

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
@objective(model, Min, sum(sum(transport_costs[string(n, "-", m)] * (SNF_t[n, m, y] + NC_t[n, m, y]) for n in nodes, m in nodes) + sum(SNF_s[n, y] + NC_s[n, y] for n in nodes) + sum(B[i, y] for i in interim_storages) + sum(SNF_t[n, hc, y] for hc in hot_cells, n in nodes) for y in years) + sum(NC_t_end[i] for i in interim_storages))

# mass balances
@constraint(model, mass_balance_SNF[n = nodes, m = nodes, y = years], sum(SNF_t[m, n, y]) - sum(SNF_t[n, m, y]) - SNF_s[n, y] + SNF_s[n, y - 1] == 0) #TODO add if y > 2030 variable
@constraint(model, mass_balance_NC[n = nodes, n ∉ hot_cells, m = nodes, y = years], sum(NC_t[m, n, y]) - sum(NC_t[n, m, y]) + NC_s[n, y] == 0) #TODO see above
# TODO can we move casks from a cisf (before moving to end storage)?
# TODO can we store SNF in a cisf?

# storage capacities cisf
@constraint(model, storage[n = nodes, y = years], SNF_s[n, y] + NC_s[n, y] <= storage_cap[n]) # are capacities fix?

# transport is possible
@constraint(model, transport_possibility[n = nodes, m = nodes, y = years], SNF_t[n, m, y] + NC_t[n, m, y] <= TRANSPORT_POSSIBLE[string(n, "-", m)] * BIG_trans)

# transport to end storage facility / clear condition
@constraint(model, transport_to_end[n = nodes], NC_t_end[n] == NC_s[n, maximum(years)])

# contraint for building cisf only once
@constraint(model, cisf_build[i = interim_storages], sum(B[i, y] for y in years) <= 1)

# dont store before you build
@constraint(model, cisf_build_before_store[i = interim_storages, y = years], sum(SNF_s[i, cur] + NC_s[i, cur] + sum(SNF_t[n, i, cur] + SNF_t[i, n, cur] + NC_t[i, n, cur] + NC_t[n, i, cur] for n in nodes) for cur in minimum(years): y) <= sum(B[i, cur] for cur in minimum(years): y) * BIG_cisf)

# cask decay
@constraint(model, cask_decay[y = (minimum(years)+1):maximum(years)], sum(SNF_t[n, hc, y] for n in nodes, hc in hot_cells) >= CASK_DECAY * sum(SNF_s[n, y - 1] for n in nodes))

# initialize snf at reactors
@constraint(model, initial_snf[n = nodes], SNF_s[n, minimum(years) - 1] == snf[n])

optimize!(model)
obj_value = objective_value(model)
println("Total costs are ", obj_value)