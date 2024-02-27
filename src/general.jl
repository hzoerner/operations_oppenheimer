### this file is supposed to contain the general model without any specifications to a plan case

using JuMP
using GLPK

interim_storages = ["Gorleben", "Finkenkrug", "Bitterfeld-Wolfen"]
hot_cells = ["Hot Cell North", "Hot Cell South"]
reactors = ["Reaktor 1", "Reaktor 2"]

nodes = [interim_storages; hot_cells; reactors]
years = 2030:2031

storage_cap = Dict(zip(nodes, [100, 25, 1000, 0, 0, 20, 20]))

BIG_cap = 1000000
BIG_cisf = BIG_cap ^2
#create model and attach solver
model = Model()
set_optimizer(model, GLPK.Optimizer)

#create all variables we need
@variable(model, SNF_t[nodes, nodes, years] >= 0)
@variable(model, SNF_s[nodes, years] >= 0)
@variable(model, NC_t[nodes, nodes, years] >= 0)
@variable(model, NC_s[nodes, years] >= 0)
@variable(model, NC_t_end[interim_storages] >= 0)
@variable(model, B[interim_storages, years], Bin)

#TODO add cost parameters to obj function!!!
@objective(model, Min, sum(sum(SNF_t[n, m, y] + NC_t[n, m, y] for n in nodes, m in nodes) + sum(SNF_s[n, y] + NC_s[n, y] for n in nodes) + sum(B[i, y] for i in interim_storages) + sum(SNF_t[n, hc, y] for hc in hot_cells, n in nodes) for y in years) + sum(NC_t_end[i] for i in interim_storages))

@constraint(model, mass_balance_SNF[n = nodes, y = years], sum(SNF_t[nodes, n, y]) - sum(SNF_t[n, nodes, y]) + SNF_s[n, y] == 0) #add if y > 2030 variable
@constraint(model, mass_balance_NC[n = nodes, n ∉ hot_cells, y = years], sum(NC_t[nodes, n, y]) - sum(NC_t[n, nodes, y]) + NC_s[n, y] == 0) #see above
@constraint(model, storage[n = nodes, y = years], SNF_s[n, y] + NC_s[n, y] <= storage_cap[n]) # are capacities fix?
@constraint(model, transport_to_end[n = nodes], NC_t_end[n] == NC_s[n, maximum(years)])

@constraint(model, cisf_build[i = interim_storages], sum(B[i, y] for y in years) <= 1)
@constraint(model, cisf_build_before_store[i = interim_storages, y = years], sum(SNF_s[i, cur] + NC_s[i, cur] + sum(SNF_t[n, i, cur] + SNF_t[i, n, cur] + NC_t[i, n, cur] + NC_t[n, i, cur] for n in nodes) for cur in minimum(years): y) <= sum(B[i, cur] for cur in minimum(years): y) * BIG_cisf)