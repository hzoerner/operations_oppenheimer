### this file is supposed to contain the general model without any specifications to a plan case

using JuMP
using GLPK

nodes = ["Gorleben", "Finkenkrug", "Bitterfeld-Wolfen"]
i = nodes

years = 2030:2100

#create model and attach solver
model = Model()
set_optimizer(model, GLPK.Optimizer)

@variable(model, SNF_t[nodes, nodes, years] >= 0)
@variable(model, SNF_s[nodes, years] >= 0)
@variable(model, NC_t[nodes, nodes, years] >= 0)
@variable(model, NC_s[nodes, years] >= 0)
@variable(model, B[i, years], Bin)

