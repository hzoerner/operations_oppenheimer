using Pkg
Pkg.add("YAML")
using YAML
Pkg.add("JuMP")
Pkg.add("HiGHS")
Pkg.add("DataFrames")
Pkg.add("CSV")
Pkg.add("XLSX")
# Pkg.add("Filesystem")
using JuMP
using DataFrames
using CSV

using HiGHS
include("utils.jl")
include("storage_functions.jl")
pwd()
path = ""

if isempty(ARGS)
    path = "/Users/simonirmer/Documents/Privat/Uni/Berlin/WS23-24/OR-INF/term_paper/operations_oppenheimer/data/ExtendedNuclearData.xlsx" #individual absolute filepath
else
    path = ARGS[1]
end

# Load parameter config
parameter_grid = YAML.load_file("/Users/simonirmer/Documents/Privat/Uni/Berlin/WS23-24/OR-INF/term_paper/operations_oppenheimer/src/parameter_config.yaml")

# initialize transport information
TRANSPORT_POSSIBLE, transport_costs = get_transport_details(path)
transport_nuclear_factor = parameter_grid["transport_costs"]["cost_factor"] * parameter_grid["transport_costs"]["nuclear_factor"]
# initialize locations
nodes = keys(get_storage_capacities(path))  # simplify function to only return nodes
interim_storages = get_cisfs(path)
hot_cells = get_hot_cells(path)
reactors = keys(get_snf_at_reactors(path))  # simplify function

# initialize CISF building options
versions = [Version(v["size"], v["capacity"], v["costs"]) for v in parameter_grid["cisf"]["versions"]]

# initialize snf at reactors
snf = get_snf_at_reactors(path)

# obsolete?
HOT_CELL_COSTS = parameter_grid["hot_cells"]["processing_cost"]

INITIAL_VOLUME = sum(values(snf))


########################### variable parameters:
# hot cells
max_hc_counts = parameter_grid["hot_cells"]["max_count"]
hc_construction_costs = parameter_grid["hot_cells"]["construction_cost"]
HC_BUILDING_COSTS = hc_construction_costs[1]
# cisf
cisf_counts = parameter_grid["cisf"]["count"]


for n_hc in max_hc_counts, n_cisf in cisf_counts

    # print("Run model with max " + n_hc + " hot cells and max " + n_cisf + " cisf.")
    #create model and attach solver
    model = Model()
    set_optimizer(model, HiGHS.Optimizer)

    #create all variables we need
    #create all variables we need
    @variable(model, SNF_t[reactors, hot_cells] >= 0)
    @variable(model, NC_t[hot_cells, interim_storages] >= 0)
    @variable(model, B[interim_storages, map(v -> v.size, versions)], Bin)
    @variable(model, HC[hot_cells], Bin)

    cost_factor = 1/1000

    @objective(model, Min,
        # transportation costs
        cost_factor * transport_nuclear_factor * (sum(transport_costs[string(r, "-", hc)] * SNF_t[r, hc] for r in reactors, hc in hot_cells) + sum(transport_costs[string(hc, "-", i)] * NC_t[hc, i] for hc in hot_cells, i in interim_storages)) +
        # hot cell repacking costs
        cost_factor * HOT_CELL_COSTS * sum(SNF_t[r, hc] for r in reactors, hc in hot_cells) + 
        # CISF construction costs
        cost_factor * sum(v.costs * B[d, v.size] for v in versions, d in interim_storages) +
        # hot cell construction costs
        cost_factor * HC_BUILDING_COSTS * sum(HC[hc] for hc in hot_cells)
    )

    ### Hot Cell constraints
    @constraint(
        model, 
        hc_production_cap[hc = hot_cells], 
        sum(SNF_t[r, hc] for r in reactors) <= HC[hc] * 2000
    )

    # # No new casks to hot cell
    # @constraint(
    #     model, 
    #     no_casks_to_hc[hc = hot_cells], 
    #     sum(NC_t[n, hc] for n in nodes) == 0)

    # # no "old" cask from hot cell
    # @constraint(
    #     model, 
    #     no_snf_from_hc[hc = hot_cells], 
    #     sum(SNF_t[hc, n] for n in nodes) == 0
    # )

    # hot cell mass balance
    @constraint(
        model, 
        mass_balance_hc[hc = hot_cells], 
        sum(SNF_t[r, hc] for r in reactors) == sum(NC_t[hc, i] for i in interim_storages)
    )

    # # everything processed in hot cells
    # @constraint(
    #     model,
    #     ensure_processing,
    #     sum(NC_t[hc, n] for hc in hot_cells, n in nodes) == INITIAL_VOLUME
    # )

    # @constraint(
    #     model,
    #     ensure_nc_quantity,
    #     sum(NC_t[n, m] for n in nodes, m in nodes) == INITIAL_VOLUME
    # )

    # other nodes   ###################################################################################
    # ISF mass balance
    # @constraint(
    #     model,
    #     mass_balance_SNF[n = nodes; n ∉ hot_cells],
    #     sum(SNF_t[m, n] for m in nodes) == sum(SNF_t[n, m] for m in nodes)
    # )

    # # no cask conversion except in hot cells
    # @constraint(
    #     model, 
    #     no_hc_from_isf[i = interim_storages], 
    #     sum(NC_t[m, i] for m in nodes) == (1-B[i]) * sum(NC_t[i, m] for m in nodes)
    # )

    # # overall mass balance
    # @constraint(
    #     model,
    #     mass_balance_isf[n = nodes; n ∉ reactors],
    #     sum(SNF_t[m, n] + NC_t[m, n] for m in nodes) == sum(SNF_t[n, m] + NC_t[n, m] for m in nodes)
    # )

    #######################################################################################################
    ### transport capacity restrictions
    # @constraint(
    #     model, 
    #     transport_possibility[n = nodes, m = nodes], 
    #     SNF_t[n, m] + NC_t[n, m] <= TRANSPORT_POSSIBLE[string(n, "-", m)] * 5000
    # )

    #NC to CISF
    @constraint(
        model, 
        target_condition, 
        sum(NC_t[hc, i] for hc in hot_cells, i in interim_storages) == INITIAL_VOLUME
    )

    # CISF related constraints
    @constraint(
        model, 
        cisf_build_before_store[i = interim_storages], 
        sum(NC_t[hc, i] for hc in hot_cells) <= sum(B[i, v.size] * v.capacity for v in versions)
    )

    if n_cisf !== nothing
        @constraint(
        model,
        number_cisf_built,
        sum(B[d, v.size] for v in versions, d in interim_storages) == n_cisf
    )
    end

    if n_hc !== nothing
        @constraint(
            model,
            number_hc_built,
            sum(HC[hc] for hc in hot_cells) == n_hc
        )
    end
    
    @constraint(
        model,
        cisf_build_only_once[d = interim_storages],
        sum(B[d, v.size] for v in versions) <= 1
    )

    # @constraint(
    #     model,
    #     no_snf_to_cisf,
    #     sum(SNF_t[r, i] for n in nodes, i in interim_storages) == 0
    # )

    # initialize snf at reactors
    @constraint(
        model, 
        initial_snf[r = reactors], 
        sum(SNF_t[r, hc] for hc in hot_cells) == get_snf_at_node(snf, r)
    )

    optimize!(model)
    # Check if a solution was found
    if termination_status(model) == MOI.OPTIMAL || termination_status(model) == MOI.LOCALLY_SOLVED
        # A solution was found, save the relevant data
        println("Optimal solution found. Saving data...")
        #save_storage(n_hc, n_cisf, nodes, SNF_s, NC_s)
        save_transport(n_hc, n_cisf, reactors, hot_cells, interim_storages, SNF_t, NC_t)
        save_cisf(n_hc, n_cisf, interim_storages, versions, B)
        save_hc(n_hc, n_cisf, hot_cells, HC)
    else
        println("No optimal solution found.")
    end
end


interim_storages