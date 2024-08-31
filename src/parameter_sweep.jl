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

# define timeframe
years = parameter_grid["time_constraints"]["min_year"]:parameter_grid["time_constraints"]["max_year"]

# initialize transport information
transport_costs, TRANSPORT_POSSIBLE = get_transport_details(path)

# initialize locations
nodes = keys(get_storage_capacities(path))  # simplify function to only return nodes
interim_storages = get_cisfs(path)
hot_cells = get_hot_cells(path)
reactors = keys(get_snf_at_reactors(path))  # simplify function

# initialize CISF building options
versions = [Version(v["size"], v["capacity"], v["costs"]) for v in parameter_grid["cisf"]["versions"]]

# initialize snf at reactors
snf = get_snf_at_reactors(path)

# other excel data
storage_cap = get_storage_capacities(path)
general = get_general_data(path)

# obsolete?
CISF_OPERATING_COSTS = general["CISF_operation_costs"]
REACTOR_OPERATING_COSTS = general["Reaktor_operation_costs"]


# obsolete?
reactor_costs = Dict()
for r in reactors
    reactor_costs[r] = REACTOR_OPERATING_COSTS
end

INITIAL_VOLUME = sum(values(snf))


########################### variable parameters:
# hot cells
max_hc_counts = parameter_grid["hot_cells"]["max_count"]
hc_construction_costs = parameter_grid["hot_cells"]["construction_cost"]
# cisf
cisf_counts = parameter_grid["cisf"]["count"]


for n_hc in max_hc_counts, hc_cost in hc_construction_costs, n_cisf in cisf_counts

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
            cost_factor * hc_cost * sum(SNF_t[n, hc, y] for n in nodes, hc in hot_cells) +
            cost_factor * CISF_OPERATING_COSTS * sum( SNF_s[i, y] + NC_s[i, y] for i in interim_storages)
            # + cost_factor * FIX_COST_RATE * sum(A[r, y] for r in reactors)
            #  + cost_factor * FIX_COST_RATE * sum(B[d, v.size, y] for v in versions, d in interim_storages) 
            for y in years) + 

        #sum(LICENSE_EXTENSION_COSTS * A[r, minimum(years) + p * MAXIMUM_RUNTIME] for p in 1:EXTENSION_PERIOD) +
        sum(
            cost_factor * sum(v.costs * B[d, v.size, y] - B[d, v.size, y - 1] for v in versions, d in interim_storages) +
            cost_factor * hc_cost * sum(HC[d, y] - HC[d, y - 1] for d in hot_cells) for y in years )
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
        sum(SNF_t[n, hc, y] for n in nodes) <= HC[hc, y] * 2000
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

    if n_hc !== nothing
        @constraint(
            model,
            hot_cell_count,
            sum(HC[hc, maximum(years)] for hc in hot_cells) <= n_hc
        )
    end

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
        SNF_t[n, m, y] + NC_t[n, m, y] <= TRANSPORT_POSSIBLE[string(n, "-", m)] * 5000
    )

    @constraint(
        model,
        countrywide_transport_cap[y=years],
        sum(SNF_t[n, m, y] + NC_t[n, m, y] for n in nodes, m in nodes) <= 5*5000
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
        SNF_s[r, y] + NC_s[r, y] <= 5000 * A[r, y]
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
        cisf_logic[d = interim_storages, v=versions],
        B[d, v.size, minimum(years)] <= B[d, v.size, minimum(years) + 1]
    )

    @constraint(
        model,
        cisf_init[d = interim_storages, v = versions],
        B[d, v.size, minimum(years)] == 0
    )

    if n_cisf !== nothing
        @constraint(
            model,
            number_cisf_built,
            sum(B[d, v.size, maximum(years)] for v in versions, d in interim_storages) <= n_cisf
        )
    end

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
    # Check if a solution was found
    if termination_status(model) == MOI.OPTIMAL || termination_status(model) == MOI.LOCALLY_SOLVED
        # A solution was found, save the relevant data
        println("Optimal solution found. Saving data...")
        save_storage(n_hc, n_cisf, hc_cost, nodes, years, SNF_s, NC_s)
        save_transport(n_hc, n_cisf, hc_cost, nodes, years, SNF_t, NC_t)
        save_cisf(n_hc, n_cisf, hc_cost, interim_storages, years, versions, B)
    else
        println("No optimal solution found.")
    end
end