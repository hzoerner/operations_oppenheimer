module ModelDefinition

using JuMP
using LinearAlgebra

export define_model

function define_model(model, nodes, interim_storages, versions, hot_cells, transport_costs, transport_nuclear_factor, HOT_CELL_COSTS, HC_BUILDING_COSTS)
    # Define variables
    @variable(model, SNF_t[nodes, nodes] >= 0)
    @variable(model, NC_t[nodes, nodes] >= 0)
    @variable(model, B[interim_storages, map(v -> v.size, versions)], Bin)
    @variable(model, HC[hot_cells], Bin)

    cost_factor = 1/1000

    # Define objective function
    @objective(model, Min,
        # transportation costs
        sum(cost_factor * transport_nuclear_factor * transport_costs[string(n, "-", m)] * (SNF_t[n, m] + NC_t[n, m]) for n in nodes, m in nodes) +
        # hot cell repacking costs
        cost_factor * HOT_CELL_COSTS * sum(SNF_t[n, hc] for n in nodes, hc in hot_cells) + 
        # CISF construction costs
        sum(v.costs * B[d, v.size] for v in versions, d in interim_storages) +
        # hot cell construction costs
        HC_BUILDING_COSTS * sum(HC[hc] for hc in hot_cells)
    )

    # Define constraints
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
    
    return model
end