using XLSX
using DataFrames

function get_snf_at_reactors(path)
    df = DataFrame(XLSX.readtable(path, "Reactors"))
    snf = Dict(df.name[i] => df.snf[i] for i in 1:nrow(df))

    return snf
end

function get_hot_cell_capacities(path)
    df = DataFrame(XLSX.readtable(path, "Hot Cells"))
    production_capacities = Dict(df.name[i] => df.production[i] for i in 1:nrow(df))

    return production_capacities
end

function get_hot_cells(path)
    df = DataFrame(XLSX.readtable(path, "Hot Cells"))
    hot_cells = df.name
    return hot_cells
end

function get_storage_capacities(path)
    storage_capacities = Dict()

    for sheet in ["Reactors", "CISF", "Hot Cells"]
        df = DataFrame(XLSX.readtable(path, sheet))
        temp_dict = Dict(df.name[i] => df.capacity[i] for i in 1:nrow(df))
        storage_capacities = merge(storage_capacities, temp_dict)
    end

    return storage_capacities
end

function get_transport_details(path)::Tuple{Dict, Dict}
    df = DataFrame(XLSX.readtable(path, "Transport"))

    costs = Dict(string(df.from[i], "-", df.to[i]) => df.costs[i] for i in 1:nrow(df))
    possible = Dict(string(df.from[i], "-", df.to[i]) => df.is_possible[i] for i in 1:nrow(df))

    return (possible, costs)
end

function get_cisf_costs(path)
    df = DataFrame(XLSX.readtable(path, "CISF"))
    costs = Dict(df.name[i] => df.costs[i] for i in 1:nrow(df))

    return costs
end

function get_cisfs(path)
    df = DataFrame(XLSX.readtable(path, "CISF"))
    return df.name
end

function get_reactor_costs(path)
    df = DataFrame(XLSX.readtable(path, "Reactors"))
    costs = Dict(df.name[i] => df.costs[i] for i in 1:nrow(df))

    return costs
end

function get_conditional_variable(node, year, f::Function, VAR)
    if f(year)
        return VAR[node, year - 1]
        else
            return 0.0
        end
end

function get_general_data(path)
    df = DataFrame(XLSX.readtable(path, "General"))
    general = Dict(df.parameter[i] => df.value[i] for i in 1:nrow(df))

    return general
end

function get_end_storage_transport_costs(path)
    df = DataFrame(XLSX.readtable(path, "End Storage"))
    costs = Dict(df.from[i] => df.costs[i] for i in 1:nrow(df))
    
    return costs
end

function get_snf_at_node(snf, node)
    if haskey(snf, node)
        return snf[node]
        end

    return 0.0
end

function get_transport_cap(n, m, cap)
    if (n,m) ∉ [("Ahaus", "Ahaus_CISF"), ("Gorleben", "Gorleben_CISF")]
        return cap
    else
        return cap * 10
    end
end

function total_transport_volume(n, m, yr, snf_t, nc_t)
    if (n,m) ∉ [("Ahaus", "Ahaus_CISF"), ("Gorleben", "Gorleben_CISF")]
        return snf_t[n, m, yr] + nc_t[n, m, yr]
    else
        return 0
    end
end

function is_transport_possible(n, m)
    if n == m
        return 0
    else
        return 1
    end
end

function get_hybrid_cap(n, y, interim_built, old_cap, new_cap)
    if B[n, y] == 0
        return old_cap[n]
    else
        return new_cap
    end
end

function get_hot_cells_and_cisf(path)
    df = DataFrame(XLSX.readtable(path, "Districts"))
    districts =df.name
    hc = Building[]
    cisf = Building[]

    for d in districts
        push!(hc, Building("hc", d))
        push!(cisf, Building("cisf", d))
    end
    
    return hc, cisf, districts
end


struct Version
    size::String
    capacity::Int32
    costs::Float32
end
