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

function get_conditional_NC_s(model, node, year, f::Function)
    if f(year)
        return model[:NC_s][node, year - 1]
        else
            return 0.0
        end
end