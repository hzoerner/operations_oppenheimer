# using Filesystem

function save_storage(n_hc, n_cisf, hc_cost, nodes, years, SNF_s, NC_s)
    # Create a directory path based on the parameter combination
    dir_path = joinpath("storage", "n_hc_$(n_hc)", "n_cisf_$(n_cisf)", "hc_cost_$(hc_cost)")
    
    # Ensure the directory exists, create it if necessary
    mkpath(dir_path)
        
    df_SNF_s = DataFrame(node = String[], year = Int[], SNF = Float32[])
    df_NC_s = DataFrame(node = String[], year = Int[], NC = Float32[])

    for n in nodes, y in years
        append!(df_SNF_s, DataFrame(node = n, year = y, SNF = round(value.(SNF_s[n,y]))))
        append!(df_NC_s, DataFrame(node = n, year = y, NC = round(value.(NC_s[n,y]))))
    end

    CSV.write(joinpath(dir_path, "snf_stored.csv"), df_SNF_s)
    CSV.write(joinpath(dir_path, "nc_stored.csv"), df_NC_s)
    
    println("Storage saved to $dir_path")
end

function save_transport(n_hc, n_cisf, hc_cost, nodes, years, SNF_t, NC_t)
    # Create a directory path based on the parameter combination
    dir_path = joinpath("storage", "n_hc_$(n_hc)", "n_cisf_$(n_cisf)", "hc_cost_$(hc_cost)")
    
    # Ensure the directory exists, create it if necessary
    mkpath(dir_path)
    
    df_SNF_t = DataFrame(from = String[], to = String[], year = Int[], SNF = Float32[])
    df_NC_t = DataFrame(from = String[], to = String[], year = Int[], NC = Float32[])

    for n in nodes, m in nodes, y in years
        append!(df_SNF_t, DataFrame(from = n, to = m, year = y, SNF = round(value.(SNF_t[n,m,y]))))
        append!(df_NC_t, DataFrame(from = n, to = m, year = y, NC = round(value.(NC_t[n,m,y]))))
    end

    CSV.write(joinpath(dir_path, "snf_shipped.csv"), df_SNF_t)
    CSV.write(joinpath(dir_path, "nc_shipped.csv"), df_NC_t)
    
    println("Transport saved to $dir_path")
end

function save_cisf(n_hc, n_cisf, hc_cost, interim_storages, years, versions, B)
    # Create a directory path based on the parameter combination
    dir_path = joinpath("storage", "n_hc_$(n_hc)", "n_cisf_$(n_cisf)", "hc_cost_$(hc_cost)")
    
    # Ensure the directory exists, create it if necessary
    mkpath(dir_path)
    
    df_Bi = DataFrame(cisf = String[], size=String[], year=Int[], build = Int[])

    for i in interim_storages, y in years, v in versions
        append!(df_Bi, DataFrame(cisf=i, size=v.size, year=y, build=round(value.(B[i,v.size,y]))))
    end

    CSV.write(joinpath(dir_path, "cisf_build.csv"), df_Bi)
    
    println("CISF saved to $joinpath(dir_path, 'cisf_build.csv')")
end

function save_hc(n_hc, n_cisf, hc_cost, interim_storages, years, versions, HC)
    # Create a directory path based on the parameter combination
    dir_path = joinpath("storage", "n_hc_$(n_hc)", "n_cisf_$(n_cisf)", "hc_cost_$(hc_cost)")
    
    # Ensure the directory exists, create it if necessary
    mkpath(dir_path)
    
    df_HCi = DataFrame(cisf = String[], size=String[], year=Int[], build = Int[])

    for i in interim_storages, y in years, v in versions
        append!(df_Bi, DataFrame(cisf=i, size=v.size, year=y, build=round(value.(HC[i,v.size,y]))))
    end

    CSV.write(joinpath(dir_path, "hc_build.csv"), df_HCi)
    
    println("CISF saved to $joinpath(dir_path, 'hc_build.csv')")
end