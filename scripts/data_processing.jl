using CSV, DataFrames
using LinearAlgebra
using Base.Threads

include("wcm.jl")

function paralell_read(dir_to_fetch::String)
    Threads.@threads for file in readdir(dir_to_fetch)
        if endswith(file,".csv")
            file_name = basename(file)
            m = load_data_matrix(joinpath(dir_to_fetch,file))            
        end
    end
end