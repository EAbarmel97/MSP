using LinearAlgebra
using DataFrames, CSV
using Plots 
using Statistics

function _load_data_matrix(file_path::String; normalize=true::Bool)::Matrix{Float64}
    df = DataFrames.DataFrame(CSV.File(file_path;header=false))
    data = Matrix{Float64}(df)
    if normalize
        data .-= mean(data, dims=1)
    end

    return data
end

function correlation_matrix(m::Matrix{Float64})::Matrix{Float64}
    return Statistics.cor(m)
end

function row_partition(r::Int64,l::Int64)::Vector{Int64}
    number_of_windows = div(r,l)

    row_partition = map(collect(1:div(r,l))) do u 
        return l*u
    end

    if r % l != 0
        push!(row_partition, r)
        return row_partition
    end

    return row_partition
end

function _windowed_correlation_matrix(m::Matrix{Float64},l::Int64)::Matrix{Float64}
    rem = size(m)[1] % l
    wcm = zeros(size(cm)...)
    if rem != 0
        @warn "one of the $(size(m)[1]) sub matrices will be windowed with $rem observations"
    end
    
    rp = row_partition(size(m)[1], l)
    for i in eachindex(rp)
        if i == 1
            cm = correlation_matrix(m[1:rp[1],:])
            wcm .+= cm
        else
            cm = correlation_matrix(m[rp[i-1]:rp[i],:])
            wcm .+= cm    
        end                 
    end
    
    return 1/div(size(m)[1],l) .* wcm
end
