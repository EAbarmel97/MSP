using LinearAlgebra
using DataFrames, CSV
using Statistics

function _load_data_matrix(file_path::String,normalize=true::Bool)::Matrix{Float64}
    df = DataFrames.DataFrame(CSV.File(file_path))
    data = Matrix{Float64}(df)
    if normalize
       data .- mean.(eachcol(data))
    end

    return data
end

function correlation_matrix(m::Matrix{Float64})::Matrix{Float54}
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
    rem = size[1] % l
    wcm = zeros(size(cm))
    if rem != 0
        @warn "one of the $(size[1]) sub matrices will be windowend with $rem observations"
    end
    
    row_partition = row_partition(size(cm)[1], l)
    for i in eachindex(row_partition)
        cm = correlation_matrix(m[row_partition[i]:row_partition[i+1],:])
        wcm .+= cm 

        if i == length(row_partition)
           cm = m[row_partition[i-1]:row_partition[i],:]
        end
        wcm .+= cm 
    end
    
    return 1/div(r,l) .* wcm
end