using Statistics

"""
    correlation_matrix(m::Matrix{Float64})::Matrix{Float64}

Compute the correlation matrix of a given matrix `m`.

# Arguments
- `m::Matrix{Float64}`: Input matrix.

# Returns
- `correlation::Matrix{Float64}`: Correlation matrix.

# Example
```julia
matrix = [1.0 2.0; 3.0 4.0]
correlation_matrix(matrix)
```
"""
function correlation_matrix(m::Matrix{Float64})::Matrix{Float64}
    return Statistics.cor(m)
end

"""
    row_partition(r::Int64, l::Int64)::Vector{Int64}

Partition the rows of a matrix into equal-sized blocks.

# Arguments
- `r::Int64`: Total number of rows.
- `l::Int64`: Size of each block.

# Returns
- `row_partition::Vector{Int64}`: Vector containing the partitioned row indices.

# Example
```julia
row_partition(10, 3)
```
"""
function row_partition(r::Int64,l::Int64)::Vector{Int64}
    row_partition = map(collect(1:div(r,l))) do u 
        return l*u
    end

    if r % l != 0
        push!(row_partition, r)
        return row_partition
    end

    return row_partition
end

"""
    windowed_correlation_matrix(m::Matrix{Float64}, l::Int64)::Matrix{Float64}

Compute the windowed correlation matrix of a matrix `m` with a window size `l`.

# Arguments
- `m::Matrix{Float64}`: Input matrix.
- `l::Int64`: Window size.

# Returns
- `wcm::Matrix{Float64}`: Windowed correlation matrix.

# Example
```julia
matrix = [1.0 2.0 3.0; 4.0 5.0 6.0; 7.0 8.0 9.0]
windowed_correlation_matrix(matrix, 2)
```
"""
function windowed_correlation_matrix(m::Matrix{Float64},l::Int64)::Matrix{Float64}
    rem = size(m)[1] % l
    wcm = zeros(size(m)[2],size(m)[2])
    if rem != 0
        @warn "one of the $(size(m)[1]) sub matrices will be windowed with $rem observations"
    end
    
    rp = row_partition(size(m)[1], l)
    for i in eachindex(rp)
        if i == 1
            cm = correlation_matrix(m[1:rp[1],:])
            wcm .+= cm
        else
            cm = correlation_matrix(m[rp[i-1]+1:rp[i],:])
            wcm .+= cm    
        end                 
    end
    
    return 1/div(size(m)[1],l) .* wcm
end