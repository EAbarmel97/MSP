using LinearAlgebra
using Base.Threads
using DataFrames
using GLM 

include("wcm.jl")

"""
    filter_singular_vals_array(m::Matrix{Float64}; atol=eps(Float64)::Float64)::Vector{Float64}

Filter the singular values of a matrix `m` based on a tolerance.

# Arguments
- `m::Matrix{Float64}`: Input matrix.
- `atol::Float64=eps(Float64)`: Absolute tolerance to filter singular values.

# Returns
- `filtered_singular_vals::Vector{Float64}`: Filtered singular values.

# Example
```julia
matrix = [1.0 2.0; 3.0 4.0]
filter_singular_vals_array(matrix)
```
"""
function filter_singular_vals_array(m::Matrix{Float64}; atol=eps(Float64)::Float64)::Vector{Float64}
    singular_vals = LinearAlgebra.svd(m).S  
    return filter(u -> u > atol, singular_vals)
end

"""
    compute_eigvals(m::Matrix{Float64}; drop_first=true::Bool)::Vector{Float64}

Compute the eigenvalues of a square matrix `m`.

# Arguments
- `m::Matrix{Float64}`: Input matrix.
- `drop_first::Bool=true`: Whether to drop the first eigenvalue.

# Returns
- `eigvals::Vector{Float64}`: Vector of eigenvalues.

# Example
```julia
matrix = [1.0 2.0; 3.0 4.0]
compute_eigvals(matrix)```
"""
function compute_eigvals(m::Matrix{Float64}; drop_first=true::Bool)::Vector{Float64}
    eigvals = Float64[]
    if drop_first 
        eigvals = abs2.(filter_singular_vals_array(m))[2:end] 
    end

    eigvals = abs2.(filter_singular_vals_array(m))
    return eigvals
end

"""
    compute_average_eigvals(m::Matrix{Float64}, l::Int64; drop_first=true::Bool)::Vector{Float64}

Compute the average eigenvalues of a matrix `m` using a windowed approach.

# Arguments
- `m::Matrix{Float64}`: Input matrix.
- `l::Int64`: Window size.
- `drop_first::Bool=true`: Whether to drop the first eigenvalue.

# Returns
- `average_eigvals::Vector{Float64}`: Vector containing the average eigenvalues.

# Example
```julia
matrix = [1.0 2.0 3.0; 4.0 5.0 6.0; 7.0 8.0 9.0]
compute_average_eigvals(matrix, 2)```
"""
function compute_average_eigvals(m::Matrix{Float64},l::Int64; drop_first=true::Bool)::Vector{Float64}
    rem = size(m)[1] % l
    spectrum = zeros(size(m)[2])
    if rem != 0
        @warn "one of the $(size(m)[1]) sub matrices will be windowed with $rem observations"
    end
    
    rp = row_partition(size(m)[1], l)
    Threads.@threads for i in eachindex(rp)
        if i == 1
            singular_vals = LinearAlgebra.svd(m[1:rp[1],:]).S
            spectrum .+= abs2.(singular_vals)
        
        else
            singular_vals = LinearAlgebra.svd(m[rp[i-1]+1:rp[i],:]).S
            spectrum .+= abs2.(singular_vals)
        end
    end
    
    if drop_first
        return 1/(div(size(m)[1],l)*(l-1)) .* spectrum[2:end]
    end

    return 1/(div(size(m)[1],l)*(l-1)) .* spectrum
end

"""
    intercept_and_exponent(x::Vector{Float64}, y::Vector{Float64})::Vector{Float64}

Compute the intercept and exponent for a linear fit between `x` and `y`.

# Arguments
- `x::Vector{Float64}`: Vector of independent variable values.
- `y::Vector{Float64}`: Vector of dependent variable values.

# Returns
- `params::Vector{Float64}`: Vector containing the intercept and exponent.

# Example
```julia
x_values = [1.0, 2.0, 3.0]
y_values = [0.1, 0.2, 0.3]
intercept_and_exponent(x_values, y_values)```
"""
function intercept_and_exponent(x::Vector{Float64},y::Vector{Float64})::Vector{Float64}
   data = DataFrame(X=x,Y=y)
   return GLM.coef(GLM.lm(@formula(Y ~ X), data))
end

"""
    compute_linear_fit_params(eigvals::Array{Float64,1})::Vector{Float64}

Compute the parameters of a linear fit for a given array of eigenvalues.

# Arguments
- `eigvals::Array{Float64,1}`: Array of eigenvalues.

# Returns
- `params::Vector{Float64}`: Vector containing the intercept and exponent of the linear fit.

# Example
```julia
eigenvalues = [0.1, 0.2, 0.3, 0.4]
compute_linear_fit_params(eigenvalues)```
"""
function compute_linear_fit_params(eigvals::Array{Float64,1})::Vector{Float64}
    log10_n = log10.(collect(Float64,1:length(eigvals)))
    log10_eigvals = log10.(eigvals)
    beta0, beta1 = intercept_and_exponent(log10_n,log10_eigvals)
    #@show [beta0, beta1]
    return [beta0, beta1]
end