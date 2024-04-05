using LinearAlgebra
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
"""
function filter_singular_vals_array(m::Matrix{Float64};atol=eps(Float64)::Float64)::Vector{Float64}
    singular_vals = svd(m).S  
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
compute_eigvals(matrix)
"""
function compute_eigvals(m::Matrix{Float64}; drop_first=true::Bool)::Vector{Float64}
    if drop_first 
        return abs2.(filter_singular_vals_array(m))[2:end]
    end
    
    return abs2.(filter_singular_vals_array(m))
end

"""
    intercept_and_exponent_from_log_psd(f::Vector{Float64}, psd::Vector{Float64})::Vector{Float64}

Compute the intercept and exponent from the logarithm of the power spectral density (PSD) and its corresponding frequencies using linear regression.

# Arguments
- `f::Vector{Float64}`: Vector of frequencies.
- `psd::Vector{Float64}`: Vector of power spectral density values.

# Returns
- `params::Vector{Float64}`: Vector containing the intercept and exponent.

# Example
```julia
frequencies = [0.1, 1.0, 10.0]
psd_values = [0.01, 0.1, 1.0]
intercept_and_exponent_from_log_psd(frequencies, psd_values
"""

function intercept_and_exponent_from_log_psd(f::Vector{Float64},psd::Vector{Float64})::Vector{Float64}
    X = hcat(ones(length(x)),x)

    return inv(X'*X)*(X'*y)
end

"""
    intercept_and_exponent_from_log_psd(f::Vector{Float64}, psd::Vector{Float64})::Vector{Float64}

Compute the intercept and exponent from the logarithm of the power spectral density (PSD) and its corresponding frequencies.

# Arguments
- `f::Vector{Float64}`: Vector of frequencies.
- `psd::Vector{Float64}`: Vector of power spectral density values.

# Returns
- `intercept_and_exponent::Vector{Float64}`: Vector containing the intercept and exponent.

# Example
```julia
frequencies = [0.1, 1.0, 10.0]
psd_values = [0.01, 0.1, 1.0]
intercept_and_exponent_from_log_psd(frequencies, psd_values)
"""
function intercept_and_exponent_from_log_psd(f::Vector{Float64},psd::Vector{Float64})::Vector{Float64}
    log10_f = log10.(f)
    log10_psd = log10.(psd)
    beta0, beta1 = intercept_and_exponent(log10_f,log10_psd)

    return [beta0,beta1]
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
compute_linear_fit_params(eigenvalues)
"""
function compute_linear_fit_params(eigvals::Array{Float64,1})::Vector{Float64}
    return intercept_and_exponent_from_log_psd(collect(Float64,1:length(eigvals)),eigvals)
end