using LinearAlgebra

function filter_singular_vals_array(m::Matrix{Float64};atol=eps(Float64)::Float64)::Vector{Float64}
    singular_vals = svd(m).S  
    return filter(u -> u > atol, singular_vals)
end

function compute_eigvals(m::Matrix{Float64}; drop_first=true::Bool)::Vector{Float64}
    if drop_first 
        return abs2.(filter_singular_vals_array(m))[2:end]
    end
    
    return abs2.(filter_singular_vals_array(m))
end

"""
    intercept_and_exponent(x::Vector{Float64},y::Vector{Float64})::Vector{Float64}

Returns an array with the parameter estimators for a linear fit
"""

function intercept_and_exponent_from_log_psd(f::Vector{Float64},psd::Vector{Float64})::Vector{Float64}
    X = hcat(ones(length(x)),x)

    return inv(X'*X)*(X'*y)
end

"""
    intercept_and_exponent_from_log_psd(f::Array{Float64,1},average_psd::Array{Float64,1})::Array{Float64,1}

Gives a 2 dimensional array containing the parameter estimators of a 2d linear fit
"""
function intercept_and_exponent_from_log_psd(f::Vector{Float64},psd::Vector{Float64})::Vector{Float64}
    log10_f = log10.(f)
    log10_psd = log10.(psd)
    beta0, beta1 = intercept_and_exponent(log10_f,log10_psd)

    return [beta0,beta1]
end

function compute_linear_fit_params(eigvals::Array{Float64,1})::Vector{Float64}
    return intercept_and_exponent_from_log_psd(collect(Float64,1:length(eigvals)),eigvals)
end