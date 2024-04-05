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

