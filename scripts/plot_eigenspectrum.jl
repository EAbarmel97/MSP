using Plots
using LaTeXStrings
using Base.Threads

include("wcm.jl")
include("eigen_analysis.jl")

function plot_eigen_spectrum end

"""
    plot_eigen_spectrum(eigvals::Vector{Float64}; dir_to_save="."::String)

Plot the eigenvalue spectrum and its linear fit on a log-log scale and save the plot to a directory.

# Arguments
- `eigvals::Vector{Float64}`: Eigenvalues to be plotted.
- `dir_to_save::String`: Directory path to save the plot. Default is the current directory.

# Example
```julia
plot_eigen_spectrum(eigenvalues)```
"""
function plot_eigen_spectrum(eigvals::Vector{Float64}; dir_to_save="."::String) 
    full_file_path = joinpath(dir_to_save,"eigen_spectrum_plot.pdf")

    #compute linear fit 
    params = compute_linear_fit_params(eigvals)
    
    #persist graph if doesn't exist
    if !isfile(full_file_path)
        # plot styling
        plt = plot(collect(1:length(eigvals)), eigvals, label=L"ev_n",xscale=:log10, yscale=:log10, lc=:blue)

        # linear fit
        x_vals = collect(1:length(eigvals))
        y_vals = exp10.(params[1] .+ params[2] .* log10.(x_vals))
        plot!(x_vals, y_vals, label="Linear Fit: beta = $(round(params[2], digits=3)), A = $(round(exp10(params[1]),digits=3))", lc=:red)

        title!("Eigen spectrum plot")
        xlabel!("rank")
        ylabel!("Eigen spectrum")
        
        #file saving
        savefig(plt, full_file_path)
    end
end

"""
    plot_eigen_spectrum(m::Matrix{Float64}, l::Int64; dir_to_save="."::String)

Plot the eigenvalue spectrum of a matrix `m` with a window size `l` and save the plot to a directory.

# Arguments
- `m::Matrix{Float64}`: Input matrix.
- `l::Int64`: Window size for computing average eigenvalues.
- `dir_to_save::String`: Directory path to save the plot. Default is the current directory.

# Example
```julia
plot_eigen_spectrum(matrix, 5)```
"""
function plot_eigen_spectrum(m::Matrix{Float64}, l::Int64; dir_to_save="."::String)
    ensamble_params = zero(2)

    rp = row_partition(size(m)[1],l)

    full_file_path = joinpath(dir_to_save,"eigen_spectra_plot.pdf")

    #compute linear fit 
    average_eig_spectrum = compute_average_eigvals(m,l)
    params = compute_linear_fit_params(average_eig_spectrum)
    
    #persist graph if doesn't exist
    if !isfile(full_file_path)
        # plot styling
        plt = plot(collect(1:length(average_eig_spectrum)), average_eig_spectrum, label=L"ev_n",xscale=:log10, yscale=:log10, lc=:blue)

        # linear fit
        x_vals = collect(1:length(average_eig_spectrum))
        y_vals = exp10.(params[1] .+ params[2] .* log10.(x_vals))
        plot!(x_vals, y_vals, label="Linear Fit: beta = $(round(params[2], digits=3)), A = $(round(exp10(params[1]),digits=3))", lc=:red)

        Threads.@threads for i in eachindex(rp)
            if i == 1
                ev = compute_eigvals(m[1:rp[1],:])
                plot!(collect(1:length(ev)), ev,xscale=:log10, yscale=:log10,alpha=0.2, lc=:grey)
            else
                ev = compute_eigvals(m[rp[i-1]+1:rp[i],:])
                plot!(collect(1:length(ev)), ev,xscale=:log10, yscale=:log10,alpha=0.2, lc=:grey) 
            end
        end
        
        # linear fit
        x_vals = collect(1:length(average_eig_spectrum))
        y_vals = exp10.(0.05 .- 2.3 .* log10.(x_vals))
        plot!(x_vals, y_vals,xscale=:log10, yscale=:log10, lc=:black)

        title!("Eigen spectrum plot")
        xlabel!("rank")
        ylabel!("Eigen spectrum")
        
        #file saving
        savefig(plt, full_file_path)
    end
end