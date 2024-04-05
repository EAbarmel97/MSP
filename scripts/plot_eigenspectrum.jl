using Plots
using LaTeXStrings

"""
    plot_eigen_spectrum(dir_to_save::String, eigvals::Vector{Float64})

Plot the eigenvalue spectrum and its linear fit on a log-log scale and save the plot to a directory.

# Arguments
- `dir_to_save::String`: Directory path to save the plot.
- `eigvals::Vector{Float64}`: Eigenvalues to be plotted.

# Example
```julia
plot_eigen_spectrum("path/to/save/directory", eigvals)
```
"""
function plot_eigen_spectrum(dir_to_save::String, eigvals::Vector{Float64})
    #build x, y axis; y being the eigenspectrum and x it's enumeration
    ploting_axes = (collect(1:length(eigvals)),eigvals)
    
    full_file_path = joinpath(dir_to_save,"eigen_spectrum_plot.pdf")

    #compute linear fit 
    params = compute_linear_fit_params(ploting_axes[2])
    
    #persist graph if doesn't exist
    if !isfile(full_file_path)
        # plot styling
        plt = plot(ploting_axes[1], ploting_axes[2], label=L"ev_n",xscale=:log10, yscale=:log10, alpha=0.2)

       # linear fit
        x_vals = collect(1:length(eigvals))
        y_vals = exp10.(params[1] .+ params[2] .* log10.(x_vals))
        plot!(x_vals, y_vals, label="Linear Fit: beta = $(round(params[2], digits=3)), A = $(round(params[1],digits=3))", lc=:red)

        title!("Eigen spectrum plot")
        xlabel!(L"n")
        ylabel!("Eigen spectrum")
        
        #file saving
        savefig(plt, full_file_path)
    end
end