using Plots
using LaTeXStrings

include("eigen_analysis.jl")

"""
    plot_eigen_spectrum(dir_to_save::String, eigvals::Vector{Float64})

Plot the eigenvalue spectrum and its linear fit on a log-log scale and save the plot to a directory.

# Arguments
- `dir_to_save::String`: Directory path to save the plot.
- `eigvals::Vector{Float64}`: Eigenvalues to be plotted.

# Example
```julia
plot_eigen_spectrum("path/to/save/directory", eigvals)
"""
function plot_eigen_spectrum(dir_to_save::String, eigvals::Vector{Float64})
    #build x, y axis; y being the eigenspectrum and x it's enumeration
    ploting_axes = (collect(1:length(eigvals),eigvals))

    #compute linear fit 
    params = compute_linear_fit_params(ploting_axes[2])
    
    #persist graph if doesn't exist
    if !isfile(full_file_path)
        #plot styling
        plt = plot(ploting_axes[1],ploting_axes[2], label=L"{Eig val}_n", legend=false, xscale=:log10, yscale=:log10,alpha=0.2)
        #linear fit
        plot!(u -> exp10(params[1] + params[2]*log10(u)),minimum(ploting_axes[1]),maximum(ploting_axes[1]), xscale=:log10,yscale=:log10,lc=:red)
        
        title!("Eigen spectrum plot")
        xlabel!(L"n")
        ylabel!("Eigen spectrum")
        
        #file saving
        savefig(plt, dir_to_save)
    end
end