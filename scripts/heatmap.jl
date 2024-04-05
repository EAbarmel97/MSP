using Plots
using LaTeXStrings

include("eigen_analysis.jl")
include("wcm.jl")

"""
    heatmap_plot(dir_to_save::String, corr_mat::Vector{Float64})

Generate a heatmap plot of a correlation matrix and save it to a directory.

# Arguments
- `dir_to_save::String`: Directory path to save the plot.
- `cm::Matrix{Float64}`: Correlation matrix.

# Example
```julia
heatmap_plot("path/to/save/directory", correlation_matrix)```
"""
function heatmap_plot(dir_to_save::String, cm::Matrix{Float64})
    full_file_path = joinpath(dir_to_save,"heatmap_plot.pdf")
    if !isfile(full_file_path)
        #plot styling
            #cols = Symbol.(names(df))
        (n,m) = size(cm)
         display(
            heatmap(cm, 
                fc = cgrad([:white,:dodgerblue4]),
                xticks = (1:m,m),
                xrot= 90,
                size= (800, 800),
                yticks = (1:m,m),
                yflip=true))
            display(
            annotate!([(j, i) for i in 1:n for j in 1:m])
        )
        
        #file saving
        savefig(full_file_path)
    end
end