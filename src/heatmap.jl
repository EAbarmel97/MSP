using Plots
using LaTeXStrings

"""
    heatmap_plot(cm::Matrix{Float64}; dir_to_save="."::String)

Generate a heatmap plot of a correlation matrix and save it to a directory.

# Arguments
- `cm::Matrix{Float64}`: Correlation matrix.
- `dir_to_save::String`: Directory path to save the plot. Default is the current directory.

# Example
```julia
heatmap_plot(correlation_matrix)```
"""
function heatmap_plot(cm::Matrix{Float64}, type::String; dir_to_save="."::String)
    full_file_path = joinpath(dir_to_save, "heatmap_plot_$type.pdf")
    if !isfile(full_file_path)
        heatmap(
            cm,
            c= :bwr,
            clim=(-1, 1),
            xticks = (1:size(cm, 2), 1:size(cm, 2)),
            xrot = 90,
            size = (500, 500),
            yticks = (1:size(cm, 1), 1:size(cm, 1)),
            yflip = true,
            aspect_ratio = 1
        )
        savefig(full_file_path)
    end
end