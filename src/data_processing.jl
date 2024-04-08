using DataFrames, CSV
include("msp.jl")
is_numeric(x) = tryparse(Float64, string(x)) !== nothing


"""
    embedd_matrix_in_row_space(m::Matrix{Float64}, row_dim::Int64)::Matrix{Float64}

Embed a matrix `m` in a row space of specified dimension.

# Arguments
- `m::Matrix{Float64}`: Input matrix.
- `row_dim::Int64`: Dimension of the row space.

# Returns
- `embedded_matrix::Matrix{Float64}`: Embedded matrix in the row space.

# Example
```julia
matrix = [1.0 2.0; 3.0 4.0]
embedd_matrix_in_row_space(matrix, 3)```
"""
function embedd_matrix_in_row_space(m::Matrix{Float64},row_dim)::Matrix{Float64}
    if size(m)[1] > row_dim
        @error "the row dimension need to be smaller"
    end
    embedded_matrix = zeros(Float64,row_dim,size(m)[2])
    
    for i in eachindex(eachrow(m))
        embedded_matrix[i,:] .= eachrow(m)[i]
    end

    return embedded_matrix    
end

"""
    load_data_matrix(file_path::String; normalize=true::Bool)::Matrix{Float64}

Load data from a CSV file into a matrix and optionally normalize it.

# Arguments
- `file_path::String`: Path to the CSV file.
- `normalize::Bool=true`: Whether to normalize the data by subtracting the mean.

# Returns
- `data::Matrix{Float64}`: Loaded data matrix.

# Example
```julia
file_path = "path/to/data.csv"
load_data_matrix(file_path)
```
"""
function load_data_matrix(file_path::String; drop_header=false,centralize::Bool=true)::Matrix{Float64}
    df = DataFrames.DataFrame(CSV.File(file_path; header=drop_header))
    data = Matrix{Float64}(df)

    if centralize
        data .-= mean(data, dims=1)
    end

    return data
end



"""
    write_csv_by_fase(dir_to_fetch::String, fase::Int64)

Combine CSV files from a directory corresponding to a specific phase and write the combined data to a new CSV file.

# Arguments
- `dir_to_fetch::String`: Directory path containing CSV files corresponding to the specified phase.
- `fase::Int64`: Phase identifier.

# Example
```julia
write_csv_by_fase("path/to/csv/files", 1)```
"""
function write_csv_by_fase(dir_to_fetch::String, fase::Int64,max_rows::Int)
    data_fases_dir = joinpath("data","fases")
    mkpath(data_fases_dir)

    files = filter(file -> endswith(file, "fase$(fase).csv"), readdir(abspath(dir_to_fetch)))
    sizes = size(load_data_matrix(joinpath(dir_to_fetch,first(files))))
    ensamble_matrix_by_fase = zeros(Float64,max_rows, sizes[2])
    
    Threads.@threads for file in files
        m = load_data_matrix(joinpath(dir_to_fetch,file))
        if size(m)[1] < max_rows
            m_embedded = embedd_matrix_in_row_space(m,max_rows)
            ensamble_matrix_by_fase = m_embedded
        else
            m_reduced = m[1:max_rows,:]
            ensamble_matrix_by_fase = m_reduced
        end        
        ensamble_matrix_by_fase .+= (1/length(files) .* ensamble_matrix_by_fase)
    end
    @show isa( ensamble_matrix_by_fase, Matrix{Float64})

    df_ensamble_by_fase = DataFrames.DataFrame(ensamble_matrix_by_fase, :auto)
    CSV.write(joinpath(data_fases_dir,"fase$(fase).csv"),df_ensamble_by_fase)
end

"""
    heatmap_and_eigenspectra_plots_by_fase(dir_to_fetch::String, fase::Int64, l::Int64)

Generate heatmap and eigenvalue spectrum plots for a specific phase's data.

# Arguments
- `dir_to_fetch::String`: Directory path containing CSV files corresponding to the specified phase.
- `fase::Int64`: Phase identifier.
- `l::Int64`: Window size for computing average eigenvalues.

# Example
```julia
heatmap_and_eigenspectra_plots_by_fase("path/to/csv/files", 1, 5)```
"""
function heatmap_and_eigenspectra_plots_by_fase(fase::Int64, l::Int64)
    #create heatmap/faseX dir. X=0,1,2
    heatmap_dir = joinpath("heatmap","fase$(fase)")
    mkpath(heatmap_dir)
    
    #crate eigspectrum/faseX dir. X=0,1,2
    eigspectrum_dir = joinpath("eigspectrum","fase$(fase)")
    mkpath(eigspectrum_dir)
    
    ensamble_matrix_by_fase = load_data_matrix("/home/enki/MSP/data/fases/fase$(fase).csv",drop_header=true)
    
    heatmap_plot(ensamble_matrix_by_fase,"fase$fase"; dir_to_save=abspath(heatmap_dir))
    plot_eigen_spectrum(ensamble_matrix_by_fase,l,"fase$fase"; dir_to_save=abspath(eigspectrum_dir))
end

"""
    heatmaps_matrix_fase_differences()

Generate heatmap plots showing the differences between matrices corresponding to different phases.

# Example
```julia
heatmaps_matrix_fase_differences()```
"""
function heatmaps_matrix_fase_differences()
   #all loaded matrices are centralized by default
   m1 = load_data_matrix("data/fases/fase0.csv")
   m2 = load_data_matrix("data/fases/fase1.csv")
   m3 = load_data_matrix("data/fases/fase2.csv")

   heatmap_plot(m2-m1,"fase21"; dir_to_save=abspath("heatmap"))
   heatmap_plot(m3-m2,"fase32"; dir_to_save=abspath("heatmap"))
   heatmap_plot(m3-m1,"fase31"; dir_to_save=abspath("heatmap"))
end