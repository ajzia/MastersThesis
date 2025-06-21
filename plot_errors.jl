using JSON 
using Plots
ENV["GKSwstype"] = "100"


"""
    join_path(file::String)

Join the given file path with the current directory.

# Parameters:
- `file::String`: The file path to be joined.

# Returns:
- `String`: The joined file path.
"""
@inline join_path(file::String) = joinpath(@__DIR__, file)


"""
    json_to_dict(file_path::String)

Convert a JSON file to a dictionary.

# Parameters:
- `file_path::String`: The path to the JSON file.
"""
function json_to_dict(file_path::String)
  file = open(file_path, "r")
  content = read(file, String)
  close(file)
  return JSON.parse(content)
end


"""
    plot_avg_errors(data::Dict{String, Any})

Plot the average errors from the given data dictionary.

# Parameters:
- `data::Dict{String, Any}`: A dictionary containing the data to be plotted.
"""
function plot_avg_errors(data::Dict{String, Any})
  algorithms::Vector{String} = collect(keys(data))
  no_nodes::Int = data["graph"]["nodes"]
  
  
  # remove StoerWagner, graph and Karger from keys
  algorithms = filter(k -> !(k in ["StoerWagner", "graph", "Karger"]), algorithms)

  # get x values from keys (Sketch_x)
  m_values::Vector{Int} = [parse(Int, split(k, "_")[2]) for k in algorithms]
  sort!(m_values)

  avg_errors::Vector{Float64} = zeros(Float64, length(m_values))
  for (i, m) in enumerate(m_values)
    avg_errors[i] = data["Sketch_$m"]["avg. error"]
  end

  # Plotting
  plot(m_values, avg_errors,
       xlabel="Rozmiar szkicu (m)",
       ylabel="Średni błąd",
       title="Średnie błędy dla egzemplarza o n = $no_nodes dla różnych rozmiarów szkicu",
       label="Średni błąd",
       marker=:circle,
       line=:solid,
       legend=:topright,
       grid=true,
       size=(700, 400),
       xticks=m_values,
       margin=5Plots.mm,
       topmargin=2Plots.mm,
       titlefontsize=10,
       guidefontsize=10
      )

  # Save the plot
  savefig(join_path("n$(no_nodes)_avg_errors_plot.svg"))
end


function main(args::Array{String})
  if length(args) != 1
    println("Usage: julia plot_sizes.jl <path_to_json_file>")
    return
  end

  json_file_path = args[1]
  if !isfile(json_file_path)
    println("File not found: $json_file_path")
    return
  end

  dict = json_to_dict(join_path(json_file_path))
  plot_avg_errors(dict)
end


if abspath(PROGRAM_FILE) == @__FILE__
  main(ARGS)
end
