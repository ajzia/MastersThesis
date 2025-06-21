include("./Sketches/Graph.jl")
using .BasicGraph
include("./Sketches/EdgeSketch.jl")
using .Sketch

include("./Algorithms/GraphMinCut.jl")
include("./Algorithms/SketchMinCut.jl")
include("./Algorithms/StoerWagner.jl")

include("./SBM.jl")

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
    run_graph(g::SimpleGraph, iter::Int, results::Dict, filename::String, weights::Vector{Float64})

Run Karger's algorithm on the given graph and store the results.

# Parameters:
- `g::SimpleGraph`: The graph on which Karger's algorithm will be run.
- `iter::Int`: The number of iterations to run Karger's algorithm.
- `results::Dict`: A dictionary to store the results of the algorithm.
- `filename::String`: The name of the file where results will be saved.
- `weights::Vector{Float64}`: A vector of weights for the edges in the graph.
"""
function run_graph(g::SimpleGraph, iter::Int, results::Dict, filename::String, weights::Vector{Float64})
  println("Min degree: $(δ(g)), Max degree: $(Δ(g))")
  adj_matrix::Matrix{Float64} = zeros(Float64, nv(g), nv(g))
  for (i, e) in enumerate(edges(g))
    adj_matrix[src(e), dst(e)] = weights[i]
    # if !g.is_directed
    adj_matrix[dst(e), src(e)] = weights[i]
    # end
  end
  graph::BasicGraph.Graph = BasicGraph.Graph(
    adj_matrix,
    false
  )
  run_graph(graph, iter, results, filename)
end


"""
    run_graph(graph::BasicGraph.Graph, iter::Int, results::Dict, filename::String)

Run Karger's algorithm on the given graph and store the results.

# Parameters:
- `graph::BasicGraph.Graph`: The graph on which Karger's algorithm will be run.
- `iter::Int`: The number of iterations to run Karger's algorithm.
- `results::Dict`: A dictionary to store the results of the algorithm.
- `filename::String`: The name of the file where results will be saved.

# Returns:
- `Vector{Float64}`: A vector containing the cuts found in each iteration.
"""
function run_graph(graph::BasicGraph.Graph, iter::Int, results::Dict, filename::String)::Vector{Float64}
  results["Karger"] = Dict{String, Union{Float64, Int}}(
    "minimum_cut" => 0.0,
    "minimum_time" => 0.0,
    "average_cut" => 0.0,
    "average_time" => 0.0,
    "failed_attempts" => 0,
    "memory" => 0
  )

  minimum_cut::Float64 = Inf
  minimum_time::Float64 = Inf
  avg_time::Float64 = 0.0
  avg_cut::Float64 = 0.0

  cut_results::Vector{Float64} = []
  failed_attempts::Int = 0
  @time for i in 1:iter
    start = time()
    cut::Float64 = GraphMinCut(graph)
    elapsed = time() - start
    push!(cut_results, cut)
    print("Iteration: $i, Cut: $cut, Elapsed: $elapsed\n")
 
    if cut == -1 || cut == 0.0
      failed_attempts += 1
      continue
    end
    
    if cut < minimum_cut
      minimum_cut = cut
    end
    if elapsed < minimum_time
      minimum_time = elapsed
    end

    avg_time += elapsed
    avg_cut += cut
  end

  
  if iter - failed_attempts == 0
    results["Karger"]["average_cut"] = 0.0
    results["Karger"]["average_time"] = 0.0
  else
    results["Karger"]["average_cut"] = avg_cut / (iter - failed_attempts)
    results["Karger"]["average_time"] = avg_time / (iter - failed_attempts)
  end
  
  results["Karger"]["minimum_cut"] = minimum_cut
  results["Karger"]["minimum_time"] = minimum_time
  results["Karger"]["failed_attempts"] = failed_attempts
  results["Karger"]["memory"] = Base.summarysize(graph)

  println(results)

  isdir("./Results") || mkdir("./Results")
  open("./Results/$(basename(filename)).json", "w") do file
    write(file, JSON.json(results))
  end

  return cut_results
end


"""
    run_sketch(g::SimpleGraph,iter::Int, results::Dict, filename::String,
               ms::Vector{Int}, weights::Vector{Float64})

Run Karger's algorithm for sketches on the given graph edges and store the results.

# Parameters:
- `g::SimpleGraph`: Generated graph to run Karger's algorithm on.
- `iter::Int`: The number of iterations to run Karger's algorithm.
- `results::Dict`: A dictionary to store the results of the algorithm.
- `filename::String`: The name of the file where results will be saved.
- `ms::Vector{Int}`: A vector of sketch sizes to be used in the algorithm.
- `weights::Vector{Float64}`: A vector of weights for the edges in the graph.
# - `save_results::Bool`: Whether to save the results to a dictionary or not.
"""
function run_sketch(
  g::SimpleGraph,
  iter::Int,
  results::Dict,
  filename::String,
  ms::Vector{Int},
  weights::Vector{Float64},
  save_results::Bool = false
)
  graph_edges::Vector{Tuple{Int, Int, Float64}} = [
    (src(e), dst(e), weights[i]) for (i, e) in enumerate(edges(g))
  ]
  return run_sketch(deepcopy(graph_edges), iter, results, filename, ms, save_results)   
end


"""
    run_sketch(filename::String, iter::Int, results::Dict, ms::Vector{Int})

Run Karger's algorithm for sketches on the given graph and store the results.

# Parameters:
- `filename::String`: The name of the file where results will be saved.
- `iter::Int`: The number of iterations to run Karger's algorithm.
- `results::Dict`: A dictionary to store the results of the algorithm.
- `ms::Vector{Int}`: A vector of sketch sizes to be used in the algorithm.
# - `save_results::Bool`: Whether to save the results to a dictionary or not.
"""
function run_sketch(filename::String, iter::Int, results::Dict, ms::Vector{Int}, save_results::Bool = false)
  edges::Vector{Tuple{Int, Int, Float64}} = []
  open(join_path(filename)) do file
    for line in eachline(file)
      if startswith(line, "%") continue end
      if startswith(line, "#") continue end
      edge = split(line, r"[\s,]+")

      i = parse(Int, edge[1]) + 1
      j = parse(Int, edge[2]) + 1

      if startswith(basename(filename), "graph_") # This is generated
        i = i - 1
        j = j - 1
      end

      w = parse(Float64, edge[3])
      push!(edges, (i, j, w))
    end
  end


  return run_sketch(edges, iter, results, basename(filename), ms, save_results)
end


"""
    run_sketch(graph_edges::Vector{Tuple{Int, Int, Float64}},
               iter::Int, results::Dict, filename::String, ms::Vector{Int})

Run Karger's algorithm for sketches on the given graph edges and store the results.

# Parameters:
- `graph_edges::Vector{Tuple{Int, Int, Float64}}`: A vector of edges in the format
                                                   (source, destination, weight).
- `iter::Int`: The number of iterations to run Karger's algorithm.
- `results::Dict`: A dictionary to store the results of the algorithm.
- `filename::String`: The name of the file where results will be saved.
- `ms::Vector{Int}`: A vector of sketch sizes to be used in the algorithm.
- `save_results::Bool`: Whether to save the results to a dictionary or not.

# Returns:
- `Dict{Int, Vector{Float64}}`: A dictionary  of cuts found in each iteration for
                                each sketch size.
"""
function run_sketch(
  graph_edges::Vector{Tuple{Int, Int, Float64}},
  iter::Int,
  results::Dict,
  filename::String,
  ms::Vector{Int},
  save_results::Bool = false
)::Dict{Int, Vector{Float64}}
  cut_results::Dict{Int, Vector{Float64}} = Dict()
  for m in ms
    cut_results[m] = Vector{Float64}()
    println("Running sketch with m=$m, iter=$iter")
    results["Sketch_$m"] = Dict{String, Union{Float64, Int}}(
      "minimum_cut" => 0.0,
      "minimum_time" => 0.0,
      "average_cut" => 0.0,
      "average_time" => 0.0,
      "failed_attempts" => 0,
      "memory" => 0
    )
    sketch::EdgeSketch = EdgeSketch(m, graph_edges)

    for node in sketch.nodes
      node.estimated_weight = EstimateNodeWeight(node, sketch.m)
    end

    minimum_cut::Float64 = Inf
    minimum_time::Float64 = Inf
    avg_time::Float64 = 0.0
    avg_cut::Float64 = 0.0

    failed_attempts::Int = 0
    for i in 1:iter
      cut::Float64 = 0.0
      if !save_results
        start = time()
        cut = SketchMinCut(sketch)
        elapsed = time() - start
      else
        while true
          start = time()
          cut = SketchMinCut(sketch)
          elapsed = time() - start
          if cut != -1 && cut != 0.0
            break
          end
        end
      end
      push!(cut_results[m], cut)
      print("Iteration: $i, Cut: $cut, Elapsed: $elapsed\n")
      if cut == -1 || cut == 0.0
        failed_attempts += 1
        continue
      end

      if cut < minimum_cut
        minimum_cut = cut
      end
 
      if elapsed < minimum_time
        minimum_time = elapsed
      end
      avg_time += elapsed
      avg_cut += cut
    end

    results["Sketch_$m"]["minimum_cut"] = minimum_cut
    if (iter - failed_attempts) == 0
      results["Sketch_$m"]["average_cut"] = 0.0
      results["Sketch_$m"]["average_time"] = 0.0
    else
      results["Sketch_$m"]["average_cut"] = avg_cut / (iter - failed_attempts)
      results["Sketch_$m"]["average_time"] = avg_time / (iter - failed_attempts)
    end
    results["Sketch_$m"]["minimum_time"] = minimum_time
    results["Sketch_$m"]["failed_attempts"] = failed_attempts
    memory::Int = 0
    for node in sketch.nodes
      memory += Base.summarysize(node.S)
      memory += Base.summarysize(node.F)
      memory += Base.summarysize(node.max)
    end
    results["Sketch_$m"]["memory"] = memory

    println(results)

    isdir("./Results") || mkdir("./Results")
    open("./Results/$(basename(filename)).json", "w") do file
      write(file, JSON.json(results))
    end
  end

  return cut_results
end


function from_file_usage()
  println("Usage: julia main.jl FROM_FILE <file> <iter> m1,m2,...,mn")
  println("file: path to the graph file")
  println("iter: number of iterations for Karger's algorithm")
  println("m1,m2,...,mn: sketch sizes (comma-separated)")
  println("Optionally you can specify the HISTOGRAM flag at the end to generate a histogram of the results.")
end


"""
    create_histogram(graph_cuts::Vector{Float64}, sketch_cuts::Vector{Float64},
                   minimum_cut::Float64, filename::String, n::Int, m::Int)

Create a histogram of the cuts found by Karger's algorithm for both graph and sketch.

# Parameters:
- `graph_cuts::Vector{Float64}`: A vector containing the cuts found in the graph.
- `sketch_cuts::Vector{Float64}`: A vector containing the cuts found in the sketch.
- `minimum_cut::Float64`: The minimum cut found by the Stoer-Wagner algorithm.
- `filename::String`: The name of the file where the histogram will be saved.
- `n::Int`: The number of nodes in the graph.
- `m::Int`: The size of the sketch used in the algorithm.
"""
function create_histogram(
  graph_cuts::Vector{Float64},
  sketch_cuts::Vector{Float64},
  minimum_cut::Float64,
  filename::String,
  n::Int,
  m::Int
)
  # Filter too big/small values
  avg_graph_cut = round(sum(graph_cuts) / length(graph_cuts); digits=2)
  min_graph_cut = round(minimum(graph_cuts); digits=2)
  avg_sketch_cut = round(sum(sketch_cuts) / length(sketch_cuts); digits=2)
  min_sketch_cut = round(minimum(sketch_cuts); digits=2)
  minimum_cut = round(minimum_cut; digits=2)
  filter!(x -> 0.8 * min_graph_cut < abs(x) < min(2 * avg_graph_cut, 2500 + avg_graph_cut), graph_cuts)
  filter!(x -> 0.8 * min_sketch_cut < abs(x) < min(2 * avg_sketch_cut, 2500 + avg_sketch_cut), sketch_cuts)
  plt = histogram(graph_cuts, label="Cięcia grafów", bins=20, alpha=0.5, color=:blue, 
                  margin=5Plots.mm, topmargin=2Plots.mm, titlefontsize=10, legend=:topright, guidefontsize=10,
                  xrotation=75, size=(800, 500))
  histogram!(sketch_cuts, label="Cięcia szkiców", bins=20, alpha=0.5, color=:red)
            
  vline!([avg_graph_cut], label="Śr. cięcie grafu ($(avg_graph_cut))", color=:blue, linestyle=:dot, linewidth=2)
  vline!([avg_sketch_cut], label="Śr. cięcie szkicu ($(avg_sketch_cut))", color=:orange, linestyle=:dot, linewidth=2)
  vline!([min_sketch_cut], label="Min. cięcie szkicu ($(min_sketch_cut))", color=:red, linestyle=:dot, linewidth=2)
  vline!([minimum_cut], label="Min. cięcie ($(minimum_cut))", color=:green, linestyle=:dot, linewidth=2)

  old_xticks_values, old_xticks_strings = xticks(plt[1])
  new_xticks_values = [minimum_cut, min_graph_cut, avg_graph_cut, min_sketch_cut, avg_sketch_cut]
  new_xticks_strings = string.(round.(new_xticks_values; digits=2))
  merged_xticks = (old_xticks_values ∪ new_xticks_values, old_xticks_strings ∪ new_xticks_strings)
  xticks!(merged_xticks)

  xlabel!("Rozmiar cięcia")
  ylabel!("Liczba wystąpień")
  title!("Histogram cięć algorytmu Kargera dla $n wierzchołków i m=$m")
  savefig(filename)
end


function main_from_file(args::Array{String})
  if length(args) < 3
    from_file_usage()
    return
  end

  filename::String = args[1]
  if !isfile(filename)
    println("File $filename does not exist.")
    return
  end
  it::Int = parse(Int, args[2])
  ms::Vector{Int} = parse.(Int, split(args[3], ","))
  make_histogram::Bool = false
  if length(args) > 3 && args[4] == "HISTOGRAM"
    make_histogram = true
  end

  println("Reading graph from $filename")
  g = GraphCreate(filename)

  degrees::Vector{Int} = [sum(g.adj[i, :] .> 0) for i in 1:g.no_nodes]

  results::Dict = Dict()
  results["graph"] = Dict(
    "nodes" => g.no_nodes,
    "edges" => g.no_edges,
    "iterations" => it,
    "min_degree" => minimum(degrees),
    "avg_degree" => sum(degrees) / length(degrees),
    "max_degree" => maximum(degrees),
  )

  println("Read graph with $(g.no_nodes) nodes and $(g.no_edges) edges")
  println("Running Karger's algorithm...")
  graph_cuts = run_graph(g, it, results, basename(filename))
  println("Running Karger's algorithm for sketches...")
  sketch_cuts = run_sketch(filename, it, results, ms, make_histogram)

  if g.no_edges > 250000 return end
  (minimum_cut, _) = StoerWagnerMinCut(g.adj)
  results["StoerWagner"] = Dict("minimum_cut" => minimum_cut)
  results["Karger"]["avg. error"] = abs(results["Karger"]["average_cut"] - minimum_cut) / minimum_cut
  results["Karger"]["min_error"] = abs(results["Karger"]["minimum_cut"] - minimum_cut) / minimum_cut

  for m in ms
    results["Sketch_$m"]["avg. error"] = abs(results["Sketch_$m"]["average_cut"] - minimum_cut) / minimum_cut
    results["Sketch_$m"]["min_error"] = abs(results["Sketch_$m"]["minimum_cut"] - minimum_cut) / minimum_cut
  end

  open("./Results/$(basename(filename))_results.json", "w") do file
    write(file, JSON.json(results, 2))
  end
  println("Results saved to ./Results/$(basename(filename))_results.json")

  if make_histogram
    println("Generating histogram...")
    for m in ms
      create_histogram(
        graph_cuts,
        sketch_cuts[m],
        minimum_cut,
        "./Results/$(basename(filename))_histogram_m$m.png",
        g.no_nodes,
        m
      )
      println("Histogram saved to ./Results/$(basename(filename))_histogram_m$m.png")
    end
  end
end


function generate_usage()
  println("Usage: julia main.jl GENERATE <n> <b> <p> <q> <iter> <m1,m2,...,mn>")
  println("n: number of nodes")
  println("b: number of blocks")
  println("p: intra-block connection probability")
  println("q: inter-block connection probability")
  println("iter: number of iterations for Karger's algorithm")
  println("m1,m2,...,mn: sketch sizes (comma-separated)")
  println("Optionally you can specify the HISTOGRAM flag at the end to generate a histogram of the results.")
end


function main_generate(args::Array{String})
  if length(args) < 6
    generate_usage()
    return
  end
  
  n = parse(Int, args[1])
  b = parse(Int, args[2])
  p = parse(Float64, args[3])
  q = parse(Float64, args[4])
  it::Int = parse(Int, args[5])
  ms::Vector{Int} = parse.(Int, split(args[6], ","))
  make_histogram::Bool = false
  if length(args) > 6 && args[7] == "HISTOGRAM"
    make_histogram = true
  end

  # generate and save the graph
  println("Generating SBM with n=$n, b=$b, p=$p, q=$q")
  (g, _) = generate_sbm_manual(n, b, p, q)
  weights::Vector{Float64} = [rand(10:20) for _ in 1:ne(g)]

  filename::String = "graph_n$(n)_e$(ne(g))_b$(b)_p$(p)_q$(q)_it$(it)_$(time())"
  isdir("./Resources") || mkdir("./Resources")
  open("./Resources/$(basename(filename))", "w") do file
    for (i, edge) in enumerate(edges(g))
      write(file, "$(src(edge)) $(dst(edge)) $(weights[i])\n")
    end
  end
  
  results::Dict = Dict()
  results["graph"] = Dict(
    "nodes" => nv(g),
    "edges" => ne(g),
    "blocks" => b,
    "p" => p,
    "q" => q,
    "iterations" => it,
    "min_degree" => δ(g),
    "avg_degree" => sum(degree(g)) / nv(g),
    "max_degree" => Δ(g)
  )

  println("Generated graph with $(nv(g)) nodes and $(ne(g)) edges")
  println("Running Karger's algorithm for sketches...")
  sketch_cuts = run_sketch(g, it, results, filename, ms, weights, make_histogram)
  println("Running Karger's algorithm...")
  graph_cuts = run_graph(g, it, results, filename, weights)


  n::Int = nv(g)
  adj_matrix::Matrix{Float64} = zeros(Float64, n, n)
  for (i, e) in enumerate(edges(g))
    adj_matrix[src(e), dst(e)] = Float64(weights[i])
    adj_matrix[dst(e), src(e)] = Float64(weights[i])
  end

  if ne(g) > 250000 return end
  (minimum_cut, _) = StoerWagnerMinCut(adj_matrix)
  results["StoerWagner"] = Dict("minimum_cut" => minimum_cut)

  results["Karger"]["avg. error"] = abs(results["Karger"]["average_cut"] - minimum_cut) / minimum_cut
  results["Karger"]["min_error"] = abs(results["Karger"]["minimum_cut"] - minimum_cut) / minimum_cut

  for m in ms
    results["Sketch_$m"]["avg. error"] = abs(results["Sketch_$m"]["average_cut"] - minimum_cut) / minimum_cut
    results["Sketch_$m"]["min_error"] = abs(results["Sketch_$m"]["minimum_cut"] - minimum_cut) / minimum_cut
  end

  open("./Results/$(basename(filename)).json", "w") do file
    write(file, JSON.json(results, 2),)
  end
  println("Results saved to ./Results/$(basename(filename)).json")

  if make_histogram
    println("Generating histogram...")
    for m in ms
      create_histogram(
        graph_cuts,
        sketch_cuts[m],
        minimum_cut,
        "./Results/$(basename(filename))_histogram_m$m.png",
        n,
        m
      )
      println("Histogram saved to ./Results/$(basename(filename))_histogram_$m.png")
    end
  end
end


if abspath(PROGRAM_FILE) == @__FILE__
  if length(ARGS) == 0
    println("No arguments provided. Please provide the necessary arguments.")
    return
  end

  if ARGS[1] == "FROM_FILE"
    main_from_file(ARGS[2:end])
  elseif ARGS[1] == "GENERATE"
    main_generate(ARGS[2:end])
  else
    println("Invalid command.")
    println("FROM_FILE mode:")
    from_file_usage()
    println("GENERATE mode:")
    generate_usage()
  end
end
