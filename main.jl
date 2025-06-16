include("./Sketches/Graph.jl")
using .BasicGraph
include("./Sketches/EdgeSketch.jl")
using .Sketch

include("./Algorithms/GraphMinCut.jl")
include("./Algorithms/SketchMinCut.jl")
include("./Algorithms/StoerWagner.jl")

include("./SBM.jl")

using JSON


@inline join_path(file::String) = joinpath(@__DIR__, file)

function run_graph(g::SimpleGraph, iter::Int, results::Dict, filename::String, weights::Vector{Float64})
  println("Min degree: $(δ(g)), Max degree: $(Δ(g))")
  graph::BasicGraph.Graph = BasicGraph.Graph(
    g.fadjlist,
    [(src(e), dst(e), weights[i]) for (i, e) in enumerate(edges(g))],
    false
  )
  run_graph(graph, iter, results, filename, weights)
end

function run_graph(graph::BasicGraph.Graph, iter::Int, results::Dict, filename::String, weights::Vector{Float64})
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

  failed_attempts::Int = 0
  @time for i in 1:iter
    start = time()
    cut::Float64 = GraphMinCut(graph)
    elapsed = time() - start
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
end

function run_sketch(g::SimpleGraph, iter::Int, results::Dict, filename::String, ms::Vector{Int}, weights::Vector{Float64})
  graph_edges::Vector{Tuple{Int, Int, Float64}} =
    [(src(e), dst(e), weights[i]) for (i, e) in enumerate(edges(g))]

  return run_sketch(graph_edges, iter, results, filename, ms)    
end

function run_sketch(g::BasicGraph.Graph, iter::Int, results::Dict, filename::String, ms::Vector{Int}, weights::Vector{Float64})
  return run_sketch(deepcopy(g.edges), iter, results, filename, ms)
end


function run_sketch(graph_edges::Vector{Tuple{Int, Int, Float64}}, iter::Int, results::Dict, filename::String, ms::Vector{Int})
  for m in ms
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
      start = time()
      cut::Float64 = SketchMinCut(sketch)
      elapsed = time() - start

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

end

function from_file_usage()
  println("Usage: julia main.jl FROM_FILE <file> <iter> m1,m2,...,mn")
  println("file: path to the graph file")
  println("iter: number of iterations for Karger's algorithm")
  println("m1,m2,...,mn: sketch sizes (comma-separated)")
end

function main_from_file(args::Array{String})
  if length(args) != 3
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

  println("Reading graph from $filename")
  g = GraphCreate(filename)
  weights::Vector{Float64} = [rand(10:20) for _ in 1:length(g.edges)]

  results::Dict = Dict()
  results["graph"] = Dict(
    "nodes" => length(g.node_ids),
    "edges" => length(g.edges),
    "iterations" => it,
  )

  println("Read graph with $(length(g.node_ids)) nodes and $(length(g.edges)) edges")
  println("Running Karger's algorithm...")
  run_graph(g, it, results, basename(filename), weights)
  println("Running Karger's algorithm for sketches...")
  run_sketch(g, it, results, basename(filename), ms, weights)

  n::Int = length(g.node_ids)
  adj_matrix::Matrix{Int} = zeros(Int, n, n)
  for (i, j, w) in g.edges
    adj_matrix[i, j] = weights[i]
    adj_matrix[j, i] = weights[i]
  end

  if length(g.edges) > 250000 return end
  (minimum_cut, _) = StoerWagnerMinCut(adj_matrix)
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
end

function generate_usage()
  println("Usage: julia main.jl GENERATE <n> <b> <p> <q> <iter> <m1,m2,...,mn>")
  println("n: number of nodes")
  println("b: number of blocks")
  println("p: intra-block connection probability")
  println("q: inter-block connection probability")
  println("iter: number of iterations for Karger's algorithm")
  println("m1,m2,...,mn: sketch sizes (comma-separated)")
end

function main_generate(args::Array{String})
  if length(args) != 6
    generate_usage()
    return
  end
  
  n = parse(Int, args[1])
  b = parse(Int, args[2])
  p = parse(Float64, args[3])
  q = parse(Float64, args[4])
  it::Int = parse(Int, args[5])
  ms::Vector{Int} = parse.(Int, split(args[6], ","))

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
  )

  println("Generated graph with $(nv(g)) nodes and $(ne(g)) edges")
  println("Running Karger's algorithm for sketches...")
  run_sketch(g, it, results, filename, ms, weights)
  println("Running Karger's algorithm...")
  run_graph(g, it, results, filename, weights)


  n::Int = nv(g)
  adj_matrix::Matrix{Int} = zeros(Int, n, n)
  for (i, e) in enumerate(edges(g))
    adj_matrix[src(e), dst(e)] = weights[i]
    adj_matrix[dst(e), src(e)] = weights[i]
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
