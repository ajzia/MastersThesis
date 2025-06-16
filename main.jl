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



function run_graph(g::SimpleGraph, iter::Int, results::Dict, filename::String, weights)
  results["Karger"] = Dict()
   graph::BasicGraph.Graph = BasicGraph.Graph(
      g.fadjlist,
      [(src(e), dst(e), weights[i]) for (i, e) in enumerate(edges(g))],
      false
    )

  minimum_cut::Float64 = Inf
  minimum_time::Float64 = Inf
  avg_time::Float64 = 0.0
  avg_cut::Float64 = 0.0

  failed_attempts::Int = 0
  @time for _ in 1:iter
    start = time()
    cut::Float64 = GraphMinCut(graph)
    elapsed = time() - start

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
  results["Karger"]["iterations"] = iter
  results["Karger"]["failed_attempts"] = failed_attempts
  results["Karger"]["memory"] = Base.summarysize(graph)

  println(results)

  isdir("./Results") || mkdir("./Results")
  open("./Results/$filename.json", "w") do file
    write(file, JSON.json(results))
  end
end



function run_sketch(g::SimpleGraph, iter::Int, results::Dict, filename::String, m::Int, weights)
  println("Running sketch with m=$m, iter=$iter")
  results["Sketch_$m"] = Dict()
  graph_edges::Vector{Tuple{Int, Int, Float64}} =
    [(src(e), dst(e), weights[i]) for (i, e) in enumerate(edges(g))]
  sketch::EdgeSketch = EdgeSketch(m, graph_edges)

  for node in sketch.nodes
    node.estimated_weight = EstimateNodeWeight(node, sketch.m)
  end

  minimum_cut::Float64 = Inf
  minimum_time::Float64 = Inf
  avg_time::Float64 = 0.0
  avg_cut::Float64 = 0.0

  failed_attempts::Int = 0
  for _ in 1:iter
    start = time()
    cut::Float64 = SketchMinCut(sketch)
    elapsed = time() - start

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
  results["Sketch_$m"]["iterations"] = iter
  results["Sketch_$m"]["failed_attempts"] = failed_attempts
  results["Sketch_$m"]["memory"] = Base.summarysize(sketch)

  println(results)

  isdir("./Results") || mkdir("./Results")
  open("./Results/$filename.json", "w") do file
    write(file, JSON.json(results))
  end
end



function main(args::Array{String})
  if length(args) != 6
    println("Usage: julia main.jl <n> <b> <p> <q> <iter> <m>")
    println("n: number of nodes")
    println("b: number of blocks")
    println("p: intra-block connection probability")
    println("q: inter-block connection probability")
    println("iter: number of iterations for Karger's algorithm")
    println("m: sketch size")
    return
  end
  
  n = parse(Int, args[1])
  b = parse(Int, args[2])
  p = parse(Float64, args[3])
  q = parse(Float64, args[4])
  it::Int = parse(Int, args[5])
  m::Int = parse(Int, args[6])


  # generate and save the graph
  println("Generating SBM with n=$n, b=$b, p=$p, q=$q")
  (g, _) = generate_sbm_manual(n, b, p, q)
  weights::Vector{Float64} = [rand(10:20) for _ in 1:ne(g)]

  filename::String = "graph_n$(n)_e$(ne(g))_b$(b)_p$(p)_q$(q)_it$(it)_m$(m)"
  isdir("./Resources") || mkdir("./Resources")
  open("./Resources/$filename", "w") do file
    for (i, edge) in enumerate(edges(g))
      write(file, "$(src(edge)) $(dst(edge)) $(weights[i])\n")
    end
  end


  println("Generated graph with $(nv(g)) nodes and $(ne(g)) edges")
  results::Dict = Dict()
  run_graph(g, it, results, filename, weights)
  run_sketch(g, it, results, filename, m, weights)


  n::Int = nv(g)
  adj_matrix::Matrix{Int} = zeros(Int, n, n)
  for (i, e) in enumerate(edges(g))
    adj_matrix[src(e), dst(e)] = weights[i]
    adj_matrix[dst(e), src(e)] = weights[i]
  end

  if ne(g) > 250000 return end
  (minimum_cut, _) = StoerWagnerMinCut(adj_matrix)
  results["StoerWagner"] = Dict("minimum_cut" => minimum_cut)
  open("./Results/$filename.json", "w") do file
    write(file, JSON.json(results))
  end
end



if abspath(PROGRAM_FILE) == @__FILE__
  main(ARGS)
end
