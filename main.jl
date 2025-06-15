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


function generate_sbm_graph(n::Int, b::Int, p::Float64, q::Float64, filename::String)
  println("Generating SBM with n=$n, b=$b, p=$p, q=$q")
  (g, _) = generate_sbm_manual(n, b, p, q)

  filename::String = "./Resources/$filename"
  open(filename, "w") do file
    for edge in edges(g)
      write(file, "$(src(edge)) $(dst(edge)) 1\n")
    end
  end

  return g
end


function run_graph(g::SimpleGraph, iter::Int, results::Dict)
  results["Karger"] = Dict()
  graph::BasicGraph.Graph =
    BasicGraph.Graph(g.fadjlist, [(src(e), dst(e), 1.) for e in edges(g)], false)

  minimum_cut::Float64 = 0.0
  minimum_time::Float64 = Inf
  avg_time::Float64 = 0.0
  avg_cut::Float64 = 0.0

  failed_attempts::Int = 0
  for i in 1:iter
    println("Iteration $i")
    start = time()
    cut::Float64 = GraphMinCut(graph)
    elapsed = time() - start

    if cut == -1
      failed_attempts += 1
      continue
    end

    if elapsed < minimum_time
      minimum_time = elapsed
      minimum_cut = cut
    end
    avg_time += elapsed
    avg_cut += cut
  end


  results["Karger"]["minimum_cut"] = minimum_cut
  
  if (iter - failed_attempts) == 0
    results["Karger"]["average_cut"] = 0.0
    results["Karger"]["average_time"] = 0.0
  else
    results["Karger"]["average_cut"] = avg_cut / (iter - failed_attempts)
    results["Karger"]["average_time"] = avg_time / (iter - failed_attempts)
  end
  
  results["Karger"]["minimum_time"] = minimum_time
  results["Karger"]["iterations"] = iter
  results["Karger"]["failed_attempts"] = failed_attempts

  println(results)

  isdir("./Results") || mkdir("./Results")
  open("./Results/$filename.json", "w") do file
    write(file, JSON.json(results))
  end
end


function run_sketch(g::SimpleGraph, iter::Int, results::Dict, m::Int)
  println("Running sketch with m=$m, iter=$iter")
  results["Sketch_$m"] = Dict()
  graph_edges::Vector{Tuple{Int, Int, Float64}} = [(src(e), dst(e), 1.0) for e in edges(g)]
  graph::EdgeSketch = EdgeSketch(m, graph_edges)


  minimum_cut::Float64 = 0.0
  minimum_time::Float64 = Inf
  avg_time::Float64 = 0.0
  avg_cut::Float64 = 0.0

  failed_attempts::Int = 0
  for i in 1:iter
    println("Iteration $i")
    start = time()
    cut::Float64 = SketchMinCut(graph)
    elapsed = time() - start

    if cut == -1
      failed_attempts += 1
      continue
    end

    if elapsed < minimum_time
      minimum_time = elapsed
      minimum_cut = cut
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
    println("m: size of the sketch")
    return
  end
  
  n = parse(Int, args[1])
  b = parse(Int, args[2])
  p = parse(Float64, args[3])
  q = parse(Float64, args[4])
  it::Int = parse(Int, args[5])
  m::Int = parse(Int, args[6])

  filename::String = "graph_n$(n)_b$(b)_p$(p)_q$(q)_it$(it)_m$(m)"

  g = generate_sbm_graph(n, b, p, q, filename)

  results::Dict = Dict()
  run_graph(g, it, results)
  run_sketch(g, it, results, m)

  (minimum_cut, _) = StoerWagnerMinCut(g)
  results["StoerWagner"] = Dict("minimum_cut" => minimum_cut)
end



if abspath(PROGRAM_FILE) == @__FILE__
  main(ARGS)
end
