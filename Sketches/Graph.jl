module BasicGraph

export Graph, GraphCreate
export addEdge!, removeEdge!
export removeNode!


struct Graph
  node_ids::Vector{Int}
  adj::Vector{Vector{Int}}
  edges::Vector{Tuple{Int, Int, Float64}}
  is_directed::Bool

  Graph(is_directed::Bool) = new([], [], [], is_directed)
  Graph(adj::Vector{Vector{Int}}, edges::Vector{Tuple{Int, Int, Float64}}, is_directed::Bool) = new(
    [i for i in 1:length(adj)],
    adj,
    edges,
    is_directed
  )
  function Graph(no_nodes::Int, edges::Vector{Tuple{Int, Int, Float64}}, is_directed::Bool)::Graph
    if no_nodes < 1
      throw(ArgumentError("Number of nodes must be at least 1"))
    end
    if length(edges) < 0
      throw(ArgumentError("Number of edges must be at least 0"))
    end
  
    adj::Vector{Vector{Int}} = [Int[] for _ in 1:no_nodes]
    for (i, j, _) in edges
      push!(adj[i], j)
      if !is_directed
        push!(adj[j], i)
      end
    end
  
    return new(
      [i for i in 1:no_nodes],
      adj,
      edges,
      is_directed
    )
  end
end # Graph



function GraphCreate(stream::String)::Graph
  graph::Graph = Graph(false) # change!

  counter = 1
  edges::Vector{Tuple{Int, Int, Float64}} = []
  no_nodes::Int = 0
  open(stream) do file
    while !eof(file)
      line = readline(file)
      if startswith(line, "%") continue end
      if startswith(line, "#") continue end
      edge = split(line, r"[\s,]+")

      i = parse(Int, edge[1]) + 1
      j = parse(Int, edge[2]) + 1
      w = parse(Float64, edge[3])
      push!(edges, (i, j, w))
      no_nodes = max(no_nodes, i, j)

      counter += 1
      if counter % 100000 == 0
        println("Processed $counter edges")
      end
    end
  end

  return Graph(no_nodes, edges, false)
end # GraphCreate


function addEdge!(graph::Graph, i::Int, j::Int, w::Float64)
  if (i, j, w) in graph.edges
    return
  end

  # check if there already is an edge with other weight
  if (i, j) in [(x, y) for (x, y, _) in graph.edges]
    return
  end

  if i ∉ graph.node_ids
    push!(graph.node_ids, i)
    push!(graph.adj, [])
  end

  if j ∉ graph.node_ids
    push!(graph.node_ids, j)
    push!(graph.adj, [])
  end

  index_i::Int = findfirst(x -> x == i, graph.node_ids)
  push!(graph.adj[index_i], j)
  push!(graph.edges, (i, j, w))

  if !graph.is_directed
    index_j::Int = findfirst(x -> x == j, graph.node_ids)
    push!(graph.adj[index_j], i)
  #   push!(graph.edges, (j, i, w))
  end
end # addEdge!


# add removing edge without weight
function removeEdge!(graph::Graph, i::Int, j::Int, w::Float64)
  if (i, j, w) ∉ graph.edges
    return
  end

  idx::Union{Int, Nothing} = findfirst(x -> x == (i, j, w), graph.edges)
  deleteat!(graph.edges, idx)

  idx = findfirst(x -> x == i, graph.node_ids)
  if !isnothing(idx)
    j_idx::Int = findfirst(x -> x == j, graph.adj[idx])
    deleteat!(graph.adj[idx], j_idx)
  end
    
  if !graph.is_directed
    idx = findfirst(x -> x == (j, i, w), graph.edges)
    if !isnothing(idx)
      deleteat!(graph.edges, idx)
    end

    idx = findfirst(x -> x == j, graph.node_ids)
    if !isnothing(idx)
      i_idx::Int = findfirst(x -> x == i, graph.adj[idx])
      deleteat!(graph.adj[idx], i_idx)
    end
  end
end # removeEdge!


function removeNode!(graph::Graph, node::Int)
  if node ∉ graph.node_ids
    return
  end

  idx::Union{Int, Nothing} = findfirst(x -> x == node, graph.node_ids)
  if !isnothing(idx)
    deleteat!(graph.node_ids, idx)
    deleteat!(graph.adj, idx)
  end

  # remove all edges that contain i
  for (i, j, w) in graph.edges
    if i == node || j == node
      removeEdge!(graph, i, j, w)
    end
  end

  for n in graph.adj
    idx = findfirst(x -> x == node, n)
    if !isnothing(idx)
      deleteat!(n, idx)
    end
  end
end # removeNode!


function printGraph(graph::Graph)::Nothing
  println("Graph:")
  println("> Directed: ", graph.is_directed)
  println("> Number of nodes: ", length(graph.node_ids))
  println("> Number of edges: ", length(graph.edges))
end # printGraph


end # BasicGraph
