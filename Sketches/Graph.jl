module BasicGraph

export Graph, GraphCreate
export addEdge!, removeEdge!
export removeNode!


mutable struct Graph
  node_ids::Vector{Int}
  no_nodes::Int
  no_edges::Int
  adj::Matrix{Float64}
  is_directed::Bool

  Graph(is_directed::Bool) = new([], 0, 0, zeros(Float64, 0, 0), is_directed)
  Graph(adj::Matrix{Float64}, is_directed::Bool) = new(
    [i for i in 1:size(adj, 1)],
    size(adj, 1),
    count(x -> x != 0.0, adj),
    adj,
    is_directed
  )
  function Graph(no_nodes::Int, edges::Vector{Tuple{Int, Int, Float64}}, is_directed::Bool)::Graph
    if no_nodes < 1
      throw(ArgumentError("Number of nodes must be at least 1"))
    end
    if length(edges) < 0
      throw(ArgumentError("Number of edges must be at least 0"))
    end

    adj::Matrix{Float64} = zeros(Float64, no_nodes, no_nodes)
    for (i, j, w) in edges
      adj[i, j] = w
      if !is_directed
        adj[j, i] = w
      end
    end
  
    return new(
      [i for i in 1:no_nodes],
      no_nodes,
      count(x -> x != 0.0, adj),
      adj,
      is_directed
    )
  end
end # Graph



function GraphCreate(stream::String)::Graph
  counter = 0
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

      if startswith(basename(stream), "graph_") # This is generated
        i = i - 1
        j = j - 1
      end

      w = parse(Float64, edge[3])
      push!(edges, (i, j, w))
      no_nodes = max(no_nodes, i, j)

      counter += 1
      if counter % 100000 == 0
        println("Processed $counter edges")
      end
    end
  end
  println("Processed total of $counter edges")
  return Graph(no_nodes, edges, false)
end # GraphCreate


function addEdge!(graph::Graph, i::Int, j::Int, w::Float64)
  
  # Add new nodes if they don't exist
  if i ∉ graph.node_ids || j ∉ graph.node_ids
    # Get current size
    old_size = size(graph.adj, 1)
    new_size = max(old_size, i, j)
    
    # Create new larger matrix
    new_adj = zeros(Float64, new_size, new_size)
    graph.no_nodes = new_size
    # Copy old values
    if old_size > 0
      new_adj[1:old_size, 1:old_size] = graph.adj
    end
    
    # Update node_ids and adjacency matrix
    graph.adj = new_adj
    for node_id in (old_size + 1):new_size
      if node_id ∉ graph.node_ids
        push!(graph.node_ids, node_id)
      end
    end
  end

  # Add the edge
  graph.adj[i, j] = w

  if !graph.is_directed
    graph.adj[j, i] = w
  end
end # addEdge!


function removeEdge!(graph::Graph, i::Int, j::Int)
  if graph.adj[i, j] == 0.0 || (graph.is_directed && graph.adj[j, i] == 0.0)
    return
  end

  # Remove edge from adjacency matrix
  graph.adj[i, j] = 0.0
  grpah.no_edges -= 1
    
  if !graph.is_directed
    graph.adj[j, i] = 0.0
    graph.no_edges -= 1
  end
end # removeEdge!


function printGraph(graph::Graph)::Nothing
  println("Graph:")
  println("> Directed: ", graph.is_directed)
  println("> Number of nodes: ", graph.no_nodes)
  println("> Number of edges: ", graph.no_edges)
end # printGraph

end # BasicGraph
