using StatsBase


"""
    contractEdge!(graph::BasicGraph.Graph, i::Int, j::Int)

Contract an edge between nodes `i` and `j` in the graph, updating the graph in place.

# Parameters:
- `graph::BasicGraph.Graph`: The graph containing the nodes and their properties.
- `i::Int`: The id of the first node.
- `j::Int`: The id of the second node.
"""
function contractEdge!(graph::BasicGraph.Graph, i::Int, j::Int)
  # check if i and j are the same
  if i == j
    return
  end

  # check if i and j are in the graph
  if i ∉ graph.node_ids || j ∉ graph.node_ids
    return
  end

  # check if edge (i, j) exists
  if graph.adj[i, j] == 0.0 || (graph.is_directed && graph.adj[j, i] == 0.0)
    return
  end

  # For each node k that was connected to j, add its edge weight to i
  for k in 1:size(graph.adj, 1)
    if k != i && k != j
      if graph.adj[j, k] != 0.0
        graph.adj[i, k] += graph.adj[j, k]
      end
      if graph.adj[k, j] != 0.0
        graph.adj[k, i] += graph.adj[k, j]
      end
    end
  end

  # Remove neighbours of j
  graph.adj[j, :] .= 0.0
  graph.adj[:, j] .= 0.0

  # Remove edge/s between i and j
  graph.adj[i, j] = 0.0

  # Remove j from node_ids
  filter!(x -> x != j, graph.node_ids)
  graph.no_nodes = length(graph.node_ids)
end # contractEdge!


"""
    GraphMinCut(graph::BasicGraph.Graph)::Float64

Find the minimum cut of a graph using the Karger's algorithm.

# Parameters:
- `graph::BasicGraph.Graph`: The graph for which the minimum cut is to be found.

# Returns:
- `Float64`: Weight of the found minimum cut.
"""
function GraphMinCut(graph::BasicGraph.Graph)::Float64
  copy_graph = deepcopy(graph)

  while copy_graph.no_nodes > 2
    if copy_graph.no_nodes % 1000 == 0
      println("Nodes left: ", copy_graph.no_nodes)
    end
   # choosing an edge to contract
    edges::Vector{Tuple{Int, Int}} = [Tuple(ij) for ij in CartesianIndices(copy_graph.adj) if copy_graph.adj[ij] != 0.0]
    weights::Vector{Float64} = [copy_graph.adj[i, j] for (i, j) in edges]
    weights /= sum(weights)
    
    (i, j) = sample(edges, Weights(weights))
    contractEdge!(copy_graph, i, j)
  end

  if all(copy_graph.adj .== 0.0)
    return -1.0
  end

  if copy_graph.is_directed
    return sum(copy_graph.adj)
  end

  # for undirected graphs, we only count the upper triangle of the adjacency matrix
  return sum([copy_graph.adj[i, j]
             for i in 1:size(copy_graph.adj, 1)
             for j in i+1:size(copy_graph.adj, 2)])
end # GraphMinCut
