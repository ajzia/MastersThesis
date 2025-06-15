using StatsBase


"""
    contractEdge!(graph::BasicGraph.Graph, i::Int, j::Int)

Contract an edge between nodes `i` and `j` in the graph,  updating the graph in place.

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
  if (i, j) ∉ [(x, y) for (x, y, _) in graph.edges]
    return
  end

  # remove any edges between i and j
  filter!(x -> x[1] != i || x[2] != j, graph.edges)
  filter!(x -> x[1] != j || x[2] != i, graph.edges)

  # remove j from i's adjacency list
  idx::Int = findfirst(x -> x == i, graph.node_ids)
  filter!(x -> x != j, graph.adj[idx])

  # remove j from node_ids and adj
  idx = findfirst(x -> x == j, graph.node_ids)
  deleteat!(graph.node_ids, idx)
  deleteat!(graph.adj, idx)

  # update adjacencies
  new_edges::Array{Int} = []
  for (k, neighbours) in enumerate(graph.adj)
    for (l, neighbour) in enumerate(neighbours)
      if neighbour == j
        graph.adj[k][l] = i
        push!(new_edges, graph.node_ids[k])
      end
    end
  end

  # append new edges to i's adjacency list
  idx = findfirst(x -> x == i, graph.node_ids)
  append!(graph.adj[idx], new_edges)

  # replace every j with i
  for (k, edge) in enumerate(graph.edges)
    (x, y, w) = edge
    if x == j
      graph.edges[k] = (i, y, w)
    end

    if y == j
      graph.edges[k] = (x, i, w)
    end
  end
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

  while length(copy_graph.node_ids) > 2
    # can happen for unconnected graphs
    if length(copy_graph.edges) == 0
      println("No edges left to contract.")
      return -1.0
    end

    # debug print
    if length(copy_graph.node_ids) % 100 == 0
      println("Nodes left: ", length(copy_graph.node_ids))
    end

    # choosing an edge to contract
    weights::Vector{Float64} = [edge[3] for edge in copy_graph.edges]
    weights /= sum(weights)
    
    (i, j, _) = sample(copy_graph.edges, Weights(weights))
    contractEdge!(copy_graph, i, j)
  end

  if isempty(copy_graph.edges)
    return -1.0
  end

  return sum(edge[3] for edge in copy_graph.edges)
end # GraphMinCut

