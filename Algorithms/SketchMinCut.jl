using StatsBase
using Random


"""
    contractEdge!(sketch::EdgeSketch, i::Int, j::Int)

Contract an edge between nodes `i` and `j` in the sketch,
merging their properties and updating the sketch in place.

# Parameters:
- `sketch::EdgeSketch`: The sketch containing the nodes and
                         their properties.
- `i::Int`: The id of the first node.
- `j::Int`: The id of the second node.
"""
function contractEdge!(sketch::EdgeSketch, i::Int, j::Int)
  # check if i and j are the same
  if i == j
    return
  end

  # check if i and j are in the sketch
  if i ∉ sketch.node_ids || j ∉ sketch.node_ids
    return
  end

  idx_i::Int = findfirst(x -> x == i, sketch.node_ids)
  idx_j::Int = findfirst(x -> x == j, sketch.node_ids)

  # create sketch (i, j) in place of the i-th node
  for k in 1:sketch.m
    if sketch.nodes[idx_i].S[k] > sketch.nodes[idx_j].S[k]
      sketch.nodes[idx_i].S[k] = sketch.nodes[idx_j].S[k]
      sketch.nodes[idx_i].F[k] = sketch.nodes[idx_j].F[k]
    end
  end


  # replace every j with i
  for (l, node) in enumerate(sketch.nodes)
    if node.id == j continue end
    for (k, edge) in enumerate(node.F)
      (u::Int, v::Int) = edge
      sketch.nodes[l].F[k] = (u == j ? i : u, v == j ? i : v)
    end
  end


  # count inside edges
  sketch.nodes[idx_i].inside_edges = count(x -> x == (i, i), sketch.nodes[idx_i].F)

  # count empty edges
  sketch.nodes[idx_i].empty_edges = count(x -> x == (0, 0), sketch.nodes[idx_i].F)

  # estimate new weight of the node
  sketch.nodes[idx_i].estimated_weight =
    EstimateNodeWeight(sketch.nodes[idx_i], sketch.m)

  # remove node j from nodes and node_ids
  deleteat!(sketch.node_ids, findfirst(x -> x == j, sketch.node_ids))
  deleteat!(sketch.nodes, findfirst(x -> x.id == j, sketch.nodes))
end # contractEdge!



"""
    GetNodeWeights(sketch::EdgeSketch)::Vector{Float64}

Get weights for each node in the sketch, normalized
to sum to 1.

# Parameters:
- `sketch::EdgeSketch`: The sketch containing the nodes
                         and their properties.

# Returns:
- `Vector{Float64}`: A vector of estimated weights for each
                     node in the sketch, normalized to sum to 1.
"""
@inline function GetNodeWeights(sketch::EdgeSketch)::Vector{Float64}
  weights::Vector{Float64} = [node.estimated_weight for node in sketch.nodes]
  return weights ./ sum(weights)
end # GetNodeWeights


"""
    EstimateNodeWeight(node::NodeSketch, m::Int)::Float64

Estimate the weight of edges in a node.

# Parameters:
- `node::NodeSketch`: The node for which the weight is to be estimated.
- `m::Int`: The total number of edges in the sketch.

# Returns:
- `Float64`: The estimated weight of the node's edges.
"""
@inline function EstimateNodeWeight(node::NodeSketch, m::Int)::Float64
  estimation::Float64 = (m - 1) / sum(s for s in node.S if s != Inf)
  estimation *= (m - node.inside_edges) / m

  return estimation
end # EstimateNodeWeight


"""
    SketchMinCut(sketch::EdgeSketch)::Float64

Estimate the minimum cut of a graph represented by a sketch using
the Karger's algorithm.

# Parameters:
- `sketch::EdgeSketch`: Sketch of a graph.

# Returns:
- `Float64`: The estimated weight of the minimum cut in the graph.
"""
function SketchMinCut(sketch::EdgeSketch)::Float64
  copy_sketch = deepcopy(sketch)

  while length(copy_sketch.node_ids) > 2
    # debug print
    if length(copy_sketch.node_ids) % 100 == 0
      println("Nodes left: ", length(copy_sketch.node_ids))
    end

    weights::Vector{Float64} = GetNodeWeights(copy_sketch)

    # Choosing first node
    i::Int = 0
    idx_i::Int = 0
    counter::Int = 1
    while true
      if counter > 10
        println("No valid node found after 10 attempts, exiting.")
        return -1.0
      end
      i = sample(MersenneTwister(), copy_sketch.node_ids, Weights(weights))
      idx_i = findfirst(x -> x == i, copy_sketch.node_ids)

      # Check if node has outside edges
      invalid_edges::Int = copy_sketch.nodes[idx_i].inside_edges + copy_sketch.nodes[idx_i].empty_edges
      if invalid_edges == copy_sketch.m
        counter += 1
        continue
      end

      break
    end

    # Choosing an edge from the selected node
    counter = 1
    j::Int = 0
    while true
      (u::Int, v::Int) = copy_sketch.nodes[idx_i].F[counter]
      if u == v || u == 0
        counter += 1
        continue
      end

      i, j = min(u, v), max(u, v)
      break
    end

    contractEdge!(copy_sketch, i, j)
  end

  estimation::Float64 = 0.0
  if copy_sketch.nodes[1] != copy_sketch.m
    estimation = (copy_sketch.m - 1) / sum(copy_sketch.nodes[1].S)
    estimation *=
      (copy_sketch.m - copy_sketch.nodes[1].inside_edges) / length(copy_sketch.m)
  else
    estimation = (copy_sketch.m - 1) / sum(copy_sketch.nodes[2].S)
    estimation *=
      (copy_sketch.m - copy_sketch.nodes[2].inside_edges) / length(copy_sketch.m)
  end

  return estimation
end # SketchMinCut
