module Sketch

using Random

export EdgeSketch, NodeSketch
export EdgeSketchCreate


mutable struct NodeSketch
  id::Int
  F::Vector{Tuple{Int, Int}}
  S::Vector{Float64}
  max::Float64
  inside_edges::Int
  empty_edges::Int
  estimated_weight::Float64

  NodeSketch(id::Int, m::Int) = new(id, [(0, 0) for _ in 1:m], [Inf for _ in 1:m],  Inf, 0, m, 0) 
end # NodeSketch


struct EdgeSketch
  m::Int # sketch size
  node_ids::Vector{Int}
  nodes::Vector{NodeSketch}

  EdgeSketch(m::Int) = new(m, [], [])
  EdgeSketch(m::Int, edges::Vector{Tuple{Int, Int, Float64}}) = begin
    sketch::EdgeSketch = EdgeSketch(m)
    for edge in edges
      EdgeSketchInit(sketch, edge)
    end
    return sketch
  end
end # EdgeSketch


function EdgeSketchCreate(m::Int, stream::String)::EdgeSketch
  sketch::EdgeSketch = EdgeSketch(m)::EdgeSketch

  # read stream of edges
  open(stream) do file
    for line in eachline(file)
      if startswith(line, "%") continue end
      if startswith(line, "#") continue end
      edge = split(line, r"[\s,]")
      i = parse(Int, edge[1])
      j = parse(Int, edge[2])

      if i == j continue end
      w = parse(Float64, edge[3])
      EdgeSketchInit(sketch, (i, j, w))
    end
  end

  return sketch
end # EdgeSketchCreate


function EdgeSketchInit(sketch::EdgeSketch, edge::Tuple{Int, Int, Float64}, is_directed::Bool=false)
  (i::Int, j::Int, _) = edge

  if i ∉ sketch.node_ids
    push!(sketch.node_ids, i)
    push!(sketch.nodes, NodeSketch(i, sketch.m))
  end

  if j ∉ sketch.node_ids && !is_directed
    push!(sketch.node_ids, j)
    push!(sketch.nodes, NodeSketch(j, sketch.m))
  end

  # update sketches
  node_i = findfirst(x -> x.id == i, sketch.nodes)
  sketch.nodes[node_i] = EdgeSketchUpdate(sketch.m, sketch.nodes[node_i], edge, is_directed)

  if !is_directed
    node_j = findfirst(x -> x.id == j, sketch.nodes)
    sketch.nodes[node_j] = EdgeSketchUpdate(sketch.m, sketch.nodes[node_j], edge, is_directed)
  end
end # EdgeSketchInit


@inline function Hash(value::String)::Float64
  # Hashing function to generate a random number
  return Float64((hash((value, 123)) % Int64(1e9))) / Float64(1e9)
end # Hash


function EdgeSketchUpdate(
  m::Int,
  node::NodeSketch,
  edge::Tuple{Int, Int, Float64},
  is_directed::Bool=false
)::NodeSketch
  # Initialization
  (i::Int, j::Int, w_ij::Float64) = edge
  P::Vector{Int} = collect(1:m)

  if !is_directed
    i, j = min(i, j), max(i, j)
  end
  
  sum::Float64 = 0.0
  update_max::Bool = false
  MAX::Float64 = node.max
  Random.seed!(hash(bitstring(i) * bitstring(j)))
  
  # Update
  @inbounds for k in 1:m
    value::String = bitstring(i) * bitstring(j) * bitstring(k)
    U::Float64 = Hash(value)
    E::Float64 = - log(U) / w_ij

    sum += E / (m - k + 1)
    if sum >= MAX break end

    r::Int = rand(k:m)
    P[k], P[r] = P[r], P[k]

    l::Int = P[k]

    if node.S[l] == MAX
      update_max = true
    end

    if sum < node.S[l]
      # Change number of empty edges
      if node.S[l] == Inf
        node.empty_edges -= 1
      end

      # Update the sketch
      node.S[l] = sum
      node.F[l] = (i, j)
    end
  end # for

  if update_max
    node.max = maximum(node.S)
  end

  return node
end # EdgeSketchUpdate

end # module
