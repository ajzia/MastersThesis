module Sketch

using Random

export EdgeSketch, NodeSketch
export EdgeSketchCreate

"""
    NodeSketch(id::Int, m::Int)

Create a new node sketch with the given ID and sketch size `m`.
"""
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


"""
    EdgeSketch(m::Int)
    EdgeSketch(m::Int, edges::Vector{Tuple{Int, Int, Float64}})

Create a new edge sketch with size `m` and initialize it with
the given edges.
"""
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


"""
    EdgeSketchCreate(m::Int, stream::String)::EdgeSketch

Create an edge sketch from a stream of edges.

# Parameters:
- `m::Int`: The size of the sketch.
- `stream::String`: The path to the file containing edges
                    in the format "i j w" or "i,j,w".

# Returns:
- `EdgeSketch`: An instance of `EdgeSketch` initialized with
                the edges from the stream.
"""
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


"""
    EdgeSketchInit(sketch::EdgeSketch, edge::Tuple{Int, Int, Float64}, is_directed::Bool=false)

Initialize the edge sketch with a new edge.

# Parameters:
- `sketch::EdgeSketch`: The edge sketch to be initialized.
- `edge::Tuple{Int, Int, Float64}`: The edge to be added.
- `is_directed::Bool`: Whether the edge is directed or not (default=false).
"""
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



"""
    Hash(value::String)::Float64

Generate a hash for a given string value, returning a
random number between 0 and 1.

# Parameters:
- `value::String`: The string value to be hashed.

# Returns:
- `Float64`: A random number between 0 and 1 generated from the hash.
"""
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
