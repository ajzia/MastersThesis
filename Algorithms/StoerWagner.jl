"""
    StoerWagnerMinCut(adj_matrix::Matrix{Int})

Finds the minimum cut of a graph using the Stoer-Wagner algorithm.

# Parameters:
- `adj_matrix::Matrix{Int}`: Adjacency matrix of the graph, where the
                             element at position `(i, j)` represents
                             the weight of the edge between nodes `i`
                             and `j`.

# Returns:
- `Tuple{Float64, Set{Int}}`: A tuple containing the weight of the
                              minimum cut and a set of nodes on one
                              side of the cut.
"""
function StoerWagnerMinCut(adj_matrix::Matrix{Int})::Tuple{Float64, Set{Int}}
  n::Int = size(adj_matrix, 1)
  nodes::Vector{Int} = collect(1:n)
  best_weight::Float64 = typemax(Float64)
  best_set::Set{Int} = Set{Int}()

  while length(nodes) > 1
    A::Vector{Int} = [nodes[1]]
    weights::Vector{Float64} = zeros(Float64, n)
    used::Vector{Bool} = falses(n)
    used[nodes[1]] = true

    @inbounds for _ in 2:length(nodes)
      for v in nodes
        if !used[v]
          weights[v] += adj_matrix[A[end], v]
        end
      end

      max_weight::Float64 = -1
      selected::Int = -1
      @inbounds for v in nodes
        if !used[v] && weights[v] > max_weight
          max_weight = weights[v]
          selected = v
        end
      end

      push!(A, selected)
      used[selected] = true
    end

    s::Int, t::Int = A[end - 1], A[end]
    cut_weight::Float64 = 0
    @inbounds for v in nodes
      cut_weight += adj_matrix[t, v]
    end

    if cut_weight < best_weight
      best_weight = cut_weight
      best_set = Set(A[1:end-1])
    end

    @inbounds for v in nodes
      if v != s && v != t
        adj_matrix[s, v] += adj_matrix[t, v]
        adj_matrix[v, s] = adj_matrix[s, v]
      end
    end

    filter!(x -> x != t, nodes)
  end

  return best_weight, best_set
end
