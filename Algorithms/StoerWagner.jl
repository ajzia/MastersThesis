# Finding opt min-cut
function StoerWagnerMinCut(adj_matrix::Matrix{Int})
  n = size(adj_matrix, 1)
  nodes = collect(1:n)
  best_cut_weight = typemax(Int)
  best_cut_set = Set{Int}()

  while length(nodes) > 1
    A = [nodes[1]]
    weights = zeros(Int, n)
    used = falses(n)
    used[nodes[1]] = true

    for _ in 2:length(nodes)
      # Find the most tightly connected unused node
      for v in nodes
        if !used[v]
          weights[v] += adj_matrix[A[end], v]
        end
      end

      # Choose the next node with max weight
      max_w = -1
      sel = -1
      for v in nodes
        if !used[v] && weights[v] > max_w
          max_w = weights[v]
          sel = v
        end
      end

      push!(A, sel)
      used[sel] = true
    end

    # Last two added nodes are s and t
    s = A[end - 1]
    t = A[end]
    cut_weight = 0
    for v in nodes
      cut_weight += adj_matrix[t, v]
    end

    if cut_weight < best_cut_weight
      best_cut_weight = cut_weight
      best_cut_set = Set(A[1:end-1])
    end

    # Merge s and t into one node
    for v in nodes
      if v != s && v != t
        adj_matrix[s, v] += adj_matrix[t, v]
        adj_matrix[v, s] = adj_matrix[s, v]
      end
    end

    # Remove t from the graph
    filter!(x -> x != t, nodes)
  end

  return best_cut_weight, best_cut_set
end
