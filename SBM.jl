using Graphs, Random


"""
    generate_sbm_manual(n::Int, b::Int, p::Float64, q::Float64)

Generates a stochastic block model (SBM) graph.

# Parameters:
- `n::Int`: The total number of nodes in the graph.
- `b::Int`: The number of blocks.
- `p::Float64`: The intra-block connection probability.
- `q::Float64`: The inter-block connection probability.

# Returns:
- `Graph`: An instance of `Graph` representing the generated SBM.
- `Vector{Int}`: The block labels for each node.
"""
function generate_sbm_manual(n::Int, b::Int, p::Float64, q::Float64)
  block_sizes = fill(n ÷ b, b)
  labels = repeat(1:b, inner=block_sizes[1])
  shuffle!(labels)

  g = SimpleGraph(n)

  probabilities = rand(n, n)

  for i in 1:n-1
    for j in i+1:n
      prob = labels[i] == labels[j] ? p : q
      probabilities[i, j] < prob && add_edge!(g, i, j)
    end
  end

  return g, labels
end
