using Graphs, Random

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
