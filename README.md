using JSON
using Random
using StatsBase
using Graphs

# Required packages

To run the code, you need to have the following packages installed:

- JSON
- Random
- StatsBase
- Graphs

# Usage

There are three ways to run the code:

1. By generating a random graph using SBM model and running Karger's algorithm on it.

```bash
julia main.jl GENERATE n b p q iter ms
```

Where:

- `n`: Number of nodes in the graph.
- `b`: Number of blocks in the graph.
- `p`: Probability of intra-block edges.
- `q`: Probability of inter-block edges.
- `iter`: Number of iterations to run Karger's algorithm.
- `ms`: A vector of integers representing sizes of the sketches in format `m1,m2,...,mk`.

Generated graphs will be saved in the `./Resources` directory, and results will be saved in a json file in the `./Results` directory.

2. By reading a graph from a file and running Karger's algorithm on it.

```bash
julia main.jl FROM_FILE path/to/file iter ms
```

Where:

- `path/to/file`: Path to the file containing the graph edges in format "i,j,w" or "i j w".
- `iter`: Number of iterations to run Karger's algorithm.
- `ms`: A vector of integers representing sizes of the sketches in format `m1,m2,...,mk`.

Results will be saved in a json file in the `./Results` directory.

Optionally you can specify the HISTOGRAM flag at the end of both commands to generate a histogram of the results.

3. You can plot average errors for different sketch sizes of the same graph using the results saved in the `./Results` directory.

```bash
julia plot_errors.jl path/to/file
```

Where:
- `path/to/file`: Path to the json file containing the results.
