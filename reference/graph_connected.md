# Check whether an undirected adjacency matrix is fully connected

Breadth-first reachability from the first node. Replaces a single-use
dependency on Claddis::is_graph_connected (the only reason hestia
imported Claddis, whose heavy dependency tree failed to build on R-devel
macOS).

## Usage

``` r
graph_connected(graph)
```

## Arguments

- graph:

  square symmetric 0/1 adjacency matrix.

## Value

TRUE if every node is reachable from the first, else FALSE.
