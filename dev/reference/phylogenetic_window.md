# Specify a phylogenetic neighbourhood

Distances use branch lengths between observed taxa, not distances
between replicated observation rows. Unused tree tips are allowed.

## Usage

``` r
phylogenetic_window(tree, taxon, radius)
```

## Arguments

- tree:

  An `ape` phylogenetic tree with finite, nonnegative branch lengths.

- taxon:

  Column containing tree tip labels, using tidy selection.

- radius:

  Inclusive patristic distance threshold.

## Value

A phylogenetic window specification for
[`neighbourhoods()`](https://gavinsimpson.github.io/naboer/dev/reference/neighbourhoods.md).
