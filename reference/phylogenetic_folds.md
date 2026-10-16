# Construct phylogenetic cross-validation folds

Supply `v` or `height` to cut complete-linkage clustering of patristic
distances. Alternatively, provide explicitly selected nonoverlapping
clades of a rooted tree. A clade is a tree node number or a character
vector of all descendant tip labels. Selected clades must cover all
observed taxa.

## Usage

``` r
phylogenetic_folds(
  data,
  tree,
  taxon,
  v = NULL,
  height = NULL,
  clades = NULL,
  .by = NULL,
  id = NULL
)
```

## Arguments

- data:

  An observation data frame.

- tree:

  An `ape` phylogenetic tree with finite, nonnegative branch lengths.

- taxon:

  Column containing tree tip labels, using tidy selection.

- v:

  Number of patristic-distance clusters.

- height:

  Maximum complete-linkage cluster height.

- clades:

  A list of nodes or complete sets of descendant tip labels.

- .by:

  Columns defining independent domains.

- id:

  Optional unique observation identifier column.
