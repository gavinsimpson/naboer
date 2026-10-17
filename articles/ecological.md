# Ecological groups, environments and phylogeny

Validation should reflect the ecological prediction task. Predicting
another observation at a sampled station differs from predicting a new
station, an unfamiliar environment or an unobserved taxon. Dependence
buffers address the relationship between held-out and training
observations; whole-group holdouts also change which observational units
are available to the model.

This article uses a small monitoring dataset with repeated observations
from six stations in two networks. All examples preserve its original
observation order.

``` r

observations <- tibble::tibble(
  observation = seq_len(24),
  network = rep(c("north", "south"), each = 12),
  station = rep(paste0("station", seq_len(6)), each = 4),
  date = rep(as.Date("2025-05-01") + 0:3, 6),
  temperature = rep(c(4, 6, 12, 14, 20, 22), each = 4),
  moisture = rep(c(90, 85, 65, 60, 35, 30), each = 4),
  habitat = rep(c("wet", "wet", "intermediate", "intermediate", "dry", "dry"), each = 4)
)
```

## Whole stations and hierarchical groups

Leave one station out to assess predictions at stations unavailable
during fitting. Select all identifying levels when station labels are
only unique within a network. The combination of the selected columns
defines a group.

``` r

new_station <- observations |>
  group_folds(group = c(network, station), id = observation)
tibble::as_tibble(new_station)
#> # A tibble: 6 × 6
#>   fold  n_assessment n_omission n_analysis n_buffer omitted_fraction
#>   <chr>        <dbl>      <dbl>      <dbl>    <dbl>            <dbl>
#> 1 1                4          4         20        0            0.167
#> 2 2                4          4         20        0            0.167
#> 3 3                4          4         20        0            0.167
#> 4 4                4          4         20        0            0.167
#> 5 5                4          4         20        0            0.167
#> 6 6                4          4         20        0            0.167
assessment_indices(new_station, folds = 1)
#> [[1]]
#> [1] 1 2 3 4
```

For V-fold CV, whole groups can be allocated together. With
`balance = TRUE`, allocation aims to balance observation counts rather
than the number of groups. The local seed makes allocation reproducible
and preserves the caller’s random number state.

``` r

station_vfold <- observations |>
  group_folds(group = c(network, station), v = 3, balance = TRUE, seed = 21)
tibble::as_tibble(station_vfold)
#> # A tibble: 3 × 6
#>   fold  n_assessment n_omission n_analysis n_buffer omitted_fraction
#>   <chr>        <dbl>      <dbl>      <dbl>    <dbl>            <dbl>
#> 1 1                8          8         16        0            0.333
#> 2 2                8          8         16        0            0.333
#> 3 3                8          8         16        0            0.333
```

Use `.by = network` when networks define independent dependence domains.
A constructor then treats a station in one network as distinct from a
station with the same label in another. V-fold allocation may score
groups from several domains in the same fold; it never splits a group.
Choose `group = network` to hold out an entire network instead.

## New environmental conditions

Environmental folds hold out clusters of observations in selected
predictor space. Numeric predictors are standardised by default before
seeded k-means clustering. Here temperature and moisture have different
scales, so standardisation gives both an equal role in the distance
calculation. Each cluster is held out as a whole, so fold sizes may
differ. `balance` applies to the existing-strata workflow below; it does
not impose balanced k-means clusters. With independent domains, `v`
requests that many clusters per domain.

``` r

new_environment <- observations |>
  environmental_folds(predictors = c(temperature, moisture), v = 3, seed = 22)
tibble::as_tibble(new_environment)
#> # A tibble: 3 × 6
#>   fold  n_assessment n_omission n_analysis n_buffer omitted_fraction
#>   <chr>        <dbl>      <dbl>      <dbl>    <dbl>            <dbl>
#> 1 1                8          8         16        0            0.333
#> 2 2                8          8         16        0            0.333
#> 3 3                8          8         16        0            0.333
```

Clustering uses the selected predictors, so its interpretation depends
on that selection. Missing or nonfinite predictors cause an error.
Predictors constant within a domain are reported and excluded; an
entirely constant domain cannot be clustered. With `.by`, each domain is
clustered separately.

Existing scientific strata can be more useful than data-derived
clusters. Supply `strata` instead of `predictors` to leave each stratum
out. An optional `v` allocates whole strata to fewer folds.

``` r

new_habitat <- observations |>
  environmental_folds(strata = habitat)
tibble::as_tibble(new_habitat)
#> # A tibble: 3 × 6
#>   fold  n_assessment n_omission n_analysis n_buffer omitted_fraction
#>   <chr>        <dbl>      <dbl>      <dbl>    <dbl>            <dbl>
#> 1 1                8          8         16        0            0.333
#> 2 2                8          8         16        0            0.333
#> 3 3                8          8         16        0            0.333
```

Environmental holdouts may involve extrapolation. Their errors answer a
different question from random V-fold errors, even when both use the
same model and data. Clustering does not automatically keep whole
stations together: observations from one station can enter different
folds if their predictors change over time. Use strata constant within
the required observational units, or add dependence buffers, when those
restrictions belong to the prediction task.

## Phylogenetic dependence and new taxa

The following examples require the optional **ape** package. Their
chunks run only when
[`requireNamespace("ape", quietly = TRUE)`](https://rdrr.io/r/base/ns-load.html)
succeeds. Install ape to execute the tree examples locally; the other
examples do not need it.

``` r

tree <- ape::read.tree(
  text = "((oak:1,beech:1):2,(pine:1,fir:1):2);"
)
taxon_data <- tibble::tibble(
  observation = seq_len(24),
  taxon = rep(c("oak", "beech", "pine", "fir"), each = 6),
  response = seq_len(24)
)

near_taxa <- taxon_data |>
  neighbourhoods(
    phylogenetic = phylogenetic_window(tree, taxon, radius = 2),
    id = observation
  )
tibble::as_tibble(near_taxa, folds = c(1, 13), membership = TRUE)
#> # A tibble: 2 × 8
#>   fold  n_assessment n_omission n_analysis n_buffer omitted_fraction assessment
#>   <chr>        <dbl>      <dbl>      <dbl>    <dbl>            <dbl> <list>    
#> 1 1                1         12         12       11              0.5 <int [1]> 
#> 2 13               1         12         12       11              0.5 <int [1]> 
#> # ℹ 1 more variable: omission <list>
```

The radius uses branch-length units and includes its boundary. In this
tree, oak and beech are two units apart, while oak and pine are six
units apart. Each focal observation omits all observations of its own
taxon and nearby taxa in its domain. Tree distances are calculated among
unique observed taxa, then mapped to observation rows; replicated taxa
do not enlarge the distance matrix. Unused tree tips are allowed, but
every observed taxon must match a tip label. Every branch needs a
finite, nonnegative length.

For new-taxon prediction, complete-linkage clustering produces blocks of
taxa. Choose either a number of clusters with `v`, or a distance
threshold with `height`. Replicated observations of a taxon stay
together. If `.by` defines independent domains, taxon holdouts apply
separately within each domain; the same taxon may remain in training in
another domain. Leave `.by` unspecified when the target is a taxon
absent from the entire training dataset.

``` r

taxon_blocks <- taxon_data |>
  phylogenetic_folds(tree, taxon, height = 2)
tibble::as_tibble(taxon_blocks)
#> # A tibble: 2 × 6
#>   fold  n_assessment n_omission n_analysis n_buffer omitted_fraction
#>   <chr>        <dbl>      <dbl>      <dbl>    <dbl>            <dbl>
#> 1 1               12         12         12        0              0.5
#> 2 2               12         12         12        0              0.5
```

Alternatively, select scientifically meaningful clades explicitly. The
tree must be rooted, and the selected nonoverlapping clades must cover
every observed taxon. Supply complete descendant tip sets or tree-node
numbers. Arbitrary tip sets that do not form a complete clade are
rejected.

``` r

new_clade <- taxon_data |>
  phylogenetic_folds(
    tree, taxon,
    clades = list(c("oak", "beech"), c("pine", "fir"))
  )
tibble::as_tibble(new_clade, membership = TRUE)
#> # A tibble: 2 × 8
#>   fold  n_assessment n_omission n_analysis n_buffer omitted_fraction assessment
#>   <chr>        <dbl>      <dbl>      <dbl>    <dbl>            <dbl> <list>    
#> 1 1               12         12         12        0              0.5 <int [12]>
#> 2 2               12         12         12        0              0.5 <int [12]>
#> # ℹ 1 more variable: omission <list>
```

## Buffered resampling and rsample

Start with a focal assessment design and add a one-day exclusion buffer
within each station. This changes omission membership while retaining
the assessment observations that contribute to the criterion. Buffers
expand from assessment observations only; they do not expand recursively
from other buffer rows.

``` r

buffered <- observations |>
  loo_folds(id = observation) |>
  buffer_folds(
    temporal = temporal_window(time_span(1, "days")),
    key = station, index = date, within = "key"
  )
tibble::as_tibble(buffered, folds = 2, membership = TRUE)
#> # A tibble: 1 × 8
#>   fold  n_assessment n_omission n_analysis n_buffer omitted_fraction assessment
#>   <chr>        <dbl>      <dbl>      <dbl>    <dbl>            <dbl> <list>    
#> 1 2                1          3         21        2            0.125 <int [1]> 
#> # ℹ 1 more variable: omission <list>
analysis_indices(buffered, folds = 2)
#> [[1]]
#>  [1]  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24
```

The next chunk requires optional **rsample** and is skipped when it is
unavailable. The exported manual resamples preserve analysis and
assessment independently. Buffer-only observations belong to neither
set. On import, omission is the complement of analysis, so those
observations remain excluded. The exported rset also saves series keys,
temporal index, independent-domain columns and identifiers.
[`from_rset()`](https://gavinsimpson.github.io/naboer/reference/as_rset.md)
restores them, so subsequent buffers respect the original domain
boundaries. For external resamples supply `key`, `index`, `.by` and `id`
as needed; explicit selections override saved metadata. Saved column
names and recovered membership values are validated during import.

``` r

resamples <- as_rset(buffered)
restored <- from_rset(resamples, id = observation)
identical(assessment_indices(restored), assessment_indices(buffered))
#> [1] TRUE
identical(omission_indices(restored), omission_indices(buffered))
#> [1] TRUE
identical(analysis_indices(restored), analysis_indices(buffered))
#> [1] TRUE
```

Ordinary rsample V-fold and repeated assessment designs can be imported.
Bootstrap analysis samples containing duplicated rows cannot be
expressed by an omission set and are rejected. All splits must refer to
identical data in the same order.

Create rsample splits on **individual observation rows**. For cubble
data, convert a nested spatial cubble with
[`cubble::face_temporal()`](https://huizezhang-sherry.github.io/cubble/reference/face.html)
before resampling, or resample the canonical table returned by
[`fold_data()`](https://gavinsimpson.github.io/naboer/reference/fold_data.md).
Importing splits whose station rows would expand into several temporal
observations is rejected: their row indices would otherwise refer to
different observations after expansion. Resampling the canonical
observation table and exporting with
[`as_rset()`](https://gavinsimpson.github.io/naboer/reference/as_rset.md)
both preserve observation order.

Conversion to rsample materialises analysis indices, which can be much
larger than compact omissions. For example, all leave-one-out analysis
sets grow quadratically with the observation count.
[`as_rset()`](https://gavinsimpson.github.io/naboer/reference/as_rset.md)
checks their estimated integer storage against
`getOption("naboer.max_bytes", 1024^3)` before expanding them. Use its
`max_bytes` argument to choose a different limit after inspecting the
design; the estimate excludes R object overhead.

Finally, export the design for mgcv with `as_nei(buffered)`, and fit
against `fold_data(buffered)`. Construct folds after preparing complete
modelling data and use `na.action = stats::na.fail` when fitting. A
valid design still requires checking whether the NCV approximation is
accurate for the intended model and omission size.
