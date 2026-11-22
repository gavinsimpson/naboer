# Ecological groups, environments and phylogeny

``` r

library("naboer")
library("tibble")
library("mgcv")
library("datasets")

have_gratia <- requireNamespace("gratia", quietly = TRUE) &&
  requireNamespace("ggplot2", quietly = TRUE)
have_ape <- requireNamespace("ape", quietly = TRUE)
have_rsample <- requireNamespace("rsample", quietly = TRUE)

if (have_gratia) {
  library("gratia")
  library("ggplot2")
}
if (have_ape) {
  library("ape")
}
if (have_rsample) {
  library("rsample")
}
```

Choose what to hold out based on where you want to make predictions. For
example, to assess predictions at a new station, fit the model without
any of that station’s observations. You can also hold out observations
from an unfamiliar environment or from taxa that the model has never
seen. When the aim is to assess a particular observation while allowing
for dependence among nearby observations, use a buffer around it. The
buffer excludes observations from fitting while leaving the set of
observations you score unchanged.

## Predict carbon dioxide uptake for a new plant

The `CO2` dataset supplied with R contains measurements from a cold
tolerance experiment on barnyard grass (*Echinochloa crus-galli*). The
researchers measured how much carbon dioxide each of 12 plants took up
at seven ambient concentrations. The plants came from Quebec and
Mississippi, and half of those from each location were chilled overnight
before the experiment. [Potvin, Lechowicz and Tardif
(1990)](https://doi.org/10.2307/1938276).

Combine origin (`Type`) and `Treatment` into one factor so the model can
fit a separate response curve for each combination.

``` r

grass <- transform(as.data.frame(CO2),
  condition = interaction(Type, Treatment)
)
```

We want to predict uptake for a new plant from one of these four
combinations. If we held out one measurement at a time, other
measurements from the same plant could help predict the one we held out.
Group by `Plant` so all seven measurements from a plant are left out and
scored together. Inspect the first of the 12 folds.

``` r

plant_folds <- group_folds(grass, group = Plant)
as_tibble(plant_folds, folds = 1)
#> # A tibble: 1 × 6
#>   fold  n_assessment n_omission n_analysis n_buffer omitted_fraction
#>   <chr>        <dbl>      <dbl>      <dbl>    <dbl>            <dbl>
#> 1 1                7          7         77        0           0.0833
```

Fit a smooth of ambient concentration (`conc`) for each combination of
origin and treatment. The separate `condition` term also lets the mean
uptake differ between combinations. Export the plant folds with
[`as_nei()`](https://gavinsimpson.github.io/naboer/reference/as_nei.md)
and pass them to [`gam()`](https://rdrr.io/pkg/mgcv/man/gam.html) so
that QNCV chooses the smoothness using approximate predictions for each
held-out plant. The resulting model uses all 12 plants.

``` r

plant_data <- fold_data(plant_folds)
plant_nei <- as_nei(plant_folds, data = plant_data)
plant_fit <- gam(uptake ~ condition + s(conc, by = condition, k = 5),
  data = plant_data, method = "QNCV", nei = plant_nei, na.action = na.fail
)
```

Use gratia’s
[`add_fitted()`](https://gavinsimpson.github.io/gratia/reference/add_fitted.html)
to evaluate the fitted curves across the observed range of
concentrations for each origin and treatment.

``` r

plant_grid <- expand.grid(
  conc = seq(min(plant_data$conc), max(plant_data$conc), length.out = 100),
  Type = levels(plant_data$Type),
  Treatment = levels(plant_data$Treatment)
)
plant_grid$condition <- with(plant_grid, interaction(Type, Treatment))
plant_curves <- add_fitted(plant_grid, plant_fit, value = "fitted_uptake")
```

Plot the mean predicted uptake for a new plant alongside the original
measurements. The curves for chilled plants are lower at both origins,
especially for plants from Mississippi.

``` r

ggplot(plant_data, aes(x = conc, colour = Treatment)) +
  geom_point(aes(y = uptake), alpha = 0.55) +
  geom_line(data = plant_curves, aes(y = fitted_uptake), linewidth = 0.9) +
  facet_wrap(vars(Type)) +
  scale_colour_manual(
    values = c(nonchilled = "#0072B2", chilled = "#D55E00"),
    labels = c(nonchilled = "Non-chilled", chilled = "Chilled")
  ) +
  labs(
    x = "Ambient carbon dioxide concentration (µL/L)",
    y = "Carbon dioxide uptake (µmol/m²/s)"
  ) +
  theme_minimal()
```

![Carbon dioxide uptake curves and observed measurements for Quebec and
Mississippi grasses. The fitted curves are lower for chilled plants,
especially those from
Mississippi.](ecological_files/figure-html/grass-curves-1.png)

## Whole stations and hierarchical groups

The remaining examples use small invented datasets to show how the folds
are constructed. In this monitoring example, six stations belong to two
networks. Each station has four daily observations and its own
temperature and moisture measurements, which are constant across the
four days.

``` r

observations <- tibble(
  observation = seq_len(24),
  network = rep(c("north", "south"), each = 12),
  station = rep(paste0("station", seq_len(6)), each = 4),
  date = rep(as.Date("2025-05-01") + 0:3, 6),
  temperature = rep(c(4, 6, 12, 14, 20, 22), each = 4),
  moisture = rep(c(90, 85, 65, 60, 35, 30), each = 4),
  habitat = rep(
    c("wet", "wet", "intermediate", "intermediate", "dry", "dry"), each = 4
  )
)
```

Use
[`group_folds()`](https://gavinsimpson.github.io/naboer/reference/group_folds.md)
to leave out one station at a time. Selecting both `network` and
`station` makes the pair the group identifier. This also works when
different networks use the same station names.

``` r

new_station <- observations |>
  group_folds(group = c(network, station), id = observation)
```

The summary shows six folds, each assessing four observations. Use
[`assessment_indices()`](https://gavinsimpson.github.io/naboer/reference/assessment_indices.md)
to see which rows contribute to the error in a fold. The indices are
positions in the original table.

``` r

as_tibble(new_station)
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

To assess multiple stations in each fold, set `v` to the number of
folds. The stations stay whole, and `balance = TRUE` aims to give the
folds similar numbers of observations. Setting `seed` makes this
allocation reproducible without changing your session’s random number
state.

``` r

station_vfold <- observations |>
  group_folds(group = c(network, station), v = 3, balance = TRUE, seed = 21)
as_tibble(station_vfold)
#> # A tibble: 3 × 6
#>   fold  n_assessment n_omission n_analysis n_buffer omitted_fraction
#>   <chr>        <dbl>      <dbl>      <dbl>    <dbl>            <dbl>
#> 1 1                8          8         16        0            0.333
#> 2 2                8          8         16        0            0.333
#> 3 3                8          8         16        0            0.333
```

Use `group = network` if the prediction target is a new network. If the
networks are independent for the purpose of defining groups or adding
buffers, use `.by = network`. This keeps group boundaries and buffers
within networks. When you request a set number of folds with `v`, a fold
can still assess whole groups from several networks.

## New environmental conditions

To assess predictions under unfamiliar conditions, use
[`environmental_folds()`](https://gavinsimpson.github.io/naboer/reference/environmental_folds.md)
to hold out observations with similar predictor values. Here `v = 3`
asks k-means to form three clusters from temperature and moisture. The
function standardises each predictor to unit standard deviation by
default, so the distance calculation uses a common scale for both. Each
cluster forms one fold, and clusters can differ in size.

``` r

new_environment <- observations |>
  environmental_folds(predictors = c(temperature, moisture), v = 3, seed = 22)
as_tibble(new_environment)
#> # A tibble: 3 × 6
#>   fold  n_assessment n_omission n_analysis n_buffer omitted_fraction
#>   <chr>        <dbl>      <dbl>      <dbl>    <dbl>            <dbl>
#> 1 1                8          8         16        0            0.333
#> 2 2                8          8         16        0            0.333
#> 3 3                8          8         16        0            0.333
```

Choose predictors that describe the environmental differences you want
to test. Supply finite numeric values for each predictor. The function
reports and drops constant predictors. It needs at least one varying
predictor and at least `v` distinct combinations of predictor values. If
you supply `.by`, these checks and the clustering apply separately
within each domain, with `v` clusters in each.

If you already have environmental categories, supply them through
`strata`. This example assesses wet, intermediate and dry habitats in
turn. To combine whole habitats into fewer folds, supply `v`. In that
case, `balance` controls whether allocation balances the number of
observations or the number of habitats.

``` r

new_habitat <- observations |>
  environmental_folds(strata = habitat)
as_tibble(new_habitat)
#> # A tibble: 3 × 6
#>   fold  n_assessment n_omission n_analysis n_buffer omitted_fraction
#>   <chr>        <dbl>      <dbl>      <dbl>    <dbl>            <dbl>
#> 1 1                8          8         16        0            0.333
#> 2 2                8          8         16        0            0.333
#> 3 3                8          8         16        0            0.333
```

These folds may require predictions outside the range of conditions seen
in training. If predictor values change over time at a station,
clustering can put its observations in different folds. When you need to
keep stations together, use categories that are constant within each
station. If you want only to exclude nearby observations from fitting,
add a buffer to the folds.

## Phylogenetic dependence and new taxa

These examples run when the optional `ape` package is installed. Use it
to make a tree of four taxa, then create a table with six observations
per taxon. All observed taxa must match tip labels in the tree, and each
branch needs a finite, non-negative length. The tree may include unused
tips.

``` r

tree <- read.tree(
  text = "((oak:1,beech:1):2,(pine:1,fir:1):2);"
)
taxon_data <- tibble(
  observation = seq_len(24),
  taxon = rep(c("oak", "beech", "pine", "fir"), each = 6)
)
```

Tree distance is the sum of branch lengths along the path between two
taxa. Oak and beech are two units apart, and oak and pine are six units
apart. To assess one observation at a time while excluding its close
relatives from training, set a radius in those same units. A radius of
two includes taxa exactly two units away.

``` r

near_taxa <- taxon_data |>
  neighbourhoods(
    phylogenetic = phylogenetic_window(tree, taxon, radius = 2),
    id = observation
  )
```

Inspect the first fold. It scores the first oak observation and leaves
all oak and beech observations out of the fit. The exclusion applies
within the same domain when you use `.by`. The function calculates
distances once between taxa and reuses them across observations of each
taxon.

``` r

assessment_indices(near_taxa, folds = 1)
#> [[1]]
#> [1] 1
omitted <- omission_indices(near_taxa, folds = 1)[[1]]
unique(taxon_data$taxon[omitted])
#> [1] "oak"   "beech"
```

To assess all observations from a group of taxa together, use
[`phylogenetic_folds()`](https://gavinsimpson.github.io/naboer/reference/phylogenetic_folds.md).
Supply `v` for a set number of groups or `height` for a maximum distance
between any pair of taxa in a group. Here `height = 2` groups oak with
beech and pine with fir. All observations of the held-out taxa are
scored.

``` r

taxon_blocks <- taxon_data |>
  phylogenetic_folds(tree, taxon, height = 2)
as_tibble(taxon_blocks)
#> # A tibble: 2 × 6
#>   fold  n_assessment n_omission n_analysis n_buffer omitted_fraction
#>   <chr>        <dbl>      <dbl>      <dbl>    <dbl>            <dbl>
#> 1 1               12         12         12        0              0.5
#> 2 2               12         12         12        0              0.5
```

If you use `.by`, a fold holds out a group within one domain, and the
same taxon can remain in training in another. Leave `.by` unspecified to
assess taxa that are absent from the entire training dataset.

You can also define the groups directly from branches of a rooted tree.
Supply the node number or the complete set of descendant tip labels for
each clade. The clades must cover every observed taxon, and each tip can
belong to only one selected clade. Here we specify the two pairs
explicitly.

``` r

new_clade <- taxon_data |>
  phylogenetic_folds(
    tree, taxon,
    clades = list(c("oak", "beech"), c("pine", "fir"))
  )
as_tibble(new_clade)
#> # A tibble: 2 × 6
#>   fold  n_assessment n_omission n_analysis n_buffer omitted_fraction
#>   <chr>        <dbl>      <dbl>      <dbl>    <dbl>            <dbl>
#> 1 1               12         12         12        0              0.5
#> 2 2               12         12         12        0              0.5
```

## Buffered resampling and rsample

Start with leave-one-out folds for the monitoring data. Each observation
is scored in its own fold.

``` r

leave_one_out <- loo_folds(observations, id = observation)
```

Add a one-day buffer around the assessment observation in each fold.
With `key = station` and `within = "key"`, the buffer stays within the
same station. Distances are measured from the assessment observation
itself.

``` r

buffered <- buffer_folds(
  leave_one_out,
  temporal = temporal_window(time_span(1, "days")),
  key = station, index = date, within = "key"
)
```

In fold 2, row 2 is still the only assessment observation. The buffer
also excludes rows 1 and 3 from fitting, so `n_buffer` is 2 and the
analysis set starts at row 4.

``` r

as_tibble(buffered, folds = 2)
#> # A tibble: 1 × 6
#>   fold  n_assessment n_omission n_analysis n_buffer omitted_fraction
#>   <chr>        <dbl>      <dbl>      <dbl>    <dbl>            <dbl>
#> 1 2                1          3         21        2            0.125
analysis_indices(buffered, folds = 2)
#> [[1]]
#>  [1]  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24
```

The next examples require the optional **rsample** package. Export with
[`as_rset()`](https://gavinsimpson.github.io/naboer/reference/as_rset.md)
to refit the model on each analysis set using rsample. The exported
resamples store the analysis and assessment sets separately. Rows
excluded by the buffer alone belong to neither set.

``` r

resamples <- as_rset(buffered)
```

Import with
[`from_rset()`](https://gavinsimpson.github.io/naboer/reference/as_rset.md)
to get a naboer design back. It treats all rows outside the analysis set
as omitted, so the buffer survives the round trip. These comparisons
check all three sets.

``` r

restored <- from_rset(resamples)
identical(assessment_indices(restored), assessment_indices(buffered))
#> [1] TRUE
identical(omission_indices(restored), omission_indices(buffered))
#> [1] TRUE
identical(analysis_indices(restored), analysis_indices(buffered))
#> [1] TRUE
```

The export also saves which columns hold series keys, time indices,
domains and observation identifiers. Import restores those selections
for use in later operations. For resamples created elsewhere, supply
`key`, `index`, `.by` and `id` as needed. Explicit selections take
precedence over saved ones.

You can import ordinary and repeated V-fold resamples. The splits must
all refer to the same data in the same row order, and each analysis set
must contain unique rows. Bootstrap samples that repeat a row cannot be
represented by this design.

Create resamples on individual observation rows. With a nested spatial
cubble, use the cubble package’s `face_temporal()` to switch to one row
per observation, or resample the observation table returned by
[`fold_data()`](https://gavinsimpson.github.io/naboer/reference/fold_data.md).
If a row represents an entire station, its index cannot identify a
single observation after expansion.

Exporting to rsample can take much more memory than storing a naboer
design. For example, leave-one-out stores almost all row indices in each
rsample analysis set.
[`as_rset()`](https://gavinsimpson.github.io/naboer/reference/as_rset.md)
estimates the memory needed for the analysis and assessment indices
before creating them. The default limit comes from
`getOption("naboer.max_bytes", 1024^3)`, and you can change it with
`max_bytes`. The estimate counts integer indices and excludes R object
overhead.

For neighbourhood cross-validation in mgcv, pass `as_nei(buffered)` as
`nei` and fit to `fold_data(buffered)`. Prepare complete modelling data
before constructing folds and use `na.action = na.fail` when fitting.
Check whether the approximation is accurate for your model and omission
size. See [Assessment, omission and prediction
targets](https://gavinsimpson.github.io/naboer/articles/naboer.md) for
an example.
