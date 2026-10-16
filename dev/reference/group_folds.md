# Hold out observational groups

Groups are identified by the combination of the selected columns. Select
all relevant levels, for example `group = c(region, site)`, for nested
groups. Groups never cross independent domains defined by `.by`.

## Usage

``` r
group_folds(
  data,
  group,
  v = NULL,
  balance = TRUE,
  seed = NULL,
  .by = NULL,
  id = NULL
)
```

## Arguments

- data:

  An observation data frame.

- group:

  Columns defining whole groups, using tidy selection.

- v:

  Number of folds. `NULL` leaves each group out in turn.

- balance:

  Balance observation counts when allocating groups to folds.

- seed:

  Optional local random seed.

- .by:

  Columns defining independent domains.

- id:

  Optional unique observation identifier column.

## Value

A `naboer_folds` object.
