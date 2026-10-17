# API and scientific review for the 0.1.0 candidate

Reviewed on 2026-10-17. This is a release candidate, not a tagged release or
CRAN submission. The review used code inspection, independent membership
calculations, deterministic integration tests and exact-refit comparisons.

## API decisions retained

- Data-first constructors use tidy selection for `id`, composite `key`,
  `index` and `.by`. Keys identify series; `.by` qualifies blocks and constrains
  generated dependence queries. It does not request a separate V-fold schedule
  per domain. Explicit custom memberships and existing omissions are retained.
- Assessment is nonempty and contained in omission; analysis is its complement.
  Repeated assessment contributes repeatedly to the criterion. Canonical data
  order determines all exported indices; identifiers provide alignment checks.
- Single folds remain useful native objects and rsample splits. `nei` export
  requires at least two folds for the tested mgcv integration.
- Native storage stays compact. Explicit accessors can expand memberships;
  `as_nei()` and `as_rset()` guard their different materialisation costs.
- `from_rset()` accepts `key`, `index` and `.by` selections and restores metadata
  retained by `as_rset()`. Explicit selections override saved metadata. Existing
  positional use of `id` is preserved.
- Construction metadata records group selectors, allocation seed/balance,
  environmental cluster count, spatial grid anchors, query CRS and engine.

## Issues corrected

| Area | Finding and resolution |
|---|---|
| Spatial radii | sf precision metadata could round candidate queries differently from exact distances and even exclude focal rows. Internal queries now consistently use stored coordinates; canonical geometries and precision metadata are unchanged. Zero-radius queries use actual geometries. |
| Spatial blocks | A fixed anchor could fail on negative coordinates, and finite hex grids could miss edge observations or change phase after subsetting. Square boundaries and an explicit regular hex lattice now have stable anchors, deterministic boundary assignment and allocation guards. |
| Cubble | Preserve the active geometry and source order. Missing coordinate columns can be carried through purely temporal designs as empty geometries; spatial queries still reject them. Composite matching is ready for metadata that supports it. |
| Calendar arithmetic | Zero shifts now preserve the exact instant at DST folds. Supplied origins use the data index's time zone for calendar arithmetic. Equivalent calendar units retain their original anchor. Non-advancing forecast steps fail clearly. |
| Forecast leakage | Per-key forecasting could retain future observations from another dependent series. It is now rejected when a declared domain contains multiple series; shared elapsed/calendar cutoffs are supported across those series. |
| Buffer overrides | Newly selected key/index columns now receive the same missing-value and selector checks as constructors. |
| Resampling alignment | Nested spatial-cubble row indices could silently become observation indices after expansion. Imports now require an observation-level row mapping. Round trips retain series/domain metadata so later buffering cannot silently cross domains. |
| Resampling memory | Analysis complements can be quadratic even for compact LOO designs. A preflight byte estimate rejects infeasible rsample conversion before expansion. |
| Phylogenies | Malformed, cyclic or disconnected trees are rejected before ape algorithms or clade traversal; branch-length and taxon checks are retained. |
| Diagnostics | Plots use `omitted` for omission-only observations. The documented `n_buffer`/`buffer_size` counts also include unavailable future and sliding-window exclusions. |

These complete unreleased functionality and are recorded once in the existing
New features entries in NEWS.

## mgcv 1.9-4 compatibility findings

The installed baseline was inspected and exercised, rather than assuming that
all structurally plausible lists were accepted:

1. `mgcv:::estimate.gam` rejects `ma` or `md` of length less than two. The bam
   path also emits a warning for an unsampled single-fold structure. Export and
   packed validation now reject single-fold lists with a compatibility error.
2. `mgcv:::bgam.fitd` expands a scalar `sample` using
   `sample(nei$ma, nei$sample)`, then treats those cumulative endpoints as fold
   indices. This fails for ordinary buffered neighbourhoods. Scalar counts are
   therefore supported only for singleton omissions; other designs require at
   least two explicit fold numbers, or no sampling specification. Export never
   draws a hidden sample or discards the full neighbourhood structure. The
   restriction is deterministic even when mgcv is not installed. Revisit it
   only after testing an upstream fix.

Gaussian NCV/QNCV and Poisson integration tests, plus sampled discrete bam,
remain part of routine testing. An independent weighted Gaussian calculation
also checks multiple assessment rows, repeated assessment and large omissions
against the full fixed-basis penalised normal equations.

## Large-omission approximation review

Run `Rscript benchmarks/approximation.R 120 tools/approximation-review.csv`.
The committed CSV records all eight radii for each family, omitted fractions,
convergence, criterion errors and maximum link-scale prediction discrepancy.
The seed is 20261017, basis size 10 and fixed smoothing parameter 0.2. All refits
converged. The same full-data basis/penalty is retained; rebuilding smooths on
each subset would compare different parameterisations.

Gaussian NCV and QNCV relative criterion errors were at most 3.3e-13. Selected
Poisson results are:

| Radius | Maximum omitted fraction | NCV relative error | QNCV relative error | Maximum absolute link error |
|---:|---:|---:|---:|---:|
| 7 | 12.5% | 0.23% | 0.47% | 0.112 |
| 15 | 25.8% | 6.04% | 15.10% | 0.272 |
| 30 | 50.8% | 0.37% | 25.37% | 0.176 |
| 48 | 80.8% | 15.36% | 68.98% | 0.805 |

These results support explicit approximation checks for representative models;
they do not establish a universally safe omitted fraction, monotonic error,
optimal smoothing-parameter accuracy or interval coverage. QNCV's finite
criterion is not an accuracy guarantee. The NCV paper's small-neighbourhood
asymptotics assume omission size does not grow with sample size. General
non-Gaussian transfer/block/forecast designs can require exact refitting.

The ecological prediction-target guidance distinguishes interpolation among
sampled units, transfer to new sites/groups/taxa, future information and
environmental extrapolation. Spatial or temporal blocking can inadvertently
create environmental extrapolation. A legal design alone does not ensure
training-model identifiability or prediction for unseen factor levels.

## Documented limits and remaining decisions

- Multi-series forecasts with common time cutoffs use elapsed/calendar spans.
  Count horizons require one series per explicitly independent domain; no
  common count alignment is invented for asynchronous dependent series.
- Forecasts retain observations in other `.by` domains, including their future
  observations. Independence does not establish information availability. For
  a global prospective cutoff, keep participating series in the same domain.
- cubble 1.0.0 accepts one key column in its public constructors. Composite keys
  work for ordinary tibbles, sf and tsibble inputs.
- Geographic grids crossing the antimeridian need an explicit appropriate
  projected CRS. Radius queries respect the current sf distance engine without
  changing global settings. Coordinate rounding must be explicit.
- Patristic distances are dense in unique observed taxa, protected by an
  allocation guard; observation replication does not enlarge the taxon matrix.
- Discrete bam NCV uses a working linear-model criterion and is costly for
  large omissions, even when packed output fits in memory.
- The LGCP paper remains background; there is no LGCP fitting/scoring support.

No unresolved implementation defect from this review remains. Application-level
choice of prediction target, window widths and model-specific approximation
checks remains the analyst's responsibility. Final maintainer release approval,
tagging and CRAN submission are separate, deliberately unperformed actions.

## Sources

- Roberts et al. (2017), *Ecography* 40:913–929,
  [doi:10.1111/ecog.02881](https://doi.org/10.1111/ecog.02881), checked against
  the supplied PDF, especially the prediction-target table and blocking discussion.
- Wood, *On Neighbourhood Cross Validation*, supplied version 4 and
  [arXiv:2404.16490](https://arxiv.org/abs/2404.16490).
- [mgcv NCV documentation](https://stat.ethz.ch/R-manual/R-devel/library/mgcv/html/NCV.html),
  [bam documentation](https://stat.ethz.ch/R-manual/R-devel/library/mgcv/html/bam.html)
  and the installed mgcv 1.9-4 R implementation.

Validation and performance results are recorded in `tools/VALIDATION.md`;
release gates are tracked in `tools/RELEASE.md`.
