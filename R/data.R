#' Daily soil moisture at Balruddery in 2019
#'
#' Daily volumetric soil water content from the COSMOS-UK monitoring station
#' at Balruddery, near Dundee. The extract contains all 365 daily records for
#' station BALRD in 2019. Contains data supplied by UK Centre for Ecology &
#' Hydrology.
#'
#' @format A data frame with 365 rows and two columns:
#'   - `date`: observation date as a `Date`.
#'   - `vwc`: volumetric soil water content, expressed as a percentage.
#' @source Stanley et al. (2026), Daily and sub-daily hydrometeorological and
#'   soil moisture data (2013–2025) (COSMOS-UK), NERC EDS Environmental
#'   Information Data Centre. [Data catalogue](https://doi.org/10.5285/dde10e2b-ee1a-4b2b-a247-176db2f7bb31).
#'   [Open Government Licence v3](https://eidc.ac.uk/licences/ogl/plain).
#' @examples
#' data("balruddery_2019", package = "naboer")
#' head(balruddery_2019)
"balruddery_2019"
