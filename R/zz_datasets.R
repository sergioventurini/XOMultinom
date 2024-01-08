#' Data: Leukaemia cases
#'
#' This is a well-known epidemiological dataset of diagnosed leukaemia cases
#' over eight counties in upstate New York. These data originated from the
#' New York State Cancer Registry, and were gathered during the 5-year period
#' 1978–1982, with a total of 584 individuals diagnosed with leukaemia over a
#' population of approximately 1 million people. The original data contain
#' spatial information about registered events split into 790 census tracts.
#'
#' @usage data(leukaemia)
#'
#' @format{
#'   A data frame with 790 observations and the following 5
#'   variables:
#'   * **`ID`** (`int`): 10 character long identification number for a cell
#'     or census district in the study area
#'   * **`x`** (`num`): x-coordinate of the geographic centroid of each cell
#'   * **`y`** (`num`): y-coordinate of the geographic centroid of each cell
#'   * **`pop`** (`int`): 1980 U.S. Census population count for each cell
#'   * **`cases`** (`num`): incident cases of leukemia (all types) occurring
#'     between 1978 and 1982 in each cell; fractional values can occur due to
#'     partially missing data
#' }
#' 
#' @references
#'   Lange, N., Ryan, L., Billard, L., Brillinger, D., Conquest, L., 
#'   Greenhouse, J. (1994), "Case Studies in Biometry",
#'   Hoboken, NJ: Wiley & Sons.
#'   
#' @source
#'   The data set has been downloaded from <https://www.stats.ox.ac.uk/pub/datasets/csb/>.
"leukaemia"
