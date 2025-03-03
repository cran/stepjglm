#'
#'
#'  Data from Injection molding experiment
#'
#' The experiment was performed to study the influence of seven controllable factors and
#' three noise factors on the mean value and the variation in the percentage of shrinkage of
#' products made by injection molding.
#'
#' @docType data
#'
#' @usage data(injection_molding)
#'
#' @format A data frame containing 32 rows and 11 variables.
#'
#'The responses were percentages of shrinkage of products made by
#'injection molding (Y).
#'
#'Controllable factors:
#'\itemize{
#' \item A: cycle time
#' \item B: mould temperature
#' \item C: cavity thickness
#' \item D: holding pressure
#' \item E: injection speed
#' \item F: holding time
#' \item G: gate size
#'}
#'
#'At each setting of the controllable factors, four
#'observations were obtained from a \eqn{2^{(3-1)}}
#'fractional factorial with three noise factors:
#'
#'\itemize{
#'
#' \item M: percentage regrind
#' \item N: moisture content
#' \item O: ambient temperature
#'}
#'
#'
#'@details
#'
#'The data set considered is well known in the literature of
#'industrial experiments and has been analyzed by several
#'authors such as Engel (1992), Engel and Huele (1996)
#'and Lee and Nelder (1998). The experiment was performed to study the influence of seven controllable factors and
#'three noise factors on the mean value and the variation in the percentage of shrinkage of
#'products made by injection molding.Noise factors are fixed during the
#'experiment but are expected to vary randomly outside the
#'experimental context.
#'
#'The aim of the experiment was to determine the process
#'parameter settings so that the shrinkage percentage was
#'close to the target value and robust against environmental
#'variations.
#'
#' @references Engel, J. (1992). Modeling variation in industrial experiments. \emph{Applied Statistics}, 41, 579-593.
#'
#' @references Engel, J. and Huele, A. F. (1996). A generalized linear modeling approach to robust Design. \emph{Technometrics}, 38, 365-373.
#'
#' @references Lee, Y. and Nelder, J.A. (1998). Generalized linear models for analysis of quality improvement experiments. \emph{The Canadian Journal of Statistics}, 26, 95-105.
#'
#'
#' @examples
#' data(injection_molding)
#' head(injection_molding)
"injection_molding"
