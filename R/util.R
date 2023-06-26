#' structure of sessionVariables:
#'
#' general settings
#' session$userData$sessionVariables$MaxProbes - maximum number of probes to use in clustering
#' session$userData$sessionVariables$P_ValMinBorder - minimum p-value to consider in clustering and to visualize
#' session$userData$sessionVariables$P_ValMaxBorder - maximum p-value to consider in clustering and to visualize

#' #' very first function
#' #' @description very first function during package load
#' #' @importFrom magrittr "%>%"
#' #' @importFrom graphics "plot.new"
#' #' @importFrom stats "as.dist"
#' #' @importFrom stats "as.formula"
#' #' @importFrom dendextend "cutree"
#' #' @importFrom foreach "%do%"
#' #' @importFrom foreach "%dopar%"
#' #' @param libname library name
#' #' @param pkgname package name
#' #' @return nothing
#' #' @keywords internal
#' ###' ##@noRd noRd
#'
#' originalWd <- NULL
#'
#' #' .onAttach()
#' .onAttach <- function(libname, pkgname) {
#'   #  globalVariables <- base::list()
#'   base::packageStartupMessage(base::paste0("start loading package "), pkgname)
#'   base::loadNamespace("PatternMatchR")
#'   base::packageStartupMessage(base::paste0("finished start loading package "), pkgname)
#' }
#'
#' .onLoad <- function(libname, pkgname) {
#'   base::print(base::paste0(Sys.time(), "set package wd"))
#'   originalWd <<- getwd()
#'   packageWd <- paste0(libname, "/", pkgname)
#'   setwd(packageWd)
#'   base::print(base::paste0(Sys.time(), " setwd():", packageWd))
#' }
#'
#' .onUnload <- function(libpath) {
#'   o <<- originalWd
#'   setwd(o)
#' }

#' starts the app
#' @description function to start the app
#' @export
PatternMatchRApp <- function() {
  shiny::shinyApp(ui, server)
}
