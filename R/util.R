#' structure of sessionVariables:
#'
#' general settings
#' session$userData$sessionVariables$MaxProbes - maximum number of probes to use in clustering
#' session$userData$sessionVariables$P_ValMinBorder - minimum p-value to consider in clustering and to visualize
#' session$userData$sessionVariables$P_ValMaxBorder - maximum p-value to consider in clustering and to visualize

#' general data
#' session$userData$globalVariables$Beta_tDF - general beta values from Illumina 450k array

#utils::globalVariables(c("globalVariables"))

#' very first function
#' @description very first function during package load
#' @importFrom magrittr "%>%"
#' @importFrom graphics "plot.new"
#' @importFrom stats "as.dist"
#' @importFrom stats "as.formula"
#' @importFrom dendextend "cutree"
#' @importFrom foreach "%do%"
#' @importFrom foreach "%dopar%"
#' @param libname library name
#' @param pkgname package name
#' @return nothing
#' @keywords internal
###' ##@noRd noRd

originalWd <- NULL

#' .onAttach()
.onAttach <- function(libname, pkgname) {
  #  globalVariables <- base::list()
  base::packageStartupMessage(base::paste0("start loading package "), pkgname)
  base::loadNamespace("PatternMatchR")
  base::packageStartupMessage(base::paste0("finished start loading package "), pkgname)
}

.onLoad <- function(libname, pkgname) {
  base::print(base::paste0(Sys.time(), "set package wd"))
  originalWd <<- getwd()
  packageWd <- paste0(libname, "/", pkgname)
  setwd(packageWd)
  base::print(base::paste0(Sys.time(), " setwd():", packageWd))
}

.onUnload <- function(libpath) {
  o <<- originalWd
  setwd(o)
}

#' starts the app
#' @description function to start the app
#' @export
PatternMatchRApp <- function() {
  #shiny::shinyApp(ui, server)
  shiny::shinyApp(generate_ui(), server)
}

loadDirLists <- function(session, input, output) {
  dfdD1 <<-
    data.table::as.data.table(base::unlist(session$userData$config$dataDir1))
  dfdD2 <<-
    data.table::as.data.table(base::unlist(session$userData$config$dataDir2))
  dfdD3 <<-
    data.table::as.data.table(base::unlist(session$userData$config$dataDir3))

  output$trait1DirList <- DT::renderDataTable(dfdD1)
  output$trait2DirList <- DT::renderDataTable(dfdD2)
  output$trait3DirList <- DT::renderDataTable(dfdD3)
  result <- list()
  result$dfdD1 <- dfdD1
  result$dfdD2 <- dfdD2
  result$dfdD3 <- dfdD3
  return(result)
}

#' loadObjects
#' loads globally needed objects (methylation matrix with beta values, annotation)
#' @param session session object
# examples loadObjects(session)
#loadObjects <- function(globalVariables) {
loadObjects <- function(session) {
  base::print(base::paste0(Sys.time(), " loading annotation."))
  annotation <- meffil::meffil.get.features("450k")
  session$userData$annotation <- annotation
  base::print(base::paste0(Sys.time(), " finished loading annotation with dim: nrow = ", nrow(annotation), ", ncol = ", ncol(annotation),"."))
  base::print(base::paste0(Sys.time(), " detecting cores."))
  numCores <- parallelly::availableCores() # parallel::detectCores()
  if (!is.null(numCores)) {
    if (numCores >= 64) {
      numCores <- as.integer(numCores / 2)
    }
    else if (numCores >= 8) {
      numCores <- numCores - 4
    }
    else {
      numCores <- numCores - 1
    }
  } else {
    numCores <- 1
  }
  if (numCores < 1) {
    numCores <- 1
  }
  session$userData$numCores <- numCores
  base::print(base::paste0(Sys.time(), " finished detecting cores."))

  #check config file for "/inst"
  betaFileName <- session$userData$config$betaFileName
  #check for existence of betaFileName
  if (!file.exists(betaFileName)) {
    base::message(base::paste0(Sys.time(), " beta file not found: ", betaFileName, ". Is your config.yml correct?"))
    #add /inst to betafilename and try again
    betaFileName <- addInstToPath(betaFileName)
    if (file.exists(betaFileName)) {
      base::message(base::paste0(Sys.time(), " beta file now found: ", betaFileName, "."))
    }
  }
  dataDirList <- as.list(session$userData$config$dataDir1)
  DirFound <- lapply(dataDirList, FUN = file.exists)
  dataDirList <- lapply(dataDirList, FUN = addInstToPath)
  session$userData$config$dataDir1 <- unlist(dataDirList)

  dataDirList <- as.list(session$userData$config$dataDir2)
  DirFound <- lapply(dataDirList, FUN = file.exists)
  dataDirList <- lapply(dataDirList, FUN = addInstToPath)
  session$userData$config$dataDir2 <- unlist(dataDirList)

  dataDirList <- as.list(session$userData$config$dataDir3)
  DirFound <- lapply(dataDirList, FUN = file.exists)
  dataDirList <- lapply(dataDirList, FUN = addInstToPath)
  session$userData$config$dataDir3 <- unlist(dataDirList)

  # read methylation data
  session$userData$BetaDF <- loadDNAm(betaFileName = betaFileName, config = session$userData$config)
  session$userData$Beta_tDF <- as.data.frame(t(session$userData$BetaDF))
  base::print(base::paste0(Sys.time(), " finished read methylation data with nrow = ", nrow(session$userData$BetaDF), " and ncol = ", ncol(session$userData$BetaDF), "."))
  return(session)
}

# setBaseSciPen<-function(){
#   options(scipen=-10, digits=5)
# }
#
# resetSciPen<-function(){
#   options(scipen=0, digits=7)
# }
#

#' loadResultFile
#' @param folder folder containing file with fileName
#' @param fileName file containing probeID, p-values, delta methylation-values and n from regression
#' @param session session object
#' @return data.frame with regression results from file
# examples loadResultFile(folder, fileName)
#loadResultFile <- function(folder, fileName, globalVariables) {
loadResultFile <- function(session, folder, fileName) {
#tbc() insert tryCatch here or options(warn=2)
options(warn = 2)
  fileName <- base::paste0(folder, fileName, ".csv")
  base::print(base::paste0(Sys.time(), " before fread()."))
  all.results <-
    data.table::fread(
      fileName,
      stringsAsFactors = FALSE,
      header = TRUE,
      sep = "\t",
      dec = ".",
      showProgress = TRUE,
      data.table = TRUE
    ) # returns data.table
  all.results <-
    data.table::setcolorder(
      all.results,
      base::c("probeID", "BETA", "SE", "P_VAL", "FDR", "DeltaMeth", "N")
    )
  all.results <- all.results[, 1:7]
  base::print(base::paste0(Sys.time(), " before merging annotation."))
  all.results <-
    data.table::merge.data.table(
      x = all.results,
      y = session$userData$annotation, #y = globalVariables$annotation,
      by.x = "probeID",
      by.y = "name",
      all.x = TRUE,
      all.y = FALSE
    )
  all.results <-
    all.results[all.results$chromosome != "chrY", ] # exclude chrY
  all.results <-
    all.results[all.results$chromosome != "chrX", ] # exclude chrX
  all.results$mLog10FDR <- base::log10(all.results$FDR) * -1
  all.results$mLog10P_VAL <- base::log10(all.results$P_VAL) * -1
  all.results <- all.results[order(all.results$FDR), ]
  # duplicated(all.results$probeID)
  rownames(all.results) <- all.results$probeID
  base::print(base::paste0(Sys.time(), " after excluding chr x and y ."))

  return(all.results)
}

#' getYoungestFile
#' gets youngest (file date and time) file from a certain folder
#' in use for checking, whether cached .rds-files have to be updated
#' @param folder folder to check in
#' @return file filename of youngest file in folder
# examples getYoungestFile(folder)
getYoungestFile <- function(folder) {
  filesInfo <- base::file.info(list.files(folder, full.names = TRUE))
  file <-
    base::rownames(filesInfo)[base::which.max(filesInfo$mtime)]
  return(file)
}

#' extractMantissaExponent
#' @param x value to extract mantissa and exponent
#' @return list with mantissa and exponent from x
# examples extractMantissaExponent(x)
# examples extractMantissaExponent(c(0, 1.234e12, 12345678901234, 123e123))
extractMantissaExponent <-
  function(x) {
    # https://r.789695.n4.nabble.com/Built-in-function-for-extracting-mantissa-and-exponent-of-a-numeric-td4670116.html
    e <- ifelse(x == 0, 0, floor(log10(x)))
    m <- x / 10^e
    list(mantissa = m, exponent = e)
  }

addInstToPath <- function (fileName) {
  if(!file.exists(fileName)) {
    base::message(base::paste0(Sys.time(), " file not found: ", fileName, ". Is your config.yml correct?"))
    fileName <- stringr::str_sub(fileName, start = 2, end = stringr::str_length(fileName))
    fileName <- paste0("./inst", fileName)
    if(file.exists(fileName)) {
      base::message(base::paste0(Sys.time(), " file now found: ", fileName, "."))
    }
    else {
      base::message(base::paste0(Sys.time(), " file furthermore not found: ", fileName, ". Is your config.yml correct?"))
    }
  }
  return(fileName)
}

#' loadDNAm
#' loads measured beta values from file config$betaFileName which was defined in configuration
#' @param betaFileName location of beta file
#' @param config config structure
#' @return data frame with measured beta values
# examples loadDNAm(betaFileName)
loadDNAm <- function(betaFileName, config) {
  base::print(base::paste0(Sys.time(), " loading configuration."))
#  config <- session$userData$config #config::get(file = "config.yml")
  base::print(base::paste0(Sys.time(), " load beta."))
  #if (TRUE) {
  if (config$debugMode == FALSE) {
    beta <-
      data.table::fread(
        betaFileName,
        stringsAsFactors = FALSE,
        header = TRUE,
        sep = "\t",
        dec = ".",
        data.table = FALSE
      )
  } else {
    beta <-
      data.table::fread(
        betaFileName,
        stringsAsFactors = FALSE,
        header = TRUE,
        sep = "\t",
        dec = ".",
        nrows = 1000,
        data.table = FALSE
      )
  }
  beta <- as.data.frame(beta)
  #    beta<-data.frame(column_to_rownames(beta, var = "PROBEID"))
  #    rownames(beta) <- beta$probeID
  rownames(beta) <- beta[, config$probeAttribut]
  beta[, config$probeAttribut] <- NULL
  return(beta)
}

#' findInFiles
#' finds strings in a bunch of files
#' @param what what (string) to find
#' @param where directory, where to search
#' @param in_files regular expression for filenames, where to search
#' @param recursive recursive search
#' @param ignore.case ignore case sensitivity
#' @return data frame with measured beta values
# examples findInFiles(what, where, in_files, recursive, ignore.case)
findInFiles <-
  function(what,
           where = ".",
           in_files = "\\.[Rr]$",
           recursive = TRUE,
           ignore.case = TRUE) {
    files <-
      list.files(
        path = where,
        pattern = in_files,
        recursive = recursive
      )
    found <- FALSE
    for (fil in files) {
      contents <- base::readLines(paste0(where, "/", fil))
      res <- grepl(what, contents, ignore.case = ignore.case)
      res <- base::which(res)
      if (base::length(res) > 0) {
        found <- TRUE
        result <- contents[res]
      }
    }
    if (!found) {
      result <- FALSE
      base::print(base::paste0(Sys.time(), " ", what, "not found in ", in_files))
    }
    return(result)
  }

#' removeAdjFromColname
#' @param colnames list with colnames probably containing 'adj' at the end of each name
#' @return colnames without 'adj' at the end
removeAdjFromColname <- function(colnames) {
  #' adj' was added to the end of the name of an analysis for a regression model with additional adjustment variables
  base::print(base::paste0(Sys.time(), " start removeAdjFromColname()"))
  result <- sub("adj$", "", colnames)
  base::print(base::paste0(Sys.time(), " end removeAdjFromColname()"))
  return(result)
}

#' is.valid
#' checks, whether an R object (x) is valid (not NULL) and exists
#' @param x object
#' @return TRUE or FALSE
#' examples is.valid(x)
is.valid <- function(x) {
  if (base::exists("x")) {
  is.null(shiny::need(x, message = FALSE))
  }
  else return(FALSE)
}

#' getEmptyPlot
#' @return empty plot
#' examples getEmptyPlot()
getEmptyPlot <- function() {
  plot.new()
}
