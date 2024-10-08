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
  base::print(base::paste0(sysTimePID(), " set package wd"))
  originalWd <<- getwd()
  packageWd <- paste0(libname, "/", pkgname)
  setwd(packageWd)
  base::print(base::paste0(sysTimePID(), " setwd():", packageWd))
}

.onUnload <- function(libpath) {
  o <- originalWd
  setwd(o)
}

#' starts the app
#' @description function to start the app
#' @export
PatternMatchRApp <- function() {
  shiny::shinyApp(generate_ui(), server)
}

loadDirLists <- function(session, input, output) {
  dfdD1 <-
    data.table::as.data.table(base::unlist(session$userData$config$dataDir1))
  dfdD2 <-
    data.table::as.data.table(base::unlist(session$userData$config$dataDir2))
  dfdD3 <-
    data.table::as.data.table(base::unlist(session$userData$config$dataDir3))

  output$trait1DirList <- DT::renderDT({DT::datatable(dfdD1, options = list(pageLength = 100))}, server = FALSE) #output$trait1DirList <- DT::renderDataTable(dfdD1)
  output$trait2DirList <- DT::renderDataTable(DT::datatable(dfdD2, options = list(pageLength = 100)), server = FALSE)
  output$trait3DirList <- DT::renderDataTable(DT::datatable(dfdD3, options = list(pageLength = 100)), server = FALSE)
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
  id <- shiny::showNotification("Loading annotation...", duration = NULL, closeButton = FALSE)
  base::on.exit(shiny::removeNotification(id), add = TRUE)
  base::print(base::paste0(sysTimePID(), " loading annotation."))
  annotation <- meffil::meffil.get.features("450k")
  rownames(annotation) <- annotation$name
  session$userData$annotation <- annotation
  base::print(base::paste0(sysTimePID(), " finished loading annotation with dim: nrow = ", nrow(annotation), ", ncol = ", ncol(annotation), "."))
  base::print(base::paste0(sysTimePID(), " detecting cores."))
  numCores <- parallelly::availableCores() # parallel::detectCores()
  nWorkers <- parallelly::availableCores(constraints = "connections")
  numCores <- base::min(numCores, nWorkers)
  if (!is.null(numCores)) {
    if (numCores >= 64) {
      #numCores <- as.integer(numCores / 2)
      numCores <- as.integer(numCores * 0.75) #take 3/4 of available cores
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
  base::print(base::paste0(sysTimePID(), " finished detecting cores."))

  #check config file for "/inst"
  betaFileName <- session$userData$config$betaFileName
  #check for existence of betaFileName
  if (!file.exists(betaFileName)) {
    base::message(base::paste0(sysTimePID(), " beta file not found: ", betaFileName, ". Is your config.yml correct?"))
    #add /inst to betafilename and try again
    betaFileName <- addInstToPath(betaFileName)
    if (file.exists(betaFileName)) {
      base::message(base::paste0(sysTimePID(), " beta file now found: ", betaFileName, "."))
    }
  }
  dataDirList <- as.list(session$userData$config$dataDir1)
  dataDirList <- lapply(dataDirList, FUN = addInstToPath)
  session$userData$config$dataDir1 <- unlist(dataDirList)

  dataDirList <- as.list(session$userData$config$dataDir2)
  dataDirList <- lapply(dataDirList, FUN = addInstToPath)
  session$userData$config$dataDir2 <- unlist(dataDirList)

  dataDirList <- as.list(session$userData$config$dataDir3)
  dataDirList <- lapply(dataDirList, FUN = addInstToPath)
  session$userData$config$dataDir3 <- unlist(dataDirList)

  # read methylation data
  session$userData$BetaDF <- loadDNAm(betaFileName = betaFileName, config = session$userData$config)
  session$userData$Beta_tDF <- as.data.frame(t(session$userData$BetaDF))
  base::print(base::paste0(sysTimePID(), " finished read methylation data with nrow = ", nrow(session$userData$BetaDF), " and ncol = ", ncol(session$userData$BetaDF), "."))

  #read base data
  baseFileName <- session$userData$config$baseFileName
  #check for existence of baseFileName
  if (!file.exists(baseFileName)) {
    base::message(base::paste0(sysTimePID(), " beta file not found: ", baseFileName, ". Is your config.yml correct?"))
    #add /inst to baseFileName and try again
    baseFileName <- addInstToPath(baseFileName)
    if (file.exists(baseFileName)) {
      base::message(base::paste0(sysTimePID(), " beta file now found: ", baseFileName, "."))
    }
  }

  session$userData$baseData <- loadBaseData(session = session, baseFileName = baseFileName)
  session$userData$baseData <- as.data.frame(session$userData$baseData)
  base::print(base::paste0(sysTimePID(), " finished read base data with nrow = ", nrow(session$userData$baseData), " and ncol = ", ncol(session$userData$baseData), "."))
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

loadBaseData <- function(session, baseFileName) {
  id <- shiny::showNotification("Loading base data...", duration = NULL, closeButton = FALSE)
  base::on.exit(shiny::removeNotification(id), add = TRUE)
  result <- NULL
  result <- data.table::fread(file=baseFileName, sep="\t", dec=".", data.table = FALSE)
  result <- base::subset(result, select=c(session$userData$config$mergeAttribut, session$userData$config$sexAttribut))
  rownames(result) <- result[, session$userData$config$mergeAttribut]
  result <- as.data.frame(result)
  return(result)
}

#' loadResultFile
#' @param folder folder containing file with fileName
#' @param fileName file containing probeID, p-values, delta methylation-values and n from regression
#' @param session session object
#' @return data.frame with regression results from file
# examples loadResultFile(folder, fileName)
loadResultFile <- function(session, folder, fileName) {
  id <- shiny::showNotification("Loading result file...", duration = NULL, closeButton = FALSE)
  base::on.exit(shiny::removeNotification(id), add = TRUE)
  options(warn = 2)
  fileName <- base::paste0(folder, fileName, ".csv")
  base::print(base::paste0(sysTimePID(), " before fread()."))
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
      base::c("probeID", "BETA", "SE", "P_VAL", "FDR", "DeltaMeth", "logFC", "N")
    )
  all.results <- all.results[, 1:8]
  base::print(base::paste0(sysTimePID(), " before merging annotation."))
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
  base::print(base::paste0(sysTimePID(), " after excluding chr x and y ."))
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
    # convert negative values
    x <- base::abs(x)
    e <- ifelse(x == 0, 0, floor(log10(x)))
    m <- x / 10^e
    list(mantissa = m, exponent = e)
  }

addInstToPath <- function(fileName) {
  if (!file.exists(fileName)) {
    base::message(base::paste0(sysTimePID(), " file not found: ", fileName, ". Is your config.yml correct?"))
    fileName <- stringr::str_sub(fileName, start = 2, end = stringr::str_length(fileName))
    fileName <- paste0("./inst", fileName)
    if (file.exists(fileName)) {
      base::message(base::paste0(sysTimePID(), " file now found: ", fileName, "."))
    }
    else {
      base::message(base::paste0(sysTimePID(), " file furthermore not found: ", fileName, ". Is your config.yml correct?"))
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
  id <- shiny::showNotification("Loading raw methylation data...", duration = NULL, closeButton = FALSE)
  base::on.exit(shiny::removeNotification(id), add = TRUE)
  waiter::waiter_show()
  waiter::transparent(alpha = 0)
  base::on.exit(waiter::waiter_hide(), add = TRUE)
  base::print(base::paste0(sysTimePID(), " loading configuration."))
  #config <- session$userData$config #config::get(file = "config.yml")
  base::print(base::paste0(sysTimePID(), " load beta."))
  if (TRUE) {
  #if (config$debugMode == FALSE) {
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
        nrows = 10000,
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
      base::print(base::paste0(sysTimePID(), " ", what, "not found in ", in_files))
    }
    return(result)
  }

#' removeAdjFromColname
#' @param colnames list with colnames probably containing 'adj' at the end of each name
#' @return colnames without 'adj' at the end
removeAdjFromColname <- function(colnames) {
  #' adj' was added to the end of the name of an analysis for a regression model with additional adjustment variables
  base::print(base::paste0(sysTimePID(), " start removeAdjFromColname()"))
  result <- sub("adj$", "", colnames)
  base::print(base::paste0(sysTimePID(), " end removeAdjFromColname()"))
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

#' addLinkToEWASDataHub
#' adds links to EWASDataHub to a data.frame as separate column
#' @param df data.frame to which links should be added
#' @param baseURL string describing link to be included
#' @param probeAttribut string describing the name of the probe variable
#' @return data.frame
#' @keywords internal
#' @noRd
# examples addLinkToEWASDataHub(data.frame, baseURL)
addLinkToEWASDataHub <- function(df, baseURL, probeAttribut) {
  #provide link to EWAS data hub
  df$EWASDataHub <- base::paste0('<a target=_blank rel="noopener noreferrer" href=', baseURL, df$probeID, '>', df[, probeAttribut], '</a>')
  return(df)
}

addLinkToEWASDataHubShort <- function(df, baseURL, probeAttribut) {
  #provide link to EWAS data hub
  df$EWASDataHubShort <- base::paste0(baseURL, df[, probeAttribut])
  return(df)
}

#' addLinkToMRCEWASCatalog
#' adds links to MRC EWAS catalog to a data.frame as separate column
#' @param df data.frame to which links should be added
#' @param baseURL string describing link to be included
#' @param probeAttribut string describing the name of the probe variable
#' @return data.frame
#' @keywords internal
#' @noRd
# examples addLinkToMRCEWASCatalog(data.frame)
addLinkToMRCEWASCatalog <- function(df, baseURL, probeAttribut) {
  #provide link to MRC EWAS catalog
  df$MRCEWASCatalog <- base::paste0('<a target=_blank rel="noopener noreferrer" href=', baseURL, df$probeID, '>', df[, probeAttribut], '</a>')
  return(df)
}

addLinkToMRCEWASCatalogShort <- function(df, baseURL, probeAttribut) {
  #provide link to MRC EWAS catalog
  df$MRCEWASCatalogShort <- base::paste0(baseURL, df[, probeAttribut])
  return(df)
}

sysTimePID <- function() {
  result <- paste0(as.character(Sys.time()), "; PID: ", as.character(Sys.getpid()))
  return(result)
}

#' plotlyPcPForDMP <- function(DMPNearRange, probe, resultDataSingleTrait, annotation) {
#'   #plotlyPcPForDMP <- function(DMPNearRange, probe, P_VAL, DeltaMeth, annotation) {
#'   #plotlyPcPForDMP <- function(globalVariables, sessionVariables, DMPNearRange) {
#'   tryCatch({
#' browser()
#'     traitName <- colnames(DMPNearRange)[3]
#'     #DMP <- sessionVariables$probe$probe
#'     #df <- sessionVariables$resultDataSingleTrait
#'     df <- resultDataSingleTrait
#'     #df <- resultDataSingleScenarioWithAnnotation(globalVariables$annotation, df)
#'     df <- resultDataSingleScenarioWithAnnotation(annotation, df)
#'     gene.symbol <- df[which(df$probeID == probe),]$gene.symbol
#'     DMPNearRange <- stats::na.omit(DMPNearRange)
#'     DMPNearRangeShort <- DMPNearRange[,4:ncol(DMPNearRange)]
#'     dimensionsList=list()
#'     P_VAL <- resultDataSingleTrait$P_VAL[resultDataSingleTrait$probeID == probe]
#'     DeltaMeth <- resultDataSingleTrait$DeltaMeth[resultDataSingleTrait$probeID == probe]
#'     for (i in 1:ncol(DMPNearRangeShort)) {
#'       lblCpG = colnames(DMPNearRangeShort)[i]
#'       lblP = signif(df[which(df$probeID == colnames(DMPNearRangeShort)[i]),]$P_VAL,3)
#'       if (shiny::isTruthy(lblP)) {
#'         lblP = paste0(",\n p: ", lblP)
#'         lblDM = paste0(",\n d: ", signif(df[which(df$probeID == colnames(DMPNearRangeShort)[i]),]$DeltaMeth,3))
#'       }
#'       else {
#'         lblP = ",\n p: n.s."
#'         lblDM = ""
#'       }
#'       lblSym = paste0(",\n sbl:", annotation[which(annotation$name == colnames(DMPNearRangeShort)[i]),]$gene.symbol)
#'       lblPos = paste0(",\n pos:", annotation[which(annotation$name == colnames(DMPNearRangeShort)[i]),]$position)
#'       label = paste0(lblCpG, lblP, lblDM, lblSym, lblPos)
#'       dimension = list(label = label, values = DMPNearRangeShort[,i],
#'                        range = c(0, 1))
#'       dimensionsList = append(dimensionsList,list(dimension))
#'     }
#'     plot <- plotly::plot_ly(data = DMPNearRange)
#'     plot <- plot %>% plotly::add_trace(type = 'parcoords',
#'                                        line = list(shape = 'spline',
#'                                                    color =  DMPNearRange[,3],
#'                                                    colorscale = 'Jet',
#'                                                    showscale = TRUE,
#'                                                    reversescale = TRUE,
#'                                                    cmin = min(DMPNearRange[,3],na.rm=TRUE),
#'                                                    cmax = max(DMPNearRange[,3],na.rm=TRUE)),
#'                                        dimensions = dimensionsList
#'     )
#'     plot <- plot %>% plotly::layout(
#'       #title = paste0(traitName, " vs. ", DMP, " gene.symbol: ", gene.symbol ," P_VAL: ", P_VAL, " DeltaMeth: ", DeltaMeth),
#'       title = paste0(traitName, " vs. ", probe, " gene.symbol: ", gene.symbol ," P_VAL: ", P_VAL, " DeltaMeth: ", DeltaMeth),
#'       xaxis = list(
#'         title = "Location",
#'         showgrid = F,
#'         tickangle = 45),
#'       yaxis = list(
#'         title = "Methylation [%]"),
#'       coloraxis = list(
#'         title = traitName)
#'     )
#'     plot <- plot %>%
#'       plotly::add_annotations(
#'         text = "Methylation [%]",
#'         x = 0,
#'         y = 0.5,
#'         yref = "paper",
#'         xref = "paper",
#'         xanchor = "center",
#'         yanchor = "center",
#'         xshift = -50,
#'         showarrow = FALSE,
#'         textangle = 270,
#'         font = list(size = 15)
#'       )
#'     plot <- plot %>%
#'       plotly::add_annotations(
#'         text = "globalArrayPosition",
#'         x = 0.5,
#'         y = 0,
#'         yref = "paper",
#'         xref = "paper",
#'         xanchor = "center",
#'         yanchor = "center",
#'         yshift = -35,
#'         showarrow = FALSE,
#'         textangle = 0,
#'         font = list(size = 15)
#'       )
#'     return (plot)
#'   }, error=function(err){
#'     print(paste0("unable to print pc plot: ", err$message, ". Maybe there was only a subset of the data loaded?"))
#'     return(empty_plot(err$message))
#'   });
#' }
#'
#' #' empty_plot
#' #' @param title title for empty plot
#' #' @return plot empty plot
#' #' @keywords internal
#' #' @noRd
#' empty_plot <- function(title = NULL){
#'   plot <- plotly::plotly_empty(type = "scatter", mode = "markers") %>%
#'     plotly::config(
#'       displayModeBar = FALSE
#'     ) %>%
#'     plotly::layout(
#'       title = list(
#'         text = title,
#'         yref = "paper",
#'         y = 0.5
#'       )
#'     )
#'   return(plot)
#' }
#'
#' resultDataSingleScenarioWithAnnotation <- function(annotation, df){
#'   #  a = reducedAnnotation(globalVariables)
#'   a = EpiVisR::reducedAnnotation(annotation)
#'   a$gene.symbolShort = stringr::str_sub(a$gene.symbol, 1, 20) #NULL
#'   a$gene.accession = NULL
#'   a$gene.region = NULL
#'   a$cpg.island.name = NULL
#'   a$relation.to.island = NULL
#'   a$snp.exclude = NULL
#'   #df = dplyr::left_join(df, a, by = c("probeID" = "name"))
#'   df <- base::merge(df, a, by.x = "probeID", by.y = "name")
#'   return (df)
#' }
