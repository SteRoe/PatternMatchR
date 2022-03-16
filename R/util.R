#' @import data.table
#' @import foreach

utils::globalVariables(c("globalVariables"))

#' very first function
#' @description very first function during package load
#' @importFrom magrittr "%>%"
#' @param libname library name
#' @param pkgname package name
#' @return nothing
#' @keywords internal
###' ##@noRd noRd
#' .onAttach()
.onAttach <- function(libname, pkgname) {
  globalVariables <- list()
  base::packageStartupMessage("start loading package")
  base::loadNamespace("PatternMatchR")
  base::packageStartupMessage("end loading package")
}

#' starts the app
#' @description function to start the app
#' @export
PatternMatchRApp <- function() {
  shiny::shinyApp(ui, server)
}

#' loadObjects
#' loads globally needed objects (methylation matrix with beta values, annotation)
#' @param globalVariables global variables
# examples loadObjects()
loadObjects <- function(globalVariables){
  print(paste0(Sys.time(), " loading annotation."))
  annotation <- meffil::meffil.get.features("450k")
  globalVariables$annotation = annotation

  print(paste0(Sys.time(), " finished loading objects."))
  print(paste0(Sys.time(), " detecting cores."))
  globalVariables$numCores = parallel::detectCores();
#  assign("numCores", numCores, envir=globalenv())
  print(paste0(Sys.time(), " finished detecting cores."))
  assign("globalVariables",globalVariables,envir=globalenv())
  return(globalVariables)
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
#' @return data.frame with regression results from file
# examples loadResultFile(folder, fileName)
loadResultFile <- function(folder, fileName) {
  fileName <- paste0(folder ,fileName, ".csv")
  print(paste0(Sys.time()," before fread()."))
  all.results <- data.table::fread(fileName, stringsAsFactors=FALSE, header=TRUE, sep="\t", showProgress=TRUE, data.table=TRUE) # returns data.table
  all.results<-setcolorder(all.results, c("probeID","BETA","SE", "P_VAL", "FDR","DeltaMeth","N"))
  all.results<-all.results[,1:7]
  print(paste0(Sys.time()," before merging annotation."))
  all.results = data.table::merge.data.table(x = all.results, y = globalVariables$annotation, by.x = "probeID", by.y = "name", all.x = TRUE, all.y = FALSE)
  all.results<-all.results[all.results$chromosome!="chrY",] #exclude chrY
  all.results<-all.results[all.results$chromosome!="chrX",] #exclude chrX
  all.results$mLog10FDR<-log10(all.results$FDR)*-1
  all.results$mLog10P_VAL = log10(all.results$P_VAL) * -1
  all.results<-all.results[order(all.results$FDR),]
  #duplicated(all.results$probeID)
  rownames(all.results)<-all.results$probeID
  print(paste0(Sys.time()," after excluding chr x and y ."))

  return(all.results)
}

#' getYoungestFile
#' gets youngest (file date and time) file from a certain folder
#' in use for checking, whether cached .rds-files have to be updated
#' @param folder folder to check in
#' @return file filename of youngest file in folder
# examples getYoungestFile(folder)
getYoungestFile <- function(folder) {
  filesInfo <- file.info(list.files(folder, full.names = T))
  file = rownames(filesInfo)[which.max(filesInfo$mtime)]
return(file)
}

#' extractMantissaExponent
#' @param x value to extract mantissa and exponent
#' @return list with mantissa and exponent from x
# examples extractMantissaExponent(x)
# examples extractMantissaExponent(c(0, 1.234e12, 12345678901234, 123e123))
extractMantissaExponent <- function(x){ #https://r.789695.n4.nabble.com/Built-in-function-for-extracting-mantissa-and-exponent-of-a-numeric-td4670116.html
  e <- ifelse(x == 0, 0, floor(log10(x)))
  m <- x/10^e
  list(mantissa = m, exponent = e)
}

#' delete.na
#' deletes rows with all NA from data.frame DF
#' @param df data.frame
#' @param n minimum number of NA
#' @return data.frame
# examples delete.na(df, 0)
delete.na <- function(df, n=0) {
  # source: https://stackoverflow.com/questions/4862178/remove-rows-with-all-or-some-nas-missing-values-in-data-frame
  df[rowSums(is.na(df)) <= n,]
}

#' loadtraitDFs
#' loads corresponding data.frames for p_val (DFP_Val), delta methylation (DFDM) and n (DFN)
#' @param traitDFs data.frames
#' @return list with three data.frames
# examples loadtraitDFs(traitDFs)
loadtraitDFs <- function(traitDFs) {
  i <- NULL
  foreach(i=1:length(traitDFs)) %do% {
    rownames(traitDFs[[i]]$P_Val)
    colnames(traitDFs[[i]]$P_Val)
    if (i==1) {
      resultDFP_Val = traitDFs[[i]]$P_Val #[[1]]
      resultDFDM = traitDFs[[i]]$DM #[[2]]
      resultDFN = traitDFs[[i]]$N #[[3]]
      resultDFP_Val$Row.names = rownames(traitDFs[[i]]$P_Val)
      resultDFDM$Row.names = rownames(traitDFs[[i]]$DM)
      resultDFN$Row.names = rownames(traitDFs[[i]]$N)
    }
    else {
      print(paste0(Sys.time(), " merge trait ", i , "."))
      traitDFs[[i]]$P_Val$Row.names = rownames(traitDFs[[i]]$P_Val)
      traitDFs[[i]]$DM$Row.names = rownames(traitDFs[[i]]$DM)
      traitDFs[[i]]$N$Row.names = rownames(traitDFs[[i]]$N)
      #Cannot merge with rownames due to:
      #A non-empty vector of column names is required for `by.x` and `by.y`.
      #resultDFP_Val = base::merge(x = resultDFP_Val, y = traitDFs[[i]]$P_Val, by.x = 0, by.y = 0, all.x = TRUE, all.y = TRUE)
      resultDFP_Val = base::merge(x = resultDFP_Val, y = traitDFs[[i]]$P_Val, by.x = "Row.names", by.y = "Row.names", all.x = TRUE, all.y = TRUE)
      resultDFDM = base::merge(x = resultDFDM, y = traitDFs[[i]]$DM, by.x = "Row.names", by.y = "Row.names", all.x = TRUE, all.y = TRUE)
      resultDFN = base::merge(x = resultDFN, y = traitDFs[[i]]$N, by.x = "Row.names", by.y = "Row.names", all.x = TRUE, all.y = TRUE)
    }
  }
  #omit Row.names
  if("Row.names" %in% colnames(resultDFP_Val)) {
    rownames(resultDFP_Val) = resultDFP_Val$Row.names
    resultDFP_Val$Row.names <- NULL
  }
  if("Row.names" %in% colnames(resultDFDM)) {
    rownames(resultDFDM) = resultDFDM$Row.names
    resultDFDM$Row.names <- NULL
  }
  if("Row.names" %in% colnames(resultDFN)) {
    rownames(resultDFN) = resultDFN$Row.names
    resultDFN$Row.names <- NULL
  }
  result = list()
  result$resultDFP_Val = resultDFP_Val
  result$resultDFDM = resultDFDM
  result$resultDFN = resultDFN
  result$rownames = rownames(resultDFP_Val)
  result$colnames = colnames(resultDFP_Val)
  return(result)
}

