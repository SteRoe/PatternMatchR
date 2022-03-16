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
  annotationRDSFileName = globalVariables$config$annotationRDSFileName
  if (utils::file_test("-f", annotationRDSFileName) == TRUE) {
    print(paste0(Sys.time(), " loading annotation RDS."))
    globalVariables$annotation = readRDS(file = annotationRDSFileName)
  }
  else {
    #load list of multimodal CpG
    multiModCpGFileName = globalVariables$config$multiModCpGFileName
    MultiModCpG<-data.table::fread(file = multiModCpGFileName, sep = "\t", dec = ".", showProgress = TRUE)
#    assign("MultiModCpG",MultiModCpG,envir=globalenv())
    print(paste0(Sys.time(), " loading beta."))
    betaFileName = globalVariables$config$betaFileName
    beta <- data.table::fread(betaFileName, stringsAsFactors = FALSE, header = TRUE, sep = "\t", showProgress = TRUE)

    #  beta<-data.frame(column_to_rownames(beta, var = "PROBEID"))
    print(paste0(Sys.time(), " assigning rownames for beta."))
    rownames(beta) <- beta$PROBEID
    #beta<-as.data.table(beta)

    print(paste0(Sys.time(), " removing multimodal CpG."))
    beta<-removeMultiModelCpGFromBeta(beta,MultiModCpG)
    rownames(beta) <- beta$PROBEID
#    assign("beta",beta,envir=globalenv())
    print(paste0(Sys.time(), " loading annotation."))
    annotation <- meffil::meffil.get.features("450k")
    annotation$relation.to.island = as.factor(annotation$relation.to.island)
    print(paste0(Sys.time(), " remove unmeasured or multimodal probeIDs from annotation."))
    annotation = annotation[which(annotation$name %in% rownames(beta)),]
    annotation = as.data.table(annotation)
    print(paste0(Sys.time(), " writing annotation."))
    saveRDS(annotation,file = annotationRDSFileName)
    #  assign("annotation",annotation,envir=globalenv())
    globalVariables$annotation = annotation
  }
  dim(globalVariables$annotation)
  print(paste0(Sys.time(), " finished loading objects."))
  print(paste0(Sys.time(), " detecting cores."))
  globalVariables$numCores = parallel::detectCores();
#  assign("numCores", numCores, envir=globalenv())
  print(paste0(Sys.time(), " finished detecting cores."))
  assign("globalVariables",globalVariables,envir=globalenv())
  return(globalVariables)
}

# loadObjects<-cmpfun(loadObjects)
#
# setBaseSciPen<-function(){
#   options(scipen=-10, digits=5)
# }
#
# resetSciPen<-function(){
#   options(scipen=0, digits=7)
# }
#
#' winsorize
#' performs winsorizing
#' @param traitDF data.frame to be used
#' @param trim value to be used for winsorizing
#' @param startVar variable to start with
#' @param endVar variable to end at
#' @return data.frame with winsorized variables
# examples winsorize(df, 0.05, 10, 20)
winsorize <- function(traitDF, trim, startVar, endVar) {
  #winsorize ScenarioDF
  i <- NULL
  foreach(i = startVar:endVar, .combine = cbind, .verbose=FALSE) %do% {
    tryCatch({
      traitDF[,i] <- psych::winsor(traitDF[,i],trim)
    }, error=function(err){
      message(paste0("Error: ", err$message))
#      i = i+1
    })
  }
  return (traitDF)
}

#winsorize<-compiler::cmpfun(winsorize)

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
#  all.results = base::merge(x = all.results, y = annotation, by.x = "probeID", by.y = "name", all.x = TRUE, all.y = FALSE)
  all.results = data.table::merge.data.table(x = all.results, y = globalVariables$annotation, by.x = "probeID", by.y = "name", all.x = TRUE, all.y = FALSE)
  all.results<-all.results[all.results$chromosome!="chrY",] #exclude chrY
  all.results<-all.results[all.results$chromosome!="chrX",] #exclude chrX
  # if (is.numeric(minN)) {
  #   all.results<-all.results[all.results$N<minN,] #exclude n<minN
  # }
  all.results$mLog10FDR<-log10(all.results$FDR)*-1
  all.results$mLog10P_VAL = log10(all.results$P_VAL) * -1
  all.results<-all.results[order(all.results$FDR),]
  #duplicated(all.results$probeID)
  rownames(all.results)<-all.results$probeID
  print(paste0(Sys.time()," after excluding chr x and y ."))
#  all.results$DeltaMeth[(all.results$BETA < 0 & all.results$DeltaMeth > 0)] <- all.results$DeltaMeth*-1
  return(all.results)
}

#loadResultFile<-compiler::cmpfun(loadResultFile)

#weiter: check, whether furthermore needed
#' removeMultiModelCpGFromBeta
#' @param df data.frame
#' @param multiModList list with known multimodal CpG to be removed
#' @return data.frame without multimodal probes
# examples removeMultiModelCpGFromBeta(data.frame, list)
removeMultiModelCpGFromBeta<-function(df, multiModList){
  #row.name to column
  df$CpGName<-row.names(df)
  #merge
  P1 <- dplyr::inner_join(df, multiModList, by = c("CpGName" = "CpG"))
  row.names(P1)<-P1$CpGName
  P1$CpGName<-NULL
  #select only CpG with NumModes=1
  #browser()
  P1<-P1[P1$NumModes<2,]
  P1$NumModes<-NULL
  P1$NormalP<-NULL
  #  structure(P1)
  #  colnames(P1)
  df<-P1
  rm(P1)
  return(df)
}

removeMultiModelCpGFromBeta<-compiler::cmpfun(removeMultiModelCpGFromBeta)

#weiter: check, whether furthermore needed
# removeOutliers<-function(probes){
# #  require(matrixStats)
#   if(nrow(probes) < ncol(probes)) warning("expecting probes are rows (long dataset)")
#   rowIQR <- matrixStats::rowIQRs(probes, na.rm = T)
#   row2575 <- matrixStats::rowQuantiles(probes, probs = c(0.25, 0.75), na.rm = T)
#   maskL <- probes < row2575[,1] - 3 * rowIQR
#   maskU <- probes > row2575[,2] + 3 * rowIQR
#   initial_NAs<-rowSums(is.na(probes))
#   probes[maskL] <- NA
#   removed_lower <- rowSums(is.na(probes))-initial_NAs
#   probes[maskU] <- NA
#   removed_upper <- rowSums(is.na(probes))-removed_lower-initial_NAs
#   N_for_probe<-rowSums(!is.na(probes))
#   Log<-data.frame(initial_NAs,removed_lower,removed_upper,N_for_probe)
#   return(list(probes, Log))
# }
# removeOutliers<-compiler::cmpfun(removeOutliers)


#proceed: check, whether furthermore needed?
#' reducedAnnotation
#' @return sorted annotation list without type, target and meth.dye
# examples reducedAnnotation()
reducedAnnotation <- function(){
    a = globalVariables$annotation
    a$type = NULL
    a$target = NULL
    a$meth.dye = NULL
    a$chromosome = as.factor(a$chromosome)
    levels(a$chromosome)[levels(a$chromosome)=="chrX"] <-"chr23" #only for sorting
    levels(a$chromosome)[levels(a$chromosome)=="chrY"] <-"chr24" #only for sorting
    a$chromosomeNum = as.factor(as.numeric(gsub("chr","",a$chromosome)))
    a = a[order(a$chromosomeNum,a$position),]
    a$globalPosition <- seq_len(nrow(a))
return (a)
}

#weiter: check, whether furthermore needed
# resultDataSingleScenarioWithAnnotation <- function(df){
#   a = reducedAnnotation()
#   a$gene.symbolShort = stringr::str_sub(a$gene.symbol, 1, 20) #NULL
#   a$gene.accession = NULL
#   a$gene.region = NULL
#   a$cpg.island.name = NULL
#   a$relation.to.island = NULL
#   a$snp.exclude = NULL
#   return (dplyr::left_join(df, a, by = c("probeID" = "name")))
# }

#weiter: check, whether furthermore needed
# resultDataSingleScenarioWithAnnotationEWAScatalogCount <- function(df){
#   return (dplyr::left_join(df, EWAScatalogCount, by = c("probeID" = "CpG")))
# }

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
      #           resultDF = base::merge(x = resultDF, y = traitDFs[[i]], by.x = "probeID", by.y = "probeID", all.x = TRUE, all.y = TRUE)
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

