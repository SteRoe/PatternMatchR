# getDistMatSPAM <- function(session, matrix) {
#   if (!is.null(session$userData$numCores)) {
#     numberCores <- session$userData$numCores - 5
#   } else {
#     numberCores <- 1
#   }
#   gc()
#   distMat <- spam::as.spam.dist(
#     parallelDist::parallelDist(
#       base::as.matrix(matrix),
#       method = "euclidean",
#       #labels = TRUE,
#       threads = numberCores
#     )
#   )
#   return(distMat)
# }

#' getDistMat
#' calculates distance matrix using parallelDist
#' @param numberCores number of CPU cores to use
#' @param matrix matrix to calculate distance matrix for
#' @return distance matrix
#' examples getDistMat(numberCores, matrix)
getDistMat <- function(numberCores, matrix) {
  #distMat <- spam::as.spam.dist(
  distMat <- stats::as.dist(
    parallelDist::parallelDist(
      base::as.matrix(matrix),
      method = "euclidean",
      labels = TRUE,
      threads = numberCores
    )
  )
  return(distMat)
}

#' getClustResFast
#' calculates hierarchical clustering using fastcluster
#' @param distanceMatrix distance matrix
#' @return hclust object
#' examples getClustResFast(distanceMatrix)
getClustResFast <- function(distanceMatrix) {
  base::print(base::paste0(sysTimePID(), " start getClustResFast()"))
  #check size of distanceMatrix
  if (is.valid(distanceMatrix)) {
    base::print(base::paste0(sysTimePID(), " clustering trait data."))
    startTime <- Sys.time()
    gc()
    ClustRes <- fastcluster::hclust(stats::as.dist(distanceMatrix), method = "ward.D2")
    endTime <- Sys.time()
    runTime <- difftime(endTime, startTime, units = "secs")
    base::print(base::paste0(sysTimePID(), " runtime for FastCluster: ", runTime))
  }
  else {
    ClustRes <- NULL
    base::print(base::paste0(sysTimePID(), " is.valid(distanceMatrix) == FALSE"))
  }
  return(ClustRes)
}

#' updateTxtClusterOut
#' generates summary text after clustering
#' @param traitReducedcombinedDFP_Val_Labels data structure with trait reduced results
#' @param minP_Val minimum p_value for model to use for clustering
#' @param maxP_Val maximum p_value for model to use for clustering
#' @param minN minimum n for model to use for clustering
#' @param sldNumClasses number of classes to use for clustering
#' @return text
#' examples updateTxtClusterOut(traitReducedcombinedDFP_Val_Labels, minP_Val, maxP_Val, minN, sldNumClasses)
updateTxtClusterOut <- function(traitReducedcombinedDFP_Val_Labels, minP_Val, maxP_Val, minN, sldNumClasses) {
  base::tryCatch(
   {
      result <- NULL
      if (is.valid(traitReducedcombinedDFP_Val_Labels)) {
        maxClasses <- length(traitReducedcombinedDFP_Val_Labels$mergedColnames)
        numRow <- nrow(traitReducedcombinedDFP_Val_Labels$dfP_Val)
        numCol <- ncol(traitReducedcombinedDFP_Val_Labels$dfP_Val)
        minClasses <- 1 #dendextend::min_depth(session$userData$sessionVariables$dendTraits)
        result <- base::paste0("finding trait clusters successful. found minClusters = ",
                            minClasses, "; maxClusters: ", maxClasses, "; Clustering made for numClasses = ", sldNumClasses, ".\n",
                            "size of result df: nrow(CpG)=", numRow, ", ncol(trait)=", numCol, ".")
      }
    },
    error = function(e) {
      message("An error occurred in updateTxtClusterOut():\n", e)
    },
    warning = function(w) {
      message("A warning occurred in updateTxtClusterOut():\n", w)
    },
    finally = {
      return(shiny::HTML(result))
    }
  )
}
