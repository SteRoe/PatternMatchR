
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
    #startTime <- Sys.time()
    gc()
    ClustRes <- fastcluster::hclust(stats::as.dist(distanceMatrix), method = "ward.D2")
    # endTime <- Sys.time()
    # runTime <- difftime(endTime, startTime, units = "secs")
    # base::tryCatch({
    #   base::print(base::paste0(sysTimePID(), " runtime for FastCluster: ", runTime))
    # })
  }
  else {
    ClustRes <- NULL
    base::message(base::paste0(sysTimePID(), " is.valid(distanceMatrix) == FALSE"))
  }
  return(ClustRes)
}

#' #' updateTxtClusterOut
#' #' generates summary text after clustering
#' #' @param traitReducedcombinedDFP_Val_Labels data structure with trait reduced results
#' #' @param minP_Val minimum p_value for model to use for clustering
#' #' @param maxP_Val maximum p_value for model to use for clustering
#' #' @param minN minimum n for model to use for clustering
#' #' @param sldNumClasses number of classes to use for clustering
#' #' @return text
#' #' examples updateTxtClusterOut(traitReducedcombinedDFP_Val_Labels, minP_Val, maxP_Val, minN, sldNumClasses)
#' updateTxtClusterOut <- function(traitReducedcombinedDFP_Val_Labels, minP_Val, maxP_Val, minN, sldNumClasses) {
#'   base::tryCatch(
#'    {
#'       result <- NULL
#'       if (is.valid(traitReducedcombinedDFP_Val_Labels)) {
#'         maxClasses <- length(traitReducedcombinedDFP_Val_Labels$mergedColnames)
#'         numRow <- nrow(traitReducedcombinedDFP_Val_Labels$dfP_Val)
#'         numCol <- ncol(traitReducedcombinedDFP_Val_Labels$dfP_Val)
#'         minClasses <- 1 #dendextend::min_depth(session$userData$sessionVariables$dendTraits)
#'         result <- base::paste0("finding trait clusters successful. found minClusters = ",
#'                             minClasses, "; maxClusters: ", maxClasses, "; Clustering made for numClasses = ", sldNumClasses, ".\n",
#'                             "size of result df: nrow(CpG)=", numRow, ", ncol(trait)=", numCol, ".")
#'       }
#'     },
#'     error = function(e) {
#'       base::message("An error occurred in updateTxtClusterOut():\n", e)
#'     },
#'     warning = function(w) {
#'       base::message("A warning occurred in updateTxtClusterOut():\n", w)
#'     },
#'     finally = {
#'       return(shiny::HTML(result))
#'     }
#'   )
#' }

#' calculateDistanceNeigboursProbes
#' calculate distance from each probe to its neigbours and gives back data frame with distance metrics
#' @param wd working directory
#' @param clustResProbes data structure with clustering result
#' @param annotation annotation of CpG (names, location etc.)
#' @param distanceToLook maximum distance to look for
#' @param numCores number of cores to use for distance calculation
#' @return data.frame with min, mean, max distance and nuber of CpG in distanceToLook window
#' examples calculateDistanceNeigboursProbes(wd, clustResProbes, annotation, distanceToLook, numCores)
calculateDistanceNeigboursProbes <- function(wd, clustResProbes, annotation, distanceToLook, numCores) {
  base::tryCatch(
    {
      base::print(base::paste0(sysTimePID(), " start calculateDistanceNeigboursProbes() with max distance (distanceToLook) = ", distanceToLook, "."))
      #get chr and location from annotation
      maxDistanceToLook <- distanceToLook
      annotation <- subset(annotation,select = c("name", "chromosome","position"))
      #merge annotation
      CpG <- as.data.frame(clustResProbes$labels[clustResProbes$order])
      colnames(CpG)[1] <- "label"
      CpG$order <- seq(1:nrow(CpG))
      distances <- base::merge(CpG, annotation, by.x = "label", by.y = "name")
      DNAdistancesUp <- base::data.frame(seq_along(distances[,2]), 2)
      DNAdistances <- base::data.frame(seq_along(distances[,2]), 5)
      #sort order given by clustering
      distances <- distances[base::order(distances$order),]
#      library(future) #we have this already in DESCRIPTION file, but without "library(future)" here, it won't work. Strange.
#      library(doFuture)
#      future::plan(strategy = future::multisession, workers = numCores)
      #calculate mean distance to distanceToLook next probes in given order, omit probes with different chr
      #      foreach::foreach(i = seq_along(distances[,2]), .combine = rbind, .verbose = TRUE) %dofuture% { #for all objects in distances
      foreach::foreach(i = seq_along(distances[,2])) %do% { #for all objects in distances
      #for(i in seq_along(distances[,2])) { #for all objects in distances
#        base::source(paste0(wd, "/R/Clustering.R")) #this is necessary for foreach %dopar% to run properly
        currentCpG <- distances[i,]
        #cut out defined number of probes (up- and downstream from CpG)
        lowerBorder <- i-distanceToLook
        if (lowerBorder < 1) {lowerBorder <- 1}#
        upperBorder <- i+distanceToLook
        if (upperBorder > nrow(distances)) {upperBorder <- nrow(distances)}
        distancesToLook <- distances [lowerBorder:upperBorder, ]
        #do not iterate over all probes, but exclude those from different chr first
        chr <- currentCpG$chr
        distancesToLook <- distancesToLook[distancesToLook$chr == chr,]
        distanceToLook <- nrow(distancesToLook)
        if (distanceToLook > 1) {
          DNAdistancesUp <- base::data.frame(seq_along(1:distanceToLook), 2)
          #upstream
          foreach::foreach(j = 1:distanceToLook) %do% { #max. distance given by distanceToLook
          #for(j in 1:distanceToLook) { #max. distance given by distanceToLook
            CpG <- distances[j,]
            if (currentCpG$label != CpG$label) {
              DNAdistancesUp[j,2] <- base::abs(currentCpG$position - CpG$position)
            }
            else {
              DNAdistancesUp[j,2] <- NA
            }
            DNAdistancesUp[j,1] <- currentCpG$label
          }
#browser() #check number of nearby CpG
          if(is.numeric(DNAdistancesUp[j,2])) {
            DNAdistances[i,2] <- base::min((DNAdistancesUp[,2]), na.rm = TRUE)
            DNAdistances[i,3] <- base::mean((DNAdistancesUp[,2]), na.rm = TRUE)
            DNAdistances[i,4] <- base::max((DNAdistancesUp[,2]), na.rm = TRUE)
            DNAdistances[i,5] <- base::length(na.omit(DNAdistancesUp[,2]))
          }
          else {
            DNAdistances[i,2] <- NA
            DNAdistances[i,3] <- NA
            DNAdistances[i,4] <- NA
            DNAdistances[i,5] <- NA
          }
          DNAdistances[i,1] <- currentCpG$label
        }
        else {
          #no near CpG on the same chromosome found
          DNAdistances[i,2] <- NA
          DNAdistances[i,3] <- NA
          DNAdistances[i,4] <- NA
          DNAdistances[i,5] <- NA
        }
        DNAdistances[i,1] <- currentCpG$label
      }
      colnames(DNAdistances) <- c("ID", "minDistance", "meanDistance", "maxDistance", "number")
      distances <- na.omit(DNAdistances)
      base::message(base::paste0(sysTimePID(), " found n = ", nrow(distances) , " neigbouring CpG with distance <=", maxDistanceToLook,""))
#tbc(): list works, but tibble not
#      DNAdistances <- tibble::rownames_to_column(as.data.frame(DNAdistances), var = "rowname")
    },
    error = function(e) {
      base::message("An error occurred in calculateDistanceNeigboursProbes():\n", e)
      browser()
    },
    warning = function(w) {
      base::message("A warning occurred in calculateDistanceNeigboursProbes():\n", w)
      browser()
    },
    finally = {
      base::print(base::paste0(sysTimePID(), " end calculateDistanceNeigboursProbes()."))
      return(DNAdistances)
    }
  )
}
