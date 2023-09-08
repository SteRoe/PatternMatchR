#' getTraitClusterMedoids
#' gets Medoids of earlier produced clusters
#' @param clustResTraits (reduced) clustering results of traits
#' @param distMatTraits (reduced) distance matrix of traits
#' @param numClusters number of clusters to produce
#' @return list with cluster medoids
#' examples getTraitClusterMedoids(clustResTraits, distMatTraits, numClusters)
getTraitClusterMedoids <- function(clustResTraits, distMatTraits, numClusters) {
  base::tryCatch(
    {
      base::print(base::paste0(sysTimePID(), " start find cluster medoids."))
      if (is.valid(clustResTraits) && numClusters > 1) {
        i <- NULL
        numClusters <- numClusters
        ClustRes <- clustResTraits
        base::print(base::paste0(sysTimePID(), " cutting tree."))
        Clusters <- cutree(ClustRes, k = numClusters)
        numClusters <- max(Clusters)
        base::print(base::paste0(sysTimePID(), " defined ", numClusters, " clusters."))
        ClusterMedoids <- base::list()
        base::print(base::paste0(sysTimePID(), " calculating clusters."))
        #iterate over each cluster:
        foreach::foreach(i = 1:numClusters) %do% {
          #select traits in cluster from DistMat
          idx <- which(Clusters == i)
          if (length(idx) > 1) {
            distMatCluster <- usedist::dist_subset(distMatTraits, idx)
            dM <- as.matrix(distMatCluster)
            CS <- colSums(dM^2)
            min(CS)
            #select trait with lowest (square) sum of distances as cluster medoid
            ClusterMedoid <- which(CS == min(CS))
            ClusterMedoid <- labels(ClusterMedoid)[[1]] #take names from ClusterMedoid
          }
          else {
            ClusterMedoid <- ClustRes$labels[[i]] #take name from ClustRes
          }
          ClusterMedoids[[i]] <- ClusterMedoid
        }
      }
      else {
        ClusterMedoids <- NULL
      }
    },
    error = function(e) {
      message("An error occurred in getTraitClusterMedoids():\n", e)
    },
    warning = function(w) {
      message("A warning occurred in getTraitClusterMedoids():\n", w)
    },
    finally = {
      return(ClusterMedoids)
      base::print(base::paste0(sysTimePID(), " end find cluster medoids."))
    }
  )
}

#' getClustersTable
#' gets a table earlier produced clusters
#' @param listClusters list with clusters
#' @param listMedoids list with medoids
#' @return data.frame with information on clusters and its medoids
#' examples getClustersTable(listClusters, listMedoids)
getClustersTable <- function(listClusters, listMedoids) {
  if (is.valid(listClusters)) {
    result <- as.data.frame(listClusters)
    colnames(result)[[1]] <- "Cluster #"
    result["Medoid name"] <- NA
    i <- NULL
    foreach::foreach(i = seq_along(listClusters)) %do% {
      result[i, 2] <- listMedoids[[result[i, 1]]]
    }
    return(result)
  }
}

#' getMedoidsTable
#' gets a table with Medoids of earlier produced clusters
#' @param listMedoids list with medoids
#' @return data.frame with information on medoids
#' examples getMedoidsTable(listMedoids)
getMedoidsTable <- function(listMedoids) {
  if (is.valid(listMedoids)) {
    result <- listMedoids
    result <- as.data.frame(result)
    colnames(result) <- NULL
    result <- t(result)
    colnames(result)[[1]] <- "Medoid"
    return(result)
  }
}

getPlot <- function(plotObject) {
  if (is.valid(plotObject)) {
    return(plot(plotObject))
  }
}

#' getDendTraits
#' gets a dendrogram with clustering results of traits
#' @param clustResTraits clustering results for traits
#' @param traitClusters list with traits for coloring
#' @return dendrogram
#' examples getDendTraits(clustResTraits, traitClusters)
getDendTraits <- function(clustResTraits, traitClusters) {
  base::tryCatch(
    {
      if (is.valid(clustResTraits) && traitClusters > 1) {
        base::print(base::paste0(sysTimePID(), " start getDendTraits()."))
        dendTraits <- stats::as.dendrogram(clustResTraits)
        dendTraits <- dendextend::color_branches(dendTraits, k = traitClusters)
        dendTraits <- dendextend::color_labels(dendTraits, k = traitClusters)
      }
      else {
        dendTraits <- getEmptyPlot()
#        dendTraits <- plot.new() #plot.window(xlim = c(0 , 1), ylim = c( 5, 10)) #NULL
      }
    },
    error = function(e) {
      message("An error occurred in getDendTraits():\n", e)
      browser()
    },
    warning = function(w) {
      message("A warning occurred in getDendTraits():\n", w)
      browser()
    },
    finally = {
      return(dendTraits)
      base::print(base::paste0(sysTimePID(), " end getDendTraits()."))
    }
  )
}

#' getTraitReducedcombinedDFP_Val_Labels
#' gets a list with data structure for regression results after reduction for traits
#' @param pReducedcombinedDFP_Val_Labels list with data structure for regression results (from p-val reduction)
#' @param ClusterMedoids medoids of clusters for name setting
#' @param keys key attributes for merging of original data
#' @return list with data structure for regression model results
#' examples getTraitReducedcombinedDFP_Val_Labels(pReducedcombinedDFP_Val_Labels, ClusterMedoids, keys)
getTraitReducedcombinedDFP_Val_Labels <- function(pReducedcombinedDFP_Val_Labels, ClusterMedoids, keys) {
  base::tryCatch(
    {
      base::print(base::paste0(sysTimePID(), " start gettraitReducedcombinedDFP_Val_Labels()."))
      traits <- getTraitSubsetcombinedDFP_Val_Labels(pReducedcombinedDFP_Val_Labels, ClusterMedoids, keys)  #select medoids from traits
    },
    error = function(e) {
      message("An error occurred in getTraitReducedcombinedDFP_Val_Labels():\n", e)
    },
    warning = function(w) {
      message("A warning occurred in getTraitReducedcombinedDFP_Val_Labels():\n", w)
    },
    finally = {
      return(traits)
      base::print(base::paste0(sysTimePID(), " end gettraitReducedcombinedDFP_Val_Labels()."))
    }
  )
}

#' getTraitSubsetcombinedDFP_Val_Labels
#' omits traits from combinedDFP_Val_Labels structure, traits to retain are defined in param traits
#' @param combinedDFP_Val_Labels structure, where traits should be omitted
#' @param traits traits to retain
#' @param keys key attributes to retain in data
#' @return combinedDFP_Val_Labels structure list()
#' examples loadResultDF(combinedDFP_Val_Labels, traits)
getTraitSubsetcombinedDFP_Val_Labels <- function(combinedDFP_Val_Labels, traits, keys) {
  base::tryCatch(
    {
      base::print(base::paste0(sysTimePID(), " start getTraitSubsetcombinedDFP_Val_Labels()."))
      if (is.numeric(traits)) {
        traitNos <- traits
      }
      else {
        traitNos <- which(combinedDFP_Val_Labels$mergedColnames %in% traits)
      }
      result <- base::list(dfP_Val = NULL, dfDM = NULL, dfN = NULL,
                           labelsDF1 = NULL, labelsDF2 = NULL, labelsDF3 = NULL,
                           mergedOriginDF = NULL, mergedColnames = NULL,
                           mergedOriginTrait = NULL, mergedDFList = NULL)
      mergedDFList <- list()

      result$dfP_Val <- combinedDFP_Val_Labels$dfP_Val[, traitNos]
      result$dfDM <- combinedDFP_Val_Labels$dfDM[, traitNos]
      result$dfN <- combinedDFP_Val_Labels$dfN[, traitNos]
      result$mergedOriginDF <- combinedDFP_Val_Labels$mergedOriginDF[traitNos]
      result$mergedColnames <- combinedDFP_Val_Labels$mergedColnames[traitNos]
      result$mergedOriginTrait <- combinedDFP_Val_Labels$mergedOriginTrait[traitNos]
      #result$mergedDFList <- combinedDFP_Val_Labels$mergedDFList[traitNos]
      #mergedDFList is more complex: it consists of three parts (red, green, blue) with PHENODF and PHENOFileName; merge them here:
      #DF1
      if (length(combinedDFP_Val_Labels$labelsDF1) > 0) {
        traitNos <-  which(combinedDFP_Val_Labels$labelsDF1 %in% result$mergedColnames)
        if (length(traitNos) > 0) {
          mergedDFL <- list(PHENODF = NULL, PHENOFileName = NULL)
          labelsDF <- combinedDFP_Val_Labels$labelsDF1[traitNos]
          result$labelsDF1 <- labelsDF
          traitVars <- which(colnames(combinedDFP_Val_Labels$mergedDFList[[1]]$PHENODF) %in% labelsDF)
          #also retain key attributes
          keyVars <- which(colnames(combinedDFP_Val_Labels$mergedDFList[[1]]$PHENODF) %in% keys)
          keyTraitVars <- c(keyVars, traitVars)

          labelsDF <- labelsDF[traitVars]
          mergedDFL$PHENODF <- combinedDFP_Val_Labels$mergedDFList[[1]]$PHENODF[, keyTraitVars]
          mergedDFL$PHENOFileName <- combinedDFP_Val_Labels$mergedDFList[[1]]$PHENOFileName
          mergedDFList[[1]] <- mergedDFL
        }
      }
      #DF2
      if (length(combinedDFP_Val_Labels$labelsDF2) > 0) {
        traitNos <-  which(combinedDFP_Val_Labels$labelsDF2 %in% result$mergedColnames)
        if (length(traitNos) > 0) {
          mergedDFL <- list(PHENODF = NULL, PHENOFileName = NULL)
          labelsDF <- combinedDFP_Val_Labels$labelsDF2[traitNos]
          result$labelsDF2 <- labelsDF #combinedDFP_Val_Labels$labelsDF2[traitNos]
          traitVars <- which(colnames(combinedDFP_Val_Labels$mergedDFList[[2]]$PHENODF) %in% labelsDF)
          keyVars <- which(colnames(combinedDFP_Val_Labels$mergedDFList[[2]]$PHENODF) %in% keys)
          keyTraitVars <- c(keyVars, traitVars)

          labelsDF <- labelsDF[traitVars]
          mergedDFL$PHENODF <- combinedDFP_Val_Labels$mergedDFList[[2]]$PHENODF[, keyTraitVars]
          mergedDFL$PHENOFileName <- combinedDFP_Val_Labels$mergedDFList[[2]]$PHENOFileName
          mergedDFList[[2]] <- mergedDFL
        }
      }
      #DF3
      if (length(combinedDFP_Val_Labels$labelsDF3) > 0) {
        traitNos <-  which(combinedDFP_Val_Labels$labelsDF3 %in% result$mergedColnames)
        if (length(traitNos) > 0) {
          mergedDFL <- list(PHENODF = NULL, PHENOFileName = NULL)
          labelsDF <- combinedDFP_Val_Labels$labelsDF3[traitNos]
          result$labelsDF3 <- labelsDF #combinedDFP_Val_Labels$labelsDF3[traitNos]
          traitVars <- which(colnames(combinedDFP_Val_Labels$mergedDFList[[3]]$PHENODF) %in% labelsDF)
          keyVars <- which(colnames(combinedDFP_Val_Labels$mergedDFList[[3]]$PHENODF) %in% keys)
          keyTraitVars <- c(keyVars, traitVars)

          labelsDF <- labelsDF[traitVars]
          mergedDFL$PHENODF <- combinedDFP_Val_Labels$mergedDFList[[3]]$PHENODF[, keyTraitVars]
          mergedDFL$PHENOFileName <- combinedDFP_Val_Labels$mergedDFList[[3]]$PHENOFileName
          mergedDFList[[3]] <- mergedDFL
        }
      }
      result$mergedDFList <- mergedDFList
      #browser() #check result
    },
    error = function(e) {
      message("An error occurred in getTraitSubsetcombinedDFP_Val_Labels():\n", e)
    },
    warning = function(w) {
      message("A warning occurred in getTraitSubsetcombinedDFP_Val_Labels():\n", w)
    },
    finally = {
      base::print(base::paste0(sysTimePID(), " size of merged data.frame: ", dim(result$dfP_Val), " ."))
      base::print(base::paste0(sysTimePID(), " end getTraitSubsetcombinedDFP_Val_Labels()."))
      return(result)
    }
  )
}

#' getplotClustergramTraitsLong
#' produces a clustergram for traits
#' @param matP_Val.t transposed matrix with result p-values
#' @param clustResTraits (reduced) clustering results for traits
#' @param traitClusters trait cluster medoids
#' @return plot with clustergram
#' examples getplotClustergramTraitsLong(matP_Val.t, clustResTraits, traitClusters)
getplotClustergramTraitsLong <- function(matP_Val.t, clustResTraits, traitClusters) {
  base::print(base::paste0(sysTimePID(), " start plotting clustergram."))
  base::tryCatch(
    {
      if (traitClusters > 1) {
        Clusters <- dendextend::cutree(clustResTraits, k = traitClusters) # Clusters <- cutree(clustResTraits, k = traitClusters)
        mycolors <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(traitClusters)
        p <- factoextra::fviz_cluster(list(data = matP_Val.t, cluster = Clusters),
                                      palette = mycolors, #palette = "jco",
                                    ggtheme = ggplot2::theme_classic())
        # p <- factoextra::fviz_mclust(list(data = matP_Val.t, cluster = Clusters),
        #                               palette = "jco",
        #                               ggtheme = ggplot2::theme_classic())
      }
      else {
        p <- getEmptyPlot()
      }
    },
    error = function(e) {
      message("An error occurred in getplotClustergramTraitsLong():\n", e)
    },
    warning = function(w) {
      message("A warning occurred in getplotClustergramTraitsLong():\n", w)
    },
    finally = {
      base::print(base::paste0(sysTimePID(), " end plotting clustergram."))
      return(p)
    }
  )
}
