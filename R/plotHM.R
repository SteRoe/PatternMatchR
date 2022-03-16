
#' getResultDfP_D_N
#' @param listOfResultDF data.frame containing a list of data.frames to be merged
#' @param P_D_N scalar value "P", "D" or "N" describing whether to merge P_VAL, DeltaMeth or N
#' @return merged data.frame
# examples getResultDfP_D_N(listDF, "P")
getResultDfP_D_N <- function(listOfResultDF,P_D_N) {
#browser()
  print(paste0(Sys.time(), " start getResultDfP_D_N()"))
  i <- NULL
  if (length(listOfResultDF) != 0) {
    foreach::foreach(i = 1:length(listOfResultDF)) %do% {
      print(paste0(Sys.time(), " processing trait ", i))
      if (P_D_N == "P") {
        DF = subset(listOfResultDF[[i]], select=c("probeID","P_VAL"))
      }
      else if (P_D_N == "D") {
        DF = subset(listOfResultDF[[i]], select=c("probeID","DeltaMeth"))
      }
      else if (P_D_N == "N") {
        DF = subset(listOfResultDF[[i]], select=c("probeID","N"))
      }
      if (i == 1) {
        merged = DF
        rownames(merged) = DF$probeID
      }
      else {
        print(paste0(Sys.time(), " merge"))
        #merge
        merged = base::merge(x = merged, y = DF, by.x = "probeID", by.y = "probeID", all.x = TRUE, all.y = TRUE)
      }
      colnames(merged)[i+1] = names(listOfResultDF)[i]
    }
  }
  rownames(merged) = merged$probeID
  merged$probeID = NULL
  print(paste0(Sys.time(), " finished getResultDfP_D_N()"))
  return(merged)
}

#' getlistOfResultsDF
#' gets list of result data.frame from folder
#' @param folder folder containing all files to read to list
#' @return list of data.frame from folder
# examples getlistOfResultsDF(folder)
getlistOfResultsDF <- function(folder) {
  #  if (isTruthy(folder)) {
  print(paste0(Sys.time()," before list.files()"))
  if (dir.exists(folder)) {
    temp <- list.files(path=folder,pattern="*.csv")
    result = list()
    for (i in 1:length(temp)) {
      trait <- stringr::str_sub(temp[i], 1, stringr::str_length(temp[i])-4)
      fileName = paste0(folder, as.character(trait),".csv")
      if (file.exists(fileName)) {
        firstlines <- utils::read.table(file = fileName, sep = "\t", header = T, nrows = 5)
        if (colnames(firstlines)[1] == "probeID") {
          if (nrow(firstlines) >= 5) {
  #          if (grepl("adj", temp[i], fixed = TRUE) == TRUE) {
              print(paste0(Sys.time()," load result file ", folder, trait))
              #read results into DF
              resultDF = loadResultFile(folder, trait)
              #omit unneccesary variables
              resultDF = resultDF[,c("probeID","P_VAL","DeltaMeth","N")]
              if (min(resultDF$P_VAL, na.rm = TRUE) < as.numeric(globalVariables$config$P_VALWarningThreshold)) {
                message (paste0(Sys.time(), " p-value below warning threshold found in ", folder))
              }
                            resultDF = list(resultDF)
              names(resultDF) = trait
              result = append(result,resultDF)
  #          }
          }
        }
      }
    }
#    saveRDS (result,file="listOfResultsDF.RDS")
    return (result)
  }
  else {
    message(paste0("folder ", folder, " does not exist."))
  }
}

#' loadFolderDFList
#' loads list of data.frame from folder
#' @param folder folder containing all files to read to list
#' @return list of data.frame from folder
# examples loadFolderDFList(folder)
loadFolderDFList <- function(folder) {
  fileNameLR = "listOfResultsDF.RDS"
  fileNameLR = paste0(folder,fileNameLR)
#browser()
  #if (file_test("-f", fileNameLR) == TRUE) {
  if (utils::file_test("-f", fileNameLR) == TRUE && getYoungestFile(folder) == fileNameLR) {
    print(paste0(Sys.time()," readRDS getlistOfResultsDF()", fileNameLR))
    listOfResultDF = readRDS (file=fileNameLR)
    getList = FALSE
    if (length(listOfResultDF) == 0) {
      getList = TRUE
    }
    else {
      getList = FALSE
    }
  }
  else {
    getList = TRUE
  }
  if (getList == TRUE) {
    #        listOfResultDF = getlistOfResultsDF(dataDir())
    print(paste0(Sys.time()," before getlistOfResultsDF()"))
    listOfResultDF = getlistOfResultsDF(folder)
    print(paste0(Sys.time()," saveRDS getlistOfResultsDF()", fileNameLR))
    saveRDS(listOfResultDF, file = fileNameLR)
  }
#browser() #listOfResultDF nicht gefunden
  return(listOfResultDF)
}

#' loadResultDF
#' loads data.frame from folder
#' @param folder folder containing all files to read to list
#' @return list with resultP_Val, resultDeltaMethylation, resultN from folder
# examples loadResultDF(folder)
loadResultDF  <- function(folder) {
  fileNameLPV = "listResultP_Val_DeltaMeth_N.RDS"
  fileNameLPV = paste0(folder,fileNameLPV)
  #if (file_test("-f", fileNameLPV) == TRUE) {
#  if (FALSE) {
  if (utils::file_test("-f", fileNameLPV) == TRUE && getYoungestFile(folder) == fileNameLPV) {
    print(paste0(Sys.time()," readRDS loadFolderDFList() ", fileNameLPV,"."))
    listResultP_Val_DeltaMeth_N = readRDS (file=fileNameLPV)
    print(paste0(Sys.time()," finished read files ", fileNameLPV,"."))
  }
  else {
    #load or generate list of DFs with results
    print(paste0(Sys.time()," before loadFolderDFList()"))
    listOfResultDF = loadFolderDFList(folder)
    if (length(listOfResultDF) != 0) {
      listResultP_Val = getResultDfP_D_N (listOfResultDF, "P")
      listResultDeltaMeth = getResultDfP_D_N (listOfResultDF, "D")
      listResultN = getResultDfP_D_N (listOfResultDF, "N")
      listResultP_Val_DeltaMeth_N=list(P_Val = listResultP_Val, DM = listResultDeltaMeth, N = listResultN)
      print(paste0(Sys.time()," saveRDS loadFolderDFList()", fileNameLPV, "."))
      saveRDS(listResultP_Val_DeltaMeth_N, file = fileNameLPV)
    }
    else {
      browser()
      message(paste0(Sys.time()," length(listOfResultDF) == 0."))
    }
  }

  if(exists("listResultP_Val_DeltaMeth_N")) {
    rn = rownames(listResultP_Val_DeltaMeth_N$P_Val)
    if (!is.data.frame(listResultP_Val_DeltaMeth_N$P_Val)) {
      listResultP_Val_DeltaMeth_N$P_Val = as.data.frame(listResultP_Val_DeltaMeth_N$P_Val)
      rownames(listResultP_Val_DeltaMeth_N$P_Val) = rn
    }
    if (!is.data.frame(listResultP_Val_DeltaMeth_N$DM)) {
      listResultP_Val_DeltaMeth_N$DM = as.data.frame(listResultP_Val_DeltaMeth_N$DM)
      rownames(listResultP_Val_DeltaMeth_N$DM) = rn
    }
    if (!is.data.frame(listResultP_Val_DeltaMeth_N$N)) {
      listResultP_Val_DeltaMeth_N$N = as.data.frame(listResultP_Val_DeltaMeth_N$N)
      rownames(listResultP_Val_DeltaMeth_N$N) = rn
    }
    rownames(listResultP_Val_DeltaMeth_N$P_Val) = rn
    rownames(listResultP_Val_DeltaMeth_N$DM) = rn
    rownames(listResultP_Val_DeltaMeth_N$N) = rn
    if (globalVariables$config$debugMode == TRUE) {
      rn = rownames(listResultP_Val_DeltaMeth_N$P_Val)[1:1000]
      listResultP_Val_DeltaMeth_N$P_Val = listResultP_Val_DeltaMeth_N$P_Val[1:1000,]
      listResultP_Val_DeltaMeth_N$DM = listResultP_Val_DeltaMeth_N$DM[1:1000,]
      listResultP_Val_DeltaMeth_N$N = listResultP_Val_DeltaMeth_N$N[1:1000,]
      rownames(listResultP_Val_DeltaMeth_N$P_Val) = rn
      rownames(listResultP_Val_DeltaMeth_N$DM) = rn
      rownames(listResultP_Val_DeltaMeth_N$N) = rn
    }
    print(paste0(Sys.time()," finished loadFolderDFList()", "."))
    return (listResultP_Val_DeltaMeth_N)
  }
}

#' getAvailNForP_VALBorder
#' counts traits with at least 2 elements > P_VAL_BORDER
#' @param DF data.frame with P_Val
#' @return data.frame with reduced data set
# examples getAvailNForP_VALBorder(data.frame)
getAvailNForP_VALBorder <- function(DF) {
  numrows = 300
  result = matrix(nrow = numrows,ncol = 2)
  for (i in 1:numrows) {
    mat = DF
    P_VAL_BORDER = 5 * 10^-i
    mat[mat > P_VAL_BORDER] = NA
    mat = delete.na(mat,ncol(mat)-1) # -1, because we need at least 2 traits to associate
    n = nrow(mat)
    if (!is.numeric(n)) break()
    if (n == 0) break()
    print(paste0(Sys.time(), " counting remaining probes at p = "), P_VAL_BORDER, " remaining n = ", n)
    result[i,1] = P_VAL_BORDER
    result[i,2] = n
  }
  colnames(result) = c("P_VAL_BORDER","Available n")
  result = result[1:i-1,]
  result = as.data.frame((result))
  result <- result[order(result[1]),]
  return (result)
}

#' getReducedP_ValMatrix
#' counts and gets back traits with at least 2 elements > P_VAL_BORDER
#' @param df data.frame with P_Val, rows are probes, cols are traits
#' @param numRows scalar with maximum allowed row number
#' @param upperP_VALborder maximum allowed P_Val
#' @param lowerP_VALborder optional minimum allowed P_Val
#' @return matrix with reduced data set
# examples getReducedP_ValMatrix(matrix)
getReducedP_Valdf <- function(df, numRows, upperP_VALborder, lowerP_VALborder) { #getReducedP_Valdf <- function(df, numRows, P_VALborder) {
    if (!is.data.frame(df)) {
#browser()
    df=as.data.frame(df)
  }
  if (missing(upperP_VALborder)) {
    upperP_VALborder = 0.05
  }
  orig_df = df
  # remove rowsProbes with all p-val > 0.05
  df[df > upperP_VALborder] = NA
#  if (!missing(lowerP_VALborder) && length(lowerP_VALborder) == 0) {
  if (!missing(lowerP_VALborder)) {
#browser()
    df[df < lowerP_VALborder] = NA
  }
  df = delete.na(df,ncol(df)-1) # -1, because we need at least 2 traits to associate
  df[is.na(df)] = 1
  result_df = orig_df[rownames(df),colnames(df)]
  print(paste0("upperP_VALborder: ", upperP_VALborder))
  print(paste0("lowerP_VALborder: ", lowerP_VALborder))
  print(paste0("nrow result matrix: ", nrow(df)))
  if (!missing(numRows)) {
    if (is.numeric(numRows)) {
#      if (numRows < nrow(mat)) {
      if (numRows < nrow(result_df)) {
#        mat = mat[1:numRows,]
        result_df = result_df[1:numRows,]
      }
    }
  }
#  return(mat)
  print(dim(result_df))
  return(result_df)
}

#' removeTraitsMinN
#' removes CpG with casecount < minN from corresponding data.frames for P_Val, for DeltaMethylation and for N
#' @param dfList list containing corresponding data.frames for P_Val, for DeltaMethylation and for N
#' @param minN minimum value for n
#' @return named list of data.frames, one df for P_Val, one for DeltaMethylation, one for N as well as labels
#' @return result$dfP_Val for p-values
#' @return result$dfDM for delta methylation values
#' @return result$dfN for n
# examples getCombinedDFP_Val_Labels(data.frame, minN)
removeTraitsMinN <- function(dfList, minN) {
  #check for minimum n in each trait
  dfN = dfList$dfN
  rn = rownames(dfN)
  dfP_Val = dfList$dfP_Val
  dfDM = dfList$dfDM
  dfN = as.data.frame(dfN)
  dfP_Val = as.data.frame(dfP_Val)
  rownames(dfP_Val) = rn
  dfDM = as.data.frame(dfDM)
  rownames(dfDM) = rn
#browser()
  traitNames = stats::na.omit(colnames(dfN)[matrixStats::colMins(as.matrix(dfN))>100]) #colnames(dfN)[na.omit(colMins(as.matrix(dfN))>as.integer(minN))]
#  rn = rownames(dfN)
  dfN =  dfN[, traitNames] #dfN[, traitNames, with = FALSE] #a = matN[,matN$Name %in% traitNames] subset(matN, colnames(matN) %in% traitNames)
  rownames(dfN) = rn

  #select the same content than in dfN
  dfP_Val = dfP_Val[rownames(dfN),colnames(dfN)] #does not work with data.table
  dfDM = dfDM[rownames(dfN),colnames(dfN)]
  dfList$dfP_Val = dfP_Val
  dfList$dfDM = dfDM
  dfList$dfN = dfN
  return(dfList)
}


#' getCombinedDFP_Val_Labels
#' merges data.frames out of sessionVariable containing p-values, delta methylation values and N by probeID
#' @param sessionVariables sessionVariables
#' @param minN minimum n
#' @return named list of data.frames, one df for merged P_Val, one for merged DeltaMethylation, one for merged N as well as labels:
#' @return result$dfP_Val for p-values
#' @return result$dfDM for delta methylation values
#' @return result$dfN for n
#' @return result$labelsDF1 for labels belonging to original df1
#' @return result$labelsDF2 for labels belonging to original df2
#' @return result$labelsDF3 for labels belonging to original df3
# examples getCombinedDFP_Val_Labels(sessionVariables)
getCombinedDFP_Val_Labels <- function(sessionVariables, minN) {
  continue = FALSE
  #merge three df
#browser()
  if (!is.null(sessionVariables$resultP_ValTrait1)) {#if (!testthat::is_null(sessionVariables$resultP_ValTrait1)) {
    dfList=list()
#browser() # check whether DF or matrix
    dfList$dfP_Val = sessionVariables$resultP_ValTrait1 #class(data.frame)
    dfList$dfDM = sessionVariables$resultDMTrait1
    dfList$dfN = sessionVariables$resultNTrait1
    dfList = removeTraitsMinN(dfList, minN)
    dfP_Val1 = dfList$dfP_Val
    dfDM1 = dfList$dfDM
    dfN1 = dfList$dfN

    if ((nrow(dfP_Val1) > 0) && (ncol(dfP_Val1) > 0)) {
      continue = TRUE
    }
    else {
      continue = FALSE
      print(paste0(Sys.time(), "nrow(DFP_Val1) or ncol(DFP_Val1) == 0"))
    }
  }
  else {
    continue = FALSE
    message(paste0(Sys.time(), "DF1 is not valid"))
  }
  if (!is.null(sessionVariables$resultP_ValTrait2)) {#(!testthat::is_null(sessionVariables$resultP_ValTrait2)) {
    dfList=list()
    dfList$dfP_Val = sessionVariables$resultP_ValTrait2
    dfList$dfDM = sessionVariables$resultDMTrait2
    dfList$dfN = sessionVariables$resultNTrait2
    dfList = removeTraitsMinN(dfList, minN)
    dfP_Val2 = dfList$dfP_Val
    dfDM2 = dfList$dfDM
    dfN2 = dfList$dfN
    if ((nrow(dfP_Val2) > 0) && (ncol(dfP_Val2) > 0)) {
      continue = TRUE
    }
    else {
      continue = FALSE
      message(paste0(Sys.time(), "nrow(DF2) or ncol(DF2) == 0"))
    }
  }
  else {
    continue = FALSE
    message(paste0(Sys.time(), "DF2 is not valid"))
  }

  if (!is.null(sessionVariables$resultP_ValTrait3)) { #(!testthat::is_null(sessionVariables$resultP_ValTrait3)) {
    dfList=list()
    dfList$dfP_Val = sessionVariables$resultP_ValTrait3
    dfList$dfDM = sessionVariables$resultDMTrait3
    dfList$dfN = sessionVariables$resultNTrait3
    dfList = removeTraitsMinN(dfList, minN)
    dfP_Val3 = dfList$dfP_Val
    dfDM3 = dfList$dfDM
    dfN3 = dfList$dfN
    if ((nrow(dfP_Val3) > 0) && (ncol(dfP_Val3) > 0)) {
      continue = TRUE
    }
    else {
      continue = FALSE
      message(paste0(Sys.time(), "nrow(DF3) or ncol(DF3) == 0"))
    }
  }
  else {
    continue = FALSE
    message(paste0(Sys.time(), "DF3 is not valid"))
  }
  if (continue == TRUE) {

    if (exists("dfP_Val1")) {
      mergedDFP_Val = as.data.frame(dfP_Val1)
      #crazy error message here: left side converting to list? due to dfP_Val1 was of class matrix
      mergedDFP_Val$Row.names = rownames(mergedDFP_Val)
      mergedDFDM = as.data.frame(dfDM1)
      mergedDFDM$Row.names = rownames(mergedDFDM)
      mergedDFN = as.data.frame(dfN1)
      mergedDFN$Row.names = rownames(dfN1)
    }
    if (exists("dfP_Val2")) {
      if (exists("mergedDFP_Val")) {
        dfP_Val2 = as.data.frame(dfP_Val2)
        dfP_Val2$Row.names = rownames(dfP_Val2)
        dfDM2 = as.data.frame(dfDM2)
        dfDM2$Row.names = rownames(dfDM2)
        dfN2 = as.data.frame(dfN2)
        dfN2$Row.names = rownames(dfN2)

        mergedDFP_Val = base::merge(mergedDFP_Val, dfP_Val2, by.x = "Row.names", by.y = "Row.names", all.x = FALSE, all.y=FALSE)
        mergedDFDM = base::merge(mergedDFDM, dfDM2, by.x = "Row.names", by.y = "Row.names", all.x = FALSE, all.y=FALSE)
        mergedDFN = base::merge(mergedDFN, dfN2, by.x = "Row.names", by.y = "Row.names", all.x = FALSE, all.y=FALSE)
        if("Row.names" %in% colnames(dfP_Val2)) {
          rownames(dfP_Val2) = dfP_Val2$Row.names
          dfP_Val2$Row.names <- NULL
        }
        if("Row.names" %in% colnames(dfDM2)) {
          rownames(dfDM2) = dfDM2$Row.names
          dfDM2$Row.names <- NULL
        }
        if("Row.names" %in% colnames(dfN2)) {
          rownames(dfN2) = dfN2$Row.names
          dfN2$Row.names <- NULL
        }
      }
      else {
        mergedDFP_Val = as.data.frame(dfP_Val2)
        mergedDFP_Val$Row.names = rownames(mergedDFP_Val)
        mergedDFDM = as.data.frame(dfDM2)
        mergedDFDM$Row.names = rownames(mergedDFDM)
        mergedDFN = as.data.frame(dfN2)
        mergedDFN$Row.names = rownames(dfN2)
      }
    }
    if (exists("dfP_Val3")) {
      if (exists("mergedDFP_Val")) {
        dfP_Val3 = as.data.frame(dfP_Val3)
        dfP_Val3$Row.names = rownames(dfP_Val3)
        dfDM3 = as.data.frame(dfDM3)
        dfDM3$Row.names = rownames(dfDM3)
        dfN3 = as.data.frame(dfN3)
        dfN3$Row.names = rownames(dfN3)

        mergedDFP_Val = base::merge(mergedDFP_Val, dfP_Val3, by.x = "Row.names", by.y = "Row.names", all.x = FALSE, all.y=FALSE)
        mergedDFDM = base::merge(mergedDFDM, dfDM3, by.x = "Row.names", by.y = "Row.names", all.x = FALSE, all.y=FALSE)
        mergedDFN = base::merge(mergedDFN, dfN3, by.x = "Row.names", by.y = "Row.names", all.x = FALSE, all.y=FALSE)
        if("Row.names" %in% colnames(dfP_Val3)) {
          rownames(dfP_Val3) = dfP_Val3$Row.names
          dfP_Val3$Row.names <- NULL
        }
        if("Row.names" %in% colnames(dfDM3)) {
          rownames(dfDM3) = dfDM3$Row.names
          dfDM3$Row.names <- NULL
        }
        if("Row.names" %in% colnames(dfN3)) {
          rownames(dfN3) = dfN3$Row.names
          dfN3$Row.names <- NULL
        }
      }
      else {
        mergedDFP_Val = as.data.frame(dfP_Val3)
        mergedDFP_Val$Row.names = rownames(mergedDFP_Val)
        mergedDFDM = as.data.frame(dfDM3)
        mergedDFDM$Row.names = rownames(mergedDFDM)
        mergedDFN = as.data.frame(dfN3)
        mergedDFN$Row.names = rownames(dfN3)
      }
    }
    if("Row.names" %in% colnames(mergedDFP_Val)) {
      rownames(mergedDFP_Val) = mergedDFP_Val$Row.names
      mergedDFP_Val$Row.names <- NULL
    }
    if("Row.names" %in% colnames(mergedDFDM)) {
      rownames(mergedDFDM) = mergedDFDM$Row.names
      mergedDFDM$Row.names <- NULL
    }
    if("Row.names" %in% colnames(mergedDFN)) {
      rownames(mergedDFN) = mergedDFN$Row.names
      mergedDFN$Row.names <- NULL
    }

    if (exists("dfP_Val1")) {
      splitPointStart = 1
      splitPointEnd = ncol(dfP_Val1)
      LabelsDF1 = colnames(mergedDFP_Val)[splitPointStart:splitPointEnd]
    }
#    splitPoint = splitPoint + 1
    if (exists("dfP_Val2")) {
      if (exists("splitPointEnd")) {
        splitPointStart = splitPointEnd + 1
      }
      else {
        splitPointStart = 1
      }
      splitPointEnd = splitPointStart + ncol(dfP_Val2) - 1
      LabelsDF2 = colnames(mergedDFP_Val)[splitPointStart:splitPointEnd]
    }
#    splitPoint = splitPoint + 1
    if (exists("dfP_Val3")) {
      if (exists("splitPointEnd")) {
        splitPointStart = splitPointEnd + 1
      }
      else {
        splitPointStart = 1
      }
      splitPointEnd = splitPointStart + ncol(dfP_Val3) - 1
      LabelsDF3 = colnames(mergedDFP_Val)[splitPointStart:splitPointEnd]
    }
    result = list()
    result$dfP_Val = mergedDFP_Val #matP_Val
    result$dfDM = mergedDFDM #matDM
    result$dfN = mergedDFN #matN
    if (exists("LabelsDF1")) {
#      result[[3]] = LabelsDF1
      result$labelsDF1 = LabelsDF1
    }
    if (exists("LabelsDF2")) {
#      result[[4]] = LabelsDF2
      result$labelsDF2 = LabelsDF2
    }
    if (exists("LabelsDF3")) {
#      result[[5]] = LabelsDF3
      result$labelsDF3 = LabelsDF3
    }
    return(result)
  }
  else {
    print("DF's could not be merged")
    return(NULL)
  }
}

#getCombinedDFP_Val_Labels <- compiler::cmpfun(getCombinedDFP_Val_Labels)

#' emptyHM
#' creates an empty heatmap
#' @return empty heatmap
# examples emptyHM()
emptyHM <- function() {
  mat = matrix(c(1,2, 3,4), nrow = 2, ncol = 2)
  ht = ComplexHeatmap::Heatmap(mat)
  if(grDevices::dev.cur() > 1) grDevices::dev.off()
  grDevices::pdf(NULL)
  grDevices::dev.off()
  ht = ComplexHeatmap::draw(ht)
}

#emptyHM <- compiler::cmpfun(emptyHM)

#' creates a regular heatmap
#' @param combinedDF_Labels list of data.frame and labels generated from function getCombinedDFP_Val_Labels()
#' @param dendProbes dendrogramProbes for probes (rows), generated externally and providing sorting information for heatmap#
#' @param dendTraits dendrogramTraits for traits (columns), generated externally and providing sorting information for heatmap
#' @return heatmap object for InteractiveComplexHeatmap::makeInteractiveComplexHeatmap
# examples combinedDFInteractiveHeatMapP_Val(combinedDF_Labels, dendProbes, dendTraits)
combinedDFInteractiveHeatMapP_Val <- function(combinedDF_Labels, dendProbes, dendTraits) {
  print(paste0(Sys.time()," start making HM"))
  mat = as.matrix(combinedDF_Labels$dfP_Val)
  matDM = as.matrix(combinedDF_Labels$dfDM)
  matN = as.matrix(combinedDF_Labels$dfN)
# ((((((((((((((((((((((((()))))))))))))))))))))))))
# reduce matrix dimensions before plotting like in
# https://gdevailly.netlify.app/post/plotting-big-matrices-in-r/
# ((((((((((((((((((((((((()))))))))))))))))))))))))
# and use rasterization like described in
# https://jokergoo.github.io/2020/06/30/rasterization-in-complexheatmap/
  print(paste0(Sys.time()," making labels"))
  labelsDF1 = combinedDF_Labels$labelsDF1
  labelsDF2 = combinedDF_Labels$labelsDF2
  labelsDF3 = combinedDF_Labels$labelsDF3
  l1 = rep("trait 1", length(labelsDF1))
  l2 = rep("trait 2", length(labelsDF2))
  l3 = rep("trait 3", length(labelsDF3))
  labels = c(l1, l2, l3)
  while (!is.null(grDevices::dev.list())) grDevices::dev.off()
  print(paste0(Sys.time()," making annotations"))
  ha = ComplexHeatmap::columnAnnotation(
    classes = labels,
    col = list(classes = c("trait 1" = "red", "trait 2" = "green", "trait 3" = "blue"))
  )
#  mat = log10(mat) * -1
# browser()
#   mat_cr = log10(mat) * -1
#   message(which(is.infinite(mat_cr), arr.ind = TRUE)) #identify infinite if given
#   mat_cr[which(is.infinite(mat_cr), arr.ind = TRUE)] = .Machine$integer.max #convert infinite to max allowed integer

  print(paste0(Sys.time()," making colors"))
  max.Col1 = 0.05
  min.Col1 = min(mat, na.rm=TRUE)
  max.Col2 = 0.05
  min.Col2 = min(mat, na.rm=TRUE)
  max.Col3 = 0.05
  min.Col3 = min(mat, na.rm=TRUE)
  col1 = circlize::colorRamp2(seq(max.Col1, min.Col1, length = 9), RColorBrewer::brewer.pal(9,"OrRd"))
  col2 = circlize::colorRamp2(seq(max.Col2, min.Col2, length = 9), RColorBrewer::brewer.pal(9,"YlGn"))
  col3 = circlize::colorRamp2(seq(max.Col3, min.Col3, length = 9), RColorBrewer::brewer.pal(9,"GnBu"))
  # col1 = colorRamp2(seq(max(mat, na.rm=TRUE), min(mat, na.rm=TRUE), length = 9), brewer.pal(9,"OrRd"))
  # col2 = colorRamp2(seq(max(mat, na.rm=TRUE), min(mat, na.rm=TRUE), length = 9), brewer.pal(9,"GnBu"))
  # col1 = colorRamp2(seq(min(mat_cr, na.rm=TRUE), max(mat_cr, na.rm=TRUE), length = 9), brewer.pal(9,"OrRd"))
  # col2 = colorRamp2(seq(min(mat_cr, na.rm=TRUE), max(mat_cr, na.rm=TRUE), length = 9), brewer.pal(9,"GnBu"))
  print(paste0(Sys.time()," making heatmap", dim(mat)))
  print(dim(mat))
#browser() #check, wheter error message occurs here
  # dim 44401 165
  # Warnung in (function (filename = "Rplot%03d.png", width = 480, height = 480,
  #                       cairo error 'invalid value (typically too big) for the size of the input (surface, pattern, etc.)'
  #                       Warnung: Error in <Anonymous>: unable to start device 'png'
  #                       [No stack trace available]
#   ht = Heatmap(mat, rect_gp = gpar(type = "none"), cluster_rows = dendProbes, cluster_columns = dendTraits,
#                top_annotation = ha,
#                layer_fun = function(j, i, x, y, w, h, fill) {
#                  l = labels[j] == "trait 1"
#                  if(any(l)) grid.rect(x[l], y[l], w[l], h[l], gp = gpar(fill = col1(pindex(mat, i[l], j[l])), col = NA))
#                  l = labels[j] == "trait 2"
#                  if(any(l)) grid.rect(x[l], y[l], w[l], h[l], gp = gpar(fill = col2(pindex(mat, i[l], j[l])), col = NA))
#                  l = labels[j] == "trait 3"
#                  if(any(l)) grid.rect(x[l], y[l], w[l], h[l], gp = gpar(fill = col3(pindex(mat, i[l], j[l])), col = NA))
# #                 grid.text(paste0("p:",sprintf("%.G", pindex(mat, i, j)),"\n",
# #                                  "d:",sprintf("%.G", pindex(matDM, i, j)),"\n",
# #                                  "n:",pindex(matN, i, j)), x, y, gp = gpar(fontsize = 8))
#                }, show_heatmap_legend = FALSE, use_raster = TRUE)
  ht = ComplexHeatmap::Heatmap(mat, rect_gp = grid::gpar(type = "none"), cluster_rows = dendProbes, cluster_columns = dendTraits,
               top_annotation = ha,
               layer_fun = function(j, i, x, y, w, h, fill) {
                 l = labels[j] == "trait 1"
                 if(any(l)) grid::grid.rect(x[l], y[l], w[l], h[l], gp = grid::gpar(fill = col1(ComplexHeatmap::pindex(mat, i[l], j[l])), col = NA))
                 l = labels[j] == "trait 2"
                 if(any(l)) grid::grid.rect(x[l], y[l], w[l], h[l], gp = grid::gpar(fill = col2(ComplexHeatmap::pindex(mat, i[l], j[l])), col = NA))
                 l = labels[j] == "trait 3"
                 if(any(l)) grid::grid.rect(x[l], y[l], w[l], h[l], gp = grid::gpar(fill = col3(ComplexHeatmap::pindex(mat, i[l], j[l])), col = NA))
               },
              show_heatmap_legend = FALSE, use_raster = TRUE)

  sub_heatmap_layer_fun = function(j, i, x, y, w, h, fill) {
    l = labels[j] == "trait 1"
    if(any(l)) grid::grid.rect(x[l], y[l], w[l], h[l], gp = grid::gpar(fill = col1(ComplexHeatmap::pindex(mat, i[l], j[l])), col = NA))
    l = labels[j] == "trait 2"
    if(any(l)) grid::grid.rect(x[l], y[l], w[l], h[l], gp = grid::gpar(fill = col2(ComplexHeatmap::pindex(mat, i[l], j[l])), col = NA))
    l = labels[j] == "trait 3"
    if(any(l)) grid::grid.rect(x[l], y[l], w[l], h[l], gp = grid::gpar(fill = col3(ComplexHeatmap::pindex(mat, i[l], j[l])), col = NA))
    grid::grid.text(paste0("p:",sprintf("%.G", ComplexHeatmap::pindex(mat, i, j)),"\n",
                     "d:",sprintf("%.G", ComplexHeatmap::pindex(matDM, i, j)),"\n",
                     "n:",ComplexHeatmap::pindex(matN, i, j)), x, y, gp = grid::gpar(fontsize = 8))
  }

  print(paste0(Sys.time()," making legend"))
  lgd = list(
    ComplexHeatmap::Legend(title = "trait 1", col_fun = col1, at = c(max.Col1, min.Col1), labels = c(max.Col1, paste0(extractMantissaExponent(min.Col1)$exponent))),
#    Legend(title = "trait 2", col_fun = col2, at = c(max.Col2, min.Col2), labels = c(max.Col2, format(min.Col2, scientific = TRUE)))
    ComplexHeatmap::Legend(title = "trait 2", col_fun = col2, at = c(max.Col2, min.Col2), labels = c(max.Col2, paste0(extractMantissaExponent(min.Col2)$exponent))),
    ComplexHeatmap::Legend(title = "trait 3", col_fun = col3, at = c(max.Col3, min.Col3), labels = c(max.Col3, paste0(extractMantissaExponent(min.Col3)$exponent)))
  )
  if(grDevices::dev.cur() > 1) grDevices::dev.off()
  grDevices::pdf(NULL)
#  message(paste0(Sys.time()," gc()"))
#  gc()
  print(paste0(Sys.time()," plotting heatmap"))
  #Error in Cairo: Failed to create Cairo backend!
#  ht_shiny(ht)
#  ht = ht_shiny(ht, heatmap_legend_list = lgd) #does not work
  tryCatch({
#with huge heatmaps, the following error occurs:
#Error in Cairo: Failed to create Cairo backend!
#  ht = draw(ht, heatmap_legend_list = lgd)
    ht = ComplexHeatmap::draw(ht, annotation_legend_list = lgd)
#    heatmap_legend_list = lgd,

  }, error=function(err){
    print(paste0(Sys.time()," Error: unable to plot HM. ", err$message))
  });

  grDevices::dev.off()
  l = list()
  l$combinedHMP_VAL = ht
  l$layer_fun = sub_heatmap_layer_fun
  print(paste0(Sys.time()," end plotting heatmap"))
  return (l) #return(ht)
}

# layer_function <- function(j, i, x, y, w, h, fill) {
# browser()
#   l = labels[j] == "trait 1"
#   if(any(l)) grid.rect(x[l], y[l], w[l], h[l], gp = gpar(fill = col1(pindex(mat, i[l], j[l])), col = NA))
#   l = labels[j] == "trait 2"
#   if(any(l)) grid.rect(x[l], y[l], w[l], h[l], gp = gpar(fill = col2(pindex(mat, i[l], j[l])), col = NA))
#   l = labels[j] == "trait 3"
#   if(any(l)) grid.rect(x[l], y[l], w[l], h[l], gp = gpar(fill = col3(pindex(mat, i[l], j[l])), col = NA))
#
#   grid.text(sprintf("%.G", pindex(mat, i, j)), x, y, gp = gpar(fontsize = 8))
# }
#extract dendrograms as described here: https://github.com/jokergoo/ComplexHeatmap/issues/136
#  r.dend <- row_dend(HM)  #If needed, extract row dendrogram
# rcl.list <- row_order(HM)  #Extract clusters (output is a list)
