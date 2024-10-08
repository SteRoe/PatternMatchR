#' reduceTraitData
#' reduces data in structure combinedDFP_Val_Labels by p-value range and minimum n; entire CpGs will be removed
#' @param combinedDFP_Val_Labels data structure which should become p_val-reduced
#' @param minP_Val minimum p-value to retain
#' @param maxP_Val maximum p-value to retain
#' @param minN minimum n for model in regression result
#' @param debugMode return smaller data structure (only session$userData$sessionVariables$debugNumber records) for faster debugging
#' @return result list()
getPReducedTraitData <- function(session, combinedDFP_Val_Labels, minP_Val, maxP_Val, minDM, maxDM, minN, maxN, debugMode) {
  id <- shiny::showNotification("Getting p reduced trait data...", duration = NULL, closeButton = FALSE)
  on.exit(shiny::removeNotification(id), add = TRUE)
  if (maxN < 1) {
    base::print(base::paste0(sysTimePID(), "Warning: maxN < 1. Please check your data.")) #that should not be the case, please check data!
    browser()
  }
  base::tryCatch(
    {
      base::print(base::paste0(sysTimePID(), " start pReduceTraitData()."))
      if (is.valid(combinedDFP_Val_Labels)) {
        result <- combinedDFP_Val_Labels
        dfP_Val <- result$dfP_Val
        dfDM <- result$dfDM
        dflogFC <- result$dflogFC
        dfN <- result$dfN
        base::print(base::paste0(sysTimePID(),
                                 " result matrix before reduce has N(row) probes=",
                                 base::nrow(dfP_Val),
                                 " and n(col) traits=",
                                 base::ncol(dfP_Val), "."))
        LabelsDF1 <- result$labelsDF1
        LabelsDF2 <- result$labelsDF2
        LabelsDF3 <- result$labelsDF3
        mergedOriginDF <- result$mergedOriginDF
        mergedColnames <- result$mergedColnames
        mergedOriginalColnames <- result$mergedOriginalColnames
        mergedOriginTrait <- result$mergedOriginTrait
        mergedDFList <- result$mergedDFList
        traitID <- result$traitID
        # omit p_values from dfs - max PVal
        #if we obtain 0 hits later, check maxP_Val and minP_Val
        cgsToRetainMaxP <- dfP_Val < maxP_Val
        #cgsToRetainMaxP <- unique(rownames(which(cgsToRetainMaxP == TRUE, arr.ind = TRUE)))
        #cgsToRetainMaxP <- unique(rownames(cgsToRetainMaxP == TRUE))
        cgsToRetainMaxP <- names(which((rowSums(cgsToRetainMaxP) > 0L) == TRUE))

        cgsToRetainMinP <- dfP_Val > minP_Val
        #cgsToRetainMinP <- unique(rownames(which(cgsToRetainMinP == TRUE, arr.ind = TRUE)))
        #cgsToRetainMinP <- unique(rownames(cgsToRetainMinP == TRUE))
        cgsToRetainMinP <- names(which((rowSums(cgsToRetainMinP) > 0L) == TRUE))
#browser() #check
        if (base::exists("cgsToRetainMaxP") && base::exists("cgsToRetainMinP")) {
          cgsToRetainP <- base::intersect(cgsToRetainMaxP, cgsToRetainMinP)
        }
        else if (base::exists("cgsToRetainMaxP")) {
          cgsToRetainP <- cgsToRetainMaxP
        }
        else if (base::exists("cgsToRetainMinP")) {
          cgsToRetainP <- cgsToRetainMinP
        }
        if (length(cgsToRetainP) == 0) {
          base::print(base::paste0(sysTimePID(), "Warning: length(cgsToRetainP) == 0. Please check your data.")) #that should not be the case, please check data!
          browser()
        }
        #take only DM outside slider defined range in cgsToRetainDM
        cgsToRetainMaxDM <- dfDM < maxDM
        #cgsToRetainMaxDM <- unique(rownames(which(cgsToRetainMaxDM == TRUE, arr.ind = TRUE)))
        #cgsToRetainMaxDM <- unique(rownames(cgsToRetainMaxDM == TRUE))
        cgsToRetainMaxDM <- names(which((rowSums(cgsToRetainMaxDM) > 0L) == TRUE))
        cgsToRetainMinDM <- dfDM > minDM
        #cgsToRetainMinDM <- unique(rownames(which(cgsToRetainMinDM == TRUE, arr.ind = TRUE)))
        #cgsToRetainMinDM <- unique(rownames(cgsToRetainMinDM == TRUE))
        cgsToRetainMinDM <- names(which((rowSums(cgsToRetainMinDM) > 0L) == TRUE))
        #cgsToRetainDM <- unique(c(cgsToRetainMaxDM, cgsToRetainMinDM))
        if (base::exists("cgsToRetainMaxDM") && base::exists("cgsToRetainMinDM")) {
          cgsToRetainDM <- base::intersect(cgsToRetainMaxDM, cgsToRetainMinDM)
        }
        else if (base::exists("cgsToRetainMaxDM")) {
          cgsToRetainDM <- cgsToRetainMaxDM
        }
        else if (base::exists("cgsToRetainMinDM")) {
          cgsToRetainDM <- cgsToRetainMinDM
        }
        if (length(cgsToRetainDM) == 0) {
          base::print(base::paste0(sysTimePID(), "Warning: length(cgsToRetainDM) == 0. Please check your data.")) #that should not be the case, please check data!
          browser()
        }
        cgsToRetainMaxN <- dfN < maxN
        #cgsToRetainMaxN <- unique(rownames(which(cgsToRetainMaxN == TRUE, arr.ind = TRUE)))
        #cgsToRetainMaxN <- unique(rownames(cgsToRetainMaxN == TRUE))
        cgsToRetainMaxN <- names(which((rowSums(cgsToRetainMaxN) > 0L) == TRUE))
        cgsToRetainMinN <- dfN > minN
        #cgsToRetainMinN <- unique(rownames(which(cgsToRetainMinN == TRUE, arr.ind = TRUE)))
        #cgsToRetainMinN <- unique(rownames(cgsToRetainMinN == TRUE))
        cgsToRetainMinN <- names(which((rowSums(cgsToRetainMinN) > 0L) == TRUE))
        if (base::exists("cgsToRetainMaxN") && base::exists("cgsToRetainMinN")) {
          cgsToRetainN <- base::intersect(cgsToRetainMaxN, cgsToRetainMinN)
        }
        else if (base::exists("cgsToRetainMaxN")) {
          cgsToRetainN <- cgsToRetainMaxN
        }
        else if (base::exists("cgsToRetainMinN")) {
          cgsToRetainN <- cgsToRetainMinN
        }
        #intersect of all three cgsToRetain here:
        cgsToRetain <- intersect(cgsToRetainP, cgsToRetainDM)
        cgsToRetain <- intersect(cgsToRetain, cgsToRetainN)
        if (!base::exists("cgsToRetain") || length(cgsToRetain) == 0) {
          base::print(base::paste0(sysTimePID(), "Warning: length(cgsToRetain) == 0. Please check your data.")) #that should not be the case, please check data!
          browser()
          base::print(base::paste0(sysTimePID(), " max p-val border too low: ",
                                   maxP_Val, "; no remaining CpG."))
        }
        dfP_Val <- dfP_Val[cgsToRetain, ]
        dfDM <- dfDM[cgsToRetain, ]
        dflogFC <- dflogFC[cgsToRetain, ]
        dfN <- dfN[cgsToRetain, ]
        base::print(base::paste0(sysTimePID(),
                                 " result matrix after reduce has N(row) probes=",
                                 base::nrow(dfP_Val),
                                 " and n(col) traits=", base::ncol(dfP_Val),
                                 "."))
        #        dfP_Val[dfP_Val > 0.05] <- NA # 1
        base::print(base::class(dfP_Val))
        base::print(Cstack_info())
        if (base::nrow(dfP_Val) >= 5) {
          if (debugMode == TRUE && base::nrow(dfP_Val) > session$userData$sessionVariables$debugNumber) {
            base::print(base::paste0(sysTimePID(), " debug mode n probes=session$userData$sessionVariables$debugNumber"))
            #dfP_Val <- dfP_Val[1:session$userData$sessionVariables$debugNumber, ]
            dfP_Val <- head(dfP_Val, session$userData$sessionVariables$debugNumber)
            dfDM <- head(dfDM, session$userData$sessionVariables$debugNumber)
            dflogFC <- head(dflogFC, session$userData$sessionVariables$debugNumber)
            dfN <- head(dfN, session$userData$sessionVariables$debugNumber)
          }
          base::print(base::paste0(
            sysTimePID(),
            " set infinite p-values to .Machine$integer.max."
          ))
          indexdfP_Val_infinite <- is.infinite(as.matrix(dfP_Val))
          dfP_Val[indexdfP_Val_infinite] <- .Machine$integer.max
          if (!base::is.data.frame(dfDM)) {
            dfDM <- base::as.data.frame(dfDM)
          }
          base::print(base::paste0(sysTimePID(), " shortening dfDM"))
          dfDM <- dfDM[rownames(dfP_Val), colnames(dfP_Val)]
          if (!base::is.data.frame(dfN)) {
            dfN <- base::as.data.frame(dfN)
          }
          base::print(base::paste0(sysTimePID(), " shortening dflogFC"))
          dflogFC <- dflogFC[rownames(dfP_Val), colnames(dfP_Val)]
          if (!base::is.data.frame(dflogFC)) {
            dflogFC <- base::as.data.frame(dflogFC)
          }
          base::print(base::paste0(sysTimePID(), " shortening dfN"))
          dfN <- dfN[rownames(dfP_Val), colnames(dfP_Val)]
          combinedDFP_Val_Labels <- base::list(dfP_Val = NULL, dfDM = NULL, dflogFC = NULL,
                                               dfN = NULL, labelsDF1 = NULL,
                                               labelsDF2 = NULL, labelsDF3 = NULL,
                                               mergedOriginDF = NULL, mergedColnames = NULL, mergedOriginalColnames = NULL,
                                               mergedOriginTrait = NULL, mergedDFList = NULL, traitID = NULL)
          combinedDFP_Val_Labels$dfP_Val <- dfP_Val
          combinedDFP_Val_Labels$dfDM <- dfDM
          combinedDFP_Val_Labels$dflogFC <- dflogFC
          combinedDFP_Val_Labels$dfN <- dfN

          combinedDFP_Val_Labels$labelsDF1 <- LabelsDF1
          combinedDFP_Val_Labels$labelsDF2 <- LabelsDF2
          combinedDFP_Val_Labels$labelsDF3 <- LabelsDF3
          combinedDFP_Val_Labels$mergedOriginDF <- mergedOriginDF
          combinedDFP_Val_Labels$mergedColnames <- mergedColnames
          combinedDFP_Val_Labels$mergedOriginalColnames <- mergedOriginalColnames
          combinedDFP_Val_Labels$mergedOriginTrait <- mergedOriginTrait
          combinedDFP_Val_Labels$mergedDFList <- mergedDFList #(consists of mergedDFList$PHENODF and mergedDFList$PHENOFileName for each source file)
          combinedDFP_Val_Labels$traitID <- traitID
        }
        else {
          base::message(base::paste0(sysTimePID(), " less than 5 probes remained. Please check your maximum and minimum p-val settings."))
        }
      }
    },
    error = function(e) {
      message("An error occurred in getPReducedTraitData():\n", e)
    },
    warning = function(w) {
      message("A warning occurred in getPReducedTraitData():\n", w)
    },
    finally = {
      #base::print(base::paste0(sysTimePID(), " size of reduced data.frame: ", dim(combinedDFP_Val_Labels$dfP_Val), " ."))
      base::print(base::paste0(sysTimePID(), " size of reduced data.frame: nrow (CpG)=",
                               nrow(combinedDFP_Val_Labels$dfP_Val), "; ncol (trait)=",
                               ncol(combinedDFP_Val_Labels$dfP_Val), " ."))
      base::print(base::paste0(sysTimePID(), " pReduceTraitData() finished."))
      return(combinedDFP_Val_Labels)
    }
  )
}

#' updateTxtpReduceOut
#' generates summary text after p-value reduction
#' @param pReducedcombinedDFP_Val_Labels data structure with trait reduced results
#' @return text
#' examples updateTxtpReduceOut(pReducedcombinedDFP_Val_Labels)
updateTxtpReduceOut <- function(pReducedcombinedDFP_Val_Labels) {
  base::tryCatch(
    {
      result <- NULL
      if (is.valid(pReducedcombinedDFP_Val_Labels)) {
        result <- base::paste0("p reduce successful. result table is: nrow (CpG): ",
                               nrow(pReducedcombinedDFP_Val_Labels$dfP_Val),
                               "; ncol (trait): ", ncol(pReducedcombinedDFP_Val_Labels$dfP_Val))
      }
      else {
        base::message(base::paste0(sysTimePID(), " is.valid(pReducedcombinedDFP_Val_Labels) == FALSE."))
      }
    },
    error = function(e) {
      base::message("An error occurred in updateTxtpReduceOut():\n", e)
    },

    warning = function(w) {
      base::message("A warning occurred in updateTxtpReduceOut():\n", w)
    },
    finally = {
      return(shiny::HTML(result))
    }
  )
}
