ReduceData_UI <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
  shiny::fluidRow(
    shiny::column(
      width = 4,
      shinyjs::disabled(shiny::sliderInput(
        ns("sldP_Val"),
        "maximum (left slider) and minimum (right slider) p-val, 5e-x",
        min = 0,
        max = 0,
        step = 0, #-1, #1
        value = c(0, 0)
      ))
    ),
    shiny::column(
      width = 4,
      shinyjs::disabled(shiny::sliderInput(
        ns("sldDM"),
        "minimum (left slider) and maximum (right slider) delta methylation",
        min = 0,
        max = 0,
        step = 0, #.01,
        value = c(0, 0)
      ))
    ),
    shiny::column(
      width = 4,
      shinyjs::disabled(shiny::sliderInput(
        ns("sldN"),
        "minimum (left slider) and maximum (right slider) n",
        min = 0,
        max = 0,
        step = 0,
        value = c(0, 0)
      ))
    )
  ),
  shinyjs::disabled(shiny::actionButton(ns("btnReduce"), label = "Step 3: Reduce data (omit CpGs) by applying thresholds for p-value, DM or n limit")),
  shiny::verbatimTextOutput(ns("txtPReduceOut"), placeholder = TRUE)
  ) #end tagList
}

ReduceData_SERVER <- function(id, session) {
  shiny::moduleServer(id, function(input, output, session) {
    base::tryCatch({
      shiny::observe({
        if (!is.valid(session$userData$sessionVariables$combinedDataStructure()$combinedDFP_Val_Labels$dfP_Val)) {
          shinyjs::disable("Count Borders")
          shinyjs::disable("btnCountProbesP_ValParallel")
          shinyjs::disable("btnCountProbesDeltaMethParallel")
          shinyjs::disable("btnCountProbesNParallel")

          shinyjs::disable("Reduce Data")
          shinyjs::disable("sldP_Val")
          shinyjs::disable("sldDM")
          shinyjs::disable("sldN")
          shinyjs::disable("btnReduce")
        }
        else {
          shinyjs::enable("Count Borders")
          shinyjs::enable("btnCountProbesP_ValParallel")
          shinyjs::enable("btnCountProbesDeltaMethParallel")
          shinyjs::enable("btnCountProbesNParallel")

          shinyjs::enable("Reduce Data")
          shinyjs::enable("sldP_Val")
          shinyjs::enable("sldDM")
          shinyjs::enable("sldN")
          shinyjs::enable("btnReduce")
        }
      })

      shiny::observe({
        result <- session$userData$sessionVariables$combinedData()
        if (is.valid(result)) {
          updateReduceDataSliders(session, result)
        }
      })

      shiny::observeEvent(input$btnReduce,
        ignoreInit = TRUE,
        {
          base::tryCatch(
            {
              if (is.valid(session$userData$sessionVariables$combinedDataStructure())) {
                base::print(base::paste0(sysTimePID(), " Step 3: start reducing data by p-value."))
                minP_Val <- 5 * 10^base::as.integer(input$sldP_Val[1]) #minP_Val <- 5 * 10^-base::as.integer(input$sldP_Val[2])
                maxP_Val <- 5 * 10^base::as.integer(input$sldP_Val[2]) #maxP_Val <- 5 * 10^-base::as.integer(input$sldP_Val[1])
                if (maxP_Val < minP_Val) { #exchange, if in wrong order
                  t <- minP_Val
                  minP_Val <- maxP_Val
                  maxP_Val <- t
                  browser() #should not happen
                }
                minDM <- input$sldDM[1]
                maxDM <- input$sldDM[2]
                minN <- base::as.integer(input$sldN[1])
                maxN <- base::as.integer(input$sldN[2])
                combinedDFP_Val_Labels <- session$userData$sessionVariables$combinedDataStructure()$combinedDFP_Val_Labels
                #if (is.valid(combinedDFP_Val_Labels)) {
                if (minP_Val != maxP_Val && minDM != maxDM && minN != maxN) {
                  result <- getPReducedTraitData(session = session,
                                                 combinedDFP_Val_Labels =
                                                   combinedDFP_Val_Labels,
                                                 minP_Val = minP_Val,
                                                 maxP_Val = maxP_Val,
                                                 minDM = minDM,
                                                 maxDM = maxDM,
                                                 minN = minN,
                                                 maxN = maxN,
                                                 debugMode = session$userData$config$debugMode)
                }
                else {
                  result <- NULL
                  base::print(base::paste0(sysTimePID(), " minP_Val == maxP_Val && minDM == maxDM && minN == maxN."))
                }
              }
              else {
                result <- NULL
                base::print(base::paste0(sysTimePID(), " is.valid(session$userData$sessionVariables$combinedDataStructure()) == FALSE."))
              }
            },
            error = function(e) {
              if (attributes(e)$class[1] != "shiny.silent.error") {
                base::message("An error occurred in shiny::observeEvent(input$btnReduce):\n", e)
              }
            },
            warning = function(w) {
              base::message("A warning occurred in shiny::observeEvent(input$btnReduce):\n", w)
            },
            finally = {
              session$userData$sessionVariables$pReducedData(result)
#              session$userData$sessionVariables$pReducedDataStructure(result)
              base::print(base::paste0(sysTimePID(), " finished reducing data by p-value."))
            }
          )
        },
        ignoreNULL = FALSE
      )

      output$txtPReduceOut <- shiny::reactive({
        return(updateTxtpReduceOut(session$userData$sessionVariables$pReducedDataStructure()$combinedDFP_Val_Labels))
      })

    },
    error = function(e) {
      if (attributes(e)$class[1] != "shiny.silent.error") {
        base::message("An error occurred in moduleServer in ReduceData_SERVER:\n", e)
        browser() #should not happen
      }
    },
    warning = function(w) {
      base::message("An error occurred in moduleServer in ReduceData_SERVER:\n", w)
      browser() #should not happen
    },
    finally = {
      base::print(base::paste0(sysTimePID(), " finished moduleServer in ReduceData_SERVER"))
    })
  } #end moduleServer
)} #end Clustering_TraitsSERVER

updateReduceDataSliders <- function(session, combinedDFP_Val_Labels) {
  DF <- combinedDFP_Val_Labels$dfP_Val
  DF <- as.matrix(DF)
  minP <- base::apply(DF, 2, FUN = function(x) {base::min(x[x > 0], na.rm = TRUE)})
  minP <- base::min(minP)
  minP <- extractMantissaExponent(minP)$exponent #base::round(extractMantissaExponent(minP)$exponent, 5)
  maxP <- base::apply(DF, 2, FUN = function(x) {base::max(x[x > 0], na.rm = TRUE)})
  maxP <- base::max(maxP)
  maxP <- extractMantissaExponent(maxP)$exponent #base::round(extractMantissaExponent(maxP)$exponent, 5)
  shiny::updateSliderInput(session = session, inputId = "sldP_Val", min = minP, max = maxP, value = c(minP, maxP))
  DF <- combinedDFP_Val_Labels$dfDM
  DF <- as.matrix(DF)
  minDM <- base::apply(DF, 2, FUN = function(x) {base::min(x[x > 0], na.rm = TRUE)})
  minDM <- base::min(minDM)
  if (minDM < 0) {
    base::message(base::paste0(sysTimePID(), " Warning: minDM < 0. Please check your data.")) #that should not be the case, please check data!
    minDM <- 0
  }
  maxDM <- base::apply(DF, 2, FUN = function(x) {base::max(x[x > 0], na.rm = TRUE)})
  maxDM <- base::max(maxDM)
  if (maxDM > 1) {
    base::message(base::paste0(sysTimePID(), "Warning: maxDM > 1. Please check your data.")) #that should not be the case, please check data!
    maxDM <- 1
  }
  #shiny::updateSliderInput(session = session, inputId = "sldDM", min = minDM, max = maxDM, value = c(minDM, maxDM), step = NULL)
  shiny::updateSliderInput(session = session, inputId = "sldDM", min = minDM, max = maxDM, value = c(minDM, maxDM), step = 0.001)
  DF <- combinedDFP_Val_Labels$dfN
  DF <- as.matrix(DF)
  minN <- base::apply(DF, 2, FUN = function(x) {base::min(as.integer(x[x > 0]), na.rm = TRUE)})
  minN <- base::min(minN)
  if (minN < 1) {
    base::message(base::paste0(sysTimePID(), "Warning: minN < 1. Please check your data.")) #that should not be the case, please check data!
    minN <- 1
  }
  if (minN != as.integer(minN)) {
    base::message(base::paste0(sysTimePID(), "Warning: minN != as.integer(minN). Please check your data.")) #that should not be the case, please check data!
    minN <- as.integer(minN)
  }
  maxN <- base::apply(DF, 2, FUN = function(x) {base::max(as.integer(x[x > 0]), na.rm = TRUE)})
  maxN <- base::max(maxN)
  if (maxN != as.integer(maxN)) {
    base::message(base::paste0(sysTimePID(), "Warning: maxN != as.integer(maxN). Please check your data.")) #that should not be the case, please check data!
    maxN <- as.integer(maxN)
  }
  if (maxN < 1) {
    base::message(base::paste0(sysTimePID(), "Warning: maxN < 1. Please check your data.")) #that should not be the case, please check data!
    browser() #should not happen
  }
  shiny::updateSliderInput(session = session, inputId = "sldN", min = minN, max = maxN, value = c(minN, maxN))
}

#' getPReducedTraitData
#' reduces data in structure combinedDFP_Val_Labels by p-value range and minimum n; entire CpGs will be removed
#' @param combinedDFP_Val_Labels data structure which should become p_val-reduced
#' @param minP_Val minimum p-value to retain
#' @param maxP_Val maximum p-value to retain
#' @param minN minimum n for model in regression result
#' @param debugMode return smaller data structure (only session$userData$sessionVariables$debugNumber records) for faster debugging
#' @return result list()
getPReducedTraitData <- function(session, combinedDFP_Val_Labels, minP_Val, maxP_Val, minDM, maxDM, minN, maxN, debugMode) {
  shinyId <- shiny::showNotification("Getting p reduced trait data...", duration = NULL, closeButton = FALSE)
  on.exit(shiny::removeNotification(shinyId), add = TRUE)
  if (maxN < 1) {
    base::print(base::paste0(sysTimePID(), "Warning: maxN < 1. Please check your data.")) #that should not be the case, please check data!
    browser() #should not happen
  }
  base::tryCatch(
    {
      base::print(base::paste0(sysTimePID(), " start pReduceTraitData()."))
      if (is.valid(combinedDFP_Val_Labels)) {
        result <- combinedDFP_Val_Labels
        dfP_Val <- result$dfP_Val
        dfDM <- result$dfDM
        dfLogFC <- result$dfLogFC
        dfN <- result$dfN
        base::print(base::paste0(sysTimePID(),
                                 " result matrix before reduce has N(row) probes=",
                                 base::nrow(dfP_Val),
                                 " and n(col) traits=",
                                 base::ncol(dfP_Val), "."))
        # LabelsDF1 <- result$labelsDF1
        # LabelsDF2 <- result$labelsDF2
        # LabelsDF3 <- result$labelsDF3
        mergedOriginDF <- result$mergedOriginDF
        mergedColnames <- result$mergedColnames
        mergedOriginalColnames <- result$mergedOriginalColnames
        mergedOriginTrait <- result$mergedOriginTrait
        mergedDFList <- result$mergedDFList
        traitID <- result$traitID
        # omit p_values from dfs - max PVal
        #if we obtain 0 hits later, check maxP_Val and minP_Val
        cgsToRetainMaxP <- dfP_Val < maxP_Val
        cgsToRetainMaxP <- names(which((rowSums(cgsToRetainMaxP) > 0L) == TRUE))

        cgsToRetainMinP <- dfP_Val > minP_Val
        cgsToRetainMinP <- names(which((rowSums(cgsToRetainMinP) > 0L) == TRUE))
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
          browser() #should not happen
        }
        #take only DM outside slider defined range in cgsToRetainDM
        cgsToRetainMaxDM <- dfDM < maxDM
        cgsToRetainMaxDM <- names(which((rowSums(cgsToRetainMaxDM) > 0L) == TRUE))
        cgsToRetainMinDM <- dfDM > minDM
        cgsToRetainMinDM <- names(which((rowSums(cgsToRetainMinDM) > 0L) == TRUE))
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
          browser() #should not happen
        }
        cgsToRetainMaxN <- dfN < maxN
        cgsToRetainMaxN <- names(which((rowSums(cgsToRetainMaxN) > 0L) == TRUE))
        cgsToRetainMinN <- dfN > minN
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
          browser() #should not happen
          base::print(base::paste0(sysTimePID(), " max p-val border too low: ",
                                   maxP_Val, "; no remaining CpG."))
        }
        dfP_Val <- dfP_Val[cgsToRetain, ]
        dfDM <- dfDM[cgsToRetain, ]
        dfLogFC <- dfLogFC[cgsToRetain, ]
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
            dfP_Val <- head(dfP_Val, session$userData$sessionVariables$debugNumber)
            dfDM <- head(dfDM, session$userData$sessionVariables$debugNumber)
            dfLogFC <- head(dfLogFC, session$userData$sessionVariables$debugNumber)
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
          base::print(base::paste0(sysTimePID(), " shortening dfLogFC"))
          dfLogFC <- dfLogFC[rownames(dfP_Val), colnames(dfP_Val)]
          if (!base::is.data.frame(dfLogFC)) {
            dfLogFC <- base::as.data.frame(dfLogFC)
          }
          base::print(base::paste0(sysTimePID(), " shortening dfN"))
          dfN <- dfN[rownames(dfP_Val), colnames(dfP_Val)]
          combinedDFP_Val_Labels <- base::list(dfP_Val = NULL, dfDM = NULL, dfLogFC = NULL,
                                               dfN = NULL, labelsDF1 = NULL,
                                               labelsDF2 = NULL, labelsDF3 = NULL,
                                               mergedOriginDF = NULL, mergedColnames = NULL, mergedOriginalColnames = NULL,
                                               mergedOriginTrait = NULL, mergedDFList = NULL, traitID = NULL)
          combinedDFP_Val_Labels$dfP_Val <- dfP_Val
          combinedDFP_Val_Labels$dfDM <- dfDM
          combinedDFP_Val_Labels$dfLogFC <- dfLogFC
          combinedDFP_Val_Labels$dfN <- dfN

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
      if (attributes(e)$class[1] != "shiny.silent.error") {
        message("An error occurred in getPReducedTraitData():\n", e)
      }
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
      if (attributes(e)$class[1] != "shiny.silent.error") {
        base::message("An error occurred in updateTxtpReduceOut():\n", e)
      }
    },

    warning = function(w) {
      base::message("A warning occurred in updateTxtpReduceOut():\n", w)
    },
    finally = {
      return(shiny::HTML(result))
      #return(shiny::renderPrint(result))
    }
  )
}
