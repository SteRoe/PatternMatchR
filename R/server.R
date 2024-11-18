server <- function(input, output, session) {
  waiter::waiter_hide() # will hide *on_load waiter
  # waiter::waiter_on_busy(
  #   html = spin_1(),
  #   color = "#333e48",
  #   logo = "",
  #   image = "",
  #   fadeout = FALSE
  # )
  reactlog::reactlog_enable()
  pid <- Sys.getpid()
  hostname <- Sys.info()["nodename"]
  output$Sys.PID <- shiny::renderText(base::paste0(hostname, ": ", pid))
  packageWd <- getwd()
  session$userData$packageWd <- packageWd
  base::print(paste0(sysTimePID(), " getwd: ", packageWd))
  base::print(paste0(sysTimePID(), " loading configuration."))
  configFileLocation <- system.file("config.yml", package = "PatternMatchR", mustWork = TRUE)
  session$userData$config <- config::get(file = configFileLocation)
  volumes <- c(Home = paste0(getwd(), sub(".", "", session$userData$config$workDir)),
               shinyFiles::getVolumes()())
  shinyFiles::shinyFileChoose(input, "file", roots = volumes, session = session)
  shinyFiles::shinyFileSave(input, "save", roots = volumes, session = session, restrictions = system.file(package = "base"))
  if (is.valid(session$userData$config$knownCpGs)) {
    #load CpG from txt file to search input
    knownCpG <- paste(unlist(data.table::fread(file = session$userData$config$knownCpGs, header = FALSE)), collapse = " ")
    shiny::updateTextInput(session, inputId = "txtSearchFullCpG", value = knownCpG)
  }
  if (is.valid(session$userData$config$knownTraits)) {
    #load Traits from txt file to search input
    knownTrait <- paste(unlist(data.table::fread(file = session$userData$config$knownTraits, header = FALSE)), collapse = " ")
    shiny::updateTextInput(session, inputId = "txtSearchFullTrait", value = knownTrait)
  }
  #base::options(spam.force64 = TRUE)

  if (session$userData$config$debugMode == TRUE) {
    shiny::updateCheckboxInput(session, "chkDebug", value = TRUE)
  }
  else {
    shiny::updateCheckboxInput(session, "chkDebug", value = FALSE)
  }

  #define sessionVariables here
  print(paste0(sysTimePID(), " defining session variables."))
  session$userData$sessionVariables <-
    shiny::reactiveValues(
      P_ValMaxBorder = double(),
      P_ValMinBorder = double(),
      MaxProbes = integer(),
      numberVariables = integer(),
      # callCounter = integer(),
      # debugNumber = integer(),
      # packageWd = character(),
      selected_row_labels = list(),
      selected_column_labels = list()
    )
  session <- loadObjects(session)

  session$userData$sessionVariables$callCounter <- 0
  session$userData$sessionVariables$debugNumber <- 1000
  session$userData$packageWd <- packageWd

  result <- loadDirLists(session = session, input = input, output = output)
  dfdD1 <- result$dfdD1
  dfdD2 <- result$dfdD2
  dfdD3 <- result$dfdD3

  session$userData$sessionVariables$resultDFListTrait1 <- shiny::reactiveVal(value = NULL, label = "resultDFListTrait1")
  session$userData$sessionVariables$resultDFListTrait2 <- shiny::reactiveVal(value = NULL, label = "resultDFListTrait2")
  session$userData$sessionVariables$resultDFListTrait3 <- shiny::reactiveVal(value = NULL, label = "resultDFListTrait3")

  session$userData$sessionVariables$combinedData <- shiny::reactiveVal(value = NULL, label = "combinedData")
  session$userData$sessionVariables$pReducedData <- shiny::reactiveVal(value = NULL, label = "pReducedData")
  session$userData$sessionVariables$traitReducedDataPVal <- shiny::reactiveVal(value = NULL, label = "traitReducedDataPVal")
  session$userData$sessionVariables$traitReducedDataLogFC <- shiny::reactiveVal(value = NULL, label = "traitReducedDataLogFC")
  #  session$userData$sessionVariables$combinedDFP_Val_Labels <- shiny::reactiveVal(value = NULL, label = "combinedDFP_Val_Labels")

  # session$userData$sessionVariables$clustResTraits <- shiny::reactiveVal(value = NULL, label = "clustResTraits")

  # session$userData$sessionVariables$traitReducedmatP_Val <- shiny::reactiveVal(value = NULL, label = "traitReducedmatP_Val")

  #new (4.2024) basic data structure:
  session$userData$sessionVariables$generalDataStructure <- shiny::reactiveVal(value = NULL, label = "generalDataStructure")

  session$userData$sessionVariables$distancesBelowThreshold <- shiny::reactiveVal(value = NULL, label = "distancesBelowThreshold")

  session$userData$sessionVariables$distNeigboursProbes10000 <- shiny::reactiveVal(value = NULL, label = "distNeigboursProbes10000")
  session$userData$sessionVariables$distNeigboursProbes1000 <- shiny::reactiveVal(value = NULL, label = "distNeigboursProbes1000")
  session$userData$sessionVariables$distNeigboursProbes100 <- shiny::reactiveVal(value = NULL, label = "distNeigboursProbes100")
  #session$userData$sessionVariables$distNeigboursProbes10 <- shiny::reactiveVal(value = NULL, label = "distNeigboursProbes10")

  session$userData$sessionVariables$selectedOriginalData <- shiny::reactiveVal(value = NULL, label = "selectedOriginalData")
  session$userData$sessionVariables$OriginalDataTraits <- shiny::reactiveVal(value = NULL, label = "OriginalDataTraits")
  session$userData$sessionVariables$selectedOriginalDataTraits <- shiny::reactiveVal(value = NULL, label = "selectedOriginalDataTraits")
  session$userData$sessionVariables$selectedOriginalDataProbes <- shiny::reactiveVal(value = NULL, label = "selectedOriginalDataProbes")
  session$userData$sessionVariables$selectedAnnotation <- shiny::reactiveVal(value = NULL, label = "selectedAnnotation")
  session$userData$sessionVariables$markingVar <- shiny::reactiveVal(value = NULL, label = "markingVar")

  session$userData$sessionVariables$selectedKey <- shiny::reactiveVal(value = NULL, label = "selectedKey")
  session$userData$sessionVariables$selectedTrait <- shiny::reactiveVal(value = NULL, label = "selectedTrait")
  session$userData$sessionVariables$selectedTraitID <- shiny::reactiveVal(value = NULL, label = "selectedTraitID")
  session$userData$sessionVariables$selectedProbe <- shiny::reactiveVal(value = NULL, label = "selectedProbe")

  shiny::updateSliderInput(session = session, inputId = "sldP_Val", min = 0, max = 0, value = c(0, 0))
  shiny::updateSliderInput(session = session, inputId = "sldDM", min = 0, max = 0, value = c(0, 0))
  shiny::updateSliderInput(session = session, inputId = "sldN", min = 0, max = 0, value = c(0, 0))

  #base::print(paste0(sysTimePID(), " starting application."))
  base::message(paste0(sysTimePID(), " starting application."))

  shinyjs::toggleClass("colRed", "red")
  shinyjs::toggleClass("colGreen", "green")
  shinyjs::toggleClass("colBlue", "blue")

################################################################################
  Clustering_SERVER("P_Val", session$userData$sessionVariables$pReducedDataStructurePVal, session$userData$sessionVariables$traitReducedDataStructurePVal, session)
  #Clustering_SERVER("P_Val", session$userData$sessionVariables$pReducedDataStructurePVal, session$userData$sessionVariables$traitReducedDataStructurePVal, session$userData$sessionVariables$distanceMultipliedTraitReducedDataStructurePVal, session)

  Clustering_SERVER("LogFC", session$userData$sessionVariables$pReducedDataStructureLogFC, session$userData$sessionVariables$traitReducedDataStructureLogFC, session)
  #Clustering_SERVER("LogFC", session$userData$sessionVariables$pReducedDataStructureLogFC, session$userData$sessionVariables$traitReducedDataStructureLogFC, session$userData$sessionVariables$distanceMultipliedTraitReducedDataStructureLogFC, session)

  Search_Full_SERVER("Search", session)
  GlobalSelection_SERVER("GlobalSelection", session)
  VolcanoPlot_SERVER("VolcanoPlot", session)
  HeatMap_SERVER("HeatMap_Full_DetailsPval", session)
  HeatMap_SERVER("HeatMap_Full_DetailsLogFC", session)
  PCPlot_SERVER("PCPlot", session)
################################################################################

  counter.invalidateLater <- local({
    static <- 0
    function() { static <<- static + 1; static }
  })

  shiny::observe({
    shiny::invalidateLater(10000, session)
    # base::print(paste0(
    #   sysTimePID(),
    #   " PatternMatchR is running in idle state."
    # ))
    a <- counter.invalidateLater()
    cat(".")
    if (a %% 10 == 0) {
      cat("*")
    }
  })

  shiny::observe({
    if (!is.valid(session$userData$sessionVariables$resultDFListTrait1()) && !is.valid(session$userData$sessionVariables$resultDFListTrait2()) && !is.valid(session$userData$sessionVariables$resultDFListTrait3())){
      shinyjs::disable("Merge")
      shinyjs::disable("btnMerge")
    }
    else {
      shinyjs::enable("Merge")
      shinyjs::enable("btnMerge")
    }
  })

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
    if (!is.valid(session$userData$sessionVariables$pReducedDataStructurePVal()$combinedDFP_Val_Labels$dfP_Val)) {
      shinyjs::disable("Omit Traits")
      shinyjs::disable("sldNumClusters")
      shinyjs::disable("sld_NumNeighbours")
    }
    else {
      shinyjs::enable("Omit Traits")
      shinyjs::enable("sldNumClusters")
      shinyjs::enable("sld_NumNeighbours")
    }
  })

  shiny::observe({
    if (!is.valid(session$userData$sessionVariables$traitReducedDataStructurePVal()$combinedDFP_Val_Labels$dfP_Val)) {
      shinyjs::disable("btnPlotCombinedHM_P_Val")
      shinyjs::disable("btnPlotCombinedCondHM_DM")
      # shinyjs::disable("numHMHSize")
      # shinyjs::disable("numHMVSize")
    }
    else {
      shinyjs::enable("btnPlotCombinedHM_P_Val")
      shinyjs::enable("btnPlotCombinedCondHM_DM")
      # shinyjs::enable("numHMHSize")
      # shinyjs::enable("numHMVSize")
    }
  })

  shiny::observe({
    if (!is.valid(session$userData$sessionVariables$traitReducedDataStructureLogFC()$combinedDFP_Val_Labels$dfLogFC)) {
      shinyjs::disable("btnPlotCombinedHM_LogFC")
    }
    else{
      shinyjs::enable("btnPlotCombinedHM_LogFC")
    }
  })

  shiny::observe({
    if (!is.valid(session$userData$sessionVariables$probeReducedDataStructure()$combinedDFP_Val_Labels$dfP_Val)) {
      shinyjs::disable("Condensed Trait-reduced Data (contains only CpG with nearby neighbours)")
      shinyjs::disable("btnPlotCombinedCondHM_P_Val")
      # shinyjs::disable("RbHighlightHM")
      shinyjs::disable("txtHMDescription_P_Val")
      # shinyjs::disable("numCondHMHSize")
      # shinyjs::disable("numCondHMVSize")
    }
    else {
      shinyjs::enable("Condensed Trait-reduced Data (contains only CpG with nearby neighbours)")
      shinyjs::enable("btnPlotCombinedCondHM_P_Val")
      # shinyjs::enable("RbHighlightHM")
      shinyjs::enable("txtHMDescription_P_Val")
      # shinyjs::enable("numCondHMHSize")
      # shinyjs::enable("numCondHMVSize")
    }
  })

  # shiny::observe({
  #   if (!is.valid(session$userData$sessionVariables$distanceMultipliedTraitReducedDataStructure()$combinedDFP_Val_Labels$dfP_Val)) {
  #     #shinyjs::disable("Full Distance Weighted Data")
  #     shinyjs::disable("btnPlotCombinedDWHM_P_Val")
  #     shinyjs::disable("numDWHMHSize")
  #     shinyjs::disable("numDWHMVSize")
  #   }
  #   else {
  #     #shinyjs::enable("Full Distance Weighted Data")
  #     shinyjs::enable("btnPlotCombinedDWHM_P_Val")
  #     shinyjs::enable("numDWHMHSize")
  #     shinyjs::enable("numDWHMVSize")
  #   }
  # })
  #
  # shiny::observe({
  #   if (!is.valid(session$userData$sessionVariables$distanceMultipliedProbeReducedDataStructure()$combinedDFP_Val_Labels$dfP_Val)) {
  #     #shinyjs::disable("Condensed Distance Weighted Data (contains only CpG with nearby neighbours)")
  #     shinyjs::disable("btnPlotCombinedCondDWHM_P_Val")
  #     shinyjs::disable("numCondDWHMHSize")
  #     shinyjs::disable("numCondDWHMVSize")
  #   }
  #   else {
  #     #shinyjs::enable("Condensed Distance Weighted Data (contains only CpG with nearby neighbours)")
  #     shinyjs::enable("btnPlotCombinedCondDWHM_P_Val")
  #     shinyjs::enable("numCondDWHMHSize")
  #     shinyjs::enable("numCondDWHMVSize")
  #   }
  # })

  session$userData$sessionVariables$numClusters <- shiny::reactive({input$sldNumClusters})
  session$userData$sessionVariables$numNeighbours <- shiny::reactive({input$sld_NumNeighbours})

  # session$userData$sessionVariables$distNeigboursProbes10 <- shiny::reactive({
  #   #calculate distance from each probe to its neigbours to build a right column in heatmap to estimate relevance of heatmap findings
  #   base::tryCatch(
  #     {
  #       result <- calculateDistanceNeigboursProbes(wd = session$userData$packageWd, clustResProbes = session$userData$sessionVariables$traitReducedDataStructurePVal()$clustResProbes, annotation = session$userData$annotation, distanceToLook = 10, numCores = session$userData$numCores)
  #     },
  #     error = function(e) {
  #       base::message("An error occurred in shiny::reactive(session$userData$sessionVariables$distNeigboursProbes10):\n", e)
  #     },
  #     warning = function(w) {
  #       base::message("A warning occurred in shiny::reactive(session$userData$sessionVariables$distNeigboursProbes10):\n", w)
  #     },
  #     finally = {
  #       base::print(base::paste0(sysTimePID(), " finished generating session$userData$sessionVariables$distNeigboursProbes10. (last step in session$userData$sessionVariables$distNeigboursProbes10 <- shiny::reactive())"))
  #       return(result)
  #     }
  #   )
  # })

  session$userData$sessionVariables$distNeigboursProbes100 <- shiny::reactive({
    #calculate distance from each probe to its neigbours to build a right column in heatmap to estimate relevance of heatmap findings
    base::tryCatch(
      {
        result <- calculateDistanceNeigboursProbes(wd = session$userData$packageWd, clustResProbes = session$userData$sessionVariables$traitReducedDataStructurePVal()$clustResProbes, annotation = session$userData$annotation, distanceToLook = 100, numCores = session$userData$numCores)
      },
      error = function(e) {
        base::message("An error occurred in shiny::reactive(session$userData$sessionVariables$distNeigboursProbes100):\n", e)
      },
      warning = function(w) {
        base::message("A warning occurred in shiny::reactive(session$userData$sessionVariables$distNeigboursProbes100):\n", w)
      },
      finally = {
        base::print(base::paste0(sysTimePID(), " finished generating session$userData$sessionVariables$distNeigboursProbes100. (last step in session$userData$sessionVariables$distNeigboursProbes100 <- shiny::reactive())"))
        return(result)
      }
    )
  })

  session$userData$sessionVariables$distNeigboursProbes1000 <- shiny::reactive({
    #calculate distance from each probe to its neigbours to build a right column in heatmap to estimate relevance of heatmap findings
    base::tryCatch(
      {
        result <- calculateDistanceNeigboursProbes(wd = session$userData$packageWd, clustResProbes = session$userData$sessionVariables$traitReducedDataStructurePVal()$clustResProbes, annotation = session$userData$annotation, distanceToLook = 1000, numCores = session$userData$numCores)
      },
      error = function(e) {
        base::message("An error occurred in shiny::reactive(session$userData$sessionVariables$distNeigboursProbes1000):\n", e)
      },
      warning = function(w) {
        base::message("A warning occurred in shiny::reactive(session$userData$sessionVariables$distNeigboursProbes1000):\n", w)
      },
      finally = {
        base::print(base::paste0(sysTimePID(), " finished generating session$userData$sessionVariables$distNeigboursProbes1000. (last step in session$userData$sessionVariables$distNeigboursProbes1000 <- shiny::reactive())"))
        return(result)
      }
    )
  })

  session$userData$sessionVariables$distNeigboursProbes10000 <- shiny::reactive({
    #calculate distance from each probe to its neigbours to build a right column in heatmap to estimate relevance of heatmap findings
    base::tryCatch(
      {
        result <- calculateDistanceNeigboursProbes(wd = session$userData$packageWd, clustResProbes = session$userData$sessionVariables$traitReducedDataStructurePVal()$clustResProbes, annotation = session$userData$annotation, distanceToLook = 10000, numCores = session$userData$numCores)
      },
      error = function(e) {
        base::message("An error occurred in shiny::reactive(session$userData$sessionVariables$distNeigboursProbes10000):\n", e)
      },
      warning = function(w) {
        base::message("A warning occurred in shiny::reactive(session$userData$sessionVariables$distNeigboursProbes10000):\n", w)
      },
      finally = {
        base::print(base::paste0(sysTimePID(), " finished generating session$userData$sessionVariables$distNeigboursProbes10000. (last step in session$userData$sessionVariables$distNeigboursProbes10000 <- shiny::reactive())"))
        return(result)
      }
    )
  })

  session$userData$sessionVariables$distancesBelowThreshold <- shiny::reactive({
    #gets CpG below a certain threshold
    base::tryCatch(
      {
        distances <- session$userData$sessionVariables$distNeigboursProbes10000()
        rownames(distances) <- distances[, 1]
        threshold <- session$userData$config$CpGDistanceThreshold
        #which distances < threshold
        result <- rownames(distances)[distances$meanDistance < threshold]
        result <- na.omit(result)
#to check, make the first elements (1,3,5) of distances orange...
        if (length(result) == 0) {
          base::message(base::paste0(sysTimePID(), "no CpG below distance threshold of ", threshold, "found."))
        }
        else {
          base::print(base::paste0(sysTimePID(), " found n = ", length(result), " CpG below distance threshold of ", threshold, "."))
        }
      },
      error = function(e) {
        base::message("An error occurred in shiny::reactive(session$userData$sessionVariables$distancesBelowThreshold):\n", e)
      },
      warning = function(w) {
        base::message("A warning occurred in shiny::reactive(session$userData$sessionVariables$distancesBelowThreshold):\n", w)
      },
      finally = {
        base::print(base::paste0(sysTimePID(), " finished generating session$userData$sessionVariables$distancesBelowThreshold. (last step in session$userData$sessionVariables$distNeigboursProbes10000 <- shiny::reactive())"))
        return(result)
      }
    )
  })

  output$txtLoadOut <- shiny::reactive({
    base::tryCatch(
      {
        base::print(base::paste0(sysTimePID(), " start generating output$txtLoadOut. (first step in output$txtLoadOut <- shiny::reactive())"))
        result <- (updateTxtLoadOut(session, session$userData$sessionVariables$resultDFListTrait1(),
                                session$userData$sessionVariables$resultDFListTrait2(),
                                session$userData$sessionVariables$resultDFListTrait3()))
      },
      error = function(e) {
        base::message("An error occurred in shiny::reactive(output$txtLoadOut):\n", e)
      },
      warning = function(w) {
        base::message("A warning occurred in shiny::reactive(output$txtLoadOut):\n", w)
      },
      finally = {
        base::print(base::paste0(sysTimePID(), " finished generating output$txtLoadOut. (last step in output$txtLoadOut <- shiny::reactive())"))
        return(result)
      }
    )
  })

  output$txtMergeOut <- shiny::reactive({
    base::tryCatch(
      {
        base::print(base::paste0(sysTimePID(), " start generating output$txtMergeOut."))
        if (is.valid(session$userData$sessionVariables$combinedDataStructure())) {
#browser() #check contents inside result and compare to txtpReduceOut
          result <- updateTxtMergeOut(session$userData$sessionVariables$combinedDataStructure())
       }
       else {
         result <- NULL
       }
      },
      error = function(e) {
        base::message("An error occurred in shiny::reactive(output$txtMergeOut):\n", e)
      },
      warning = function(w) {
        base::message("A warning occurred in shiny::reactive(output$txtMergeOut):\n", w)
      },
      finally = {
        base::print(base::paste0(sysTimePID(), " finished generating output$txtMergeOut."))
        return(result)
      }
    )
  })

  output$txtPReduceOut <- shiny::reactive({
    base::tryCatch(
      {
        base::print(base::paste0(sysTimePID(), " start generating output$txtPReduceOut."))
        if (is.valid(session$userData$sessionVariables$pReducedDataStructurePVal()$combinedDFP_Val_Labels$dfP_Val)) {
          result <- updateTxtpReduceOut(session$userData$sessionVariables$pReducedDataStructurePVal()$combinedDFP_Val_Labels)
        }
        else {
          result <- NULL
        }
      },
      error = function(e) {
        base::message("An error occurred in shiny::reactive(output$txtPReduceOut):\n", e)
      },
      warning = function(w) {
        base::message("A warning occurred in shiny::reactive(output$txtPReduceOut):\n", w)
      },
      finally = {
        base::print(base::paste0(sysTimePID(), " finished generating output$txtPReduceOut."))
        return(result)
      }
    )
  })

  # output$txtCondOut <- shiny::reactive({
  #   base::tryCatch(
  #     {
  #       base::print(base::paste0(sysTimePID(), " start generating output$txtCondOut"))
  #       if (is.valid(session$userData$sessionVariables$probeReducedDataStructure()$combinedDFP_Val_Labels)) {
  #         result <- updateTxtOmitTraitsOut(session$userData$sessionVariables$probeReducedDataStructure()$combinedDFP_Val_Labels)
  #       }
  #       else {
  #         result <- NULL
  #       }
  #     },
  #     error = function(e) {
  #       base::message("An error occurred in shiny::reactive(output$txtCondOut):\n", e)
  #     },
  #     warning = function(w) {
  #       base::message("A warning occurred in shiny::reactive(output$txtCondOut):\n", w)
  #     },
  #     finally = {
  #       base::print(base::paste0(sysTimePID(), " finished generating output$txtCondOut"))
  #       return(result)
  #     }
  #   )
  # })

  shiny::observeEvent(input$save,
    ignoreInit = TRUE,
    {
      if (is.integer(input$save)) {
        #                        cat("No file have been selected for save.")
      } else {
        result <- shinyFiles::parseSavePath(volumes, input$save)
        filePath <- as.character(result$datapath)
        cat(paste0(filePath, " has been selected."))
        # insert
        base::saveRDS(file = filePath, session$userData$sessionVariables)
        base::print(base::paste0(sysTimePID(), " session data has been saved to ", filePath))
      }
    },
    ignoreNULL = FALSE
  )

  shiny::observeEvent(input$btnBrowser,
    ignoreInit = TRUE,
    {
      browser() #for debugging purposes
    },
    ignoreNULL = FALSE
  )

  shiny::observeEvent(input$btnLoadDir1,
    ignoreInit = TRUE,
    {
      base::tryCatch(
        {
          base::print(base::paste0(sysTimePID(), " start loading trait 1 folders."))
          if (base::is.numeric(input$trait1DirList_rows_selected)) {
            traitDirList <-
              base::as.list(dfdD1[input$trait1DirList_rows_selected, ]) #base::as.list(dfdD1[input$trait1DirList_rows_selected, ][[1]])
            base::print(base::paste0(sysTimePID(), " selected folders: ", as.character(traitDirList)))
            session$userData$sessionVariables$resultDFListTrait1(loadDir(session = session, traitDirList = traitDirList))
          }
          else {
            session$userData$sessionVariables$resultDFListTrait1(NULL)
            base::message(base::paste0(sysTimePID(), " no entries selected from trait1 folders."))
          }
        },
        error = function(e) {
          base::message("An error occurred in shiny::observeEvent(input$btnLoadDir1):\n", e)
        },
        warning = function(w) {
          base::message("A warning occurred in shiny::observeEvent(input$btnLoadDir1):\n", w)
        },
        finally = {
          base::print(base::paste0(sysTimePID(), " end loadDir1."))
        }
      )
    },
    ignoreNULL = FALSE
  )

  shiny::observeEvent(input$btnLoadDir2,
    ignoreInit = TRUE,
    {
      base::tryCatch(
        {
          base::print(base::paste0(sysTimePID(), " start loading trait 2 folders."))
          if (base::is.numeric(input$trait2DirList_rows_selected)) {
            traitDirList <-
              base::as.list(dfdD2[input$trait2DirList_rows_selected, ]) #base::as.list(dfdD2[input$trait2DirList_rows_selected, ][[1]])
            session$userData$sessionVariables$resultDFListTrait2(loadDir(session = session, traitDirList = traitDirList))
          }
          else {
            session$userData$sessionVariables$resultDFListTrait2(NULL)
            base::message(base::paste0(sysTimePID(), " no entries selected from trait2 folders."))
          }
        },
        error = function(e) {
          base::message("An error occurred in shiny::observeEvent(input$btnLoadDir2):\n", e)
        },
        warning = function(w) {
          base::message("A warning occurred in shiny::observeEvent(input$btnLoadDir2):\n", w)
        },
        finally = {
          base::print(base::paste0(sysTimePID(), " end loadDir2."))
        }
      )
    },
    ignoreNULL = FALSE
  )

  shiny::observeEvent(input$btnLoadDir3,
    ignoreInit = TRUE,
    {
      base::tryCatch(
        {
          base::print(base::paste0(sysTimePID(), " start loading trait3 folders."))
          if (base::is.numeric(input$trait3DirList_rows_selected)) {
            traitDirList <-
              base::as.list(dfdD3[input$trait3DirList_rows_selected, ]) #base::as.list(dfdD3[input$trait3DirList_rows_selected, ][[1]])
            session$userData$sessionVariables$resultDFListTrait3(loadDir(session = session, traitDirList = traitDirList))
          } else {
            session$userData$sessionVariables$resultDFListTrait3(NULL)
            base::message(base::paste0(sysTimePID(), " no entries selected from trait3 folders."))
          }
        },
        error = function(e) {
          base::message("An error occurred in shiny::observeEvent(input$btnLoadDir3):\n", e)
        },
        warning = function(w) {
          base::message("A warning occurred in shiny::observeEvent(input$btnLoadDir3):\n", w)
        },
        finally = {
          base::print(base::paste0(sysTimePID(), " end loadDir3."))
        }
      )
    },
    ignoreNULL = FALSE
  )

  shiny::observeEvent(input$btnLoadDirAll,
    ignoreInit = TRUE,
    {
      base::tryCatch(
        {
          base::print(base::paste0(sysTimePID(), " Step1: start loading all folders. (first step in shiny::observeEvent(btnLoadDirAll))"))
          if (base::is.numeric(input$trait1DirList_rows_selected)) {
            traitDirList <-
              base::as.list(dfdD1[input$trait1DirList_rows_selected, ])
            base::print(base::paste0(sysTimePID(), " traitDirList1: ", as.character(traitDirList)))
            session$userData$sessionVariables$resultDFListTrait1(NULL)
            session$userData$sessionVariables$resultDFListTrait1(loadDir(session = session, traitDirList = traitDirList))
          } else {
            session$userData$sessionVariables$resultDFListTrait1(NULL)
            base::message(base::paste0(sysTimePID(), " no entries selected from trait1 folders."))
          }
          if (base::is.numeric(input$trait2DirList_rows_selected)) {
            traitDirList <-
              base::as.list(dfdD2[input$trait2DirList_rows_selected, ])
            base::print(base::paste0(sysTimePID(), " traitDirList2: ", as.character(traitDirList)))
            session$userData$sessionVariables$resultDFListTrait2(NULL)
            session$userData$sessionVariables$resultDFListTrait2(loadDir(session = session, traitDirList = traitDirList))
          } else {
            session$userData$sessionVariables$resultDFListTrait2(NULL)
            base::message(base::paste0(sysTimePID(), " no entries selected from trait2 folders."))
          }
          if (base::is.numeric(input$trait3DirList_rows_selected)) {
            traitDirList <-
              base::as.list(dfdD3[input$trait3DirList_rows_selected, ])
            base::print(base::paste0(sysTimePID(), " traitDirList3: ", as.character(traitDirList)))
            session$userData$sessionVariables$resultDFListTrait3(NULL)
            session$userData$sessionVariables$resultDFListTrait3(loadDir(session = session, traitDirList = traitDirList))
          } else {
            session$userData$sessionVariables$resultDFListTrait3(NULL)
            base::message(base::paste0(sysTimePID(), " no entries selected from trait3 folders."))
          }
        },
        error = function(e) {
          base::message("An error occurred in shiny::observeEvent(input$btnLoadDirAll):\n", e)
        },
        warning = function(w) {
          base::message("A warning occurred in shiny::observeEvent(input$btnLoadDirAll):\n", w)
        },
        finally = {
          base::print(base::paste0(sysTimePID(), " finished loading all folders. (last step in shiny::observeEvent(btnLoadDirAll))"))
        }
      )
    },
    ignoreNULL = FALSE
  )

  shiny::observeEvent(input$btnDebug,
    ignoreInit = TRUE,
    {
      base::tryCatch(
        {
          minN <- base::as.integer(input$txtCases)
          # output$plotDebug1 <-
          #   shiny::renderPlot(session$userData$sessionVariables$clustResTraits())
        },
        error = function(e) {
          base::message("An error occurred in shiny::observeEvent(input$btnDebug):\n", e)
        },
        warning = function(w) {
          base::message("A warning occurred shiny::observeEvent(input$btnDebug):\n", w)
        },
        finally = {
          base::print(base::paste0(sysTimePID(), " end debug test."))
        }
      )
    },
    ignoreNULL = FALSE
  )

  shiny::observeEvent(input$btnMerge,
    ignoreInit = TRUE,
    {
      base::tryCatch(
        {
          base::print(base::paste0(sysTimePID(), " Step 2: start merging data. (first step in shiny::observeEvent(btnMerge))"))
          minN <- base::as.integer(input$txtCases)
          #if (session$userData$sessionVariables$LoadInitialized == TRUE) {
          if (is.valid(session$userData$sessionVariables$resultDFListTrait1()) || is.valid(session$userData$sessionVariables$resultDFListTrait2())  || is.valid(session$userData$sessionVariables$resultDFListTrait3())) {
            result <- mergeDFP_Val_Labels(session$userData$sessionVariables$resultDFListTrait1(),
                                                          session$userData$sessionVariables$resultDFListTrait2(),
                                                          session$userData$sessionVariables$resultDFListTrait3(),
                                                          minN)
          }
          else {
            base::print(base::paste0(sysTimePID(), " is.valid(session$userData$sessionVariables$resultDFListTrait1()) || is.valid(session$userData$sessionVariables$resultDFListTrait2())  || is.valid(session$userData$sessionVariables$resultDFListTrait3()) == FALSE."))
            result <- NULL
          }
#browser() #check: here we have dflogFC
          updateReduceDataSliders(session, result)
          #updateSliders(session, session$userData$sessionVariables$combinedDataStructure()$combinedDFP_Val_Labels)
          session$userData$sessionVariables$combinedData(result)
        },
        error = function(e) {
          base::message("An error occurred in shiny::observeEvent(input$btnMerge):\n", e)
        },
        warning = function(w) {
          base::message("A warning occurred in shiny::observeEvent(input$btnMerge):\n", w)
        },
        finally = {
          base::print(base::paste0(sysTimePID(), " end merging data. (last step in shiny::observeEvent(btnMerge))"))
        }
      )
    },
    ignoreNULL = FALSE
  )

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
          #browser() #check for dflogFC and for return result in the end
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
              if (is.valid(result$dfP_Val)) {
                maxTraits <- ncol(result$dfP_Val)
                value <- input$sldNumClusters #maxTraits
                if (value == 0) {value <- maxTraits}
                shinyjs::enable("sldNumClusters")
                shiny::updateSliderInput(session = session, inputId = "sldNumClusters", max = maxTraits, min = 1, value = value, step = 1)
#                shinyjs::disable("sldNumClusters")
              }
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
          base::message("An error occurred in shiny::observeEvent(input$btnReduce):\n", e)
        },
        warning = function(w) {
          base::message("A warning occurred in shiny::observeEvent(input$btnReduce):\n", w)
        },
        finally = {
          session$userData$sessionVariables$pReducedData(result)
          base::print(base::paste0(sysTimePID(), " finished reducing data by p-value."))
        }
      )
    },
    ignoreNULL = FALSE
  )

  shiny::observeEvent(input$btnCountProbesP_ValParallel,
    ignoreInit = TRUE,
    {
      base::tryCatch(
        {
          base::print(base::paste0(sysTimePID(), " start counting probes p-value."))
          minN <- base::as.integer(input$txtCases)
          if (is.valid(session$userData$sessionVariables$combinedDataStructure()$combinedDFP_Val_Labels)) {
            P_VALNTable <-
              getAvailNForP_VALBorderParallel(session = session, wd = session$userData$packageWd, numCores = session$userData$numCores, DF = session$userData$sessionVariables$combinedDataStructure()$combinedDFP_Val_Labels$dfP_Val)
            output$DTP_VALborder <- DT::renderDataTable(P_VALNTable)
            #insert scatterplot with table results
            plot <- plotly::plot_ly(x = P_VALNTable$P_VAL_BORDER, y = P_VALNTable$'Available CpG' , type = "scatter", mode = "lines+markers", name = "scatterP_ValBorder") %>%
              plotly::layout(xaxis = list(title = "p-val", type = "log")) %>%
              plotly::layout(yaxis = list(title = "n"))
            output$plotDendrogramP_VALborder <- plotly::renderPlotly(plot)
            base::print(base::paste0(sysTimePID(), " finished counting probes p-value."))
          }
          else {
            base::print(base::paste0(sysTimePID(), " is.valid(session$userData$sessionVariables$combinedDataStructure()$combinedDFP_Val_Labels) == FALSE."))
          }
        },
        error = function(e) {
          base::message("An error occurred in shiny::observeEvent(input$btnCountProbesP_ValParallel):\n", e)
        },
        warning = function(w) {
          base::message("A warning occurred in shiny::observeEvent(input$btnCountProbesP_ValParallel):\n", w)
        },
        finally = {
          base::print(base::paste0(sysTimePID(), " end counting probes p-value."))
        }
      )
    },
    ignoreNULL = FALSE
  )

  shiny::observeEvent(input$btnCountProbesDeltaMethParallel,
    ignoreInit = TRUE,
    {
      base::tryCatch(
        {
          base::print(base::paste0(sysTimePID(), " start counting probes delta methylation."))
          minN <- base::as.integer(input$txtCases)
          if (is.valid(session$userData$sessionVariables$combinedDataStructure()$combinedDFP_Val_Labels)) {
            DMNTable <-
              getAvailNForDMBorderParallel(session = session, wd = session$userData$packageWd, numCores = session$userData$numCores, DF = session$userData$sessionVariables$combinedDataStructure()$combinedDFP_Val_Labels$dfDM)
            output$DTDMborder <- DT::renderDataTable(DMNTable)
            plot <- plotly::plot_ly(x = DMNTable$DM_BORDER, y = DMNTable$'Available CpG', type = "scatter", mode = "lines+markers", name = "scatterDeltaMethylationBorder") %>%
              plotly::layout(xaxis = list(title = "DeltaMethylation", type = "linear")) %>%
              plotly::layout(yaxis = list(title = "n"))
            output$plotDendrogramDMborder <- plotly::renderPlotly(plot)
          }
          else {
            base::print(base::paste0(sysTimePID(), " is.valid(session$userData$sessionVariables$combinedDataStructure()$combinedDFP_Val_Labels) == FALSE."))
          }
        },
        error = function(e) {
          base::message("An error occurred in shiny::observeEvent(input$btnCountProbesDeltaMethParallel):\n", e)
        },
        warning = function(w) {
          base::message("A warning occurred in shiny::observeEvent(input$btnCountProbesDeltaMethParallel):\n", w)
        },
        finally = {
          base::print(base::paste0(sysTimePID(), " end counting probes delta methylation."))
        }
      )
    },
    ignoreNULL = FALSE
  )

  shiny::observeEvent(input$btnCountProbesNParallel,
    ignoreInit = TRUE,
    {
      base::tryCatch(
        {
          base::print(base::paste0(sysTimePID(), " start counting probes n."))
          minN <- base::as.integer(input$txtCases)
          if (is.valid(session$userData$sessionVariables$combinedDataStructure()$combinedDFP_Val_Labels)) {
            NNTable <-
              getAvailNForNBorderParallel(session = session, wd = session$userData$packageWd, numCores = session$userData$numCores, DF = session$userData$sessionVariables$combinedDataStructure()$combinedDFP_Val_Labels$dfN)
            output$DTNborder <- DT::renderDataTable(NNTable)
            plot <- plotly::plot_ly(x = NNTable$N_BORDER, y = NNTable$'Available CpG', type = "scatter", mode = "lines+markers", name = "scatterNBorder") %>%
              plotly::layout(xaxis = list(title = "n", type = "linear")) %>%
              plotly::layout(yaxis = list(title = "n"))
            output$plotDendrogramNborder <- plotly::renderPlotly(plot)
          }
          else {
            base::print(base::paste0(sysTimePID(), " is.valid(session$userData$sessionVariables$combinedDataStructure()$combinedDFP_Val_Labels) == FALSE."))
          }
        },
        error = function(e) {
          base::message("An error occurred in shiny::observeEvent(input$btnCountProbesNParallel):\n", e)
        },
        warning = function(w) {
          base::message("A warning occurred in shiny::observeEvent(input$btnCountProbesNParallel):\n", w)
        },
        finally = {
          base::print(base::paste0(sysTimePID(), " finished counting probes n."))
        }
      )
    },
    ignoreNULL = FALSE
  )

  histP_Val <- shiny::reactive({
    base::tryCatch({
      base::print(base::paste0(sysTimePID(), " start render plotly histP_Val()."))
      h <- sort(as.numeric(unlist(session$userData$sessionVariables$traitReducedDataStructurePVal()$combinedDFP_Val_Labels$dfP_Val)))
      result <- plotly::plot_ly(x = h, type = "histogram", name = "histP_Val")
    },
    error = function(e) {
      base::message("An error occurred in histP_Val <- shiny::reactive():\n", e)
      browser() #should not happen
    },
    warning = function(w) {
      base::message("A warning occurred in histP_Val <- shiny::reactive():\n", w)
    },
    finally = {
      base::print(base::paste0(sysTimePID(), " finished render plotly histP_Val()."))
      return(result)
    })
  })

  histLogFC <- shiny::reactive({
    base::tryCatch({
      base::print(base::paste0(sysTimePID(), " start render plotly histLogFC()."))
      h <- sort(as.numeric(unlist(session$userData$sessionVariables$traitReducedDataStructureLogFC()$combinedDFP_Val_Labels$dfP_Val)))
      result <- plotly::plot_ly(x = h, type = "histogram", name = "histLogFC")
    },
    error = function(e) {
      base::message("An error occurred in histLogFC <- shiny::reactive():\n", e)
      browser() #should not happen
    },
    warning = function(w) {
      base::message("A warning occurred in histLogFC <- shiny::reactive():\n", w)
    },
    finally = {
      base::print(base::paste0(sysTimePID(), " finished render plotly histLogFC()."))
      return(result)
    })
  })

  # fullDWHistP_Val <- shiny::reactive({
  #   base::tryCatch({
  #     base::print(base::paste0(sysTimePID(), " start render plotly fullDWHistP_Val()."))
  #     P_Val <- sort(as.numeric(unlist(session$userData$sessionVariables$distanceMultipliedTraitReducedDataStructurePVal()$combinedDFP_Val_Labels$dfP_Val)))
  #     result <- plotly::plot_ly(x = P_Val, type = "histogram", name = "fullDWHistP_Val")
  #   },
  #   error = function(e) {
  #     base::message("An error occurred in fullDWHistP_Val <- shiny::reactive():\n", e)
  #     browser() #should not happen
  #   },
  #   warning = function(w) {
  #     base::message("A warning occurred in fullDWHistP_Val <- shiny::reactive():\n", w)
  #   },
  #   finally = {
  #     base::print(base::paste0(sysTimePID(), " finished render plotly fullDWHistP_Val()."))
  #     return(result)
  #   })
  # })

  condHistP_Val <- shiny::reactive({
    base::tryCatch({
      base::print(base::paste0(sysTimePID(), " start render plotly condhistP_Val()."))
      P_Val <- sort(as.numeric(unlist(session$userData$sessionVariables$probeReducedDataStructure()$combinedDFP_Val_Labels$dfP_Val)))
      result <- plotly::plot_ly(x = P_Val, type = "histogram", name = "histP_Val")
    },
    error = function(e) {
      base::message("An error occurred in condHistP_Val <- shiny::reactive():\n", e)
      browser() #should not happen
    },
    warning = function(w) {
      base::message("A warning occurred in condHistP_Val <- shiny::reactive():\n", w)
    },
    finally = {
      base::print(base::paste0(sysTimePID(), " finished render plotly condHistP_Val()."))
      return(result)
    })
  })

  # condDWHistP_Val <- shiny::reactive({
  #   base::tryCatch({
  #     base::print(base::paste0(sysTimePID(), " start render plotly condDWhistP_Val()."))
  #     P_Val <- sort(as.numeric(unlist(session$userData$sessionVariables$distanceMultipliedProbeReducedDataStructure()$combinedDFP_Val_Labels$dfP_Val)))
  #     result <- plotly::plot_ly(x = P_Val, type = "histogram", name = "histP_Val")
  #   },
  #   error = function(e) {
  #     base::message("An error occurred in condDWHistP_Val <- shiny::reactive():\n", e)
  #     browser() #should not happen
  #   },
  #   warning = function(w) {
  #     base::message("A warning occurred in condDWHistP_Val <- shiny::reactive():\n", w)
  #   },
  #   finally = {
  #     base::print(base::paste0(sysTimePID(), " finished render plotly condDWHistP_Val()."))
  #     return(result)
  #   })
  # })

  shiny::observeEvent(input$btnExportCombinedData,
    ignoreInit = TRUE,
    {
      base::tryCatch(
        {
          base::print(base::paste0(sysTimePID(), " start export combined data. (first step in shiny::observeEvent(btnExportCombinedData))"))
          fileNameCombinedHM <- base::paste0("CombinedHM.RDS")
          fileName <-
            base::paste0(session$userData$config$workDir, fileNameCombinedHM) #base::paste0(globalVariables$config$workDir, fileNameCombinedHM)
          base::print(base::paste0(sysTimePID(), " start exporting session data to ", fileName, "."))
          base::saveRDS(file = fileName, session$userData)
          base::print(base::paste0(sysTimePID(), " end exporting session data to ", fileName, "."))
        },
        error = function(e) {
          base::message("An error occurred in shiny::observeEvent(input$btnExportCombinedData):\n", e)
        },
        warning = function(w) {
          base::message("A warning occurred in shiny::observeEvent(input$btnExportCombinedData):\n", w)
        },
        finally = {
          base::print(base::paste0(sysTimePID(), " finished export combined data. (last step in shiny::observeEvent(btnExportCombinedData))"))
        }
      )
    },
    ignoreNULL = FALSE
  )

  shiny::observeEvent(input$btnImportCombinedData,
    ignoreInit = TRUE,
    {
      base::tryCatch(
        {
          base::print(base::paste0(sysTimePID(), " start import combined data. (first step in shiny::observeEvent(btnImportCombinedData))"))
          fileNameCombinedHM <- base::paste0("CombinedHM.RDS")
          fileName <-
            base::paste0(session$userData$config$workDir, fileNameCombinedHM) #base::paste0(globalVariables$config$workDir, fileNameCombinedHM)
          base::print(base::paste0(sysTimePID(), " start importing data from ", fileName, "."))
          if (utils::file_test("-f", fileName) == TRUE) {

            session$userData <-
              base::readRDS(file = fileName)
            base::print(base::paste0(sysTimePID(), " end reading data"))
          }
          base::print(base::paste0(sysTimePID(), " end importing session data from ", fileName, "."))
        },
        error = function(e) {
          base::message("An error occurred in shiny::observeEvent(input$btnImportCombinedData):\n", e)
        },
        warning = function(w) {
          base::message("A warning occurred in shiny::observeEvent(input$btnImportCombinedData):\n", w)
        },
        finally = {
          base::print(base::paste0(sysTimePID(), " finished import combined data. (last step in shiny::observeEvent(btnImportCombinedData))"))
        }
      )
    },
    ignoreNULL = FALSE
  )

  #("Full Data")
  output$traitReducedDTProbesPval <- DT::renderDataTable(as.data.frame(traitReducedDTProbesPval()),
                                                     options = list(pageLength = 1000, info = FALSE,
                                                                    lengthMenu = list(c(100, 1000, -1), c("100", "1000", "All"))))
  output$traitReducedDTProbesLogFC <- DT::renderDataTable(as.data.frame(traitReducedDTProbesLogFC()),
                                                         options = list(pageLength = 1000, info = FALSE,
                                                                        lengthMenu = list(c(100, 1000, -1), c("100", "1000", "All"))))

  output$traitReducedPlotDendrogramProbesPval <- plotly::renderPlotly(ggdendro::ggdendrogram(session$userData$sessionVariables$traitReducedDataStructurePVal()$clustResProbes,
                                                                                         rotate = TRUE, theme_dendro = FALSE))
  output$traitReducedPlotDendrogramProbesLogFC <- plotly::renderPlotly(ggdendro::ggdendrogram(session$userData$sessionVariables$traitReducedDataStructureLogFC()$clustResProbes,
                                                                                            rotate = TRUE, theme_dendro = FALSE))

  output$traitReducedDTTraitsPval <- DT::renderDataTable(as.data.frame(DTTraitsPval()),
                                                     options = list(pageLength = 1000, info = FALSE,
                                                                    lengthMenu = list(c(100, 1000, -1), c("100", "1000", "All"))))
  output$traitReducedDTTraitsLogFC <- DT::renderDataTable(as.data.frame(DTTraitsLogFC()),
                                                          options = list(pageLength = 1000, info = FALSE,
                                                                         lengthMenu = list(c(100, 1000, -1), c("100", "1000", "All"))))

  output$traitReducedPlotDendrogramTraitsPval <- plotly::renderPlotly(ggdendro::ggdendrogram(session$userData$sessionVariables$traitReducedDataStructurePVal()$clustResTraits,
                                                                                         rotate = TRUE, theme_dendro = FALSE))
  output$traitReducedPlotDendrogramTraitsLogFC <- plotly::renderPlotly(ggdendro::ggdendrogram(session$userData$sessionVariables$traitReducedDataStructureLogFC()$clustResTraits,
                                                                                              rotate = TRUE, theme_dendro = FALSE))

  output$traitReducedHistP_Val <- plotly::renderPlotly(histP_Val())
  output$traitReducedHistLogFC <- plotly::renderPlotly(histLogFC())

  #("Condensed Distance Weighted Data (contains only CpG with nearby neighbours)")
  # output$condDWDTP_VAL <- DT::renderDataTable(as.data.frame(session$userData$sessionVariables$distanceMultipliedProbeReducedDataStructure()$combinedDFP_Val_Labels$dfP_Val_w_number))
  # output$condDWDTDM <- DT::renderDataTable(as.data.frame(session$userData$sessionVariables$distanceMultipliedProbeReducedDataStructure()$combinedDFP_Val_Labels$dfDM_w_number))
  # output$condDWDTN <- DT::renderDataTable(as.data.frame(session$userData$sessionVariables$distanceMultipliedProbeReducedDataStructure()$combinedDFP_Val_Labels$dfN_w_number))
  # output$condDWDTProbes <- DT::renderDataTable(as.data.frame(condDWDTProbes()),
  #                                              options = list(pageLength = 1000, info = FALSE,
  #                                                             lengthMenu = list(c(100, 1000, -1), c("100", "1000", "All"))))
  # output$condDWPlotDendrogramProbes <- plotly::renderPlotly(ggdendro::ggdendrogram(session$userData$sessionVariables$distanceMultipliedProbeReducedDataStructure()$clustResProbes,
  #                                                                                 rotate = TRUE, theme_dendro = FALSE))
  # output$condDWDTTraits <- DT::renderDataTable(as.data.frame(condDWDTTraits()),
  #                                              options = list(pageLength = 1000, info = FALSE,
  #                                                             lengthMenu = list(c(100, 1000, -1), c("100", "1000", "All"))))
  # output$condDWPlotDendrogramTraits <- plotly::renderPlotly(ggdendro::ggdendrogram(session$userData$sessionVariables$distanceMultipliedProbeReducedDataStructure()$clustResTraits,
  #                                                                                  rotate = TRUE, theme_dendro = FALSE))
  # output$condDWHistP_Val <- plotly::renderPlotly(condDWHistP_Val())
  #
  # #("Full Distance Weighted Data")
  # output$fullDWDTP_VALPval <- DT::renderDataTable(as.data.frame(session$userData$sessionVariables$distanceMultipliedTraitReducedDataStructurePVal()$combinedDFP_Val_Labels$dfP_Val_w_number))
  # output$fullDWDTDMPval <- DT::renderDataTable(as.data.frame(session$userData$sessionVariables$distanceMultipliedTraitReducedDataStructurePVal()$combinedDFP_Val_Labels$dfDM_w_number))
  # output$fullDWDTLogFCPval <- DT::renderDataTable(as.data.frame(session$userData$sessionVariables$distanceMultipliedTraitReducedDataStructurePVal()$combinedDFP_Val_Labels$dfLogFC_w_number))
  # output$fullDWDTNPval <- DT::renderDataTable(as.data.frame(session$userData$sessionVariables$distanceMultipliedTraitReducedDataStructurePVal()$combinedDFP_Val_Labels$dfN_w_number))
  # output$fullDWDTProbesPval <- DT::renderDataTable(as.data.frame(fullDWDTProbes()),
  #                                              options = list(pageLength = 1000, info = FALSE,
  #                                                             lengthMenu = list(c(100, 1000, -1), c("100", "1000", "All"))))
  # output$fullDWPlotDendrogramProbesPval <- plotly::renderPlotly(ggdendro::ggdendrogram(session$userData$sessionVariables$distanceMultipliedTraitReducedDataStructurePVal()$clustResProbes,
  #                                                                                 rotate = TRUE, theme_dendro = FALSE))
  # output$fullDWDTTraitsPval <- DT::renderDataTable(as.data.frame(fullDWDTTraits()),
  #                                              options = list(pageLength = 1000, info = FALSE,
  #                                                             lengthMenu = list(c(100, 1000, -1), c("100", "1000", "All"))))
  # output$fullDWPlotDendrogramTraitsPval <- plotly::renderPlotly(ggdendro::ggdendrogram(session$userData$sessionVariables$distanceMultipliedTraitReducedDataStructurePVal()$clustResTraits,
  #                                                                                        rotate = TRUE, theme_dendro = FALSE))
  # output$fullDWHistP_Val <- plotly::renderPlotly(fullDWHistP_Val())

  #("Condensed Data (contains only CpG with nearby neighbours)")
  output$condDTP_VAL <- DT::renderDataTable(as.data.frame(session$userData$sessionVariables$probeReducedDataStructure()$combinedDFP_Val_Labels$dfP_Val_w_number))
  output$condDTDM <- DT::renderDataTable(as.data.frame(session$userData$sessionVariables$probeReducedDataStructure()$combinedDFP_Val_Labels$dfDM_w_number))
  output$condDTLogFC <- DT::renderDataTable(as.data.frame(session$userData$sessionVariables$probeReducedDataStructure()$combinedDFP_Val_Labels$dfLogFC_w_number))

  output$condDTVolcano <- DT::renderDataTable(as.data.frame(session$userData$sessionVariables$probeReducedDataStructure()$dfVolcano))
  output$condPlotVolcano <- plotly::renderPlotly(probeReducedVolcano())

  output$condDTN <- DT::renderDataTable(as.data.frame(session$userData$sessionVariables$probeReducedDataStructure()$combinedDFP_Val_Labels$dfN_w_number))
  output$condDTProbes <- DT::renderDataTable(as.data.frame(condDTProbes()),
                                     options = list(pageLength = 1000, info = FALSE,
                                                    lengthMenu = list(c(100, 1000, -1), c("100", "1000", "All"))))
  output$condPlotDendrogramProbes <- plotly::renderPlotly(ggdendro::ggdendrogram(session$userData$sessionVariables$probeReducedDataStructure()$clustResProbes,
                                                                          rotate = TRUE, theme_dendro = FALSE))
  output$condDTTraits <- DT::renderDataTable(as.data.frame(condDTTraits()),
                                      options = list(pageLength = 1000, info = FALSE,
                                                     lengthMenu = list(c(100, 1000, -1), c("100", "1000", "All"))))

  output$condPlotDendrogramTraits <- plotly::renderPlotly(ggdendro::ggdendrogram(session$userData$sessionVariables$probeReducedDataStructure()$clustResTraits,
                                                                          rotate = TRUE, theme_dendro = FALSE))
  output$condHistP_Val <- plotly::renderPlotly(condHistP_Val())

  #this data structure holds everything (as a named list), that is needed for working with trait reduced (by selecting a subset of traits) HM
  session$userData$sessionVariables$combinedDataStructure <- shiny::reactive({
    id <- shiny::showNotification("Creating combined reduced data structure...", duration = NULL, closeButton = FALSE)
    on.exit(shiny::removeNotification(id), add = TRUE)
     base::tryCatch(
       {
         result <- base::list(combinedDFP_Val_Labels = session$userData$sessionVariables$combinedData()
         )
 #browser() #check for missings (should be there) and negative values (too) -> checked, right
 #View(result$combinedDFP_Val_Labels$dfLogFC)
      },
      error = function(e) {
        base::message("An error occurred in shiny::reactive(session$userData$sessionVariables$combinedDataStructure):\n", e)
        browser() #should not happen
      },
      warning = function(w) {
        base::message("A warning occurred in shiny::reactive(session$userData$sessionVariables$combinedDataStructure):\n", w)
        browser() #should not happen
      },
      finally = {
        base::print(base::paste0(sysTimePID(), " finished generating session$userData$sessionVariables$combinedDataStructure"))
        return(result)
      }
    )
  })

  session$userData$sessionVariables$pReducedDataStructurePVal <- shiny::reactive({
    id <- shiny::showNotification("Creating p reduced data structure (p-val)...", duration = NULL, closeButton = FALSE)
    on.exit(shiny::removeNotification(id), add = TRUE)
    base::tryCatch(
      {
        if (is.valid(session$userData$sessionVariables$pReducedData())) {
          result <- base::list(combinedDFP_Val_Labels = session$userData$sessionVariables$pReducedData(),
                           matP_Val.t = NULL,
                           distMatTraits = NULL,
                           clustResTraits = NULL,
                           traitClusters = NULL,
                           traitClusterMedoids = NULL,
                           traitDendrogram = NULL,
                           traitClustergram = NULL,
                           distMatProbes = NULL,
                           clustResProbes = NULL,
                           dendProbes = NULL
          )
#browser() #check for missings (should be there) and negative values (too) -> checked, it's fine now
#View(result$combinedDFP_Val_Labels$dfLogFC)
#View(result$combinedDFP_Val_Labels$dfP_Val)
          result$matP_Val.t <- t(as.matrix(result$combinedDFP_Val_Labels$dfP_Val))
          numberCores <- session$userData$numCores
          base::print(paste0(sysTimePID(), " (pReducedDataStructurePVal) before distance matrix (p-val) for n(reduced traits) = ", base::nrow(result$matP_Val.t), " (takes some time). Using n(cores) = ", numberCores, "."))
          if (is.valid(result$matP_Val.t)) {
            result$distMatTraits <- getDistMat(numberCores = numberCores, matrix = result$matP_Val.t)
          }
          else {
            result$distMatTraits <- NULL
          }
          base::print(paste0(sysTimePID(), " (pReducedDataStructurePVal) after distance matrix (p-val) for reduced traits."))
          base::print(paste0(sysTimePID(), " (pReducedDataStructurePVal) before clustering for traits.", nrow(result$matP_Val.t)))
          if (is.valid(result$distMatTraits)) {
            result$clustResTraits <- getClustResFast(result$distMatTraits)
          }
          else {
            result$clustResTraits <- NULL
          }
          base::print(paste0(sysTimePID(), " after clustering results for traits."))
        }
        else {
          base::print(paste0(sysTimePID(), " is.valid(session$userData$sessionVariables$pReducedData()) == FALSE"))
          result <- NULL
        }
      },
      error = function(e) {
        base::message("An error occurred in shiny::reactive(session$userData$sessionVariables$pReducedDataStructurePVal):\n", e)
        browser() #should not happen
      },
      warning = function(w) {
        base::message("A warning occurred in shiny::reactive(session$userData$sessionVariables$pReducedDataStructurePVal):\n", w)
        browser() #should not happen
      },
      finally = {
        base::print(base::paste0(sysTimePID(), " finished generating session$userData$sessionVariables$pReducedDataStructurePVal"))
        return(result)
      }
    )
  })

  session$userData$sessionVariables$pReducedDataStructureLogFC <- shiny::reactive({
    id <- shiny::showNotification("Creating p reduced data structure (log(FC))...", duration = NULL, closeButton = FALSE)
    on.exit(shiny::removeNotification(id), add = TRUE)
    base::tryCatch(
      {
        if (is.valid(session$userData$sessionVariables$pReducedData())) {
          result <- base::list(combinedDFP_Val_Labels = session$userData$sessionVariables$pReducedData(),
                               #matP_Val.t = NULL,
                               matLogFC.t = NULL,
                               distMatTraits = NULL,
                               clustResTraits = NULL,
                               traitClusters = NULL,
                               traitClusterMedoids = NULL,
                               traitDendrogram = NULL,
                               traitClustergram = NULL,
                               distMatProbes = NULL,
                               clustResProbes = NULL,
                               dendProbes = NULL
          )
          #check for missings (should be there) and negative values (too) -> checked, it's fine now
          matLogFC <- na.omit(result$combinedDFP_Val_Labels$dfLogFC) #omit empty rows
          result$matLogFC.t <- t(as.matrix(matLogFC))
          numberCores <- session$userData$numCores
          base::print(paste0(sysTimePID(), " (pReducedDataStructureLogFC) before distance matrix (log(FC)) for n(reduced traits) = ", base::nrow(result$matP_Val.t), " (takes some time). Using n(cores) = ", numberCores, "."))
          if (is.valid(result$matLogFC.t)) {
            result$distMatTraits <- getDistMat(numberCores = numberCores, matrix = result$matLogFC.t)
          }
          else {
            result$distMatTraits <- NULL
          }
          base::print(paste0(sysTimePID(), " (pReducedDataStructureLogFC) after distance matrix (log(FC)) for reduced traits."))
          base::print(paste0(sysTimePID(), " (pReducedDataStructureLogFC) before clustering for traits.", nrow(result$matP_Val.t)))
          if (is.valid(result$distMatTraits)) {
            result$clustResTraits <- getClustResFast(result$distMatTraits)
          }
          else {
            result$clustResTraits <- NULL
          }
          base::print(paste0(sysTimePID(), " after clustering results for traits."))
        }
        else {
          base::print(paste0(sysTimePID(), " is.valid(session$userData$sessionVariables$pReducedData()) == FALSE"))
          result <- NULL
        }
      },
      error = function(e) {
        base::message("An error occurred in shiny::reactive(session$userData$sessionVariables$pReducedDataStructureLogFC):\n", e)
        browser() #should not happen
      },
      warning = function(w) {
        base::message("A warning occurred in shiny::reactive(session$userData$sessionVariables$pReducedDataStructureLogFC):\n", w)
        browser() #should not happen
      },
      finally = {
        base::print(base::paste0(sysTimePID(), " finished generating session$userData$sessionVariables$pReducedDataStructureLogFC"))
        return(result)
      }
    )
  })

  session$userData$sessionVariables$traitReducedDataStructurePVal <- shiny::reactive({
    id <- shiny::showNotification("Creating trait reduced data structure (p-val)...", duration = NULL, closeButton = FALSE)
    on.exit(shiny::removeNotification(id), add = TRUE)
    base::tryCatch(
      {
#browser() #check for error "attempt to apply non-function"
        if (is.valid(session$userData$sessionVariables$pReducedDataStructurePVal()$combinedDFP_Val_Labels$dfP_Val)) {
          result <- base::list(combinedDFP_Val_Labels = NULL,
                               #matP_Val = session$userData$sessionVariables$matP_Val(), #is part of combinedDFP_Val_Labels()
                               matP_Val.t = NULL,
                               distMatTraits = NULL,
                               clustResTraits = NULL,
                               traitClusters = NULL,
                               traitClusterMedoids = NULL,
                               traitDendrogram = NULL,
                               traitClustergram = NULL,
                               distMatProbes = NULL,
                               clustResProbes = NULL,
                               probeDendrogram = NULL,
                               DNAdistances = NULL,
                               dfVolcano = NULL,
                               dfKeyShadow = NULL
          )
          if (is.valid(session$userData$sessionVariables$traitReducedDataPVal())) {
            result$combinedDFP_Val_Labels <- session$userData$sessionVariables$traitReducedDataPVal()
            #check whats going on with dfVolcano... everything is fine, except too few cases from debug mode
            result$matP_Val.t <- t(as.matrix(result$combinedDFP_Val_Labels$dfP_Val))
            numberCores <- session$userData$numCores
            base::print(paste0(sysTimePID(), " (traitReducedDataStructure) before distance matrix for n(reduced traits) = ", base::nrow(result$matP_Val.t), " (takes some time). Using n(cores) = ", numberCores, "."))
            if (is.valid(result$matP_Val.t)) {
              result$distMatTraits <- getDistMat(numberCores = numberCores, matrix = result$matP_Val.t)
            }
            else {
              result$distMatTraits <- NULL
            }
            base::print(paste0(sysTimePID(), " (traitReducedDataStructure) after distance matrix for reduced traits."))
            #for unknown reason getClustResFast() crashes, if executed without Sys.sleep in advance...
            Sys.sleep(1)
            base::print(paste0(sysTimePID(), " (traitReducedDataStructure) before clustering for traits.", base::nrow(result$matP_Val.t)))
            if (is.valid(result$distMatTraits)) {
              result$clustResTraits <- getClustResFast(result$distMatTraits)
            }
            else {
              result$clustResTraits <- NULL
            }
            base::print(paste0(sysTimePID(), " (traitReducedDataStructure) after clustering results for traits."))
            base::print(base::paste0(sysTimePID(), " (traitReducedDataStructure) start generating traitClusters."))
            numClusters <- length(result$clustResTraits$order)
            if (is.valid(result$clustResTraits) && numClusters > 1) {
              result$traitClusters <- cutree(result$clustResTraits,
                               k = numClusters)
            }
            else {
              result$traitClusters <- NULL
            }
            base::print(base::paste0(sysTimePID(), " (traitReducedDataStructure) start generating traitClusterMedoids."))
            if (is.valid(result$clustResTraits) && is.valid(result$distMatTraits)) {
              result$traitClusterMedoids <- getTraitClusterMedoids(clustResTraits = result$clustResTraits,
                                             distMatTraits = result$distMatTraits,
                                             numClusters = numClusters)
            }
            else {
              result$traitClusterMedoids <- NULL
            }

            base::print(base::paste0(sysTimePID(), " (traitReducedDataStructure) start generating traitDendrogram."))
            #do this only, if a dataset is already loaded
            if (is.valid(result$clustResTraits)) {
              result$traitDendrogram <- getDendTraits(clustResTraits = result$clustResTraits, traitClusters = numClusters)
              base::print(base::paste0(sysTimePID(), " (traitReducedDataStructure) after making traitDendrogram."))
            }
            else {
              result$traitDendrogram <- NULL
            }

            base::print(base::paste0(sysTimePID(), " (traitReducedDataStructure) start generating traitClustergram."))
            if (is.valid(result$matP_Val.t) && is.valid(result$clustResTraits)) {
              result$traitClustergram <- getplotClustergramTraitsLong(mat.t = result$matP_Val.t,
                                                   clustResTraits = result$clustResTraits,
                                                   traitClusters = numClusters)
            }
            else {
              result$traitClustergram <- NULL
            }

            base::print(base::paste0(sysTimePID(), " (traitReducedDataStructure) start generating distMatProbes."))
            dfP_Val <- result$combinedDFP_Val_Labels$dfP_Val
            if (is.valid(dfP_Val)) {
              dfP_Val[dfP_Val > 0.05] <- NA # 1
              base::print(base::paste0(sysTimePID(), " (traitReducedDataStructure) calculating distance matrix with rows= ", nrow(dfP_Val), " cols= ", ncol(dfP_Val)))
              base::print(base::class(dfP_Val))
              base::print(base::paste0(sysTimePID(), " set missing p-values to 1."))
              dfP_Val[base::is.na(dfP_Val)] <- 1 # set missing P_VAL to 1
              base::print(Cstack_info())
              if (base::nrow(dfP_Val) >= 5) {
                numberCores <- session$userData$numCores
                # clustering for rows
                base::print(base::paste0(sysTimePID(), " (traitReducedDataStructure) before distance matrix for n(probes) = ", base::nrow(dfP_Val), " (takes some time). Using n(cores) = ", numberCores, "."))
                gc()
                result$distMatProbes <- getDistMat(numberCores = numberCores, matrix = dfP_Val)
                base::print(base::paste0(sysTimePID(), " (traitReducedDataStructure) after distance matrix for probes.", base::nrow(dfP_Val)))
              }
            }
            else {
              result$distMatProbes <- NULL
              browser() # should not happen
            }
            base::print(base::paste0(sysTimePID(), " (traitReducedDataStructure) start clustResProbes."))
            distMat <- result$distMatProbes
            if (is.valid(distMat)) {
              result$clustResProbes <- getClustResFast(distMat)
            }
            else {
              result$clustResProbes <- NULL
              browser() # should not happen
            }
            if (is.valid(result$clustResTraits)) {
              # add "number" and reorder columns; order comes from result$clustResTraits
              result$combinedDFP_Val_Labels$dfP_Val <- result$combinedDFP_Val_Labels$dfP_Val[, result$clustResTraits$order]
              result$combinedDFP_Val_Labels$dfDM <- result$combinedDFP_Val_Labels$dfP_Val[, result$clustResTraits$order]
              result$combinedDFP_Val_Labels$dfN <- result$combinedDFP_Val_Labels$dfN[, result$clustResTraits$order]
              result$combinedDFP_Val_Labels$dfLogFC <- result$combinedDFP_Val_Labels$dfLogFC[, result$clustResTraits$order]
            }
            else {
              browser() #should not happen
            }

            if (is.valid(result$clustResProbes)) {
              result$probeDendrogram <- stats::as.dendrogram(result$clustResProbes)
              #reorder rows; order comes from result$clustResProbes
              result$combinedDFP_Val_Labels$dfP_Val <- result$combinedDFP_Val_Labels$dfP_Val[result$clustResProbes$order, ]
              result$combinedDFP_Val_Labels$dfDM <- result$combinedDFP_Val_Labels$dfP_Val[result$clustResProbes$order, ]
              result$combinedDFP_Val_Labels$dfN <- result$combinedDFP_Val_Labels$dfN[result$clustResProbes$order, ]
              result$combinedDFP_Val_Labels$dfLogFC <- result$combinedDFP_Val_Labels$dfLogFC[result$clustResProbes$order, ]
            }
            else {
              result$probeDendrogram <- NULL
              browser() # should not happen
            }

            result$combinedDFP_Val_Labels$dfP_Val_w_number <- result$combinedDFP_Val_Labels$dfP_Val
            nprobes <- nrow(result$combinedDFP_Val_Labels$dfP_Val_w_number)
            result$combinedDFP_Val_Labels$dfP_Val_w_number$number <- seq(1:nprobes)
            col_order <- c("number", colnames(result$combinedDFP_Val_Labels$dfP_Val_w_number))
            result$combinedDFP_Val_Labels$dfP_Val_w_number <- result$combinedDFP_Val_Labels$dfP_Val_w_number[, col_order]
            result$combinedDFP_Val_Labels$dfP_Val_w_number <- result$combinedDFP_Val_Labels$dfP_Val_w_number[ , -which(colnames(result$combinedDFP_Val_Labels$dfP_Val_w_number) %in% "number.1")]

            result$combinedDFP_Val_Labels$dfDM_w_number <- result$combinedDFP_Val_Labels$dfDM
            nprobes <- nrow(result$combinedDFP_Val_Labels$dfDM_w_number)
            result$combinedDFP_Val_Labels$dfDM_w_number$number <- seq(1:nprobes)
            col_order <- c("number", colnames(result$combinedDFP_Val_Labels$dfDM_w_number))
            result$combinedDFP_Val_Labels$dfDM_w_number <- result$combinedDFP_Val_Labels$dfDM_w_number[, col_order]

            result$combinedDFP_Val_Labels$dfN_w_number <- result$combinedDFP_Val_Labels$dfN
            nprobes <- nrow(result$combinedDFP_Val_Labels$dfN_w_number)
            result$combinedDFP_Val_Labels$dfN_w_number$number <- seq(1:nprobes)
            col_order <- c("number", colnames(result$combinedDFP_Val_Labels$dfN_w_number))
            result$combinedDFP_Val_Labels$dfN_w_number <- result$combinedDFP_Val_Labels$dfN_w_number[, col_order]
            result$combinedDFP_Val_Labels$dfN_w_number <- result$combinedDFP_Val_Labels$dfN_w_number[ , -which(colnames(result$combinedDFP_Val_Labels$dfN_w_number) %in% "number.1")]

            result$combinedDFP_Val_Labels$dfLogFC_w_number <- result$combinedDFP_Val_Labels$dfLogFC
            nprobes <- nrow(result$combinedDFP_Val_Labels$dfLogFC_w_number)
            result$combinedDFP_Val_Labels$dfLogFC_w_number$number <- seq(1:nprobes)
            col_order <- c("number", colnames(result$combinedDFP_Val_Labels$dfLogFC_w_number))
            result$combinedDFP_Val_Labels$dfLogFC_w_number <- result$combinedDFP_Val_Labels$dfLogFC_w_number[, col_order]
            result$combinedDFP_Val_Labels$dfLogFC_w_number <- result$combinedDFP_Val_Labels$dfLogFC_w_number[ , -which(colnames(result$combinedDFP_Val_Labels$dfLogFC_w_number) %in% "number.1")]

            Distance <- input$sld_NumNeighbours
            DNAdistances <- calculateDistanceNeigboursProbes(wd = session$userData$packageWd, clustResProbes = result$clustResProbes, annotation = session$userData$annotation, distanceToLook = Distance, numCores = session$userData$numCores)
            result$DNAdistances <- DNAdistances

            dfLogFC <- result$combinedDFP_Val_Labels$dfLogFC
            if (is.valid(dfLogFC)) {
              dfP_Val <- result$combinedDFP_Val_Labels$dfP_Val
              if (is.valid(dfP_Val)) {
                #take everything into one table with columns p-val and logFC ...
                dfLogFC$probe <- row.names(dfLogFC)
                dfLogFC <- tidyr::pivot_longer(dfLogFC, cols  = -probe, names_to = c("trait"))
                colnames(dfLogFC)[3] <- "LogFC"
                dfP_Val$probe <- row.names(dfP_Val)
                dfP_Val <- tidyr::pivot_longer(dfP_Val, cols  = -probe, names_to = c("trait"))
                colnames(dfP_Val)[3] <- "P_Val"
                #dfVolcano <- base::merge(dfP_Val, dfLogFC, by.x = c("probe","trait"), by.y = c("probe","trait"), all.x = FALSE, all.y = FALSE)
                dfVolcano <- base::merge(dfP_Val, dfLogFC, by.x = c("probe","trait"), by.y = c("probe","trait"), all.x = TRUE, all.y = FALSE)
                #merge chr and position
                annotation <- base::subset(session$userData$annotation, select = c("name", "chromosome", "position", "gene.symbol"))
                dfVolcano <- base::merge(dfVolcano, annotation, by.x = "probe", by.y = "name", all.x = FALSE, all.y = FALSE)

                #add distances to dfVolcano
                DNAdistances <- result$DNAdistances
                row.names(DNAdistances) <- DNAdistances$ID
                DNAdistances$cg <- DNAdistances$ID #row.names(DNAdistances)
                DNAdistancesNumber <- DNAdistances[,c("cg", "number")]
                DNAdistancesNumber <- tidyr::pivot_longer(DNAdistancesNumber, cols  = -cg, names_to = c("number"))
                DNAdistancesNumber <- DNAdistancesNumber[,c("cg", "value")]
                colnames(DNAdistancesNumber)[2] <- "DistanceNumber"
                DNAdistancesMin <- DNAdistances[,c("cg", "minDistance")]
                DNAdistancesMin <- tidyr::pivot_longer(DNAdistancesMin, cols  = -cg, names_to = c("min"))
                DNAdistancesMin <- DNAdistancesMin[,c("cg", "value")]
                colnames(DNAdistancesMin)[2] <- "DistanceMin"

                DNAdistancesMean <- DNAdistances[,c("cg", "meanDistance")]
                DNAdistancesMean <- tidyr::pivot_longer(DNAdistancesMean, cols  = -cg, names_to = c("mean"))
                DNAdistancesMean <- DNAdistancesMean[,c("cg", "value")]
                colnames(DNAdistancesMean)[2] <- "DistanceMean"

                DNAdistancesMax <- DNAdistances[,c("cg", "maxDistance")]
                DNAdistancesMax <- tidyr::pivot_longer(DNAdistancesMax, cols  = -cg, names_to = c("max"))
                DNAdistancesMax <- DNAdistancesMax[,c("cg", "value")]
                colnames(DNAdistancesMax)[2] <- "DistanceMax"

                dfVolcano <- base::merge(dfVolcano, DNAdistancesNumber, by.x = c("probe"), by.y = c("cg"), all.x = TRUE, all.y = FALSE)
                dfVolcano <- base::merge(dfVolcano, DNAdistancesMin, by.x = c("probe"), by.y = c("cg"), all.x = TRUE, all.y = FALSE)
                dfVolcano <- base::merge(dfVolcano, DNAdistancesMean, by.x = c("probe"), by.y = c("cg"), all.x = TRUE, all.y = FALSE)
                dfVolcano <- base::merge(dfVolcano, DNAdistancesMax, by.x = c("probe"), by.y = c("cg"), all.x = TRUE, all.y = FALSE)
                #create shadow table for key assignment
                #add probe trait and traitsource to key
                originTrait <- result$combinedDFP_Val_Labels$mergedOriginTrait
                originTrait <- rep(originTrait, nprobes)
                keys <- seq(1:nrow(dfVolcano))
                dfKeyShadow <- base::data.frame(key = keys)
                dfKeyShadow$probe <- dfVolcano$probe
                traitLabels <- session$userData$sessionVariables$traitReducedDataPVal()$mergedOriginalColnames
                dfKeyShadow$trait <- traitLabels #dfVolcano$trait
                traitID <- session$userData$sessionVariables$traitReducedDataPVal()$traitID
                dfKeyShadow$traitID <- traitID
                dfKeyShadow$traitSource <- originTrait
                dfVolcano$key <- dfKeyShadow$key
                rownames(dfVolcano) <- dfVolcano$key
                rownames(dfKeyShadow) <- dfKeyShadow$key
                #sort by p-val and logFC
                dfVolcano <- dfVolcano[base::order(dfVolcano$P_Val, dfVolcano$LogFC, decreasing = c(FALSE, TRUE), na.last = c(TRUE,TRUE)),]
                #take only first 2^16 entities to be able to plot volcano plot using plotly
                if(nrow(dfVolcano)>2^16) {
                  dfVolcano <- dfVolcano[1:2^16-1,]
                }
                result$dfVolcano <- dfVolcano
                result$dfKeyShadow <- dfKeyShadow
              }
              else {
                result$dfVolcano <- NULL
                browser() # should not happen
              }
            }
            else {
              result$dfVolcano <- NULL
              browser() # should not happen
            }

          }
        }
      },
      error = function(e) {
        base::message("An error occurred in shiny::reactive(session$userData$sessionVariables$traitReducedDataStructurePVal):\n", e)
        browser() #should not happen
      },
      warning = function(w) {
        base::message("A warning occurred in shiny::reactive(session$userData$sessionVariables$traitReducedDataStructurePVal):\n", w)
        browser() #should not happen
      },
      finally = {
        base::print(base::paste0(sysTimePID(), " finished generating session$userData$sessionVariables$traitReducedDataStructurePVal."))
        return(result)
      }
    )
  })

  session$userData$sessionVariables$traitReducedDataStructureLogFC <- shiny::reactive({
    id <- shiny::showNotification("Creating trait reduced data structure (log(FC))...", duration = NULL, closeButton = FALSE)
    on.exit(shiny::removeNotification(id), add = TRUE)
    base::tryCatch(
      {
        if (is.valid(session$userData$sessionVariables$pReducedDataStructurePVal()$combinedDFP_Val_Labels$dfLogFC)) {
          result <- base::list(combinedDFP_Val_Labels = NULL,
                               #matP_Val = session$userData$sessionVariables$matP_Val(), #is part of combinedDFP_Val_Labels()
                               matP_Val.t = NULL,
                               matLogFC.t = NULL,
                               distMatTraits = NULL,
                               clustResTraits = NULL,
                               traitClusters = NULL,
                               traitClusterMedoids = NULL,
                               traitDendrogram = NULL,
                               traitClustergram = NULL,
                               distMatProbes = NULL,
                               clustResProbes = NULL,
                               probeDendrogram = NULL,
                               DNAdistances = NULL,
                               dfVolcano = NULL,
                               dfKeyShadow = NULL
          )
          if (is.valid(session$userData$sessionVariables$traitReducedDataLogFC())) {
            result$combinedDFP_Val_Labels <- session$userData$sessionVariables$traitReducedDataLogFC()
            result$combinedDFP_Val_Labels$dfLogFC <- na.omit(result$combinedDFP_Val_Labels$dfLogFC)
            #check dim(result$combinedDFP_Val_Labels)
            matLogFC <- na.omit(result$combinedDFP_Val_Labels$dfLogFC) #omit empty rows
            result$matLogFC.t <- t(as.matrix(matLogFC))
            numberCores <- session$userData$numCores
            base::print(paste0(sysTimePID(), " (pReducedDataStructurePVal) before distance matrix (log(FC)) for n(reduced traits) = ", base::nrow(result$matP_Val.t), " (takes some time). Using n(cores) = ", numberCores, "."))
            if (is.valid(result$matLogFC.t)) {
              result$distMatTraits <- getDistMat(numberCores = numberCores, matrix = result$matLogFC.t)
            }
            else {
              result$distMatTraits <- NULL
            }
            base::print(paste0(sysTimePID(), " (traitReducedDataStructure) after distance matrix for reduced traits."))
            #for unknown reason getClustResFast() crashes, if executed without Sys.sleep in advance...
            Sys.sleep(1)
            base::print(paste0(sysTimePID(), " (traitReducedDataStructure) before clustering for traits.", base::nrow(result$matP_Val.t)))
            if (is.valid(result$distMatTraits)) {
              result$clustResTraits <- getClustResFast(result$distMatTraits)
            }
            else {
              result$clustResTraits <- NULL
            }
            base::print(paste0(sysTimePID(), " (traitReducedDataStructure (log(FC))) after clustering results for traits."))
            base::print(base::paste0(sysTimePID(), " (traitReducedDataStructure (log(FC))) start generating traitClusters."))
            numClusters <- length(result$clustResTraits$order)
            if (is.valid(result$clustResTraits) && numClusters > 1) {
              result$traitClusters <- cutree(result$clustResTraits,
                                             k = numClusters)
            }
            else {
              result$traitClusters <- NULL
            }
            base::print(base::paste0(sysTimePID(), " (traitReducedDataStructure (log(FC))) start generating traitClusterMedoids."))
            if (is.valid(result$clustResTraits) && is.valid(result$distMatTraits)) {
              result$traitClusterMedoids <- getTraitClusterMedoids(clustResTraits = result$clustResTraits,
                                                                   distMatTraits = result$distMatTraits,
                                                                   numClusters = numClusters)
            }
            else {
              result$traitClusterMedoids <- NULL
            }

            base::print(base::paste0(sysTimePID(), " (traitReducedDataStructure (log(FC))) start generating traitDendrogram."))
            #do this only, if a dataset is already loaded
            if (is.valid(result$clustResTraits)) {
              result$traitDendrogram <- getDendTraits(clustResTraits = result$clustResTraits, traitClusters = numClusters)
              base::print(base::paste0(sysTimePID(), " (traitReducedDataStructure (log(FC))) after making traitDendrogram."))
            }
            else {
              result$traitDendrogram <- NULL
            }

            base::print(base::paste0(sysTimePID(), " (traitReducedDataStructure (log(FC))) start generating traitClustergram."))
            if (is.valid(result$matLogFC.t) && is.valid(result$clustResTraits)) {
              result$traitClustergram <- getplotClustergramTraitsLong(mat.t = result$matLogFC.t,
                                                                      clustResTraits = result$clustResTraits,
                                                                      traitClusters = numClusters)
            }
            else {
              result$traitClustergram <- NULL
            }
            base::print(base::paste0(sysTimePID(), " (traitReducedDataStructure (log(FC))) start generating distMatProbes."))

            dfLogFC <- na.omit(result$combinedDFP_Val_Labels$dfLogFC) #omit empty rows
            if (is.valid(dfLogFC)) {
              # dflogFC[dflogFC > 0.05] <- NA # 1
              base::print(base::paste0(sysTimePID(), " (traitReducedDataStructure (log(FC))) calculating distance matrix with rows= ", nrow(dfLogFC), " cols= ", ncol(dfLogFC)))
              base::print(base::class(dfLogFC))
              base::print(base::paste0(sysTimePID(), " set missing p-values to 1."))
              # dflogFC[base::is.na(dflogFC)] <- 1 # set missing P_VAL to 1
              base::print(Cstack_info())
              if (base::nrow(dfLogFC) >= 5) {
                numberCores <- session$userData$numCores
                # clustering for rows
                base::print(base::paste0(sysTimePID(), " (traitReducedDataStructure (log(FC))) before distance matrix for n(probes) = ", base::nrow(dfLogFC), " (takes some time). Using n(cores) = ", numberCores, "."))
                gc()
                dfLogFC <- as.matrix(dfLogFC)
                result$distMatProbes <- getDistMat(numberCores = numberCores, matrix = dfLogFC)
                base::print(base::paste0(sysTimePID(), " (traitReducedDataStructure) after distance matrix for probes.", base::nrow(dfLogFC)))
              }
            }
            else {
              result$distMatProbes <- NULL
              browser() # should not happen
            }
            base::print(base::paste0(sysTimePID(), " (traitReducedDataStructure (log(FC)) start clustResProbes."))
            distMat <- result$distMatProbes
            if (is.valid(distMat)) {
              result$clustResProbes <- getClustResFast(distMat)
            }
            else {
              result$clustResProbes <- NULL
              browser() # should not happen
            }

            if (is.valid(result$clustResTraits)) {
              # add "number" and reorder columns; order comes from result$clustResTraits
              result$combinedDFP_Val_Labels$dfP_Val <- result$combinedDFP_Val_Labels$dfP_Val[, result$clustResTraits$order]
              result$combinedDFP_Val_Labels$dfDM <- result$combinedDFP_Val_Labels$dfP_Val[, result$clustResTraits$order]
              result$combinedDFP_Val_Labels$dfN <- result$combinedDFP_Val_Labels$dfN[, result$clustResTraits$order]
              result$combinedDFP_Val_Labels$dfLogFC <- result$combinedDFP_Val_Labels$dfLogFC[, result$clustResTraits$order]
            }
            else {
              browser() #should not happen
            }

            if (is.valid(result$clustResProbes)) {
              result$probeDendrogram <- stats::as.dendrogram(result$clustResProbes)
              #reorder rows; order comes from result$clustResProbes
              result$combinedDFP_Val_Labels$dfP_Val <- result$combinedDFP_Val_Labels$dfP_Val[result$clustResProbes$order, ]
              result$combinedDFP_Val_Labels$dfDM <- result$combinedDFP_Val_Labels$dfP_Val[result$clustResProbes$order, ]
              result$combinedDFP_Val_Labels$dfN <- result$combinedDFP_Val_Labels$dfN[result$clustResProbes$order, ]
              result$combinedDFP_Val_Labels$dfLogFC <- result$combinedDFP_Val_Labels$dfLogFC[result$clustResProbes$order, ]
            }
            else {
              result$probeDendrogram <- NULL
              browser() # should not happen
            }
            result$combinedDFP_Val_Labels$dfP_Val_w_number <- result$combinedDFP_Val_Labels$dfP_Val
            nprobes <- nrow(result$combinedDFP_Val_Labels$dfP_Val_w_number)
            result$combinedDFP_Val_Labels$dfP_Val_w_number$number <- seq(1:nprobes)
            col_order <- c("number", colnames(result$combinedDFP_Val_Labels$dfP_Val_w_number))
            result$combinedDFP_Val_Labels$dfP_Val_w_number <- result$combinedDFP_Val_Labels$dfP_Val_w_number[, col_order]
            result$combinedDFP_Val_Labels$dfP_Val_w_number <- result$combinedDFP_Val_Labels$dfP_Val_w_number[ , -which(colnames(result$combinedDFP_Val_Labels$dfP_Val_w_number) %in% "number.1")]

            result$combinedDFP_Val_Labels$dfDM_w_number <- result$combinedDFP_Val_Labels$dfDM
            nprobes <- nrow(result$combinedDFP_Val_Labels$dfDM_w_number)
            result$combinedDFP_Val_Labels$dfDM_w_number$number <- seq(1:nprobes)
            col_order <- c("number", colnames(result$combinedDFP_Val_Labels$dfDM_w_number))
            result$combinedDFP_Val_Labels$dfDM_w_number <- result$combinedDFP_Val_Labels$dfDM_w_number[, col_order]

            result$combinedDFP_Val_Labels$dfN_w_number <- result$combinedDFP_Val_Labels$dfN
            nprobes <- nrow(result$combinedDFP_Val_Labels$dfN_w_number)
            result$combinedDFP_Val_Labels$dfN_w_number$number <- seq(1:nprobes)
            col_order <- c("number", colnames(result$combinedDFP_Val_Labels$dfN_w_number))
            result$combinedDFP_Val_Labels$dfN_w_number <- result$combinedDFP_Val_Labels$dfN_w_number[, col_order]
            result$combinedDFP_Val_Labels$dfN_w_number <- result$combinedDFP_Val_Labels$dfN_w_number[ , -which(colnames(result$combinedDFP_Val_Labels$dfN_w_number) %in% "number.1")]

            result$combinedDFP_Val_Labels$dfLogFC_w_number <- result$combinedDFP_Val_Labels$dfLogFC
            nprobes <- nrow(result$combinedDFP_Val_Labels$dfLogFC_w_number)
            result$combinedDFP_Val_Labels$dfLogFC_w_number$number <- seq(1:nprobes)
            col_order <- c("number", colnames(result$combinedDFP_Val_Labels$dfLogFC_w_number))
            result$combinedDFP_Val_Labels$dfLogFC_w_number <- result$combinedDFP_Val_Labels$dfLogFC_w_number[, col_order]
            result$combinedDFP_Val_Labels$dfLogFC_w_number <- result$combinedDFP_Val_Labels$dfLogFC_w_number[ , -which(colnames(result$combinedDFP_Val_Labels$dfLogFC_w_number) %in% "number.1")]

            Distance <- input$sld_NumNeighbours
            DNAdistances <- calculateDistanceNeigboursProbes(wd = session$userData$packageWd, clustResProbes = result$clustResProbes, annotation = session$userData$annotation, distanceToLook = Distance, numCores = session$userData$numCores)
            result$DNAdistances <- DNAdistances
            dfLogFC <- result$combinedDFP_Val_Labels$dfLogFC
            if (is.valid(dfLogFC)) {
              dfP_Val <- result$combinedDFP_Val_Labels$dfP_Val
              if (is.valid(dfP_Val)) {
                #take everything into one table with columns p-val and logFC ...
                dfLogFC$probe <- row.names(dfLogFC)
                dfLogFC <- tidyr::pivot_longer(dfLogFC, cols  = -probe, names_to = c("trait"))
                colnames(dfLogFC)[3] <- "LogFC"
                dfP_Val$probe <- row.names(dfP_Val)
                dfP_Val <- tidyr::pivot_longer(dfP_Val, cols  = -probe, names_to = c("trait"))
                colnames(dfP_Val)[3] <- "P_Val"
                dfVolcano <- base::merge(dfP_Val, dfLogFC, by.x = c("probe","trait"), by.y = c("probe","trait"), all.x = TRUE, all.y = FALSE)
                #merge chr and position
                annotation <- base::subset(session$userData$annotation, select = c("name", "chromosome", "position", "gene.symbol"))
                dfVolcano <- base::merge(dfVolcano, annotation, by.x = "probe", by.y = "name", all.x = FALSE, all.y = FALSE)
                #add distances to dfVolcano
                DNAdistances <- result$DNAdistances
                row.names(DNAdistances) <- DNAdistances$ID
                DNAdistances$cg <- DNAdistances$ID #row.names(DNAdistances)
                DNAdistancesNumber <- DNAdistances[,c("cg", "number")]
                DNAdistancesNumber <- tidyr::pivot_longer(DNAdistancesNumber, cols  = -cg, names_to = c("number"))
                DNAdistancesNumber <- DNAdistancesNumber[,c("cg", "value")]
                colnames(DNAdistancesNumber)[2] <- "DistanceNumber"
                DNAdistancesMin <- DNAdistances[,c("cg", "minDistance")]
                DNAdistancesMin <- tidyr::pivot_longer(DNAdistancesMin, cols  = -cg, names_to = c("min"))
                DNAdistancesMin <- DNAdistancesMin[,c("cg", "value")]
                colnames(DNAdistancesMin)[2] <- "DistanceMin"

                DNAdistancesMean <- DNAdistances[,c("cg", "meanDistance")]
                DNAdistancesMean <- tidyr::pivot_longer(DNAdistancesMean, cols  = -cg, names_to = c("mean"))
                DNAdistancesMean <- DNAdistancesMean[,c("cg", "value")]
                colnames(DNAdistancesMean)[2] <- "DistanceMean"

                DNAdistancesMax <- DNAdistances[,c("cg", "maxDistance")]
                DNAdistancesMax <- tidyr::pivot_longer(DNAdistancesMax, cols  = -cg, names_to = c("max"))
                DNAdistancesMax <- DNAdistancesMax[,c("cg", "value")]
                colnames(DNAdistancesMax)[2] <- "DistanceMax"

                dfVolcano <- base::merge(dfVolcano, DNAdistancesNumber, by.x = c("probe"), by.y = c("cg"), all.x = TRUE, all.y = FALSE)
                dfVolcano <- base::merge(dfVolcano, DNAdistancesMin, by.x = c("probe"), by.y = c("cg"), all.x = TRUE, all.y = FALSE)
                dfVolcano <- base::merge(dfVolcano, DNAdistancesMean, by.x = c("probe"), by.y = c("cg"), all.x = TRUE, all.y = FALSE)
                dfVolcano <- base::merge(dfVolcano, DNAdistancesMax, by.x = c("probe"), by.y = c("cg"), all.x = TRUE, all.y = FALSE)
                #create shadow table for key assignment
                #add probe trait and traitsource to key
                originTrait <- result$combinedDFP_Val_Labels$mergedOriginTrait
                originTrait <- rep(originTrait, nprobes)
                keys <- seq(1:nrow(dfVolcano))
                dfKeyShadow <- base::data.frame(key = keys)
                dfKeyShadow$probe <- dfVolcano$probe
                traitLabels <- result$combinedDFP_Val_Labels$mergedOriginalColnames
                dfKeyShadow$trait <- traitLabels #dfVolcano$trait
                traitID <- result$combinedDFP_Val_Labels$traitID
                dfKeyShadow$traitID <- traitID
                dfKeyShadow$traitSource <- originTrait #this one works for clustering by p-value (ln 1580), but not for clustering by log(FC), originTrait passt, dfVolcano ist zu klein 75268...
                dfVolcano$key <- dfKeyShadow$key
                rownames(dfVolcano) <- dfVolcano$key
                rownames(dfKeyShadow) <- dfKeyShadow$key
                #sort by p-val and logFC
                dfVolcano <- dfVolcano[base::order(dfVolcano$P_Val, dfVolcano$LogFC, decreasing = c(FALSE, TRUE), na.last = c(TRUE,TRUE)),]
                #take only first 2^16 entities to be able to plot volcano plot using plotly
                if(nrow(dfVolcano)>2^16) {
                  dfVolcano <- dfVolcano[1:2^16-1,]
                  dfKeyShadow <- dfKeyShadow[1:2^16-1,]
                }
                result$dfVolcano <- dfVolcano
                result$dfKeyShadow <- dfKeyShadow
              }
              else {
                result$dfVolcano <- NULL
                browser() # should not happen
              }
            }
            else {
              result$dfVolcano <- NULL
              browser() # should not happen
            }

          }
        }
      },
      error = function(e) {
        base::message("An error occurred in shiny::reactive(session$userData$sessionVariables$traitReducedDataStructureLogFC):\n", e)
        browser() #should not happen
      },
      warning = function(w) {
        base::message("A warning occurred in shiny::reactive(session$userData$sessionVariables$traitReducedDataStructureLogFC):\n", w)
        browser() #should not happen
      },
      finally = {
        base::print(base::paste0(sysTimePID(), " finished generating session$userData$sessionVariables$traitReducedDataStructureLogFC."))
        return(result)
      }
    )
  })

  # session$userData$sessionVariables$probeReducedDataStructure contains probes that are within a defined range (DNAdistances) around each selected probe
  session$userData$sessionVariables$probeReducedDataStructure <- shiny::reactive({
    id <- shiny::showNotification("Creating probe reduced data structure...", duration = NULL, closeButton = FALSE)
    on.exit(shiny::removeNotification(id), add = TRUE)
    base::tryCatch(
      {
        if (is.valid(session$userData$sessionVariables$traitReducedDataStructurePVal()$combinedDFP_Val_Labels$dfP_Val) || is.valid(session$userData$sessionVariables$traitReducedDataStructureLogFC()$combinedDFP_Val_Labels$dfP_Val)) {
          if (is.valid(session$userData$sessionVariables$traitReducedDataStructurePVal()$combinedDFP_Val_Labels$dfP_Val)) {
            traitReducedDataStructure <- session$userData$sessionVariables$traitReducedDataStructurePVal()
          }
          else if (is.valid(session$userData$sessionVariables$traitReducedDataStructureLogFC()$combinedDFP_Val_Labels$dfP_Val)) {
            traitReducedDataStructure <- session$userData$sessionVariables$traitReducedDataStructureLogFC()
          }
          else {
            browser() #should not happen
          }
          result <- base::list(combinedDFP_Val_Labels = traitReducedDataStructure$combinedDFP_Val_Labels,
                              matP_Val.t = NULL,
                              distMatTraitsP_Val = NULL,
                              clustResTraits = NULL,
                              traitClusters = NULL,
                              traitClusterMedoids = NULL,
                              traitDendrogram = NULL,
                              traitClustergram = NULL,
                              distMatProbes = NULL,
                              clustResProbes = NULL,
                              probeDendrogram = NULL,
                              DNAdistances = NULL,
                              dfVolcano = NULL,
                              dfKeyShadow = NULL
          )
          DNAdistances <- traitReducedDataStructure$DNAdistances
          DNAdistances <- na.omit(DNAdistances)
          result$DNAdistances <- DNAdistances
          dfP_Val <- result$combinedDFP_Val_Labels$dfP_Val
          dfP_Val <- dfP_Val[which(rownames(dfP_Val) %in% DNAdistances$ID), ]
          result$combinedDFP_Val_Labels$dfP_Val <- dfP_Val

          rm(dfP_Val)
          dfDM <- result$combinedDFP_Val_Labels$dfDM
          dfDM <- dfDM[which(rownames(dfDM) %in% DNAdistances$ID), ]
          result$combinedDFP_Val_Labels$dfDM <- dfDM
          rm(dfDM)
          dfLogFC <- result$combinedDFP_Val_Labels$dfLogFC
          #check for missings (should be there) and negative values (too), yes we have negatives and missings here
          dfLogFC <- dfLogFC[which(rownames(dfLogFC) %in% DNAdistances$ID), ]
          result$combinedDFP_Val_Labels$dfLogFC <- dfLogFC
          rm(dfLogFC)
          dfN <- result$combinedDFP_Val_Labels$dfN
          dfN <- dfN[which(rownames(dfN) %in% DNAdistances$ID), ]
          result$combinedDFP_Val_Labels$dfN <- dfN
          rm(dfN)
          # distance matrix for traits based on p-val
          result$matP_Val.t <- t(as.matrix(result$combinedDFP_Val_Labels$dfP_Val))
          numberCores <- session$userData$numCores
          base::print(paste0(sysTimePID(), " (traitReducedDataStructure) before distance matrix for n(reduced traits) = ", base::nrow(result$matP_Val.t), " (takes some time). Using n(cores) = ", numberCores, "."))
          if (is.valid(result$matP_Val.t)) {
            result$distMatTraitsP_Val <- getDistMat(numberCores = numberCores, matrix = result$matP_Val.t)
          }
          else {
            result$distMatTraitsP_Val <- NULL
          }
          # identical (result$distMatTraits, session$userData$sessionVariables$traitReducedDataStructure()$distMatTraits) #they are not identical, but similar
          base::print(paste0(sysTimePID(), " (traitReducedDataStructure) after distance matrix for reduced traits."))
          #for unknown reason getClustResFast() crashes, if executed without Sys.sleep in advance...
          Sys.sleep(1)
          base::print(paste0(sysTimePID(), " (traitReducedDataStructure) before clustering for traits.", base::nrow(result$matP_Val.t)))
          if (is.valid(result$distMatTraitsP_Val)) {
            result$clustResTraits <- getClustResFast(result$distMatTraitsP_Val)
          }
          else {
            result$clustResTraits <- NULL
          }
          base::print(paste0(sysTimePID(), " (traitReducedDataStructure) after clustering results for traits."))
          # identical(result$clustResTraits, session$userData$sessionVariables$traitReducedDataStructure()$clustResTraits) #they are not identical, but similar
          base::print(base::paste0(sysTimePID(), " (traitReducedDataStructure) start generating traitClusters."))
          numClusters <- length(result$clustResTraits$order)
          if (is.valid(result$clustResTraits) && numClusters > 1) {
            result$traitClusters <- cutree(result$clustResTraits,
                                           k = numClusters)
          }
          else {
            result$traitClusters <- NULL
          }
          base::print(base::paste0(sysTimePID(), " (traitReducedDataStructure) start generating traitClusterMedoids."))
          if (is.valid(result$clustResTraits) && is.valid(result$distMatTraitsP_Val)) {
            result$traitClusterMedoids <- getTraitClusterMedoids(clustResTraits = result$clustResTraits,
                                                                 distMatTraits = result$distMatTraitsP_Val,
                                                                 numClusters = numClusters)
          }
          else {
            result$traitClusterMedoids <- NULL
          }

          base::print(base::paste0(sysTimePID(), " (traitReducedDataStructure) start generating traitDendrogram."))
          #do this only, if a dataset is already loaded
          if (is.valid(result$clustResTraits)) {
            result$traitDendrogram <- getDendTraits(clustResTraits = result$clustResTraits, traitClusters = numClusters)
            base::print(base::paste0(sysTimePID(), " (traitReducedDataStructure) after making traitDendrogram."))
          }
          else {
            result$traitDendrogram <- NULL
          }

          base::print(base::paste0(sysTimePID(), " (traitReducedDataStructure) start generating traitClustergram."))
          if (is.valid(result$matP_Val.t) && is.valid(result$clustResTraits)) {
            result$traitClustergram <- getplotClustergramTraitsLong(mat.t = result$matP_Val.t,
                                                                    clustResTraits = result$clustResTraits,
                                                                    traitClusters = numClusters)
          }
          else {
            result$traitClustergram <- NULL
          }
          # add "number" and reorder columns
          result$combinedDFP_Val_Labels$dfP_Val_w_number <- result$combinedDFP_Val_Labels$dfP_Val
          nprobes <- nrow(result$combinedDFP_Val_Labels$dfP_Val_w_number)
          result$combinedDFP_Val_Labels$dfP_Val_w_number$number <- seq(1:nprobes)
          col_order <- c("number", colnames(result$combinedDFP_Val_Labels$dfP_Val_w_number))
          result$combinedDFP_Val_Labels$dfP_Val_w_number <- result$combinedDFP_Val_Labels$dfP_Val_w_number[, col_order]
          result$combinedDFP_Val_Labels$dfP_Val_w_number <- result$combinedDFP_Val_Labels$dfP_Val_w_number[ , -which(colnames(result$combinedDFP_Val_Labels$dfP_Val_w_number) %in% "number.1")]

          result$combinedDFP_Val_Labels$dfDM_w_number <- result$combinedDFP_Val_Labels$dfDM
          nprobes <- nrow(result$combinedDFP_Val_Labels$dfDM_w_number)
          result$combinedDFP_Val_Labels$dfDM_w_number$number <- seq(1:nprobes)
          col_order <- c("number", colnames(result$combinedDFP_Val_Labels$dfDM_w_number))
          result$combinedDFP_Val_Labels$dfDM_w_number <- result$combinedDFP_Val_Labels$dfDM_w_number[, col_order]

          result$combinedDFP_Val_Labels$dfLogFC_w_number <- result$combinedDFP_Val_Labels$dfLogFC
          nprobes <- nrow(result$combinedDFP_Val_Labels$dfLogFC_w_number)
          result$combinedDFP_Val_Labels$dfLogFC_w_number$number <- seq(1:nprobes) #this line crashes...
          col_order <- c("number", colnames(result$combinedDFP_Val_Labels$dfLogFC_w_number))
          result$combinedDFP_Val_Labels$dfLogFC_w_number <- result$combinedDFP_Val_Labels$dfLogFC_w_number[, col_order]
          result$combinedDFP_Val_Labels$dfLogFC_w_number <- result$combinedDFP_Val_Labels$dfLogFC_w_number[ , -which(colnames(result$combinedDFP_Val_Labels$dfLogFC_w_number) %in% "number.1")]

          result$combinedDFP_Val_Labels$dfN_w_number <- result$combinedDFP_Val_Labels$dfN
          nprobes <- nrow(result$combinedDFP_Val_Labels$dfN_w_number)
          result$combinedDFP_Val_Labels$dfN_w_number$number <- seq(1:nprobes)
          col_order <- c("number", colnames(result$combinedDFP_Val_Labels$dfN_w_number))
          result$combinedDFP_Val_Labels$dfN_w_number <- result$combinedDFP_Val_Labels$dfN_w_number[, col_order]
          result$combinedDFP_Val_Labels$dfN_w_number <- result$combinedDFP_Val_Labels$dfN_w_number[ , -which(colnames(result$combinedDFP_Val_Labels$dfN_w_number) %in% "number.1")]

          result$dfVolcano <- traitReducedDataStructure$dfVolcano
          result$dfKeyShadow <- traitReducedDataStructure$dfKeyShadow
        }
      },
      error = function(e) {
        base::message("An error occurred in shiny::reactive(session$userData$sessionVariables$probeReducedDataStructure):\n", e)
        browser() #should not happen
      },
      warning = function(w) {
        base::message("A warning occurred in shiny::reactive(session$userData$sessionVariables$probeReducedDataStructure):\n", w)
        browser() #should not happen
      },
      finally = {
        base::print(base::paste0(sysTimePID(), " finished generating session$userData$sessionVariables$probeReducedGeneralDataStructure.\n"))
        return(result)
      }
    )
  })

#   session$userData$sessionVariables$distanceMultipliedTraitReducedDataStructurePVal <- shiny::reactive({
#     id <- shiny::showNotification("Creating distance multiplied trait reduced data structure p-val...", duration = NULL, closeButton = FALSE)
#     on.exit(shiny::removeNotification(id), add = TRUE)
# browser()
#     base::tryCatch(
#       {
#         if (is.valid(session$userData$sessionVariables$traitReducedDataStructurePVal()$combinedDFP_Val_Labels$dfP_Val)) {
#           traitReducedDataStructure <- session$userData$sessionVariables$traitReducedDataStructurePVal()
#           result <- base::list(combinedDFP_Val_Labels = traitReducedDataStructure$combinedDFP_Val_Labels,
#                               matP_Val.t = NULL,
#                               distMatTraits = NULL,
#                               clustResTraits = NULL,
#                               traitClusters = NULL,
#                               traitClusterMedoids = NULL,
#                               traitDendrogram = NULL,
#                               traitClustergram = NULL,
#                               distMatProbes = NULL,
#                               clustResProbes = NULL,
#                               probeDendrogram = NULL,
#                               DNAdistances = NULL
#           )
# #          Distance <- input$sld_NumNeighbours
# #          DNAdistances <- calculateDistanceNeigboursProbes(wd = session$userData$packageWd, clustResProbes = session$userData$sessionVariables$traitReducedDataStructure()$clustResProbes, annotation = session$userData$annotation, distanceToLook = Distance, numCores = session$userData$numCores)
#           DNAdistances <- session$userData$sessionVariables$traitReducedDataStructurePVal()$DNAdistances
#           if(is.valid(DNAdistances)) {
#             result$DNAdistances <- DNAdistances
#             dfP_Val <- result$combinedDFP_Val_Labels$dfP_Val
#             dfP_Val <- dfP_Val[which(rownames(dfP_Val) %in% DNAdistances$ID), ]
#             if(nrow(dfP_Val)>0) {
#               rownames <- rownames(dfP_Val)
#               result$combinedDFP_Val_Labels$dfP_Val <- dfP_Val
#               matP_Val <- base::as.matrix(dfP_Val)
#               dt <- data.table::data.table(matP_Val)
#               # build set union of vec and dt
#               DNAdistances <- DNAdistances[which(DNAdistances$ID %in% rownames(dfP_Val)), ]
#               vec <- DNAdistances$meanDistance # take means# or
#               # vec <- DNAdistances$number # take numbers
#               # normalize vec to -1...1
#               vec <-  scales::rescale(vec, to = c(-1, 1))
#               # invert vec, so that small distances become large multiplies
#               vec <- vec * -1
#
#               result$combinedDFP_Val_Labels$dfP_Val <- data.table::data.table(t(t(dt) * vec))
#               rownames(result$combinedDFP_Val_Labels$dfP_Val) <- rownames
#               rm(dt)
#             }
#             else {
#               base::print(base::paste0(sysTimePID(), " no set union between DNAdistances and dfP_Val in session$userData$sessionVariables$distanceMultipliedTraitReducedDataStructurePVal\n"))
#               result$combinedDFP_Val_Labels$dfP_Val <- result$combinedDFP_Val_Labels$dfP_Val
#             }
# # browser()
# #             #build subset
# #             rowlabels <- rownames(result$combinedDFP_Val_Labels$dfP_Val)
# #             columnlabels <- colnames(result$combinedDFP_Val_Labels$dfP_Val) #tbc() apply for all the other df's; notice, that if no results inside DW results, then there is no union set between neighbouring cpg and prior selection (distance to look too small)
# #             dfDM <- session$userData$sessionVariables$traitReducedDataStructurePVal()$combinedDFP_Val_Labels$dfDM
# #             dfDM <- dfDM[which(rownames(dfDM) %in% rowlabels), which(colnames(dfDM) %in% columnlabels)]
# #             result$combinedDFP_Val_Labels$dfDM <- dfDM
# #             result$combinedDFP_Val_Labels$dfDM
# #             dfN <- session$userData$sessionVariables$traitReducedDataStructurePVal()$combinedDFP_Val_Labels$dfN
# #             dfN <- dfN[which(rownames(dfN) %in% rowlabels), which(colnames(dfN) %in% columnlabels)]
# #             result$combinedDFP_Val_Labels$dfN <- dfN
#
#             # calculate transposed matP_Val and clustering
#             matP_Val <- na.omit(result$combinedDFP_Val_Labels$dfP_Val) #omit empty rows
#             result$matP_Val.t <- t(as.matrix(matP_Val))
#             numberCores <- session$userData$numCores
#             base::print(paste0(sysTimePID(), " (distanceMultipliedTraitReducedDataStructurePVal) before distance matrix (log(FC)) for n(reduced traits) = ", base::nrow(result$matP_Val.t), " (takes some time). Using n(cores) = ", numberCores, "."))
#             if (!is.null(result$matP_Val.t)) {
#               result$distMatTraits <- getDistMat(numberCores = numberCores, matrix = result$matP_Val.t)
#             }
#             else {
#               result$distMatTraits <- NULL
#             }
#             base::print(paste0(sysTimePID(), " (distanceMultipliedTraitReducedDataStructurePVal) after distance matrix (log(FC)) for reduced traits."))
#
#             #calculate result$traitClusters and beyond from above clustering results
#             #for unknown reason getClustResFast() crashes, if executed without Sys.sleep in advance...
#             Sys.sleep(1)
#             base::print(paste0(sysTimePID(), " (distanceMultipliedTraitReducedDataStructurePVal) before clustering for traits.", base::nrow(result$matP_Val.t)))
# browser() #tbc() check is.valid -> yes, it has too many missings, if neighbours window is too narrow ; therefore no clustering
#             if (is.valid(result$distMatTraits)) {
#               result$clustResTraits <- getClustResFast(result$distMatTraits)
#             }
#             else {
#               result$clustResTraits <- NULL
#             }
#             base::print(paste0(sysTimePID(), " (distanceMultipliedTraitReducedDataStructurePVal) after clustering results for traits."))
#             base::print(base::paste0(sysTimePID(), " (distanceMultipliedTraitReducedDataStructurePVal) start generating traitClusters."))
#             numClusters <- length(result$clustResTraits$order)
#             if (is.valid(result$clustResTraits) && numClusters > 1) {
#               result$traitClusters <- cutree(result$clustResTraits,
#                                              k = numClusters)
#             }
#             else {
#               result$traitClusters <- NULL
#             }
#             base::print(base::paste0(sysTimePID(), " (distanceMultipliedTraitReducedDataStructurePVal) start generating traitClusterMedoids."))
#             if (is.valid(result$clustResTraits) && is.valid(result$distMatTraits)) {
#               result$traitClusterMedoids <- getTraitClusterMedoids(clustResTraits = result$clustResTraits,
#                                                                    distMatTraits = result$distMatTraits,
#                                                                    numClusters = numClusters)
#             }
#             else {
#               result$traitClusterMedoids <- NULL
#             }
#
#             base::print(base::paste0(sysTimePID(), " (distanceMultipliedTraitReducedDataStructurePVal) start generating traitDendrogram."))
#             #do this only, if a dataset is already loaded
#             if (is.valid(result$clustResTraits)) {
#               result$traitDendrogram <- getDendTraits(clustResTraits = result$clustResTraits, traitClusters = numClusters)
#               base::print(base::paste0(sysTimePID(), " (distanceMultipliedTraitReducedDataStructurePVal) after making traitDendrogram."))
#             }
#             else {
#               result$traitDendrogram <- NULL
#             }
#
#             base::print(base::paste0(sysTimePID(), " (distanceMultipliedTraitReducedDataStructurePVal) start generating traitClustergram."))
#             if (is.valid(result$matP_Val.t) && is.valid(result$clustResTraits)) {
#               result$traitClustergram <- getplotClustergramTraitsLong(mat.t = result$matP_Val.t,
#                                                                       clustResTraits = result$clustResTraits,
#                                                                       traitClusters = numClusters)
#             }
#             else {
#               result$traitClustergram <- NULL
#             }
#
#             result$combinedDFP_Val_Labels$dfP_Val_w_number <- result$combinedDFP_Val_Labels$dfP_Val
#             nprobes <- nrow(result$combinedDFP_Val_Labels$dfP_Val_w_number)
#             result$combinedDFP_Val_Labels$dfP_Val_w_number$number <- seq(1:nprobes)
#             col_order <- c("number", colnames(result$combinedDFP_Val_Labels$dfP_Val_w_number))
#             result$combinedDFP_Val_Labels$dfP_Val_w_number <- result$combinedDFP_Val_Labels$dfP_Val_w_number[, col_order]
#             result$combinedDFP_Val_Labels$dfP_Val_w_number <- result$combinedDFP_Val_Labels$dfP_Val_w_number[ , -which(colnames(result$combinedDFP_Val_Labels$dfP_Val_w_number) %in% "number.1")]
#
#             result$combinedDFP_Val_Labels$dfDM_w_number <- result$combinedDFP_Val_Labels$dfDM
#             nprobes <- nrow(result$combinedDFP_Val_Labels$dfDM_w_number)
#             result$combinedDFP_Val_Labels$dfDM_w_number$number <- seq(1:nprobes)
#             col_order <- c("number", colnames(result$combinedDFP_Val_Labels$dfDM_w_number))
#             result$combinedDFP_Val_Labels$dfDM_w_number <- result$combinedDFP_Val_Labels$dfDM_w_number[, col_order]
#             #              result$dfDM_w_number <- result$dfDM[ , -which(colnames(result$dfDM_w_number) %in% "number.1")]
#
#             result$combinedDFP_Val_Labels$dfLogFC_w_number <- result$combinedDFP_Val_Labels$dfLogFC
#             nprobes <- nrow(result$combinedDFP_Val_Labels$dfLogFC_w_number)
#             result$combinedDFP_Val_Labels$dfLogFC_w_number$number <- seq(1:nprobes)
#             col_order <- c("number", colnames(result$combinedDFP_Val_Labels$dfLogFC_w_number))
#             result$combinedDFP_Val_Labels$dfLogFC_w_number <- result$combinedDFP_Val_Labels$dfLogFC_w_number[, col_order]
#             result$combinedDFP_Val_Labels$dfLogFC_w_number <- result$combinedDFP_Val_Labels$dfLogFC_w_number[ , -which(colnames(result$combinedDFP_Val_Labels$dfLogFC_w_number) %in% "number.1")]
#
#             result$combinedDFP_Val_Labels$dfN_w_number <- result$combinedDFP_Val_Labels$dfN
#             nprobes <- nrow(result$combinedDFP_Val_Labels$dfN_w_number)
#             result$combinedDFP_Val_Labels$dfN_w_number$number <- seq(1:nprobes)
#             col_order <- c("number", colnames(result$combinedDFP_Val_Labels$dfN_w_number))
#             result$combinedDFP_Val_Labels$dfN_w_number <- result$combinedDFP_Val_Labels$dfN_w_number[, col_order]
#             result$combinedDFP_Val_Labels$dfN_w_number <- result$combinedDFP_Val_Labels$dfN_w_number[ , -which(colnames(result$combinedDFP_Val_Labels$dfN_w_number) %in% "number.1")]
#
#             base::print(base::paste0(sysTimePID(), " (distanceMultipliedTraitReducedDataStructurePVal) start generating distMatProbes."))
#             dfP_Val <- result$combinedDFP_Val_Labels$dfP_Val
#             if (is.valid(dfP_Val)) {
#               dfP_Val[dfP_Val > 0.05] <- NA # 1
#               base::print(base::paste0(sysTimePID(), " (distanceMultipliedTraitReducedDataStructurePVal) calculating distance matrix with rows= ", nrow(dfP_Val), " cols= ", ncol(dfP_Val)))
#               base::print(base::class(dfP_Val))
#               base::print(base::paste0(sysTimePID(), " set missing p-values to 1."))
#               dfP_Val[base::is.na(dfP_Val)] <- 1 # set missing P_VAL to 1
#               base::print(Cstack_info())
#               if (base::nrow(dfP_Val) >= 5) {
#                 numberCores <- session$userData$numCores
#                 # clustering for rows
#                 base::print(base::paste0(sysTimePID(), " (distanceMultipliedTraitReducedDataStructurePVal) before distance matrix for n(probes) = ", base::nrow(dfP_Val), " (takes some time). Using n(cores) = ", numberCores, "."))
#                 gc()
#                 result$distMatProbes <- getDistMat(numberCores = numberCores, matrix = dfP_Val)
#                 base::print(base::paste0(sysTimePID(), " (distanceMultipliedTraitReducedDataStructurePVal) after distance matrix for probes.", base::nrow(dfP_Val)))
#               }
#             }
#             else {
#               result$distMatProbes <- NULL
#               browser() # should not happen
#             }
#             base::print(base::paste0(sysTimePID(), " (distanceMultipliedTraitReducedDataStructurePVal) start clustResProbes."))
#             distMat <- result$distMatProbes
#             if (is.valid(distMat)) {
#               result$clustResProbes <- getClustResFast(distMat)
#             }
#             else {
#               result$clustResProbes <- NULL
#               browser() # should not happen
#             }
#             if (is.valid(result$clustResTraits)) {
#               # add "number" and reorder columns; order comes from result$clustResTraits
#               result$combinedDFP_Val_Labels$dfP_Val <- result$combinedDFP_Val_Labels$dfP_Val[, result$clustResTraits$order]
#               result$combinedDFP_Val_Labels$dfDM <- result$combinedDFP_Val_Labels$dfP_Val[, result$clustResTraits$order]
#               result$combinedDFP_Val_Labels$dfN <- result$combinedDFP_Val_Labels$dfN[, result$clustResTraits$order]
#               result$combinedDFP_Val_Labels$dfLogFC <- result$combinedDFP_Val_Labels$dfLogFC[, result$clustResTraits$order]
#             }
#             else {
#               browser() #should not happen
#             }
#
#             if (is.valid(result$clustResProbes)) {
#               result$probeDendrogram <- stats::as.dendrogram(result$clustResProbes)
#               #reorder rows; order comes from result$clustResProbes
#               result$combinedDFP_Val_Labels$dfP_Val <- result$combinedDFP_Val_Labels$dfP_Val[result$clustResProbes$order, ]
#               result$combinedDFP_Val_Labels$dfDM <- result$combinedDFP_Val_Labels$dfP_Val[result$clustResProbes$order, ]
#               result$combinedDFP_Val_Labels$dfN <- result$combinedDFP_Val_Labels$dfN[result$clustResProbes$order, ]
#               result$combinedDFP_Val_Labels$dfLogFC <- result$combinedDFP_Val_Labels$dfLogFC[result$clustResProbes$order, ]
#             }
#             else {
#               result$probeDendrogram <- NULL
#               browser() # should not happen
#             }
#             # result$distMatProbes <- session$userData$sessionVariables$traitReducedDataStructurePVal()$distMatProbes
#             # result$clustResProbes <- session$userData$sessionVariables$traitReducedDataStructurePVal()$clustResProbes
#             # result$probeDendrogram <- session$userData$sessionVariables$traitReducedDataStructurePVal()$probeDendrogram
#           }
#           else {
#             result <- NULL
#           }
#         }
#       },
#       error = function(e) {
#         base::message("An error occurred in shiny::reactive(session$userData$sessionVariables$distanceMultipliedTraitReducedDataStructurePVal):\n", e)
#         browser() #should not happen
#       },
#       warning = function(w) {
#         base::message("A warning occurred in shiny::reactive(session$userData$sessionVariables$distanceMultipliedTraitReducedDataStructurePVal):\n", w)
#         browser() #should not happen
#       },
#       finally = {
#         base::print(base::paste0(sysTimePID(), " finished generating session$userData$sessionVariables$distanceMultipliedTraitReducedDataStructurePVal\n"))
#         return(result)
#       }
#     )
#   })
#
#   session$userData$sessionVariables$distanceMultipliedTraitReducedDataStructureLogFC <- shiny::reactive({
#     id <- shiny::showNotification("Creating distance multiplied trait reduced data structure log(FC)...", duration = NULL, closeButton = FALSE)
#     on.exit(shiny::removeNotification(id), add = TRUE)
#     browser()
#     base::tryCatch(
#       {
#         if (is.valid(session$userData$sessionVariables$traitReducedDataStructureLogFC()$combinedDFP_Val_Labels$dfP_Val)) {
#           traitReducedDataStructure <- session$userData$sessionVariables$traitReducedDataStructurePVal()
#           result <- base::list(combinedDFP_Val_Labels = traitReducedDataStructure$combinedDFP_Val_Labels,
#                                matP_Val.t = NULL,
#                                distMatTraits = NULL,
#                                clustResTraits = NULL,
#                                traitClusters = NULL,
#                                traitClusterMedoids = NULL,
#                                traitDendrogram = NULL,
#                                traitClustergram = NULL,
#                                distMatProbes = NULL,
#                                clustResProbes = NULL,
#                                probeDendrogram = NULL,
#                                DNAdistances = NULL
#           )
#           #          Distance <- input$sld_NumNeighbours
#           #          DNAdistances <- calculateDistanceNeigboursProbes(wd = session$userData$packageWd, clustResProbes = session$userData$sessionVariables$traitReducedDataStructure()$clustResProbes, annotation = session$userData$annotation, distanceToLook = Distance, numCores = session$userData$numCores)
#           DNAdistances <- session$userData$sessionVariables$traitReducedDataStructurePVal()$DNAdistances
#           if(is.valid(DNAdistances)) {
#             result$DNAdistances <- DNAdistances
#             dfLogFC <- result$combinedDFP_Val_Labels$dfLogFC
#             dfLogFC <- dfLogFC[which(rownames(dfLogFC) %in% DNAdistances$ID), ]
#             if(nrow(dfLogFC)>0) {
#               rownames <- rownames(dfLogFC)
#               matLogFC <- base::as.matrix(dfLogFC)
#               dt <- data.table::data.table(matLogFC)
#               # build set union of vec and dt
#               DNAdistances <- DNAdistances[which(DNAdistances$ID %in% rownames(dfLogFC)), ]
#               vec <- DNAdistances$meanDistance # take means# or
#               # vec <- DNAdistances$number # take numbers
#               # normalize vec to -1...1
#               vec <-  scales::rescale(vec, to = c(-1, 1))
#               # invert vec, so that small distances become large multiplies
#               vec <- vec * -1
#
#               result$combinedDFP_Val_Labels$dfLogFC <- data.table::data.table(t(t(dt) * vec))
#               rownames(result$combinedDFP_Val_Labels$dfLogFC) <- rownames
#               rm(dt)
#             }
#             else {
#               base::print(base::paste0(sysTimePID(), " no set union between DNAdistances and dfP_Val in session$userData$sessionVariables$distanceMultipliedTraitReducedDataStructurePVal\n"))
#               result$combinedDFP_Val_Labels$dfLogFC <- result$combinedDFP_Val_Labels$dfLogFC
#             }
#             # calculate transposed LogFC
#             matLogFC <- na.omit(result$combinedDFP_Val_Labels$dfLogFC) #omit empty rows
#             result$matLogFC.t <- t(as.matrix(matLogFC))
#             numberCores <- session$userData$numCores
#             base::print(paste0(sysTimePID(), " (pReducedDataStructureLogFC) before distance matrix (log(FC)) for n(reduced traits) = ", base::nrow(result$matP_Val.t), " (takes some time). Using n(cores) = ", numberCores, "."))
#             if (!is.null(result$matLogFC.t)) {
#               result$distMatTraits <- getDistMat(numberCores = numberCores, matrix = result$matLogFC.t)
#             }
#             else {
#               result$distMatTraits <- NULL
#             }
#             base::print(paste0(sysTimePID(), " (distanceMultipliedTraitReducedDataStructureLogFC) after distance matrix (log(FC)) for reduced traits."))
#             #calculate result$traitClusters and beyond from above clustering results
#             #for unknown reason getClustResFast() crashes, if executed without Sys.sleep in advance...
#             Sys.sleep(1)
#             base::print(paste0(sysTimePID(), " (distanceMultipliedTraitReducedDataStructureLogFC) before clustering for traits.", base::nrow(result$matP_Val.t)))
# browser() #tbc() check is.valid
#             if (is.valid(result$distMatTraits)) {
#               result$clustResTraits <- getClustResFast(result$distMatTraits)
#             }
#             else {
#               result$clustResTraits <- NULL
#             }
#             base::print(paste0(sysTimePID(), " (distanceMultipliedTraitReducedDataStructureLogFC) after clustering results for traits."))
#             base::print(base::paste0(sysTimePID(), " (distanceMultipliedTraitReducedDataStructureLogFC) start generating traitClusters."))
#             numClusters <- length(result$clustResTraits$order)
#             if (is.valid(result$clustResTraits) && numClusters > 1) {
#               result$traitClusters <- cutree(result$clustResTraits,
#                                              k = numClusters)
#             }
#             else {
#               result$traitClusters <- NULL
#             }
#             base::print(base::paste0(sysTimePID(), " (distanceMultipliedTraitReducedDataStructureLogFC) start generating traitClusterMedoids."))
#             if (is.valid(result$clustResTraits) && is.valid(result$distMatTraits)) {
#               result$traitClusterMedoids <- getTraitClusterMedoids(clustResTraits = result$clustResTraits,
#                                                                    distMatTraits = result$distMatTraits,
#                                                                    numClusters = numClusters)
#             }
#             else {
#               result$traitClusterMedoids <- NULL
#             }
#
#             base::print(base::paste0(sysTimePID(), " (distanceMultipliedTraitReducedDataStructureLogFC) start generating traitDendrogram."))
#             #do this only, if a dataset is already loaded
#             if (is.valid(result$clustResTraits)) {
#               result$traitDendrogram <- getDendTraits(clustResTraits = result$clustResTraits, traitClusters = numClusters)
#               base::print(base::paste0(sysTimePID(), " (distanceMultipliedTraitReducedDataStructureLogFC) after making traitDendrogram."))
#             }
#             else {
#               result$traitDendrogram <- NULL
#             }
#
#             base::print(base::paste0(sysTimePID(), " (distanceMultipliedTraitReducedDataStructureLogFC) start generating traitClustergram."))
#             if (is.valid(result$matP_Val.t) && is.valid(result$clustResTraits)) {
#               result$traitClustergram <- getplotClustergramTraitsLong(mat.t = result$matP_Val.t,
#                                                                       clustResTraits = result$clustResTraits,
#                                                                       traitClusters = numClusters)
#             }
#             else {
#               result$traitClustergram <- NULL
#             }
#             # result$traitClusters <- session$userData$sessionVariables$traitReducedDataStructurePVal()$traitClusters
#             # result$traitClusterMedoids <- session$userData$sessionVariables$traitReducedDataStructurePVal()$traitClusterMedoids
#             # result$traitDendrogram <- session$userData$sessionVariables$traitReducedDataStructurePVal()$traitDendrogram
#             # result$traitClustergram <- session$userData$sessionVariables$traitReducedDataStructurePVal()$traitClustergram
#
#             result$combinedDFP_Val_Labels$dfP_Val_w_number <- result$combinedDFP_Val_Labels$dfP_Val
#             nprobes <- nrow(result$combinedDFP_Val_Labels$dfP_Val_w_number)
#             result$combinedDFP_Val_Labels$dfP_Val_w_number$number <- seq(1:nprobes)
#             col_order <- c("number", colnames(result$combinedDFP_Val_Labels$dfP_Val_w_number))
#             result$combinedDFP_Val_Labels$dfP_Val_w_number <- result$combinedDFP_Val_Labels$dfP_Val_w_number[, col_order]
#             result$combinedDFP_Val_Labels$dfP_Val_w_number <- result$combinedDFP_Val_Labels$dfP_Val_w_number[ , -which(colnames(result$combinedDFP_Val_Labels$dfP_Val_w_number) %in% "number.1")]
#
#             result$combinedDFP_Val_Labels$dfDM_w_number <- result$combinedDFP_Val_Labels$dfDM
#             nprobes <- nrow(result$combinedDFP_Val_Labels$dfDM_w_number)
#             result$combinedDFP_Val_Labels$dfDM_w_number$number <- seq(1:nprobes)
#             col_order <- c("number", colnames(result$combinedDFP_Val_Labels$dfDM_w_number))
#             result$combinedDFP_Val_Labels$dfDM_w_number <- result$combinedDFP_Val_Labels$dfDM_w_number[, col_order]
#             #              result$dfDM_w_number <- result$dfDM[ , -which(colnames(result$dfDM_w_number) %in% "number.1")]
#
#             result$combinedDFP_Val_Labels$dfLogFC_w_number <- result$combinedDFP_Val_Labels$dfLogFC
#             nprobes <- nrow(result$combinedDFP_Val_Labels$dfLogFC_w_number)
#             result$combinedDFP_Val_Labels$dfLogFC_w_number$number <- seq(1:nprobes)
#             col_order <- c("number", colnames(result$combinedDFP_Val_Labels$dfLogFC_w_number))
#             result$combinedDFP_Val_Labels$dfLogFC_w_number <- result$combinedDFP_Val_Labels$dfLogFC_w_number[, col_order]
#             result$combinedDFP_Val_Labels$dfLogFC_w_number <- result$combinedDFP_Val_Labels$dfLogFC_w_number[ , -which(colnames(result$combinedDFP_Val_Labels$dfLogFC_w_number) %in% "number.1")]
#
#             result$combinedDFP_Val_Labels$dfN_w_number <- result$combinedDFP_Val_Labels$dfN
#             nprobes <- nrow(result$combinedDFP_Val_Labels$dfN_w_number)
#             result$combinedDFP_Val_Labels$dfN_w_number$number <- seq(1:nprobes)
#             col_order <- c("number", colnames(result$combinedDFP_Val_Labels$dfN_w_number))
#             result$combinedDFP_Val_Labels$dfN_w_number <- result$combinedDFP_Val_Labels$dfN_w_number[, col_order]
#             result$combinedDFP_Val_Labels$dfN_w_number <- result$combinedDFP_Val_Labels$dfN_w_number[ , -which(colnames(result$combinedDFP_Val_Labels$dfN_w_number) %in% "number.1")]
#
#             base::print(base::paste0(sysTimePID(), " (distanceMultipliedTraitReducedDataStructureLogFC) start generating distMatProbes."))
#             dfP_Val <- result$combinedDFP_Val_Labels$dfP_Val
#             if (is.valid(dfP_Val)) {
#               dfP_Val[dfP_Val > 0.05] <- NA # 1
#               base::print(base::paste0(sysTimePID(), " (distanceMultipliedTraitReducedDataStructureLogFC) calculating distance matrix with rows= ", nrow(dfP_Val), " cols= ", ncol(dfP_Val)))
#               base::print(base::class(dfP_Val))
#               base::print(base::paste0(sysTimePID(), " set missing p-values to 1."))
#               dfP_Val[base::is.na(dfP_Val)] <- 1 # set missing P_VAL to 1
#               base::print(Cstack_info())
#               if (base::nrow(dfP_Val) >= 5) {
#                 numberCores <- session$userData$numCores
#                 # clustering for rows
#                 base::print(base::paste0(sysTimePID(), " (distanceMultipliedTraitReducedDataStructureLogFC) before distance matrix for n(probes) = ", base::nrow(dfP_Val), " (takes some time). Using n(cores) = ", numberCores, "."))
#                 gc()
#                 result$distMatProbes <- getDistMat(numberCores = numberCores, matrix = dfP_Val)
#                 base::print(base::paste0(sysTimePID(), " (distanceMultipliedTraitReducedDataStructurePVal) after distance matrix for probes.", base::nrow(dfP_Val)))
#               }
#             }
#             else {
#               result$distMatProbes <- NULL
#               browser() # should not happen
#             }
#             base::print(base::paste0(sysTimePID(), " (distanceMultipliedTraitReducedDataStructureLogFC) start clustResProbes."))
#             distMat <- result$distMatProbes
#             if (is.valid(distMat)) {
#               result$clustResProbes <- getClustResFast(distMat)
#             }
#             else {
#               result$clustResProbes <- NULL
#               browser() # should not happen
#             }
#             if (is.valid(result$clustResTraits)) {
#               # add "number" and reorder columns; order comes from result$clustResTraits
#               result$combinedDFP_Val_Labels$dfP_Val <- result$combinedDFP_Val_Labels$dfP_Val[, result$clustResTraits$order]
#               result$combinedDFP_Val_Labels$dfDM <- result$combinedDFP_Val_Labels$dfP_Val[, result$clustResTraits$order]
#               result$combinedDFP_Val_Labels$dfN <- result$combinedDFP_Val_Labels$dfN[, result$clustResTraits$order]
#               result$combinedDFP_Val_Labels$dfLogFC <- result$combinedDFP_Val_Labels$dfLogFC[, result$clustResTraits$order]
#             }
#             else {
#               browser() #should not happen; happens mostly in debug mode, if there are too few valid cases
#             }
#
#             if (is.valid(result$clustResProbes)) {
#               result$probeDendrogram <- stats::as.dendrogram(result$clustResProbes)
#               #reorder rows; order comes from result$clustResProbes
#               result$combinedDFP_Val_Labels$dfP_Val <- result$combinedDFP_Val_Labels$dfP_Val[result$clustResProbes$order, ]
#               result$combinedDFP_Val_Labels$dfDM <- result$combinedDFP_Val_Labels$dfP_Val[result$clustResProbes$order, ]
#               result$combinedDFP_Val_Labels$dfN <- result$combinedDFP_Val_Labels$dfN[result$clustResProbes$order, ]
#               result$combinedDFP_Val_Labels$dfLogFC <- result$combinedDFP_Val_Labels$dfLogFC[result$clustResProbes$order, ]
#             }
#             else {
#               result$probeDendrogram <- NULL
#               browser() # should not happen
#             }
#
#             # result$distMatProbes <- session$userData$sessionVariables$traitReducedDataStructurePVal()$distMatProbes
#             # result$clustResProbes <- session$userData$sessionVariables$traitReducedDataStructurePVal()$clustResProbes
#             # result$probeDendrogram <- session$userData$sessionVariables$traitReducedDataStructurePVal()$probeDendrogram
#           }
#           else {
#             result <- NULL
#           }
#         }
#       },
#       error = function(e) {
#         base::message("An error occurred in shiny::reactive(session$userData$sessionVariables$distanceMultipliedTraitReducedDataStructureLogFC):\n", e)
#         browser() #should not happen
#       },
#       warning = function(w) {
#         base::message("A warning occurred in shiny::reactive(session$userData$sessionVariables$distanceMultipliedTraitReducedDataStructureLogFC):\n", w)
#         browser() #should not happen
#       },
#       finally = {
#         base::print(base::paste0(sysTimePID(), " finished generating session$userData$sessionVariables$distanceMultipliedTraitReducedDataStructureLogFC\n"))
#         return(result)
#       }
#     )
#   })
#
#   session$userData$sessionVariables$distanceMultipliedProbeReducedDataStructure <- shiny::reactive({
#     id <- shiny::showNotification("Creating distance multiplied probe reduced data structure...", duration = NULL, closeButton = FALSE)
#     on.exit(shiny::removeNotification(id), add = TRUE)
# browser()
#     base::tryCatch(
#       {
#         if (is.valid(session$userData$sessionVariables$probeReducedDataStructure()$combinedDFP_Val_Labels$dfP_Val)) {
#           probeReducedDataStructure <- session$userData$sessionVariables$probeReducedDataStructure()
#           result <- base::list(combinedDFP_Val_Labels = probeReducedDataStructure$combinedDFP_Val_Labels,
#                                matP_Val.t = NULL,
#                                distMatTraits = NULL,
#                                clustResTraits = NULL,
#                                traitClusters = NULL,
#                                traitClusterMedoids = NULL,
#                                traitDendrogram = NULL,
#                                traitClustergram = NULL,
#                                distMatProbes = NULL,
#                                clustResProbes = NULL,
#                                probeDendrogram = NULL,
#                                DNAdistances = NULL
#           )
#           DNAdistances <- probeReducedDataStructure$DNAdistances
#           if(is.valid(DNAdistances)) {
#             DNAdistances <- na.omit(DNAdistances)
#             result$DNAdistances <- DNAdistances
#             dfP_Val <- result$combinedDFP_Val_Labels$dfP_Val
#             dfP_Val <- dfP_Val[which(rownames(dfP_Val) %in% DNAdistances$ID), ]
#             if (nrow(dfP_Val)>0) {
#               result$combinedDFP_Val_Labels$dfP_Val <- dfP_Val
#               matP_Val <- base::as.matrix(dfP_Val)
#               dt <- data.table::data.table(matP_Val)
#               # build set union of vec and dt
#               DNAdistances <- DNAdistances[which(DNAdistances$ID %in% rownames(dfP_Val)), ]
#               vec <- DNAdistances$meanDistance # take means# or
#               # vec <- DNAdistances$number # take numbers
#               # normalize vec to -1...1
#               vec <-  scales::rescale(vec, to = c(-1, 1))
#               # invert vec, so that small distances become large multiplies
#               vec <- vec * -1
#               result$combinedDFP_Val_Labels$dfP_Val <- data.table::data.table(t(t(dt) * vec))
#               rm(dt)
#             }
#             else {
#               base::print(base::paste0(sysTimePID(), " no set union between DNAdistances and dfP_Val in session$userData$sessionVariables$distanceMultipliedProbeReducedDataStructure\n"))
#               result$combinedDFP_Val_Labels$dfP_Val <- result$combinedDFP_Val_Labels$dfP_Val
#             }
#             result$combinedDFP_Val_Labels$dfDM <- session$userData$sessionVariables$probeReducedDataStructure()$combinedDFP_Val_Labels$dfDM
#             result$combinedDFP_Val_Labels$dfN <- session$userData$sessionVariables$probeReducedDataStructure()$combinedDFP_Val_Labels$dfN
#             result$matP_Val.t <- session$userData$sessionVariables$probeReducedDataStructure()$combinedDFP_Val_Labels$matP_Val.t
#             result$distMatTraits <- session$userData$sessionVariables$probeReducedDataStructure()$distMatTraits
#             result$clustResTraits <- session$userData$sessionVariables$probeReducedDataStructure()$clustResTraits
#             result$traitClusters <- session$userData$sessionVariables$probeReducedDataStructure()$traitClusters
#             result$traitClusterMedoids <- session$userData$sessionVariables$probeReducedDataStructure()$traitClusterMedoids
#             result$traitDendrogram <- session$userData$sessionVariables$probeReducedDataStructure()$traitDendrogram
#             result$traitClustergram <- session$userData$sessionVariables$probeReducedDataStructure()$traitClustergram
#
#             result$combinedDFP_Val_Labels$dfP_Val_w_number <- result$combinedDFP_Val_Labels$dfP_Val
#             nprobes <- nrow(result$combinedDFP_Val_Labels$dfP_Val_w_number)
#             result$combinedDFP_Val_Labels$dfP_Val_w_number$number <- seq(1:nprobes)
#             col_order <- c("number", colnames(result$combinedDFP_Val_Labels$dfP_Val_w_number))
#             result$combinedDFP_Val_Labels$dfP_Val_w_number <- result$combinedDFP_Val_Labels$dfP_Val_w_number[, col_order]
#             result$combinedDFP_Val_Labels$dfP_Val_w_number <- result$combinedDFP_Val_Labels$dfP_Val_w_number[ , -which(colnames(result$combinedDFP_Val_Labels$dfP_Val_w_number) %in% "number.1")]
#
#             result$combinedDFP_Val_Labels$dfDM_w_number <- result$combinedDFP_Val_Labels$dfDM
#             nprobes <- nrow(result$combinedDFP_Val_Labels$dfDM_w_number)
#             result$combinedDFP_Val_Labels$dfDM_w_number$number <- seq(1:nprobes)
#             col_order <- c("number", colnames(result$combinedDFP_Val_Labels$dfDM_w_number))
#             result$combinedDFP_Val_Labels$dfDM_w_number <- result$combinedDFP_Val_Labels$dfDM_w_number[, col_order]
#             #              result$dfDM_w_number <- result$dfDM[ , -which(colnames(result$dfDM_w_number) %in% "number.1")]
#
#             result$combinedDFP_Val_Labels$dfLogFC_w_number <- result$combinedDFP_Val_Labels$dfLogFC
#             nprobes <- nrow(result$combinedDFP_Val_Labels$dfLogFC_w_number)
#             result$combinedDFP_Val_Labels$dfLogFC_w_number$number <- seq(1:nprobes) #this line crashes...
#             col_order <- c("number", colnames(result$combinedDFP_Val_Labels$dfLogFC_w_number))
#             result$combinedDFP_Val_Labels$dfLogFC_w_number <- result$combinedDFP_Val_Labels$dfLogFC_w_number[, col_order]
#             result$combinedDFP_Val_Labels$dfLogFC_w_number <- result$combinedDFP_Val_Labels$dfLogFC_w_number[ , -which(colnames(result$combinedDFP_Val_Labels$dfLogFC_w_number) %in% "number.1")]
#
#             result$combinedDFP_Val_Labels$dfN_w_number <- result$combinedDFP_Val_Labels$dfN
#             nprobes <- nrow(result$combinedDFP_Val_Labels$dfN_w_number)
#             result$combinedDFP_Val_Labels$dfN_w_number$number <- seq(1:nprobes)
#             col_order <- c("number", colnames(result$combinedDFP_Val_Labels$dfN_w_number))
#             result$combinedDFP_Val_Labels$dfN_w_number <- result$combinedDFP_Val_Labels$dfN_w_number[, col_order]
#             result$combinedDFP_Val_Labels$dfN_w_number <- result$combinedDFP_Val_Labels$dfN_w_number[ , -which(colnames(result$combinedDFP_Val_Labels$dfN_w_number) %in% "number.1")]
#
#                       result$distMatProbes <- session$userData$sessionVariables$probeReducedDataStructure()$distMatProbes
#             result$clustResProbes <- session$userData$sessionVariables$probeReducedDataStructure()$clustResProbes
#             result$probeDendrogram <- session$userData$sessionVariables$probeReducedDataStructure()$probeDendrogram
#           }
#           else {
#             result <- NULL
#           }
#         }
#       },
#       error = function(e) {
#         base::message("An error occurred in shiny::reactive(session$userData$sessionVariables$distanceMultipliedProbeReducedDataStructure):\n", e)
#         browser() #should not happen
#       },
#       warning = function(w) {
#         base::message("A warning occurred in shiny::reactive(session$userData$sessionVariables$distanceMultipliedProbeReducedDataStructure):\n", w)
#         browser() #should not happen
#       },
#       finally = {
#         base::print(base::paste0(sysTimePID(), " finished generating session$userData$sessionVariables$distanceMultipliedProbeReducedDataStructure\n"))
#         return(result)
#       }
#     )
#   })

  traitReducedDTProbesPval <- shiny::reactive({
    id <- shiny::showNotification("Creating trait reduced data table probes (p-val)...", duration = NULL, closeButton = FALSE)
    on.exit(shiny::removeNotification(id), add = TRUE)
    base::tryCatch(
      {
        base::print(base::paste0(sysTimePID(), " start generating DTProbes p-val."))
        base::print(base::paste0(sysTimePID(), " before making probe table."))
        if (!is.null(session$userData$sessionVariables$traitReducedDataStructurePVal()$clustResProbes)) {
          dendProbes <- session$userData$sessionVariables$traitReducedDataStructurePVal()$probeDendrogram
          listProbes <- (base::labels(dendProbes))
          DTProbes <-
            base::data.frame(row.names = seq_along(listProbes))
          DTProbes$probeID <- listProbes
          DTProbes$order <- base::seq_len(base::nrow(DTProbes))
          rownames(DTProbes) <- DTProbes$probeID
          # add annotation
          DTProbes <-
            base::merge(
              x = DTProbes,
              y = session$userData$annotation,
              by.x = "probeID",
              by.y = "name",
              all.x = TRUE,
              all.y = FALSE
            )
          # sort
          DTProbes <- DTProbes[base::order(DTProbes$order), ]
          rownames(DTProbes) <- DTProbes$probeID
#           DTProbes <- addLinkToEWASDataHubShort(DTProbes, session$userData$config$baseURL_EWASDataHub, session$userData$config$probeAttribut)
#           DTProbes <- addLinkToMRCEWASCatalogShort(DTProbes, session$userData$config$baseURL_MRCEWASCatalog, session$userData$config$probeAttribut)
#
#           DTProbes <- addLinkToEWASDataHub(DTProbes, session$userData$config$baseURL_EWASDataHub, session$userData$config$probeAttribut)
#           DTProbes <- addLinkToMRCEWASCatalog(DTProbes, session$userData$config$baseURL_MRCEWASCatalog, session$userData$config$probeAttribut)
          DTProbes$probeID <- NULL
          result <- DTProbes
          base::print(base::paste0(sysTimePID(), " after making probe table."))
        }
        else {
          result <- NULL
        }
      },
      error = function(e) {
        base::message("An error occurred in shiny::reactive(traitReducedDTProbesPval):\n", e)
        browser() #should not happen
      },
      warning = function(w) {
        base::message("A warning occurred in shiny::reactive(traitReducedDTProbesPval):\n", w)
        browser() #should not happen
      },
      finally = {
        base::print(base::paste0(sysTimePID(), " finished generating traitReducedDTProbesPval."))
        return(result)
      }
    )
  })

  traitReducedDTProbesLogFC <- shiny::reactive({
    id <- shiny::showNotification("Creating trait reduced data table probes (log(FC))...", duration = NULL, closeButton = FALSE)
    on.exit(shiny::removeNotification(id), add = TRUE)
    base::tryCatch(
      {
        base::print(base::paste0(sysTimePID(), " start generating DTProbes log(FC)."))
        base::print(base::paste0(sysTimePID(), " before making probe table."))
        if (!is.null(session$userData$sessionVariables$traitReducedDataStructureLogFC()$clustResProbes)) {
          dendProbes <- session$userData$sessionVariables$traitReducedDataStructureLogFC()$probeDendrogram
          listProbes <- (base::labels(dendProbes))
          DTProbes <-
            base::data.frame(row.names = seq_along(listProbes))
          DTProbes$probeID <- listProbes
          DTProbes$order <- base::seq_len(base::nrow(DTProbes))
          rownames(DTProbes) <- DTProbes$probeID
          # add annotation
          DTProbes <-
            base::merge(
              x = DTProbes,
              y = session$userData$annotation,
              by.x = "probeID",
              by.y = "name",
              all.x = TRUE,
              all.y = FALSE
            )
          # sort
          DTProbes <- DTProbes[base::order(DTProbes$order), ]
          rownames(DTProbes) <- DTProbes$probeID
          #           DTProbes <- addLinkToEWASDataHubShort(DTProbes, session$userData$config$baseURL_EWASDataHub, session$userData$config$probeAttribut)
          #           DTProbes <- addLinkToMRCEWASCatalogShort(DTProbes, session$userData$config$baseURL_MRCEWASCatalog, session$userData$config$probeAttribut)
          #
          #           DTProbes <- addLinkToEWASDataHub(DTProbes, session$userData$config$baseURL_EWASDataHub, session$userData$config$probeAttribut)
          #           DTProbes <- addLinkToMRCEWASCatalog(DTProbes, session$userData$config$baseURL_MRCEWASCatalog, session$userData$config$probeAttribut)
          DTProbes$probeID <- NULL
          result <- DTProbes
          base::print(base::paste0(sysTimePID(), " after making probe table."))
        }
        else {
          result <- NULL
        }
      },
      error = function(e) {
        base::message("An error occurred in shiny::reactive(traitReducedDTProbesLogFC):\n", e)
        browser() #should not happen
      },
      warning = function(w) {
        base::message("A warning occurred in shiny::reactive(traitReducedDTProbesLogFC):\n", w)
        browser() #should not happen
      },
      finally = {
        base::print(base::paste0(sysTimePID(), " finished generating traitReducedDTProbesLogFC."))
        return(result)
      }
    )
  })

  # fullDWDTProbes <- shiny::reactive({
  #   base::tryCatch(
  #     {
  #       base::print(base::paste0(sysTimePID(), " start generating fullDWDTProbes"))
  #       base::print(base::paste0(sysTimePID(), " before making probe table."))
  #       if (!is.null(session$userData$sessionVariables$distanceMultipliedTraitReducedDataStructurePVal()$clustResProbes)) {
  #         dendProbes <- session$userData$sessionVariables$distanceMultipliedTraitReducedDataStructurePVal()$probeDendrogram
  #         listProbes <- (base::labels(dendProbes))
  #         DTProbes <-
  #           base::data.frame(row.names = seq_along(listProbes))
  #         DTProbes$probeID <- listProbes
  #         DTProbes$order <- base::seq_len(base::nrow(DTProbes))
  #         rownames(DTProbes) <- DTProbes$probeID
  #         # add annotation
  #         DTProbes <-
  #           base::merge(
  #             x = DTProbes,
  #             y = session$userData$annotation,
  #             by.x = "probeID",
  #             by.y = "name",
  #             all.x = TRUE,
  #             all.y = FALSE
  #           )
  #         # sort
  #         DTProbes <- DTProbes[base::order(DTProbes$order), ]
  #         rownames(DTProbes) <- DTProbes$probeID
  #         result <- DTProbes
  #         base::print(base::paste0(sysTimePID(), " after making probe table."))
  #       }
  #       else {
  #         result <- NULL
  #       }
  #     },
  #     error = function(e) {
  #       base::message("An error occurred in shiny::reactive(fullDWDTProbes):\n", e)
  #     },
  #     warning = function(w) {
  #       base::message("A warning occurred in shiny::reactive(fullDWDTProbes):\n", w)
  #     },
  #     finally = {
  #       base::print(base::paste0(sysTimePID(), " finished generating fullDWDTProbes"))
  #       return(result)
  #     }
  #   )
  # })
  #
  # DWDTProbes <- shiny::reactive({
  #   base::tryCatch(
  #     {
  #       base::print(base::paste0(sysTimePID(), " start generating DWDTProbes"))
  #       base::print(base::paste0(sysTimePID(), " before making probe table."))
  #       if (!is.null(session$userData$sessionVariables$probeReducedDataStructure()$clustResProbes)) {
  #         dendProbes <- session$userData$sessionVariables$probeReducedDataStructure()$probeDendrogram
  #         listProbes <- (base::labels(dendProbes)) # base::as.numeric
  #         DTProbes <-
  #           base::data.frame(row.names = seq_along(listProbes))
  #         DTProbes$probeID <- listProbes
  #         DTProbes$order <- base::seq_len(base::nrow(DTProbes))
  #         rownames(DTProbes) <- DTProbes$probeID
  #         # add annotation
  #         DTProbes <-
  #           base::merge(
  #             x = DTProbes,
  #             y = session$userData$annotation, #y = globalVariables$annotation,
  #             by.x = "probeID",
  #             by.y = "name",
  #             all.x = TRUE,
  #             all.y = FALSE
  #           )
  #         # sort
  #         DTProbes <- DTProbes[base::order(DTProbes$order), ]
  #         rownames(DTProbes) <- DTProbes$probeID
  #         result <- DTProbes
  #         base::print(base::paste0(sysTimePID(), " after making probe table."))
  #       }
  #       else {
  #         result <- NULL
  #       }
  #     },
  #     error = function(e) {
  #       base::message("An error occurred in shiny::reactive(DWDTProbes):\n", e)
  #     },
  #     warning = function(w) {
  #       base::message("A warning occurred in shiny::reactive(DWDTProbes):\n", w)
  #     },
  #     finally = {
  #       base::print(base::paste0(sysTimePID(), " finished generating DWDTProbes"))
  #       return(result)
  #     }
  #   )
  # })

  condDTProbes <- shiny::reactive({
    base::tryCatch(
      {
        base::print(base::paste0(sysTimePID(), " start generating condDTProbes"))
        base::print(base::paste0(sysTimePID(), " before making probe table."))
        if (!is.null(session$userData$sessionVariables$probeReducedDataStructure()$clustResProbes)) {
          dendProbes <- session$userData$sessionVariables$probeReducedDataStructure()$probeDendrogram
          listProbes <- (base::labels(dendProbes)) # base::as.numeric
          DTProbes <-
            base::data.frame(row.names = seq_along(listProbes))
          DTProbes$probeID <- listProbes
          DTProbes$order <- base::seq_len(base::nrow(DTProbes))
          rownames(DTProbes) <- DTProbes$probeID
          # add annotation
          DTProbes <-
            base::merge(
              x = DTProbes,
              y = session$userData$annotation, #y = globalVariables$annotation,
              by.x = "probeID",
              by.y = "name",
              all.x = TRUE,
              all.y = FALSE
            )
          # sort
          DTProbes <- DTProbes[base::order(DTProbes$order), ]
          rownames(DTProbes) <- DTProbes$probeID
          DTProbes$probeID <- NULL
          result <- DTProbes
          base::print(base::paste0(sysTimePID(), " after making probe table."))
        }
        else {
          result <- NULL
        }
      },
      error = function(e) {
        base::message("An error occurred in shiny::reactive(condDTProbes):\n", e)
      },
      warning = function(w) {
        base::message("A warning occurred in shiny::reactive(condDTProbes):\n", w)
      },
      finally = {
        base::print(base::paste0(sysTimePID(), " finished generating condDTProbes"))
        return(result)
      }
    )
  })

  # condDWDTProbes <- shiny::reactive({
  #   base::tryCatch(
  #     {
  #       base::print(base::paste0(sysTimePID(), " start generating condDWDTProbes"))
  #       base::print(base::paste0(sysTimePID(), " before making probe table."))
  #       if (!is.null(session$userData$sessionVariables$distanceMultipliedProbeReducedDataStructure()$clustResProbes)) {
  #         dendProbes <- session$userData$sessionVariables$distanceMultipliedProbeReducedDataStructure()$probeDendrogram
  #         listProbes <- (base::labels(dendProbes)) # base::as.numeric
  #         DTProbes <-
  #           base::data.frame(row.names = seq_along(listProbes))
  #         DTProbes$probeID <- listProbes
  #         DTProbes$order <- base::seq_len(base::nrow(DTProbes))
  #         rownames(DTProbes) <- DTProbes$probeID
  #         # add annotation
  #         DTProbes <-
  #           base::merge(
  #             x = DTProbes,
  #             y = session$userData$annotation, #y = globalVariables$annotation,
  #             by.x = "probeID",
  #             by.y = "name",
  #             all.x = TRUE,
  #             all.y = FALSE
  #           )
  #         # sort
  #         DTProbes <- DTProbes[base::order(DTProbes$order), ]
  #         rownames(DTProbes) <- DTProbes$probeID
  #         result <- DTProbes
  #         base::print(base::paste0(sysTimePID(), " after making probe table."))
  #       }
  #       else {
  #         result <- NULL
  #       }
  #     },
  #     error = function(e) {
  #       base::message("An error occurred in shiny::reactive(condDWDTProbes):\n", e)
  #     },
  #     warning = function(w) {
  #       base::message("A warning occurred in shiny::reactive(condDWDTProbes):\n", w)
  #     },
  #     finally = {
  #       base::print(base::paste0(sysTimePID(), " finished generating condDWDTProbes"))
  #       return(result)
  #     }
  #   )
  # })

  DTTraitsPval <- shiny::reactive({
    id <- shiny::showNotification("Creating data table traits (p-val)...", duration = NULL, closeButton = FALSE)
    on.exit(shiny::removeNotification(id), add = TRUE)
    base::tryCatch(
      {
        base::print(base::paste0(sysTimePID(), " start generating DTTraitsPval"))
        base::print(base::paste0(sysTimePID(), " before making traits table."))
        if (!is.null(session$userData$sessionVariables$traitReducedDataStructurePVal()$clustResTraits)) {
          listTraits <- ggdendro::dendro_data(session$userData$sessionVariables$traitReducedDataStructurePVal()$clustResTraits, type = "rectangle")$labels$label

          base::print(base::paste0(sysTimePID(), " before rendering dendrogram tables traits"))
          DTTraits <-
            base::data.frame(row.names = seq_along(listTraits))
          DTTraits$Name <- listTraits
          DTTraits$order <- base::seq_len(base::nrow(DTTraits))
          rownames(DTTraits) <- DTTraits$Name
          DTTraits <- DTTraits[order(DTTraits$order), ]
          base::print(base::paste0(sysTimePID(), " after making traits table."))
          result <- DTTraits
        }
        else {
          result <- NULL
        }
      },
      error = function(e) {
        base::message("An error occurred in shiny::reactive(DTTraitsPval):\n", e)
      },
      warning = function(w) {
        base::message("A warning occurred in shiny::reactive(DTTraitsPval):\n", w)
      },
      finally = {
        base::print(base::paste0(sysTimePID(), " finished generating DTTraitsPval"))
        return(result)
      }
    )
  })

  DTTraitsLogFC  <- shiny::reactive({
    id <- shiny::showNotification("Creating data table traits (log(FC))...", duration = NULL, closeButton = FALSE)
    on.exit(shiny::removeNotification(id), add = TRUE)
    base::tryCatch(
      {
        base::print(base::paste0(sysTimePID(), " start generating DTTraitsLogFC"))
        base::print(base::paste0(sysTimePID(), " before making traits table."))
        if (!is.null(session$userData$sessionVariables$traitReducedDataStructureLogFC()$clustResTraits)) {
          listTraits <- ggdendro::dendro_data(session$userData$sessionVariables$traitReducedDataStructureLogFC()$clustResTraits, type = "rectangle")$labels$label

          base::print(base::paste0(sysTimePID(), " before rendering dendrogram tables traits"))
          DTTraits <-
            base::data.frame(row.names = seq_along(listTraits))
          DTTraits$Name <- listTraits
          DTTraits$order <- base::seq_len(base::nrow(DTTraits))
          rownames(DTTraits) <- DTTraits$Name
          DTTraits <- DTTraits[order(DTTraits$order), ]
          base::print(base::paste0(sysTimePID(), " after making traits table."))
          result <- DTTraits
        }
        else {
          result <- NULL
        }
      },
      error = function(e) {
        base::message("An error occurred in shiny::reactive(DTTraitsLogFC):\n", e)
      },
      warning = function(w) {
        base::message("A warning occurred in shiny::reactive(DTTraitsLogFC):\n", w)
      },
      finally = {
        base::print(base::paste0(sysTimePID(), " finished generating DTTraitsLogFC"))
        return(result)
      }
    )
  })

  # fullDWDTTraits <- shiny::reactive({
  #   base::tryCatch(
  #     {
  #       base::print(base::paste0(sysTimePID(), " start generating fullDWDTTraits"))
  #       base::print(base::paste0(sysTimePID(), " before making traits table."))
  #       if (!is.null(session$userData$sessionVariables$distanceMultipliedTraitReducedDataStructurePVal()$clustResTraits)) {
  #         listTraits <- ggdendro::dendro_data(session$userData$sessionVariables$distanceMultipliedTraitReducedDataStructurePVal()$clustResTraits, type = "rectangle")$labels$label
  #
  #         base::print(base::paste0(sysTimePID(), " before rendering dendrogram tables traits"))
  #         DTTraits <-
  #           base::data.frame(row.names = seq_along(listTraits))
  #         DTTraits$Name <- listTraits
  #         DTTraits$order <- base::seq_len(base::nrow(DTTraits))
  #         rownames(DTTraits) <- DTTraits$Name
  #         DTTraits <- DTTraits[order(DTTraits$order), ]
  #         base::print(base::paste0(sysTimePID(), " after making traits table."))
  #         result <- DTTraits
  #       }
  #       else {
  #         result <- NULL
  #       }
  #     },
  #     error = function(e) {
  #       base::message("An error occurred in shiny::reactive(fullDWDTTraits):\n", e)
  #     },
  #     warning = function(w) {
  #       base::message("A warning occurred in shiny::reactive(fullDWDTTraits):\n", w)
  #     },
  #     finally = {
  #       base::print(base::paste0(sysTimePID(), " finished generating fullDWDTTraits"))
  #       return(result)
  #     }
  #   )
  # })

  condDTTraits <- shiny::reactive({
    base::tryCatch(
      {
        base::print(base::paste0(sysTimePID(), " start generating condDTTraits"))
        base::print(base::paste0(sysTimePID(), " before making traits table."))
        if (!is.null(session$userData$sessionVariables$probeReducedDataStructure()$clustResTraits)) {
          listTraits <- ggdendro::dendro_data(session$userData$sessionVariables$probeReducedDataStructure()$clustResTraits, type = "rectangle")$labels$label

          base::print(base::paste0(sysTimePID(), " before rendering dendrogram tables traits"))
          DTTraits <-
            base::data.frame(row.names = seq_along(listTraits))
          DTTraits$Name <- listTraits
          DTTraits$order <- base::seq_len(base::nrow(DTTraits))
          rownames(DTTraits) <- DTTraits$Name
          DTTraits <- DTTraits[order(DTTraits$order), ]
          base::print(base::paste0(sysTimePID(), " after making traits table."))
          result <- DTTraits
        }
        else {
          result <- NULL
        }
      },
      error = function(e) {
        base::message("An error occurred in shiny::reactive(condDTTraits):\n", e)
      },
      warning = function(w) {
        base::message("A warning occurred in shiny::reactive(condDTTraits):\n", w)
      },
      finally = {
        base::print(base::paste0(sysTimePID(), " finished generating condDTTraits"))
        return(result)
      }
    )
  })

  # condDWDTTraits <- shiny::reactive({
  #   base::tryCatch(
  #     {
  #       base::print(base::paste0(sysTimePID(), " start generating condDWDTTraits"))
  #       base::print(base::paste0(sysTimePID(), " before making traits table."))
  #       if (!is.null(session$userData$sessionVariables$distanceMultipliedProbeReducedDataStructure()$clustResTraits)) {
  #         listTraits <- ggdendro::dendro_data(session$userData$sessionVariables$distanceMultipliedProbeReducedDataStructure()$clustResTraits, type = "rectangle")$labels$label
  #
  #         base::print(base::paste0(sysTimePID(), " before rendering dendrogram tables traits"))
  #         DTTraits <-
  #           base::data.frame(row.names = seq_along(listTraits))
  #         DTTraits$Name <- listTraits
  #         DTTraits$order <- base::seq_len(base::nrow(DTTraits))
  #         rownames(DTTraits) <- DTTraits$Name
  #         DTTraits <- DTTraits[order(DTTraits$order), ]
  #         base::print(base::paste0(sysTimePID(), " after making traits table."))
  #         result <- DTTraits
  #       }
  #       else {
  #         result <- NULL
  #       }
  #     },
  #     error = function(e) {
  #       base::message("An error occurred in shiny::reactive(condDWDTTraits):\n", e)
  #     },
  #     warning = function(w) {
  #       base::message("A warning occurred in shiny::reactive(condDWDTTraits):\n", w)
  #     },
  #     finally = {
  #       base::print(base::paste0(sysTimePID(), " finished generating condDWDTTraits"))
  #       return(result)
  #     }
  #   )
  # })

  observeEvent(input$keypressed,
   {
     #catch ESC key to prevent unwanted stopping from ESC key press
     if(input$keypressed==27) {
        # browser()
        # stopApp()
      }
   })

  shiny::observeEvent(input$btnPlotCombinedHM_P_Val,
    ignoreInit = TRUE,
    {
      base::tryCatch(
        {
          base::print(base::paste0(sysTimePID(), " Step 5a: start plotting heatmap for P_Val. (first step in shiny::observeEvent(input$btnPlotCombinedHM_P_Val))"))
          plotCombinedHM_P_Val(input = input, output = output, session = session)
          #          plotHMDNADistances(input = input, output = output, session = session)
          session$userData$sessionVariables$callCounter <- session$userData$sessionVariables$callCounter + 1
        },
        error = function(e) {
          base::message("An error occurred in shiny::observeEvent(input$plotCombinedHM_P_Val):\n", e)
        },
        warning = function(w) {
          base::message("A warning occurred in shiny::observeEvent(input$plotCombinedHM_P_Val):\n", w)
        },
        finally = {
          base::print(base::paste0(sysTimePID(), " finished plotting heatmap for P_Val. (last step in shiny::observeEvent(input$btnPlotCombinedHM_P_Val))"))
        }
      )
    },
    ignoreNULL = FALSE
  )

  shiny::observeEvent(input$btnPlotCombinedHM_LogFC,
    ignoreInit = TRUE,
    {
      base::tryCatch(
        {
          base::print(base::paste0(sysTimePID(), " Step 5b: start plotting heatmap for log(FC). (first step in shiny::observeEvent(input$btnPlotCombinedHM_LogFC))"))
          plotCombinedHM_LogFC(input = input, output = output, session = session)
          session$userData$sessionVariables$callCounter <- session$userData$sessionVariables$callCounter + 1
        },
        error = function(e) {
          base::message("An error occurred in shiny::observeEvent(input$btnPlotCombinedHM_LogFC):\n", e)
        },
        warning = function(w) {
          base::message("A warning occurred in shiny::observeEvent(input$btnPlotCombinedHM_LogFC):\n", w)
        },
        finally = {
          base::print(base::paste0(sysTimePID(), " finished plotting heatmap for log(FC). (last step in shiny::observeEvent(input$btnPlotCombinedHM_LogFC))"))
        }
      )
    },
    ignoreNULL = FALSE
  )

  shiny::observeEvent(input$btnPlotCombinedCondHM_DM,
    ignoreInit = TRUE,
    {
      base::tryCatch(
        {
          base::print(base::paste0(sysTimePID(), " Step 6b: start plotting heatmap for Delta Methylation (logFC). (first step in shiny::observeEvent(input$btnPlotCombinedCondHM_DM))"))
          plotCombinedHM_DMLogFC(input = input, output = output, session = session)
        },
        error = function(e) {
          base::message("An error occurred in shiny::observeEvent(input$btnPlotCombinedCondHM_DM):\n", e)
        },
        warning = function(w) {
          base::message("A warning occurred in shiny::observeEvent(input$btnPlotCombinedCondHM_DM):\n", w)
        },
        finally = {
          base::print(base::paste0(sysTimePID(), " finished plotting heatmap for Delta Methylation (log(FC)). (last step in shiny::observeEvent(input$btnPlotCombinedCondHM_DM))"))
        }
      )
    },
    ignoreNULL = FALSE
  )

  shiny::observeEvent(input$btnPlotCombinedDWHM_P_Val,
    ignoreInit = TRUE,
    {
      base::tryCatch(
        {
          base::print(base::paste0(sysTimePID(), " Step 6: start plotting heatmap for distance weighted DM. (first step in shiny::observeEvent(input$btnPlotCombinedDWHM_P_Val))"))
          plotCombinedDWHM_P_Val(input = input, output = output, session = session)
          #session$userData$sessionVariables$callCounter <- session$userData$sessionVariables$callCounter + 1
        },
        error = function(e) {
          base::message("An error occurred in shiny::observeEvent(input$btnPlotCombinedDWHM_P_Val):\n", e)
        },
        warning = function(w) {
          base::message("A warning occurred in shiny::observeEvent(input$btnPlotCombinedDWHM_P_Val):\n", w)
        },
        finally = {
          base::print(base::paste0(sysTimePID(), " finished plotting heatmap for P_Val. (last step in shiny::observeEvent(input$btnPlotCombinedDWHM_P_Val))"))
        }
      )
    },
    ignoreNULL = FALSE
  )

  shiny::observeEvent(input$btnPlotCombinedCondHM_P_Val,
                      ignoreInit = TRUE,
                      {
                        base::tryCatch(
                          {
                            base::print(base::paste0(sysTimePID(), " Step 6: start plotting condensed heatmap for distance weighted DM. (first step in shiny::observeEvent(input$btnPlotCombinedCondHM_P_Val))"))
                            plotCombinedCondHM_P_Val(input = input, output = output, session = session)
                          },
                          error = function(e) {
                            base::message("An error occurred in shiny::observeEvent(input$btnPlotCombinedCondHM_P_Val):\n", e)
                          },
                          warning = function(w) {
                            base::message("A warning occurred in shiny::observeEvent(input$btnPlotCombinedCondHM_P_Val):\n", w)
                          },
                          finally = {
                            base::print(base::paste0(sysTimePID(), " finished plotting heatmap for P_Val. (last step in shiny::observeEvent(input$btnPlotCombinedCondHM_P_Val))"))
                          }
                        )
                      },
                      ignoreNULL = FALSE
  )

  shiny::observeEvent(input$btnPlotCombinedCondDWHM_P_Val,
                      ignoreInit = TRUE,
                      {
                        base::tryCatch(
                          {
                            base::print(base::paste0(sysTimePID(), " Step 6: start plotting condensed heatmap for distance weighted DM. (first step in shiny::observeEvent(input$btnPlotCombinedCondDWHM_P_Val))"))
                            plotCombinedCondDWHM_P_Val(input = input, output = output, session = session)
                          },
                          error = function(e) {
                            base::message("An error occurred in shiny::observeEvent(input$btnPlotCombinedCondDWHM_P_Val):\n", e)
                          },
                          warning = function(w) {
                            base::message("A warning occurred in shiny::observeEvent(input$btnPlotCombinedCondDWHM_P_Val):\n", w)
                          },
                          finally = {
                            base::print(base::paste0(sysTimePID(), " finished plotting heatmap for P_Val. (last step in shiny::observeEvent(input$btnPlotCombinedCondDWHM_P_Val))"))
                          }
                        )
                      },
                      ignoreNULL = FALSE
  )

  shiny::observeEvent(input$chkDebug,
    ignoreInit = TRUE, #FALSE,
    {
      base::tryCatch(
        {
          if (input$chkDebug == TRUE) {
            session$userData$config$debugMode <- TRUE
            base::print(base::paste0(sysTimePID(), " set debugMode = TRUE."))
          }
          else {
            session$userData$config$debugMode <- FALSE
            base::print(base::paste0(sysTimePID(), " set debugMode = FALSE."))
          }
          result <- loadDirLists(session = session, input = input, output = output)
          dfdD1 <- result$dfdD1
          dfdD2 <- result$dfdD2
          dfdD3 <- result$dfdD3
        },
        error = function(e) {
          base::message("An error occurred in shiny::observeEvent(input$chkDebug):\n", e)
        },
        warning = function(w) {
          base::message("A warning occurred in shiny::observeEvent(input$chkDebug):\n", w)
        },
        finally = {
        }
      )
    },
    ignoreNULL = FALSE
  )

}

