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
  configFileLocation <- base::options()$configFileLocation
  if (missing(configFileLocation)) {
    configFileLocation <- system.file("config.yml", package = "PatternMatchR", mustWork = TRUE)
  }
  else {
    configFileLocation <- paste0(configFileLocation ,"/config.yml")
  }
  session$userData$config <- config::get(file = configFileLocation)
  volumes <- c(Home = paste0(getwd(), sub(".", "", session$userData$config$workDir)),
               shinyFiles::getVolumes()())
  shinyFiles::shinyFileChoose(input, "file", roots = volumes, session = session)
  shinyFiles::shinyFileSave(input, "save", roots = volumes, session = session, restrictions = system.file(package = "base"))
  if (is.valid(session$userData$config$knownCpGs)) {
    #load CpG from txt file to search input
    knownCpG <- paste(unlist(data.table::fread(file = session$userData$config$knownCpGs, header = FALSE)), collapse = " ")
    shiny::updateTextInput(session, inputId = "txtSearchFullCpGPVal", value = knownCpG)
  }
  if (is.valid(session$userData$config$knownTraits)) {
    #load Traits from txt file to search input
    knownTrait <- paste(unlist(data.table::fread(file = session$userData$config$knownTraits, header = FALSE)), collapse = " ")
    shiny::updateTextInput(session, inputId = "txtSearchFullTraitPVal", value = knownTrait)
  }

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
  session$userData$sessionVariables$pReducedDataStructure <- shiny::reactiveVal(value = NULL, label = "pReducedDataStructure")
  session$userData$sessionVariables$pReducedData <- shiny::reactiveVal(value = NULL, label = "pReducedData")
  session$userData$sessionVariables$traitReducedDataStructurePVal <- shiny::reactiveVal(value = NULL, label = "traitReducedDataStructurePVal")
  session$userData$sessionVariables$traitReducedDataStructureLogFC <- shiny::reactiveVal(value = NULL, label = "traitReducedDataStructureLogFC")
  session$userData$sessionVariables$traitReducedDataPVal <- shiny::reactiveVal(value = NULL, label = "traitReducedDataPVal")
  session$userData$sessionVariables$traitReducedDataLogFC <- shiny::reactiveVal(value = NULL, label = "traitReducedDataLogFC")
  session$userData$sessionVariables$probeReducedDataStructurePVal <- shiny::reactiveVal(value = NULL, label = "probeReducedDataStructurePVal")
  session$userData$sessionVariables$probeReducedDataStructureLogFC <- shiny::reactiveVal(value = NULL, label = "probeReducedDataStructureLogFC")
  session$userData$sessionVariables$probeReducedDataPVal <- shiny::reactiveVal(value = NULL, label = "probeReducedDataPVal")
  session$userData$sessionVariables$probeReducedDataLogFC <- shiny::reactiveVal(value = NULL, label = "probeReducedDataLogFC")
  session$userData$sessionVariables$probeReducedDataWOGapPVal <- shiny::reactiveVal(value = NULL, label = "probeReducedDataWOGapPVal")
  session$userData$sessionVariables$probeReducedDataWOGapLogFC <- shiny::reactiveVal(value = NULL, label = "probeReducedDataWOGapLogFC")
  session$userData$sessionVariables$probeReducedDataStructureWOGapPVal <- shiny::reactiveVal(value = NULL, label = "probeReducedDataStructurePVal")
  session$userData$sessionVariables$probeReducedDataStructureWOGapLogFC <- shiny::reactiveVal(value = NULL, label = "probeReducedDataStructureLogFC")

  session$userData$sessionVariables$generalDataStructure <- shiny::reactiveVal(value = NULL, label = "generalDataStructure")

  session$userData$sessionVariables$distancesBelowThreshold <- shiny::reactiveVal(value = NULL, label = "distancesBelowThreshold")

  session$userData$sessionVariables$distNeigboursProbes10000 <- shiny::reactiveVal(value = NULL, label = "distNeigboursProbes10000")
  session$userData$sessionVariables$distNeigboursProbes1000 <- shiny::reactiveVal(value = NULL, label = "distNeigboursProbes1000")
  session$userData$sessionVariables$distNeigboursProbes100 <- shiny::reactiveVal(value = NULL, label = "distNeigboursProbes100")

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

####Start module server integration
  ReduceData_SERVER("Reduce", session)
  Clustering_SERVER("Traits", session)
  ClusteringTraits_SERVER("PVal", session$userData$sessionVariables$pReducedDataStructure, session$userData$sessionVariables$traitReducedDataStructurePVal, session)
  ClusteringTraits_SERVER("LogFC", session$userData$sessionVariables$pReducedDataStructure, session$userData$sessionVariables$traitReducedDataStructureLogFC, session)

  ClusteringProbesGeneral_SERVER("Probes", session)
  ClusteringProbes_SERVER("PVal", session$userData$sessionVariables$traitReducedDataStructurePVal, session$userData$sessionVariables$probeReducedDataStructurePVal, session$userData$sessionVariables$probeReducedDataStructureWOGapPVal, session)
  ClusteringProbes_SERVER("LogFC", session$userData$sessionVariables$traitReducedDataStructureLogFC, session$userData$sessionVariables$probeReducedDataStructureLogFC, session$userData$sessionVariables$probeReducedDataStructureWOGapLogFC, session)

  HeatMap_SERVER("PVal", session)
  HeatMap_SERVER("PValWOGap", session)
  HeatMap_SERVER("LogFC", session)
  HeatMap_SERVER("LogFCWOGap", session)

  Search_Full_SERVER("Search", session)
  GlobalSelection_SERVER("GlobalSelection", session)
  VolcanoPlot_SERVER("VolcanoPlot", session)
  PCPlot_SERVER("PCPlot", session)
####End module server integration

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
    if (!is.valid(session$userData$sessionVariables$probeReducedDataStructurePVal()$combinedDFP_Val_Labels$dfP_Val)) {
      shinyjs::disable("btnPlotCombinedHM_P_Val")
    }
    else {
      shinyjs::enable("btnPlotCombinedHM_P_Val")
    }
  })

  shiny::observe({
    if (!is.valid(session$userData$sessionVariables$probeReducedDataStructureWOGapPVal()$combinedDFP_Val_Labels$dfP_Val)) {
      shinyjs::disable("btnPlotCombinedHM_P_ValWOGap")
    }
    else {
      shinyjs::enable("btnPlotCombinedHM_P_ValWOGap")
    }
  })

  shiny::observe({
    if (!is.valid(session$userData$sessionVariables$probeReducedDataStructureLogFC()$combinedDFP_Val_Labels$dfLogFC)) {
      shinyjs::disable("btnPlotCombinedHM_LogFC")
    }
    else{
      shinyjs::enable("btnPlotCombinedHM_LogFC")
    }
  })

  shiny::observe({
    if (!is.valid(session$userData$sessionVariables$probeReducedDataStructureWOGapLogFC()$combinedDFP_Val_Labels$dfLogFC)) {
      shinyjs::disable("btnPlotCombinedHM_LogFCWOGap")
    }
    else{
      shinyjs::enable("btnPlotCombinedHM_LogFCWOGap")
    }
  })

  session$userData$sessionVariables$distNeigboursProbes100 <- shiny::reactive({
    #calculate distance from each probe to its neigbours to build a right column in heatmap to estimate relevance of heatmap findings
    base::tryCatch(
      {
        result <- calculateDistanceNeigboursProbes(wd = session$userData$packageWd, clustResProbes = session$userData$sessionVariables$traitReducedDataStructurePVal()$clustResProbes, annotation = session$userData$annotation, distanceToLook = 100, numCores = session$userData$numCores)
      },
      error = function(e) {
        if (attributes(e)$class[1] != "shiny.silent.error") {
          base::message("An error occurred in shiny::reactive(session$userData$sessionVariables$distNeigboursProbes100):\n", e)
        }
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
        if (attributes(e)$class[1] != "shiny.silent.error") {
          base::message("An error occurred in shiny::reactive(session$userData$sessionVariables$distNeigboursProbes1000):\n", e)
        }
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
        if (attributes(e)$class[1] != "shiny.silent.error") {
          base::message("An error occurred in shiny::reactive(session$userData$sessionVariables$distNeigboursProbes10000):\n", e)
        }
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
        if (attributes(e)$class[1] != "shiny.silent.error") {
          base::message("An error occurred in shiny::reactive(session$userData$sessionVariables$distancesBelowThreshold):\n", e)
        }
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
        if (attributes(e)$class[1] != "shiny.silent.error") {
          base::message("An error occurred in shiny::reactive(output$txtLoadOut):\n", e)
        }
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
          result <- updateTxtMergeOut(session$userData$sessionVariables$combinedDataStructure())
       }
       else {
         result <- NULL
       }
      },
      error = function(e) {
        if (attributes(e)$class[1] != "shiny.silent.error") {
          base::message("An error occurred in shiny::reactive(output$txtMergeOut):\n", e)
        }
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
          if (attributes(e)$class[1] != "shiny.silent.error") {
            base::message("An error occurred in shiny::observeEvent(input$btnLoadDir1):\n", e)
          }
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
          if (attributes(e)$class[1] != "shiny.silent.error") {
            base::message("An error occurred in shiny::observeEvent(input$btnLoadDir2):\n", e)
          }
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
          if (attributes(e)$class[1] != "shiny.silent.error") {
            base::message("An error occurred in shiny::observeEvent(input$btnLoadDir3):\n", e)
          }
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
#browser() #check for rownames in session$userData$sessionVariables$resultDFListTrait1()$PHENODF -> we don't have them
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
          if (attributes(e)$class[1] != "shiny.silent.error") {
            base::message("An error occurred in shiny::observeEvent(input$btnLoadDirAll):\n", e)
          }
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
          if (attributes(e)$class[1] != "shiny.silent.error") {
            base::message("An error occurred in shiny::observeEvent(input$btnDebug):\n", e)
          }
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
            #change names in dendrogram and heatmap to more readable ones
            prefix <- dplyr::case_match(result$mergedOriginTrait, "1" ~ "red", "2" ~ "green", "3" ~ "blue")
            cn <- paste0(prefix, "_", result$mergedOriginalColnames)
            colnames(result$dfP_Val) <- cn
            colnames(result$dfDM) <- cn
            colnames(result$dfN) <- cn
            colnames(result$dfLogFC) <- cn
          }
          else {
            base::print(base::paste0(sysTimePID(), " is.valid(session$userData$sessionVariables$resultDFListTrait1()) || is.valid(session$userData$sessionVariables$resultDFListTrait2())  || is.valid(session$userData$sessionVariables$resultDFListTrait3()) == FALSE."))
            result <- NULL
          }
#          updateReduceDataSliders(session, result)
          session$userData$sessionVariables$combinedData(result)
        },
        error = function(e) {
          if (attributes(e)$class[1] != "shiny.silent.error") {
            base::message("An error occurred in shiny::observeEvent(input$btnMerge):\n", e)
          }
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
          if (attributes(e)$class[1] != "shiny.silent.error") {
            base::message("An error occurred in shiny::observeEvent(input$btnCountProbesP_ValParallel):\n", e)
          }
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
          if (attributes(e)$class[1] != "shiny.silent.error") {
            base::message("An error occurred in shiny::observeEvent(input$btnCountProbesDeltaMethParallel):\n", e)
          }
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
          if (attributes(e)$class[1] != "shiny.silent.error") {
            base::message("An error occurred in shiny::observeEvent(input$btnCountProbesNParallel):\n", e)
          }
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
          if (attributes(e)$class[1] != "shiny.silent.error") {
            base::message("An error occurred in shiny::observeEvent(input$btnExportCombinedData):\n", e)
          }
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
          if (attributes(e)$class[1] != "shiny.silent.error") {
            base::message("An error occurred in shiny::observeEvent(input$btnImportCombinedData):\n", e)
          }
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
  # output$condDTP_VAL <- DT::renderDataTable(as.data.frame(session$userData$sessionVariables$probeReducedDataStructure()$combinedDFP_Val_Labels$dfP_Val_w_number))
  # output$condDTDM <- DT::renderDataTable(as.data.frame(session$userData$sessionVariables$probeReducedDataStructure()$combinedDFP_Val_Labels$dfDM_w_number))
  # output$condDTLogFC <- DT::renderDataTable(as.data.frame(session$userData$sessionVariables$probeReducedDataStructure()$combinedDFP_Val_Labels$dfLogFC_w_number))
  #
  # output$condDTVolcano <- DT::renderDataTable(as.data.frame(session$userData$sessionVariables$probeReducedDataStructure()$dfVolcano))
  # output$condPlotVolcano <- plotly::renderPlotly(probeReducedVolcano())
  #
  # output$condDTN <- DT::renderDataTable(as.data.frame(session$userData$sessionVariables$probeReducedDataStructure()$combinedDFP_Val_Labels$dfN_w_number))
  # output$condDTProbes <- DT::renderDataTable(as.data.frame(condDTProbes()),
  #                                    options = list(pageLength = 1000, info = FALSE,
  #                                                   lengthMenu = list(c(100, 1000, -1), c("100", "1000", "All"))))
  # output$condPlotDendrogramProbes <- plotly::renderPlotly(ggdendro::ggdendrogram(session$userData$sessionVariables$probeReducedDataStructure()$clustResProbes,
  #                                                                         rotate = TRUE, theme_dendro = FALSE))
  # output$condDTTraits <- DT::renderDataTable(as.data.frame(condDTTraits()),
  #                                     options = list(pageLength = 1000, info = FALSE,
  #                                                    lengthMenu = list(c(100, 1000, -1), c("100", "1000", "All"))))
  #
  # output$condPlotDendrogramTraits <- plotly::renderPlotly(ggdendro::ggdendrogram(session$userData$sessionVariables$probeReducedDataStructure()$clustResTraits,
  #                                                                         rotate = TRUE, theme_dendro = FALSE))
  # output$condHistP_Val <- plotly::renderPlotly(condHistP_Val())

  #this data structure holds everything (as a named list), that is needed for working with trait reduced (by selecting a subset of traits) HM
  session$userData$sessionVariables$combinedDataStructure <- shiny::reactive({
    shinyId <- shiny::showNotification("Creating combined reduced data structure...", duration = NULL, closeButton = FALSE)
    on.exit(shiny::removeNotification(shinyId), add = TRUE)
     base::tryCatch(
       {
         result <- base::list(combinedDFP_Val_Labels = session$userData$sessionVariables$combinedData()
         )
      },
      error = function(e) {
        if (attributes(e)$class[1] != "shiny.silent.error") {
          base::message("An error occurred in shiny::reactive(session$userData$sessionVariables$combinedDataStructure):\n", e)
          browser() #should not happen
        }
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

  session$userData$sessionVariables$pReducedDataStructure <- shiny::reactive({
    shinyId <- shiny::showNotification("Creating p reduced data structure (p-val)...", duration = NULL, closeButton = FALSE)
    on.exit(shiny::removeNotification(shinyId), add = TRUE)
    base::tryCatch(
      {
        if (is.valid(session$userData$sessionVariables$pReducedData())) {
          result$combinedDFP_Val_Labels <- session$userData$sessionVariables$pReducedData()
          result <- base::list(combinedDFP_Val_Labels = session$userData$sessionVariables$pReducedData(),
                           matP_Val.t = NULL
                           # distMatTraits = NULL,
                           # clustResTraits = NULL,
                           # traitClusters = NULL,
                           # traitClusterMedoids = NULL,
                           # traitDendrogram = NULL,
                           # traitClustergram = NULL,
                           # distMatProbes = NULL,
                           # clustResProbes = NULL,
                           # dendProbes = NULL
          )
        }
        else {
          base::print(paste0(sysTimePID(), " is.valid(session$userData$sessionVariables$pReducedData()) == FALSE"))
          result <- NULL
        }
      },
      error = function(e) {
        if (attributes(e)$class[1] != "shiny.silent.error") {
          base::message("An error occurred in shiny::reactive(session$userData$sessionVariables$pReducedDataStructure):\n", e)
          browser() #should not happen
        }
      },
      warning = function(w) {
        base::message("A warning occurred in shiny::reactive(session$userData$sessionVariables$pReducedDataStructure):\n", w)
        browser() #should not happen
      },
      finally = {
        base::print(base::paste0(sysTimePID(), " finished generating session$userData$sessionVariables$pReducedDataStructure"))
        return(result)
      }
    )
  })

  session$userData$sessionVariables$traitReducedDataStructurePVal <- shiny::reactive({
    shinyId <- shiny::showNotification("Creating trait reduced data structure (p-val)...", duration = NULL, closeButton = FALSE)
    on.exit(shiny::removeNotification(shinyId), add = TRUE)
    return(session$userData$sessionVariables$traitReducedDataPVal())
  })

  session$userData$sessionVariables$traitReducedDataStructureLogFC <- shiny::reactive({
    shinyId <- shiny::showNotification("Creating trait reduced data structure (log(FC))...", duration = NULL, closeButton = FALSE)
    on.exit(shiny::removeNotification(shinyId), add = TRUE)
    return(session$userData$sessionVariables$traitReducedDataLogFC())
  })

  session$userData$sessionVariables$probeReducedDataStructurePVal <- shiny::reactive({
    shinyId <- shiny::showNotification("Creating probe reduced data structure (p-val)...", duration = NULL, closeButton = FALSE)
    on.exit(shiny::removeNotification(shinyId), add = TRUE)
    return(session$userData$sessionVariables$probeReducedDataPVal())
  })

  session$userData$sessionVariables$probeReducedDataStructureWOGapPVal <- shiny::reactive({
    shinyId <- shiny::showNotification("Creating probe reduced data structure (p-val) w/o gap...", duration = NULL, closeButton = FALSE)
    on.exit(shiny::removeNotification(shinyId), add = TRUE)
    return(session$userData$sessionVariables$probeReducedDataWOGapPVal())
  })

  session$userData$sessionVariables$probeReducedDataStructureLogFC <- shiny::reactive({
    shinyId <- shiny::showNotification("Creating probe reduced data structure (log(FC))...", duration = NULL, closeButton = FALSE)
    on.exit(shiny::removeNotification(shinyId), add = TRUE)
    return(session$userData$sessionVariables$probeReducedDataLogFC())
  })

  session$userData$sessionVariables$probeReducedDataStructureWOGapLogFC <- shiny::reactive({
    shinyId <- shiny::showNotification("Creating probe reduced data structure (log(FC)) w/o gap...", duration = NULL, closeButton = FALSE)
    on.exit(shiny::removeNotification(shinyId), add = TRUE)
    return(session$userData$sessionVariables$probeReducedDataWOGapLogFC())
  })

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
          base::print(base::paste0(sysTimePID(), " Step 6a: start plotting heatmap for P_Val. (first step in shiny::observeEvent(input$btnPlotCombinedHM_P_Val))"))
          plotCombinedHM(id = "PVal", input = input, output = output, session = session)
          #          plotHMDNADistances(input = input, output = output, session = session)
          session$userData$sessionVariables$callCounter <- session$userData$sessionVariables$callCounter + 1
        },
        error = function(e) {
          if (attributes(e)$class[1] != "shiny.silent.error") {
            base::message("An error occurred in shiny::observeEvent(input$plotCombinedHM_P_Val):\n", e)
          }
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

  shiny::observeEvent(input$btnPlotCombinedHM_P_ValWOGap,
    ignoreInit = TRUE,
    {
      base::tryCatch(
        {
          base::print(base::paste0(sysTimePID(), " Step 6b: start plotting heatmap for P_Val w/o gaps. (first step in shiny::observeEvent(input$btnPlotCombinedHM_P_ValWOGap))"))
          plotCombinedHM(id = "PValWOGap", input = input, output = output, session = session)
          session$userData$sessionVariables$callCounter <- session$userData$sessionVariables$callCounter + 1
        },
        error = function(e) {
          if (attributes(e)$class[1] != "shiny.silent.error") {
            base::message("An error occurred in shiny::observeEvent(input$btnPlotCombinedHM_P_ValWOGap):\n", e)
          }
        },
        warning = function(w) {
          base::message("A warning occurred in shiny::observeEvent(input$btnPlotCombinedHM_P_ValWOGap):\n", w)
        },
        finally = {
          base::print(base::paste0(sysTimePID(), " finished plotting heatmap for P_Val. (last step in shiny::observeEvent(input$btnPlotCombinedHM_P_ValWOGap))"))
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
          base::print(base::paste0(sysTimePID(), " Step 6c: start plotting heatmap for log(FC). (first step in shiny::observeEvent(input$btnPlotCombinedHM_LogFC))"))
          plotCombinedHM(id = "LogFC", input = input, output = output, session = session)
          session$userData$sessionVariables$callCounter <- session$userData$sessionVariables$callCounter + 1
        },
        error = function(e) {
          if (attributes(e)$class[1] != "shiny.silent.error") {
            base::message("An error occurred in shiny::observeEvent(input$btnPlotCombinedHM_LogFC):\n", e)
          }
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

  shiny::observeEvent(input$btnPlotCombinedHM_LogFCWOGap,
    ignoreInit = TRUE,
    {
      base::tryCatch(
        {
          base::print(base::paste0(sysTimePID(), " Step 6d: start plotting heatmap for log(FC) w/o gaps. (first step in shiny::observeEvent(input$btnPlotCombinedHM_LogFCWOGap))"))
          plotCombinedHM(id = "LogFCWOGap", input = input, output = output, session = session)
          session$userData$sessionVariables$callCounter <- session$userData$sessionVariables$callCounter + 1
        },
        error = function(e) {
          if (attributes(e)$class[1] != "shiny.silent.error") {
            base::message("An error occurred in shiny::observeEvent(input$btnPlotCombinedHM_LogFCWOGap):\n", e)
          }
        },
        warning = function(w) {
          base::message("A warning occurred in shiny::observeEvent(input$btnPlotCombinedHM_LogFCWOGap):\n", w)
        },
        finally = {
          base::print(base::paste0(sysTimePID(), " finished plotting heatmap for log(FC). (last step in shiny::observeEvent(input$btnPlotCombinedHM_LogFCWOGap))"))
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
          if (attributes(e)$class[1] != "shiny.silent.error") {
            base::message("An error occurred in shiny::observeEvent(input$chkDebug):\n", e)
          }
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

