server <- function(input, output, session) {
  #draw empty HM, without the github version won't work. Why?
  InteractiveComplexHeatmap::makeInteractiveComplexHeatmap(
    input = input,
    output = output,
    session = session,
    ht_list = emptyHM(),
    heatmap_id = "heatmap_2",
    show_layer_fun = FALSE,
    click_action = NULL,
    brush_action = NULL,
    hover_action = NULL
  )
  #define sessionVariables here
  reactlog::reactlog_enable()
  pid <- Sys.getpid()
  hostname <- Sys.info()["nodename"]
  output$Sys.PID <- shiny::renderText(base::paste0(hostname, ": ", pid))
  packageWd <<- getwd()
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
    shiny::updateTextInput(session, inputId = "txtSearchCpG", value = knownCpG)
  }
  if (is.valid(session$userData$config$knownTraits)) {
    #load Traits from txt file to search input
    knownTrait <- paste(unlist(data.table::fread(file = session$userData$config$knownTraits, header = FALSE)), collapse = " ")
    shiny::updateTextInput(session, inputId = "txtSearchTrait", value = knownTrait)
  }
  #base::options(spam.force64 = TRUE)

  if (session$userData$config$debugMode == TRUE) {
    shiny::updateCheckboxInput(session, "chkDebug", value = TRUE)
  }
  else {
    shiny::updateCheckboxInput(session, "chkDebug", value = FALSE)
  }

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
  session$userData$sessionVariables$traitReducedData <- shiny::reactiveVal(value = NULL, label = "traitReducedData")
  #  session$userData$sessionVariables$combinedDFP_Val_Labels <- shiny::reactiveVal(value = NULL, label = "combinedDFP_Val_Labels")

  # session$userData$sessionVariables$clustResTraits <- shiny::reactiveVal(value = NULL, label = "clustResTraits")

  # session$userData$sessionVariables$traitReducedmatP_Val <- shiny::reactiveVal(value = NULL, label = "traitReducedmatP_Val")

  #new (4.2024) basic data structure:
  session$userData$sessionVariables$generalDataStructure <- shiny::reactiveVal(value = NULL, label = "generalDataStructure")

  session$userData$sessionVariables$distancesBelowThreshold <- shiny::reactiveVal(value = NULL, label = "distancesBelowThreshold")

  session$userData$sessionVariables$distNeigboursProbes10000 <- shiny::reactiveVal(value = NULL, label = "distNeigboursProbes10000")
  session$userData$sessionVariables$distNeigboursProbes1000 <- shiny::reactiveVal(value = NULL, label = "distNeigboursProbes1000")
  session$userData$sessionVariables$distNeigboursProbes100 <- shiny::reactiveVal(value = NULL, label = "distNeigboursProbes100")
  session$userData$sessionVariables$distNeigboursProbes10 <- shiny::reactiveVal(value = NULL, label = "distNeigboursProbes10")

  session$userData$sessionVariables$selectedOriginalData <- shiny::reactiveVal(value = NULL, label = "selectedOriginalData")
  session$userData$sessionVariables$selectedOriginalDataTraits <- shiny::reactiveVal(value = NULL, label = "selectedOriginalDataTraits")
  session$userData$sessionVariables$selectedOriginalDataProbes <- shiny::reactiveVal(value = NULL, label = "selectedOriginalDataProbes")
  session$userData$sessionVariables$selectedCpG <- shiny::reactiveVal(value = NULL, label = "selectedCpG")
  session$userData$sessionVariables$selectedTrait <- shiny::reactiveVal(value = NULL, label = "selectedTrait")
  session$userData$sessionVariables$selectedTraitReducedcombinedDFP_Val_Labels <- shiny::reactiveVal(value = NULL, label = "selectedTraitReducedcombinedDFP_Val_Labels")

#  session$userData$sessionVariables$SPLOM <- FALSE
  session$userData$sessionVariables$markingVar <- shiny::reactiveVal(value = NULL, label = "markingVar")

  shiny::updateSliderInput(session = session, inputId = "sldP_Val", min = 0, max = 0, value = c(0, 0))
  shiny::updateSliderInput(session = session, inputId = "sldDM", min = 0, max = 0, value = c(0, 0))
  shiny::updateSliderInput(session = session, inputId = "sldN", min = 0, max = 0, value = c(0, 0))

  #base::print(paste0(sysTimePID(), " starting application."))
  base::message(paste0(sysTimePID(), " starting application."))

  shinyjs::toggleClass("colRed", "red")
  shinyjs::toggleClass("colGreen", "green")
  shinyjs::toggleClass("colBlue", "blue")

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

  session$userData$sessionVariables$distNeigboursProbes10 <- shiny::reactive({
    #calculate distance from each probe to its neigbours to build a right column in heatmap to estimate relevance of heatmap findings
    base::tryCatch(
      {
        result <- calculateDistanceNeigboursProbes(wd = session$userData$packageWd, clustResProbes = session$userData$sessionVariables$traitReducedDataStructure()$clustResProbes, annotation = session$userData$annotation, distanceToLook = 10, numCores = session$userData$numCores)
      },
      error = function(e) {
        base::message("An error occurred in shiny::reactive(session$userData$sessionVariables$distNeigboursProbes10):\n", e)
      },
      warning = function(w) {
        base::message("A warning occurred in shiny::reactive(session$userData$sessionVariables$distNeigboursProbes10):\n", w)
      },
      finally = {
        base::print(base::paste0(sysTimePID(), " finished generating session$userData$sessionVariables$distNeigboursProbes10. (last step in session$userData$sessionVariables$distNeigboursProbes10 <- shiny::reactive())"))
        return(result)
      }
    )
  })

  session$userData$sessionVariables$distNeigboursProbes100 <- shiny::reactive({
    #calculate distance from each probe to its neigbours to build a right column in heatmap to estimate relevance of heatmap findings
    base::tryCatch(
      {
        result <- calculateDistanceNeigboursProbes(wd = session$userData$packageWd, clustResProbes = session$userData$sessionVariables$traitReducedDataStructure()$clustResProbes, annotation = session$userData$annotation, distanceToLook = 100, numCores = session$userData$numCores)
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
        result <- calculateDistanceNeigboursProbes(wd = session$userData$packageWd, clustResProbes = session$userData$sessionVariables$traitReducedDataStructure()$clustResProbes, annotation = session$userData$annotation, distanceToLook = 1000, numCores = session$userData$numCores)
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
        result <- calculateDistanceNeigboursProbes(wd = session$userData$packageWd, clustResProbes = session$userData$sessionVariables$traitReducedDataStructure()$clustResProbes, annotation = session$userData$annotation, distanceToLook = 10000, numCores = session$userData$numCores)
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
        rownames(distances) <- distances[,1]
#        distances <- distances[,3]
#        distances <- na.omit(distances)
        threshold <- session$userData$config$CpGDistanceThreshold
        #which distances < threshold
        result <- rownames(distances)[distances$meanDistance < threshold]
        result <- na.omit(result)
#to check, make the first elements (1,3,5) of distances orange...
#        browser()
#        result <- rownames(distances)[c(1,3,5)]
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

  session$userData$sessionVariables$markingVar <- shiny::reactive({

    return(input$markingVar)
  })

  output$SPLOM <- plotly::renderPlotly({
    if (!is.null(session$userData$sessionVariables$selectedOriginalData())) {
      height <- input$numSPLOMVSize
      width <- input$numSPLOMHSize
      base::print(base::paste0(sysTimePID(), " number of traits and probes in SPLOM (columns in selectedDF): ", ncol(session$userData$sessionVariables$selectedOriginalData()))) #thats sum of probes and traits
      base::print(base::paste0(sysTimePID(), " number of cases in SPLOM (rows in selectedDF): ", nrow(session$userData$sessionVariables$selectedOriginalData()))) #thats number of cases in data set
      base::print(base::paste0(sysTimePID(), " number of traits in SPLOM (selectedTraits): ", ncol(session$userData$sessionVariables$selectedOriginalDataTraits())))
      base::print(base::paste0(sysTimePID(), " number of probes in SPLOM: (selectedProbes)", ncol(session$userData$sessionVariables$selectedOriginalDataProbes())))
      fig <- getSPLOM(session$userData$sessionVariables$selectedOriginalData(), XVars = colnames(session$userData$sessionVariables$selectedOriginalDataTraits()), YVars = colnames(session$userData$sessionVariables$selectedOriginalDataProbes()), markingVar = session$userData$sessionVariables$markingVar(), height = height, width = width)
      return(fig)
    }
  })

  output$SPLOMTrait <- plotly::renderPlotly({
    if (!is.null(session$userData$sessionVariables$selectedOriginalDataTraits())) {
      height <- input$numSPLOMVSize
      width <- input$numSPLOMHSize
      fig <- getSPLOM(session$userData$sessionVariables$selectedOriginalDataTraits(), XVars = session$userData$sessionVariables$selectedOriginalDataTraits(), YVars = session$userData$sessionVariables$selectedOriginalDataTraits(), markingVar = session$userData$sessionVariables$markingVar(), height = height, width = width)
      return(fig)
    }
  })

  output$SPLOMProbe <- plotly::renderPlotly({
    if (!is.null(session$userData$sessionVariables$selectedOriginalDataProbes())) {
      height <- input$numSPLOMVSize
      width <- input$numSPLOMHSize
      fig <- getSPLOM(session$userData$sessionVariables$selectedOriginalDataProbes(), XVars = session$userData$sessionVariables$selectedOriginalDataProbes(), YVars = session$userData$sessionVariables$selectedOriginalDataProbes(), markingVar = session$userData$sessionVariables$markingVar(), height = height, width = width)
      return(fig)
    }
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
#        if (is.valid(session$userData$sessionVariables$combinedDataStructure())) {
          result <- updateTxtMergeOut(session$userData$sessionVariables$combinedDataStructure())
#        }
#        else {
#          result <- NULL
#        }
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

  # counter.output.txtPReduceOut <- local({
  #   static <- 0
  #   function() { static <<- static + 1; static }
  # })

  output$txtPReduceOut <- shiny::reactive({
    base::tryCatch(
      {
        base::print(base::paste0(sysTimePID(), " start generating output$txtPReduceOut."))
        #if (is.valid(session$userData$sessionVariables$pReducedData())) {
        if (is.valid(session$userData$sessionVariables$pReducedDataStructure()$combinedDFP_Val_Labels$dfP_Val)) {
#        maxTraits <- ncol(session$userData$sessionVariables$pReducedDataStructure()$combinedDFP_Val_Labels$dfP_Val)
#        value <- maxTraits
        #if (is.valid(session$userData$sessionVariables$traitReducedmatP_Val())) {
        # if (is.valid(session$userData$sessionVariables$traitReducedDataStructure()$combinedDFP_Val_Labels$dfP_Val)) {
        #   value <- ncol(session$userData$sessionVariables$traitReducedDataStructure()$combinedDFP_Val_Labels$dfP_Val)
        # }
        # else {
        #   value <- maxTraits
        # }
#        shiny::updateSliderInput(session = session, inputId = "sldNumClusters", max = maxTraits, min = 1, value = value, step = 1)
        result <- updateTxtpReduceOut(session$userData$sessionVariables$pReducedDataStructure()$combinedDFP_Val_Labels)
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
#        counter.output.txtPReduceOut()
        base::print(base::paste0(sysTimePID(), " finished generating output$txtPReduceOut."))
        return(result)
      }
    )
  })

  output$txtOmitOut <- shiny::reactive({
    base::tryCatch(
      {
        base::print(base::paste0(sysTimePID(), " start generating output$txtOmitOut."))
        if (is.valid(session$userData$sessionVariables$traitReducedDataStructure()$combinedDFP_Val_Labels)) {
          result <- updateTxtOmitTraitsOut(session$userData$sessionVariables$traitReducedDataStructure()$combinedDFP_Val_Labels)
        }
        else {
          result <- NULL
        }
      },
      error = function(e) {
        base::message("An error occurred in shiny::reactive(output$txtOmitOut):\n", e)
      },
      warning = function(w) {
        base::message("A warning occurred in shiny::reactive(output$txtOmitOut):\n", w)
      },
      finally = {
        base::print(base::paste0(sysTimePID(), " finished generating output$txtOmitOut."))
        return(result)
      }
    )
  })

  output$plotDendrogramTraitsLong <- shiny::renderPlot(getPlot(session$userData$sessionVariables$traitReducedDataStructure()$traitDendrogram))

  output$plotClustergramTraitsLong <- shiny::renderPlot(getPlot(session$userData$sessionVariables$traitReducedDataStructure()$traitClustergram))

  output$DTTraitsMedoids <- DT::renderDataTable(as.data.frame(getMedoidsTable(session$userData$sessionVariables$traitReducedDataStructure()$traitClusterMedoids)))

  output$DTTraitsClusters <- DT::renderDataTable(as.data.frame(getClustersTable(session$userData$sessionVariables$traitReducedDataStructure()$traitClusters,
                                                                                session$userData$sessionVariables$traitReducedDataStructure()$traitClusterMedoids)))
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

  shiny::observeEvent(input$numHMHSize,
    ignoreInit = TRUE,
    {
#      session$userData$sessionVariables$markingVar <- input$markingVar
      # redraw HM
      #browser()
      height <- input$numHMVSize
      width <- input$numHMHSize
      if (session$userData$sessionVariables$callCounter > 1) {
browser() #check, whether this is called initially and why plotCombinedHM_P_Val is called twice
      }
      plotCombinedHM_P_Val(input = input, output = output, session = session)
    },
    ignoreNULL = FALSE
  )
  shiny::observeEvent(input$numHMVSize,
    ignoreInit = TRUE,
    {
#      session$userData$sessionVariables$markingVar <- input$markingVar
      # redraw HM
      #browser()
      height <- input$numHMVSize
      width <- input$numHMHSize
      if (session$userData$sessionVariables$callCounter > 1) {
browser() #check, whether this is called initially and why plotCombinedHM_P_Val is called twice
      }
      plotCombinedHM_P_Val(input = input, output = output, session = session)
    },
    ignoreNULL = FALSE
  )

  shiny::observeEvent(input$btnBrowser,
    ignoreInit = TRUE,
    {
      browser()
    },
    ignoreNULL = FALSE
  )
  shiny::observeEvent(input$btnLoadDir1,
    ignoreInit = TRUE,
    {
      base::tryCatch(
        {
          base::print(base::paste0(sysTimePID(), " start loading trait 1 folders."))
          base::print(base::paste0(sysTimePID(), " creating empty heatmap during load process."))
          combinedHMP_VAL <- emptyHM()
          InteractiveComplexHeatmap::makeInteractiveComplexHeatmap(
            input = input,
            output = output,
            session = session,
            ht_list = combinedHMP_VAL,
            heatmap_id = "Heatmap_P_Val",
            show_layer_fun = FALSE,
            click_action = NULL,
            brush_action = NULL,
            hover_action = NULL
            )
          base::print(base::paste0(sysTimePID(), " before is.numeric()."))
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
          base::print(base::paste0(sysTimePID(), " creating empty heatmap during load process."))
          combinedHMP_VAL <- emptyHM()
          InteractiveComplexHeatmap::makeInteractiveComplexHeatmap(
            input = input,
            output = output,
            session = session,
            ht_list = combinedHMP_VAL,
            heatmap_id = "Heatmap_P_Val",
            show_layer_fun = FALSE,
            click_action = NULL,
            brush_action = NULL,
            hover_action = NULL
            )
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
          base::print(base::paste0(sysTimePID(), " creating empty heatmap during load process."))
          combinedHMP_VAL <- emptyHM()
          InteractiveComplexHeatmap::makeInteractiveComplexHeatmap(
            input = input,
            output = output,
            session = session,
            ht_list = combinedHMP_VAL,
            heatmap_id = "Heatmap_P_Val",
            show_layer_fun = FALSE,
            click_action = NULL,
            brush_action = NULL,
            hover_action = NULL
            )
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
          base::print(base::paste0(sysTimePID(), " creating empty heatmap."))
          combinedHMP_VAL <- emptyHM()
          InteractiveComplexHeatmap::makeInteractiveComplexHeatmap(
            input = input,
            output = output,
            session = session,
            ht_list = combinedHMP_VAL,
            heatmap_id = "Heatmap_P_Val",
            show_layer_fun = FALSE,
            click_action = NULL,
            brush_action = NULL,
            hover_action = NULL
            )
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
          base::print(base::paste0(sysTimePID(), " creating empty heatmap."))
          # combinedHMP_VAL <- emptyHM()
          # InteractiveComplexHeatmap::makeInteractiveComplexHeatmap(input, output, session, combinedHMP_VAL, "Heatmap_P_Val")
          minN <- base::as.integer(input$txtCases)
          #if (session$userData$sessionVariables$LoadInitialized == TRUE) {
          if (is.valid(session$userData$sessionVariables$resultDFListTrait1()) || is.valid(session$userData$sessionVariables$resultDFListTrait2())  || is.valid(session$userData$sessionVariables$resultDFListTrait3())) {
            combinedDFP_Val_Labels <- mergeDFP_Val_Labels(session$userData$sessionVariables$resultDFListTrait1(),
                                                          session$userData$sessionVariables$resultDFListTrait2(),
                                                          session$userData$sessionVariables$resultDFListTrait3(),
                                                          minN)
          }
          else {
            base::print(base::paste0(sysTimePID(), " is.valid(session$userData$sessionVariables$resultDFListTrait1()) || is.valid(session$userData$sessionVariables$resultDFListTrait2())  || is.valid(session$userData$sessionVariables$resultDFListTrait3()) == FALSE."))
            result <- NULL
          }
          session$userData$sessionVariables$combinedData(combinedDFP_Val_Labels)
          updateSliders(session, session$userData$sessionVariables$combinedDataStructure()$combinedDFP_Val_Labels)
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
          base::print(base::paste0(sysTimePID(), " Step 3: start reducing data by p-value."))
          base::print(base::paste0(sysTimePID(), " creating empty heatmap."))
          minP_Val <- 5 * 10^base::as.integer(input$sldP_Val[1]) #minP_Val <- 5 * 10^-base::as.integer(input$sldP_Val[2])
          maxP_Val <- 5 * 10^base::as.integer(input$sldP_Val[2]) #maxP_Val <- 5 * 10^-base::as.integer(input$sldP_Val[1])
          if (maxP_Val < minP_Val) { #exchange, if in wrong order
            t <- minP_Val
            minP_Val <- maxP_Val
            maxP_Val <- t
            browser()
          }
          minDM <- input$sldDM[1]
          maxDM <- input$sldDM[2]
          minN <- base::as.integer(input$sldN[1])
          maxN <- base::as.integer(input$sldN[2])
          combinedDFP_Val_Labels <- session$userData$sessionVariables$combinedDataStructure()$combinedDFP_Val_Labels
          if (is.valid(combinedDFP_Val_Labels)) {
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
                shiny::updateSliderInput(session = session, inputId = "sldNumClusters", max = maxTraits, min = 1, value = value, step = 1)
              }
            }
            else {
              result <- NULL
              base::print(base::paste0(sysTimePID(), " minP_Val == maxP_Val && minDM == maxDM && minN == maxN."))
            }
          }
          else {
            result <- NULL
            base::print(base::paste0(sysTimePID(), " is.valid(combinedDFP_Val_Labels) == FALSE."))
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

  shiny::observeEvent(input$btnOmitTraits,
    ignoreInit = TRUE,
    {
      base::tryCatch(
        {
#browser()
          if (is.valid(session$userData$sessionVariables$pReducedDataStructure()$combinedDFP_Val_Labels$dfP_Val)) {
            if (is.valid(input$sldNumClusters)) {
              result <- session$userData$sessionVariables$pReducedDataStructure()
              base::print(base::paste0(sysTimePID(), " (btnOmitTraits) start generating traitClusters."))
              if (is.valid(result$clustResTraits) && input$sldNumClusters > 1) {
                traitClusters <- cutree(result$clustResTraits,
                                               k = input$sldNumClusters)
              }
              else {
                traitClusters <- NULL
              }
              base::print(base::paste0(sysTimePID(), " (btnOmitTraits) start generating traitClusterMedoids."))
              if (is.valid(result$clustResTraits) && is.valid(result$distMatTraits)) {
                traitClusterMedoids <- getTraitClusterMedoids(clustResTraits = result$clustResTraits,
                                                                     distMatTraits = result$distMatTraits,
                                                                     numClusters = input$sldNumClusters)
              }
              else {
                traitClusterMedoids <- NULL
              }
            }
          }
          if (is.valid(session$userData$sessionVariables$pReducedDataStructure()$combinedDFP_Val_Labels$dfP_Val) &&
              is.valid(traitClusterMedoids)) {
            keys <- session$userData$config$keyAttributes
            result <- getTraitReducedData(session$userData$sessionVariables$pReducedDataStructure()$combinedDFP_Val_Labels,
                                          traitClusterMedoids, keys)
          }
          else{
            result <- NULL
          }
        },
        error = function(e) {
          base::message("An error occurred in shiny::observeEvent(input$btnOmitTraits):\n", e)
        },
        warning = function(w) {
          base::message("A warning occurred in shiny::observeEvent(input$btnOmitTraits):\n", w)
        },
        finally = {
          session$userData$sessionVariables$traitReducedData(result)
          base::print(base::paste0(sysTimePID(), " finished omitting traits."))
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
            #output$plotDendrogramP_VALborder <- plotly::renderPlotly(plotly::plot_ly(x = P_Val, type = "histogram", name = "histP_ValBorder"))
            plot <- plotly::plot_ly(x = P_VALNTable$P_VAL_BORDER, y = P_VALNTable$'Available CpG' , type = "scatter", mode = "lines+markers", name = "scatterP_ValBorder") %>%
              plotly::layout(xaxis = list(title = 'p-val', type = "log")) %>%
              plotly::layout(yaxis = list(title = 'n' ))
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
              plotly::layout(xaxis = list(title = 'DeltaMethylation', type = "linear")) %>%
              plotly::layout(yaxis = list(title = 'n' ))
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
            plot <- plotly::plot_ly(x = NNTable$N_BORDER, y = NNTable$'Available CpG' , type = "scatter", mode = "lines+markers", name = "scatterNBorder") %>%
              plotly::layout(xaxis = list(title = 'n', type = "linear")) %>%
              plotly::layout(yaxis = list(title = 'n' ))
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

  shiny::observeEvent(input$btnCountP_ValProbes,
    ignoreInit = TRUE,
    {
      base::tryCatch(
        {
          base::print(base::paste0(sysTimePID(), " start counting probes."))
          minN <- base::as.integer(input$txtCases)
          #if (session$userData$sessionVariables$LoadInitialized == TRUE) {
          if (is.valid(session$userData$sessionVariables$resultDFListTrait1()) || is.valid(session$userData$sessionVariables$resultDFListTrait2())  || is.valid(session$userData$sessionVariables$resultDFListTrait3())) {
            combinedDFP_Val_Labels <- mergeDFP_Val_Labels(session$userData$sessionVariables$resultDFListTrait1(),
                                                          session$userData$sessionVariables$resultDFListTrait2(),
                                                          session$userData$sessionVariables$resultDFListTrait3(), minN)
          }
          else {
            base::print(base::paste0(sysTimePID(), " is.valid(session$userData$sessionVariables$resultDFListTrait1()) || is.valid(session$userData$sessionVariables$resultDFListTrait2())  || is.valid(session$userData$sessionVariables$resultDFListTrait3()) == FALSE."))
            result <- NULL
          }
          P_VALNTable <-
            getAvailNForP_VALBorder(combinedDFP_Val_Labels$dfP_Val)
          output$DTP_VALborder <- DT::renderDataTable(P_VALNTable)
          base::print(base::paste0(sysTimePID(), " finished counting probes."))
        },
        error = function(e) {
          base::message("An error occurred in shiny::observeEvent(input$btnCountP_ValProbes):\n", e)
        },
        warning = function(w) {
          base::message("A warning occurred in shiny::observeEvent(input$btnCountP_ValProbes):\n", w)
        },
        finally = {
          base::print(base::paste0(sysTimePID(), " finished counting probes."))
        }
      )
    },
    ignoreNULL = FALSE
  )

  output$DTSelectedP_Val <- DT::renderDataTable(as.data.frame(session$userData$sessionVariables$selectedTraitReducedcombinedDFP_Val_Labels()$dfP_Val), escape = FALSE, extensions = c("Buttons"), options = list(dom = "Bfrtip", buttons = c("csv"), pageLength = 10000))  #TBC()

  output$DTSelectedTrait <- DT::renderDataTable(as.data.frame(session$userData$sessionVariables$selectedTrait()), escape = FALSE, extensions = c("Buttons"), options = list(dom = "Bfrtip", buttons = c("csv"), pageLength = 10000))

  output$selectHistDM <- plotly::renderPlotly(selectHistDM())

  output$histP_Val <- plotly::renderPlotly(histP_Val())

  output$histMinDistance10 <- plotly::renderPlotly(histMinDistance10())
  output$histMeanDistance10 <- plotly::renderPlotly(histMeanDistance10())
  output$histMaxDistance10 <- plotly::renderPlotly(histMaxDistance10())

  output$histMinDistance100 <- plotly::renderPlotly(histMinDistance100())
  output$histMeanDistance100 <- plotly::renderPlotly(histMeanDistance100())
  output$histMaxDistance100 <- plotly::renderPlotly(histMaxDistance100())

  output$histMinDistance1000 <- plotly::renderPlotly(histMinDistance1000())
  output$histMeanDistance1000 <- plotly::renderPlotly(histMeanDistance1000())
  output$histMaxDistance1000 <- plotly::renderPlotly(histMaxDistance1000())

  output$histMinDistance10000 <- plotly::renderPlotly(histMinDistance10000())
  output$histMeanDistance10000 <- plotly::renderPlotly(histMeanDistance10000())
  output$histMaxDistance10000 <- plotly::renderPlotly(histMaxDistance10000())

  output$DTDistance10 <- DT::renderDataTable(DTDistance10(), escape = FALSE, extensions = c("Buttons"), options = list(dom = "Bfrtip", buttons = c("csv"), pageLength = 10000))
  output$DTDistance10reduced <- DT::renderDataTable(DTDistance10reduced(), escape = FALSE, extensions = c("Buttons"), options = list(dom = "Bfrtip", buttons = c("csv"), pageLength = 10000))

  output$DTDistance100 <- DT::renderDataTable(DTDistance100(), escape = FALSE, extensions = c("Buttons"), options = list(dom = "Bfrtip", buttons = c("csv"), pageLength = 10000))
  output$DTDistance100reduced <- DT::renderDataTable(DTDistance100reduced(), escape = FALSE, extensions = c("Buttons"), options = list(dom = "Bfrtip", buttons = c("csv"), pageLength = 10000))

  output$DTDistance1000 <- DT::renderDataTable(DTDistance1000(), escape = FALSE, extensions = c("Buttons"), options = list(dom = "Bfrtip", buttons = c("csv"), pageLength = 10000))
  output$DTDistance1000reduced <- DT::renderDataTable(DTDistance1000reduced(), escape = FALSE, extensions = c("Buttons"), options = list(dom = "Bfrtip", buttons = c("csv"), pageLength = 10000))

  output$DTDistance10000 <- DT::renderDataTable(DTDistance10000(), escape = FALSE, extensions = c("Buttons"), options = list(dom = "Bfrtip", buttons = c("csv"), pageLength = 10000))
  output$DTDistance10000reduced <- DT::renderDataTable(DTDistance10000reduced(), escape = FALSE, extensions = c("Buttons"), options = list(dom = "Bfrtip", buttons = c("csv"), pageLength = 10000))

  selectHistDM <- shiny::reactive({
    base::tryCatch({
      base::print(base::paste0(sysTimePID(), " start render plotly selectHistDM(). (first step in renderPlotly(selectHistDM()))"))
      if(is.valid(session$userData$sessionVariables$selectedTraitReducedcombinedDFP_Val_Labels()$dfDM)) {
        DM <- sort(as.numeric(unlist(session$userData$sessionVariables$selectedTraitReducedcombinedDFP_Val_Labels()$dfDM)))
        result <- plotly::plot_ly(x = DM, type = "histogram", name = "selectHistDM")
      }
      else {
        result <- NULL
      }
    },
    error = function(e) {
      base::message("An error occurred in selectHistDM <- shiny::reactive():\n", e)
      browser()
    },
    warning = function(w) {
      base::message("A warning occurred in selectHistDM <- shiny::reactive():\n", w)
    },
    finally = {
      base::print(base::paste0(sysTimePID(), " finished render plotly selectHistDM(). (last step in renderPlotly(selectHistDM()))"))
      return(result)
    })
  })

  histMinDistance10 <- shiny::reactive({
    base::tryCatch({
      base::print(base::paste0(sysTimePID(), " start render plotly histMeanDistance10(). (first step in renderPlotly(histMeanDistance10()))"))
      MinDistance <- session$userData$sessionVariables$distNeigboursProbes10()[,2]
      result <- plotly::plot_ly(x = MinDistance, type = "histogram", name = "histMinDistance10")
    },
    error = function(e) {
      base::message("An error occurred in histMinDistance10 <- shiny::reactive():\n", e)
      browser()
    },
    warning = function(w) {
      base::message("A warning occurred in histMinDistance10 <- shiny::reactive():\n", w)
    },
    finally = {
      base::print(base::paste0(sysTimePID(), " finished render plotly histMinDistance10(). (last step in renderPlotly(histMinDistance10()))"))
      return(result)
    })
  })

  histMeanDistance10 <- shiny::reactive({
    base::tryCatch({
      base::print(base::paste0(sysTimePID(), " start render plotly histMeanDistance10(). (first step in renderPlotly(histMeanDistance10()))"))
      MeanDistance <- session$userData$sessionVariables$distNeigboursProbes10()[,3]
      result <- plotly::plot_ly(x = MeanDistance, type = "histogram", name = "histMeanDistance10")
    },
    error = function(e) {
      base::message("An error occurred in histMeanDistance10 <- shiny::reactive():\n", e)
      browser()
    },
    warning = function(w) {
      base::message("A warning occurred in histMeanDistance10 <- shiny::reactive():\n", w)
    },
    finally = {
      base::print(base::paste0(sysTimePID(), " finished render plotly histMeanDistance10(). (last step in renderPlotly(histMeanDistance10()))"))
      return(result)
    })
  })

  histMaxDistance10 <- shiny::reactive({
    base::tryCatch({
      base::print(base::paste0(sysTimePID(), " start render plotly histMaxDistance10(). (first step in renderPlotly(histMaxDistance10()))"))
      MaxDistance <- session$userData$sessionVariables$distNeigboursProbes10()[,4]
      result <- plotly::plot_ly(x = MaxDistance, type = "histogram", name = "histMaxDistance10")
    },
    error = function(e) {
      base::message("An error occurred in histMaxDistance10 <- shiny::reactive():\n", e)
      browser()
    },
    warning = function(w) {
      base::message("A warning occurred in histMaxDistance10 <- shiny::reactive():\n", w)
    },
    finally = {
      base::print(base::paste0(sysTimePID(), " finished render plotly histMaxDistance10(). (last step in renderPlotly(histMaxDistance10()))"))
      return(result)
    })
  })

  histMinDistance100 <- shiny::reactive({
    base::tryCatch({
      base::print(base::paste0(sysTimePID(), " start render plotly histMinDistance100(). (first step in renderPlotly(histMinDistance100()))"))
      MinDistance <- session$userData$sessionVariables$distNeigboursProbes100()[,2]
      result <- plotly::plot_ly(x = MinDistance, type = "histogram", name = "histMinDistance100")
    },
    error = function(e) {
      base::message("An error occurred in histMinDistance100 <- shiny::reactive():\n", e)
      browser()
    },
    warning = function(w) {
      base::message("A warning occurred in histMinDistance100 <- shiny::reactive():\n", w)
    },
    finally = {
      base::print(base::paste0(sysTimePID(), " finished render plotly histMinDistance100(). (last step in renderPlotly(histMinDistance100()))"))
      return(result)
    })
  })

  histMeanDistance100 <- shiny::reactive({
    base::tryCatch({
      base::print(base::paste0(sysTimePID(), " start render plotly histMeanDistance100(). (first step in renderPlotly(histMeanDistance100()))"))
      MeanDistance <- session$userData$sessionVariables$distNeigboursProbes100()[,3]
      result <- plotly::plot_ly(x = MeanDistance, type = "histogram", name = "histMeanDistance100")
    },
    error = function(e) {
      base::message("An error occurred in histMeanDistance100 <- shiny::reactive():\n", e)
      browser()
    },
    warning = function(w) {
      base::message("A warning occurred in histMeanDistance100 <- shiny::reactive():\n", w)
    },
    finally = {
      base::print(base::paste0(sysTimePID(), " finished render plotly histMeanDistance100(). (last step in renderPlotly(histMeanDistance100()))"))
      return(result)
    })
  })

  histMaxDistance100 <- shiny::reactive({
    base::tryCatch({
      base::print(base::paste0(sysTimePID(), " start render plotly histMaxDistance100(). (first step in renderPlotly(histMaxDistance100()))"))
      MaxDistance <- session$userData$sessionVariables$distNeigboursProbes100()[,4]
      result <- plotly::plot_ly(x = MaxDistance, type = "histogram", name = "histMaxDistance100")
    },
    error = function(e) {
      base::message("An error occurred in histMaxDistance100 <- shiny::reactive():\n", e)
      browser()
    },
    warning = function(w) {
      base::message("A warning occurred in histMaxDistance100 <- shiny::reactive():\n", w)
    },
    finally = {
      base::print(base::paste0(sysTimePID(), " finished render plotly histMaxDistance100(). (last step in renderPlotly(histMaxDistance100()))"))
      return(result)
    })
  })

  histMinDistance1000 <- shiny::reactive({
    base::tryCatch({
      base::print(base::paste0(sysTimePID(), " start render plotly histMinDistance1000(). (first step in renderPlotly(histMinDistance1000()))"))
      MinDistance <- session$userData$sessionVariables$distNeigboursProbes1000()[,2]
      result <- plotly::plot_ly(x = MinDistance, type = "histogram", name = "histMinDistance1000")
    },
    error = function(e) {
      base::message("An error occurred in histMinDistance1000 <- shiny::reactive():\n", e)
      browser()
    },
    warning = function(w) {
      base::message("A warning occurred in histMinDistance1000 <- shiny::reactive():\n", w)
    },
    finally = {
      base::print(base::paste0(sysTimePID(), " finished render plotly histMinDistance1000(). (last step in renderPlotly(histMinDistance1000()))"))
      return(result)
    })
  })

  histMeanDistance1000 <- shiny::reactive({
    base::tryCatch({
      base::print(base::paste0(sysTimePID(), " start render plotly histMeanDistance1000(). (first step in renderPlotly(histMeanDistance1000()))"))
      MeanDistance <- session$userData$sessionVariables$distNeigboursProbes1000()[,3]
      result <- plotly::plot_ly(x = MeanDistance, type = "histogram", name = "histMeanDistance1000")
    },
    error = function(e) {
      base::message("An error occurred in histMeanDistance1000 <- shiny::reactive():\n", e)
      browser()
    },
    warning = function(w) {
      base::message("A warning occurred in histMeanDistance1000 <- shiny::reactive():\n", w)
    },
    finally = {
      base::print(base::paste0(sysTimePID(), " finished render plotly histMeanDistance1000(). (last step in renderPlotly(histMeanDistance1000()))"))
      return(result)
    })
  })

  histMaxDistance1000 <- shiny::reactive({
    base::tryCatch({
      base::print(base::paste0(sysTimePID(), " start render plotly histMaxDistance1000(). (first step in renderPlotly(histMaxDistance1000()))"))
      MaxDistance <- session$userData$sessionVariables$distNeigboursProbes1000()[,4]
      result <- plotly::plot_ly(x = MaxDistance, type = "histogram", name = "histMaxDistance1000")
    },
    error = function(e) {
      base::message("An error occurred in histMaxDistance1000 <- shiny::reactive():\n", e)
      browser()
    },
    warning = function(w) {
      base::message("A warning occurred in histMaxDistance1000 <- shiny::reactive():\n", w)
    },
    finally = {
      base::print(base::paste0(sysTimePID(), " finished render plotly histMaxDistance1000(). (last step in renderPlotly(histMaxDistance1000()))"))
      return(result)
    })
  })

  histMinDistance10000 <- shiny::reactive({
    base::tryCatch({
      base::print(base::paste0(sysTimePID(), " start render plotly histMinDistance10000(). (first step in renderPlotly(histMinDistance10000()))"))
      MinDistance <- session$userData$sessionVariables$distNeigboursProbes10000()[,2]
      result <- plotly::plot_ly(x = MinDistance, type = "histogram", name = "histMinDistance10000")
    },
    error = function(e) {
      base::message("An error occurred in histMinDistance10000 <- shiny::reactive():\n", e)
      browser()
    },
    warning = function(w) {
      base::message("A warning occurred in histMinDistance10000 <- shiny::reactive():\n", w)
    },
    finally = {
      base::print(base::paste0(sysTimePID(), " finished render plotly histMinDistance10000(). (last step in renderPlotly(histMinDistance10000()))"))
      return(result)
    })
  })

  histMeanDistance10000 <- shiny::reactive({
    base::tryCatch({
      base::print(base::paste0(sysTimePID(), " start render plotly histMeanDistance10000(). (first step in renderPlotly(histMeanDistance10000()))"))
      MeanDistance <- session$userData$sessionVariables$distNeigboursProbes10000()[,3]
      result <- plotly::plot_ly(x = MeanDistance, type = "histogram", name = "histMeanDistance10000")
    },
    error = function(e) {
      base::message("An error occurred in histMeanDistance10000 <- shiny::reactive():\n", e)
      browser()
    },
    warning = function(w) {
      base::message("A warning occurred in histMeanDistance10000 <- shiny::reactive():\n", w)
    },
    finally = {
      base::print(base::paste0(sysTimePID(), " finished render plotly histMeanDistance10000(). (last step in renderPlotly(histMeanDistance10000()))"))
      return(result)
    })
  })

  histMaxDistance10000 <- shiny::reactive({
    base::tryCatch({
      base::print(base::paste0(sysTimePID(), " start render plotly histMaxDistance10000(). (first step in renderPlotly(histMaxDistance10000()))"))
      MaxDistance <- session$userData$sessionVariables$distNeigboursProbes10000()[,4]
      result <- plotly::plot_ly(x = MaxDistance, type = "histogram", name = "histMaxDistance10000")
    },
    error = function(e) {
      base::message("An error occurred in histMaxDistance10000 <- shiny::reactive():\n", e)
      browser()
    },
    warning = function(w) {
      base::message("A warning occurred in histMaxDistance10000 <- shiny::reactive():\n", w)
    },
    finally = {
      base::print(base::paste0(sysTimePID(), " finished render plotly histMaxDistance10000(). (last step in renderPlotly(histMaxDistance10000()))"))
      return(result)
    })
  })

  DTDistance10 <- shiny::reactive({
    base::tryCatch({
      base::print(base::paste0(sysTimePID(), " start render DT DTMeanDistance10()."))
      DNAdistances <- session$userData$sessionVariables$distNeigboursProbes10()
      result <- DNAdistances
    },
    error = function(e) {
      base::message("An error occurred in DTMeanDistance10 <- shiny::reactive():\n", e)
      browser()
    },
    warning = function(w) {
      base::message("A warning occurred in DTMeanDistance10 <- shiny::reactive():\n", w)
    },
    finally = {
      base::print(base::paste0(sysTimePID(), " finished render DT DTMeanDistance10()"))
      return(result)
    })
  })

  DTDistance10reduced <- shiny::reactive({
    base::tryCatch({
      base::print(base::paste0(sysTimePID(), " start render DT DTMeanDistance10()."))
      DNAdistances <- session$userData$sessionVariables$distNeigboursProbes10()
      DNAdistances <- na.omit(DNAdistances)
      result <- DNAdistances
    },
    error = function(e) {
      base::message("An error occurred in DTMeanDistance10 <- shiny::reactive():\n", e)
      browser()
      a <- object.size(session$userData$sessionVariables$traitReducedmatP_Val())
    },
    warning = function(w) {
      base::message("A warning occurred in DTMeanDistance10 <- shiny::reactive():\n", w)
    },
    finally = {
      base::print(base::paste0(sysTimePID(), " finished render DT DTMeanDistance10()"))
      return(result)
    })
  })

  DTDistance100 <- shiny::reactive({
    base::tryCatch({
      base::print(base::paste0(sysTimePID(), " start render DT DTMeanDistance100()."))
      DNAdistances <- session$userData$sessionVariables$distNeigboursProbes100()
      result <- DNAdistances
    },
    error = function(e) {
      base::message("An error occurred in DTMeanDistance100 <- shiny::reactive():\n", e)
      browser()
    },
    warning = function(w) {
      base::message("A warning occurred in DTMeanDistance100 <- shiny::reactive():\n", w)
    },
    finally = {
      base::print(base::paste0(sysTimePID(), " finished render DT DTMeanDistance100()"))
      return(result)
    })
  })

  DTDistance100reduced <- shiny::reactive({
    base::tryCatch({
      base::print(base::paste0(sysTimePID(), " start render DT DTMeanDistance100()."))
      DNAdistances <- session$userData$sessionVariables$distNeigboursProbes100()
      DNAdistances <- na.omit(DNAdistances)
      result <- DNAdistances
    },
    error = function(e) {
      base::message("An error occurred in DTMeanDistance100 <- shiny::reactive():\n", e)
      browser()
    },
    warning = function(w) {
      base::message("A warning occurred in DTMeanDistance100 <- shiny::reactive():\n", w)
    },
    finally = {
      base::print(base::paste0(sysTimePID(), " finished render DT DTMeanDistance100()"))
      return(result)
    })
  })

  DTDistance1000 <- shiny::reactive({
    base::tryCatch({
      base::print(base::paste0(sysTimePID(), " start render DT DTMeanDistance1000()."))
      DNAdistances <- session$userData$sessionVariables$distNeigboursProbes1000()
      result <- DNAdistances
    },
    error = function(e) {
      base::message("An error occurred in DTMeanDistance1000 <- shiny::reactive():\n", e)
      browser()
    },
    warning = function(w) {
      base::message("A warning occurred in DTMeanDistance1000 <- shiny::reactive():\n", w)
    },
    finally = {
      base::print(base::paste0(sysTimePID(), " finished render DT DTMeanDistance1000()"))
      return(result)
    })
  })

  DTDistance1000reduced <- shiny::reactive({
    base::tryCatch({
      base::print(base::paste0(sysTimePID(), " start render DT DTMeanDistance1000()."))
      DNAdistances <- session$userData$sessionVariables$distNeigboursProbes1000()
      DNAdistances <- na.omit(DNAdistances)
      result <- DNAdistances
    },
    error = function(e) {
      base::message("An error occurred in DTMeanDistance1000 <- shiny::reactive():\n", e)
      browser()
    },
    warning = function(w) {
      base::message("A warning occurred in DTMeanDistance1000 <- shiny::reactive():\n", w)
    },
    finally = {
      base::print(base::paste0(sysTimePID(), " finished render DT DTMeanDistance1000()"))
      return(result)
    })
  })

  DTDistance10000 <- shiny::reactive({
    base::tryCatch({
      base::print(base::paste0(sysTimePID(), " start render DT DTMeanDistance10000()."))
      DNAdistances <- session$userData$sessionVariables$distNeigboursProbes10000()
      result <- DNAdistances
    },
    error = function(e) {
      base::message("An error occurred in DTMeanDistance10000 <- shiny::reactive():\n", e)
      browser()

    },
    warning = function(w) {
      base::message("A warning occurred in DTMeanDistance10000 <- shiny::reactive():\n", w)
    },
    finally = {
      base::print(base::paste0(sysTimePID(), " finished render DT DTMeanDistance10000()"))
      return(result)
    })
  })

  DTDistance10000reduced <- shiny::reactive({
    base::tryCatch({
      base::print(base::paste0(sysTimePID(), " start render DT DTMeanDistance10000()."))
      DNAdistances <- session$userData$sessionVariables$distNeigboursProbes10000()
      DNAdistances <- na.omit(DNAdistances)
      result <- DNAdistances
    },
    error = function(e) {
      base::message("An error occurred in DTMeanDistance10000 <- shiny::reactive():\n", e)
      browser()
    },
    warning = function(w) {
      base::message("A warning occurred in DTMeanDistance10000 <- shiny::reactive():\n", w)
    },
    finally = {
      base::print(base::paste0(sysTimePID(), " finished render DT DTMeanDistance10000()"))
      return(result)
    })
  })

  histP_Val <- shiny::reactive({
    base::tryCatch({
      base::print(base::paste0(sysTimePID(), " start render plotly histP_Val()."))
      #P_Val <- sort(as.numeric(unlist(session$userData$sessionVariables$traitReducedDataStructure()$matP_Val())))
      P_Val <- sort(as.numeric(unlist(session$userData$sessionVariables$traitReducedDataStructure()$combinedDFP_Val_Labels$dfP_Val)))
      result <- plotly::plot_ly(x = P_Val, type = "histogram", name = "histP_Val")
    },
    error = function(e) {
      base::message("An error occurred in histP_Val <- shiny::reactive():\n", e)
      browser()
    },
    warning = function(w) {
      base::message("A warning occurred in histP_Val <- shiny::reactive():\n", w)
    },
    finally = {
      base::print(base::paste0(sysTimePID(), " finished render plotly histP_Val()."))
      return(result)
    })
  })

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

  output$DTP_VAL <- DT::renderDataTable(as.data.frame(session$userData$sessionVariables$traitReducedDataStructure()$combinedDFP_Val_Labels$dfP_Val_w_number))

  output$DTDM <- DT::renderDataTable(as.data.frame(session$userData$sessionVariables$traitReducedDataStructure()$combinedDFP_Val_Labels$dfDM_w_number))

  output$DTN <- DT::renderDataTable(as.data.frame(session$userData$sessionVariables$traitReducedDataStructure()$combinedDFP_Val_Labels$dfN_w_number))

  #this data structure holds everything (as a named list), that is needed for working wih trait reduced (by selecting a subset of traits) HM
  session$userData$sessionVariables$combinedDataStructure <- shiny::reactive({

     base::tryCatch(
       {
         result <- base::list(combinedDFP_Val_Labels = session$userData$sessionVariables$combinedData()
         )
         #numberCores <- session$userData$numCores
      },
      error = function(e) {
        base::message("An error occurred in shiny::reactive(session$userData$sessionVariables$combinedDataStructure):\n", e)
      },
      warning = function(w) {
        base::message("A warning occurred in shiny::reactive(session$userData$sessionVariables$combinedDataStructure):\n", w)
      },
      finally = {
        base::print(base::paste0(sysTimePID(), " finished generating session$userData$sessionVariables$combinedDataStructure"))
        return(result)
      }
    )
  })

  session$userData$sessionVariables$pReducedDataStructure <- shiny::reactive({
    base::tryCatch(
      {
        if (is.valid(session$userData$sessionVariables$combinedDataStructure()$combinedDFP_Val_Labels$dfP_Val)) {
          result <- base::list(combinedDFP_Val_Labels = session$userData$sessionVariables$pReducedData(),
                           #                            #matP_Val = session$userData$sessionVariables$matP_Val(), #is part of combinedDFP_Val_Labels()
                           matP_Val.t = NULL,
                           distMatTraits = NULL,
                           clustResTraits = NULL, #session$userData$sessionVariables$clustResTraits(),
                           traitClusters = NULL, #session$userData$sessionVariables$traitClusters(),
                           traitClusterMedoids = NULL, #session$userData$sessionVariables$traitClusterMedoids(),
                           traitDendrogram = NULL, #session$userData$sessionVariables$traitDendrogram(),
                           traitClustergram = NULL, #session$userData$sessionVariables$traitClustergram(),
                           distMatProbes = NULL, #session$userData$sessionVariables$distMatProbes(),
                           clustResProbes = NULL, #session$userData$sessionVariables$clustResProbes(),
                           dendProbes = NULL #session$userData$sessionVariables$dendProbes()
          )
          result$matP_Val.t <- t(as.matrix(result$combinedDFP_Val_Labels$dfP_Val))
          numberCores <- session$userData$numCores
          base::print(paste0(sysTimePID(), " (pReducedDataStructure) before distance matrix for n(reduced traits) = ", base::nrow(result$matP_Val.t), " (takes some time). Using n(cores) = ", numberCores, "."))
          if (is.valid(result$matP_Val.t)) {
            result$distMatTraits <- getDistMat(numberCores = numberCores, matrix = result$matP_Val.t)
          }
          else {
            result$distMatTraits <- NULL
          }
          base::print(paste0(sysTimePID(), " (pReducedDataStructure) after distance matrix for reduced traits."))
          base::print(paste0(sysTimePID(), " (pReducedDataStructure) before clustering for traits.", nrow(result$matP_Val.t)))
          if (is.valid(result$distMatTraits)) {
            result$clustResTraits <- getClustResFast(result$distMatTraits)
          }
          else {
            result$clustResTraits <- NULL
          }
          base::print(paste0(sysTimePID(), " after clustering results for traits."))

          # if (is.valid(input$sldNumClusters)) {
          #   base::print(base::paste0(sysTimePID(), " (pReducedDataStructure) start generating traitClusters."))
          #   if (is.valid(result$clustResTraits) && input$sldNumClusters > 1) {
          #     result$traitClusters <- cutree(result$clustResTraits,
          #                                    k = input$sldNumClusters)
          #   }
          #   else {
          #     result$traitClusters <- NULL
          #   }
          #   base::print(base::paste0(sysTimePID(), " (pReducedDataStructure) start generating traitClusterMedoids."))
          #   if (is.valid(result$clustResTraits) && is.valid(result$distMatTraits)) {
          #     result$traitClusterMedoids <- getTraitClusterMedoids(clustResTraits = result$clustResTraits,
          #                                                          distMatTraits = result$distMatTraits,
          #                                                          numClusters = input$sldNumClusters)
          #   }
          #   else {
          #     result$traitClusterMedoids <- NULL
          #   }
          # }

          # base::print(base::paste0(sysTimePID(), " (pReducedDataStructure) start generating traitDendrogram."))
          # #do this only, if a dataset is already loaded
          #
          # if (is.valid(result$clustResTraits)) {
          #   shiny::reactive({
          #     result$traitDendrogram <- getDendTraits(clustResTraits = result$clustResTraits, traitClusters = input$sldNumClusters)
          #   })
          #   base::print(base::paste0(sysTimePID(), " (pReducedDataStructure) after making traitDendrogram."))
          # }
          # else {
          #   result$traitDendrogram <- NULL
          # }
          #
          # base::print(base::paste0(sysTimePID(), " (pReducedDataStructure) start generating traitClustergram."))
          # if (is.valid(result$matP_Val.t) && is.valid(result$clustResTraits)) {
          #   result$traitClustergram <- getplotClustergramTraitsLong(matP_Val.t = result$matP_Val.t,
          #                                                           clustResTraits = result$clustResTraits,
          #                                                           traitClusters = input$sldNumClusters)
          # }
          # else {
          #   result$traitClustergram <- NULL
          # }
        }
      },
      error = function(e) {
        base::message("An error occurred in shiny::reactive(session$userData$sessionVariables$pReducedDataStructure):\n", e)
      },
      warning = function(w) {
        base::message("A warning occurred in shiny::reactive(session$userData$sessionVariables$pReducedDataStructure):\n", w)
      },
      finally = {
        base::print(base::paste0(sysTimePID(), " finished generating session$userData$sessionVariables$pReducedDataStructure"))
        return(result)
      }
    )
  })

  session$userData$sessionVariables$traitReducedDataStructure <- shiny::reactive({
    base::tryCatch(
      {
        if (is.valid(session$userData$sessionVariables$pReducedDataStructure()$combinedDFP_Val_Labels$dfP_Val)) {
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
                               probeDendrogram = NULL
          )
          if (is.valid(session$userData$sessionVariables$traitReducedData())) {
            result$combinedDFP_Val_Labels <- session$userData$sessionVariables$traitReducedData()
            result$matP_Val.t = t(as.matrix(result$combinedDFP_Val_Labels$dfP_Val))
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
              result$traitClustergram <- getplotClustergramTraitsLong(matP_Val.t = result$matP_Val.t,
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
            result$combinedDFP_Val_Labels$dfDM_w_number$number <- seq(1:nprobes)
            col_order <- c("number", colnames(result$combinedDFP_Val_Labels$dfDM_w_number))
            result$combinedDFP_Val_Labels$dfDM_w_number <- result$combinedDFP_Val_Labels$dfDM_w_number[, col_order]
            #              result$dfDM_w_number <- result$dfDM[ , -which(colnames(result$dfDM_w_number) %in% "number.1")]
            result$combinedDFP_Val_Labels$dfN_w_number <- result$combinedDFP_Val_Labels$dfN
            result$combinedDFP_Val_Labels$dfN_w_number$number <- seq(1:nprobes)
            col_order <- c("number", colnames(result$combinedDFP_Val_Labels$dfN_w_number))
            result$combinedDFP_Val_Labels$dfN_w_number <- result$combinedDFP_Val_Labels$dfN_w_number[, col_order]
            result$combinedDFP_Val_Labels$dfN_w_number <- result$combinedDFP_Val_Labels$dfN_w_number[ , -which(colnames(result$combinedDFP_Val_Labels$dfN_w_number) %in% "number.1")]

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
            }
            base::print(base::paste0(sysTimePID(), " (traitReducedDataStructure) start clustResProbes."))
            distMat <- result$distMatProbes
            if (is.valid(distMat)) {
              result$clustResProbes <- getClustResFast(distMat)
            }
            else {
              result$clustResProbes <- NULL
            }
            if (is.valid(result$clustResProbes)) {
              result$probeDendrogram <- stats::as.dendrogram(result$clustResProbes)
            }
            else {
              result$probeDendrogram <- NULL
            }
          }
          ###tbc()
        }
      },
      error = function(e) {
browser()
        base::message("An error occurred in shiny::reactive(session$userData$sessionVariables$traitReducedDataStructure):\n", e)
      },
      warning = function(w) {
        base::message("A warning occurred in shiny::reactive(session$userData$sessionVariables$traitReducedDataStructure):\n", w)
      },
      finally = {
        base::print(base::paste0(sysTimePID(), " finished generating session$userData$sessionVariables$traitReducedDataStructure."))
        return(result)
      }
    )
  })

#   session$userData$sessionVariables$probeReducedDataStructure <- shiny::reactive({
#     base::tryCatch(
#       {
#         traitReducedDataStructure <- session$userData$sessionVariables$traitReducedDataStructure()
#         result <- base::list(combinedDFP_Val_Labels = traitReducedDataStructure$combinedDFP_Val_Labels,
#                              matP_Val.t = NULL, #traitReducedDataStructure$matP_Val.t,
#                              distMatTraits = NULL, #traitReducedDataStructure$distMatTraits,
#                              clustResTraits = NULL, #traitReducedDataStructure$clustResTraits,
#                              traitClusters = NULL, #traitReducedDataStructure$traitClusters,
#                              traitClusterMedoids = NULL, #traitReducedDataStructure$traitClusterMedoids,
#                              traitDendrogram = NULL, #traitReducedDataStructure$traitDendrogram,
#                              traitClustergram = NULL, #traitReducedDataStructure$traitClustergram
#                              distMatProbes = NULL,
#                              clustResProbes = NULL,
#                              probeDendrogram = NULL
#                              )
#         DNAdistances <- session$userData$sessionVariables$distNeigboursProbes10()
#         DNAdistances <- calculateDistanceNeigboursProbes(wd = session$userData$packageWd, clustResProbes = session$userData$sessionVariables$traitReducedDataStructure()$clustResProbes, annotation = session$userData$annotation, distanceToLook = 10, numCores = session$userData$numCores)
#         DNAdistances <- na.omit(DNAdistances)
#         dfP_Val <- result$combinedDFP_Val_Labels$dfP_Val
#         dfP_Val <- dfP_Val[which(rownames(dfP_Val) %in% DNAdistances$ID),]
#         result$combinedDFP_Val_Labels$dfP_Val <- dfP_Val
#         rm(dfP_Val)
#         dfDM <- result$combinedDFP_Val_Labels$dfDM
#         dfDM <- dfDM[which(rownames(dfDM) %in% DNAdistances$ID),]
#         result$combinedDFP_Val_Labels$dfDM <- dfDM
#         rm(dfDM)
#         dfN <- result$combinedDFP_Val_Labels$dfN
#         dfN <- dfN[which(rownames(dfN) %in% DNAdistances$ID),]
#         result$combinedDFP_Val_Labels$dfN <- dfN
#         rm(dfN)
# browser()
#         #do this only, if session$userData$sessionVariables$matP_Val() is.valid
#         if (is.valid(result$combinedDFP_Val_Labels$dfP_Val)) {
#           dfP_Val <- result$combinedDFP_Val_Labels$dfP_Val
#           dfP_Val[base::is.na(dfP_Val)] <- 1 # set missing P_VAL to 1
#           result$matP_Val.t <- t(dfP_Val)
#         }
#         else {
#           base::print(paste0(sysTimePID(), " is.valid(result$combinedDFP_Val_Labels$dfP_Val) == FALSE."))
# browser()
#           result$matP_Val.t <- NULL
#         }
#
#         base::print(base::paste0(sysTimePID(), " start generating session$userData$sessionVariables$probeReducedGeneralDataStructure$distMatTraits."))
#         numberCores <- session$userData$numCores
#         if (is.valid(result$combinedDFP_Val_Labels$dfP_Val)) {
#           base::print(paste0(sysTimePID(), " before distance matrix for n(traits) = ", base::nrow(result$matP_Val.t),  " (takes some time). Using n(cores) = ", numberCores, "."))
#           result$distMatTraits <- getDistMat(numberCores = numberCores, matrix = result$matP_Val.t)
#           base::print(paste0(sysTimePID(), " after distance matrix for traits."))
#         }
#         else {
#           base::print(base::paste0(sysTimePID(), " is.valid(result$combinedDFP_Val_Labels$dfP_Val) == FALSE"))
# browser()
#           result$distMatTraits <- NULL
#         }
#         base::print(base::paste0(sysTimePID(), " start generating session$userData$sessionVariables$probeReducedGeneralDataStructure$clustResTraits."))
#         if(is.valid(result$distMatTraits)) {
#           base::print(paste0(sysTimePID(), " before clustering for traits.", nrow(result$matP_Val.t)))
#           result$clustResTraits <- getClustResFast(result$distMatTraits)
# #          base::print(paste0(sysTimePID(), " after clustering results for traits.")) #, nrow(result$matP_Val.t)))
#         }
#         else {
#           base::print(base::paste0(sysTimePID(), " is.valid(result$distMatTraits) == FALSE"))
# browser()
#           result$clustResTraits <- NULL
#         }
#
#         base::print(base::paste0(sysTimePID(), " start generating session$userData$sessionVariables$probeReducedGeneralDataStructure$traitClusters."))
#         if (is.valid(result$clustResTraits) && input$sldNumClusters > 1) {
#           result$traitClusters <- cutree(result$clustResTraits,
#                            k = input$sldNumClusters)
#         }
#         else {
#           base::print(base::paste0(sysTimePID(), " is.valid(result$clustResTraits) == FALSE && input$sldNumClusters > 1"))
# browser()
#           result$traitClusters <- NULL
#         }
#         base::print(paste0(sysTimePID(), " after generating session$userData$sessionVariables$probeReducedGeneralDataStructure$traitClusters."))
#
#         base::print(base::paste0(sysTimePID(), " start generating session$userData$sessionVariables$probeReducedGeneralDataStructure$traitClusterMedoids."))
#         result$traitClusterMedoids <- getTraitClusterMedoids(clustResTraits = result$clustResTraits,
#                                          distMatTraits = result$distMatTraits,
#                                          numClusters = input$sldNumClusters)
#         base::print(paste0(sysTimePID(), " after generating session$userData$sessionVariables$probeReducedGeneralDataStructure$traitClusterMedoids."))
#
#         if (is.valid(result$clustResTraits)) {
#           result$traitDendrogram <- getDendTraits(clustResTraits = result$clustResTraits, traitClusters = input$sldNumClusters)
#           base::print(base::paste0(sysTimePID(), " after making traitDendrogram."))
#         }
#         else {
#           base::print(base::paste0(sysTimePID(), " is.valid(result$clustResTraits == FALSE."))
# browser()
#           result$traitDendrogram <- NULL
#         }
#
#         base::print(base::paste0(sysTimePID(), " start session$userData$sessionVariables$probeReducedGeneralDataStructure$traitClustergram."))
#         if (is.valid(result$matP_Val.t) && is.valid(result$clustResTraits)) {
#           result$traitClustergram <- getplotClustergramTraitsLong(matP_Val.t = result$matP_Val.t,
#                                                clustResTraits = result$clustResTraits,
#                                                traitClusters = input$sldNumClusters)
#         }
#         else {
#           base::print(base::paste0(sysTimePID(), " (is.valid(result$matP_Val.t) && is.valid(result$clustResTraits)) == FALSE."))
# browser()
#           result$traitClustergram <- NULL
#         }
#         base::print(base::paste0(sysTimePID(), " start generating session$userData$sessionVariables$probeReducedGeneralDataStructure$distMatProbes."))
#         dfP_Val <- result$combinedDFP_Val_Labels$dfP_Val #session$userData$sessionVariables$traitReducedmatP_Val()
#         if (is.valid(dfP_Val)) {
#           dfP_Val[dfP_Val > 0.05] <- NA # 1
#           base::print(base::paste0(sysTimePID(), " calculating distance matrix with rows= ", nrow(dfP_Val), " cols= ", ncol(dfP_Val)))
#           base::print(base::class(dfP_Val))
#           base::print(base::paste0(sysTimePID(), " set missing p-values to 1."))
#           dfP_Val[base::is.na(dfP_Val)] <- 1 # set missing P_VAL to 1
#           base::print(Cstack_info())
#           if (base::nrow(dfP_Val) >= 5) {
#             numberCores <- session$userData$numCores
#             # clustering for rows
#             base::print(base::paste0(sysTimePID(), " before distance matrix for n(probes) = ", base::nrow(dfP_Val), " (takes some time). Using n(cores) = ", numberCores, "."))
#             gc()
#             result$distMatProbes <- getDistMat(numberCores = numberCores, matrix = dfP_Val)
#             base::print(base::paste0(sysTimePID(), " after distance matrix for probes.", base::nrow(dfP_Val)))
#           } else {
#             base::message(base::paste0(sysTimePID(), " less than 5 probes remained."))
#             result$distMatProbes <- NULL
#           }
#           rm(dfP_Val)
#         }
#         else {
# browser()
#           result$distMatProbes <- NULL
#         }
#
#         base::print(base::paste0(sysTimePID(), " start generating session$userData$sessionVariables$probeReducedGeneralDataStructure$clustResProbes."))
#         if (is.valid(result$distMatProbes)) {
#           result$clustResProbes <- getClustResFast(result$distMatProbes)
# #          base::print(base::paste0(sysTimePID(), " finished generating session$userData$sessionVariables$probeReducedGeneralDataStructure$clustResProbes."))
#         }
#         else {
#           base::print(base::paste0(sysTimePID(), " is.valid(result$distMatProbes) == FALSE."))
# browser()
#           result$clustResProbes <- NULL
#         }
#
#         base::print(base::paste0(sysTimePID(), " start generating session$userData$sessionVariables$probeReducedGeneralDataStructure$dendProbes."))
#         if(is.valid(result$clustResProbes)){
#           result$probeDendrogram <- stats::as.dendrogram(result$clustResProbes)
#           length(unlist(result))
#           base::print(base::paste0(sysTimePID(), " after making dendrogram for probes"))
#         }
#         else {
#           base::print(base::paste0(sysTimePID(), " is.valid(result$clustResProbes) == FALSE."))
# browser()
#           result$probeDendrogram = NULL
#         }
#       },
#       error = function(e) {
#         base::message("An error occurred in shiny::reactive(session$userData$sessionVariables$probeReducedGeneralDataStructure):\n", e)
#       },
#       warning = function(w) {
#         base::message("A warning occurred in shiny::reactive(session$userData$sessionVariables$probeReducedGeneralDataStructure):\n", w)
#       },
#       finally = {
#         base::print(base::paste0(sysTimePID(), " finished generating session$userData$sessionVariables$probeReducedGeneralDataStructure. (last step in session$userData$sessionVariables$probeReducedGeneralDataStructure <- shiny::reactive())"))
#         return(result)
#       }
#     )
#   })

  DTProbes <- shiny::reactive({
    base::tryCatch(
      {
        base::print(base::paste0(sysTimePID(), " start generating DTProbes."))
        base::print(base::paste0(sysTimePID(), " before making probe table."))
        if (!is.null(session$userData$sessionVariables$traitReducedDataStructure()$clustResProbes)) {
          dendProbes <- session$userData$sessionVariables$traitReducedDataStructure()$probeDendrogram
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
          result <- DTProbes
          base::print(base::paste0(sysTimePID(), " after making probe table."))
        }
        else {
          result <- NULL
        }
      },
      error = function(e) {
        base::message("An error occurred in shiny::reactive(DTProbes):\n", e)
      },
      warning = function(w) {
        base::message("A warning occurred in shiny::reactive(DTProbes):\n", w)
      },
      finally = {
        base::print(base::paste0(sysTimePID(), " finished generating DTProbes."))
        return(result)
      }
    )
  })

  output$DTProbes <- DT::renderDataTable(as.data.frame(DTProbes()),
                                         options = list(pageLength = 1000, info = FALSE,
                                                                                   lengthMenu = list(c(100, 1000, -1), c("100", "1000", "All")) ))
  DTTraits <- shiny::reactive({
    base::tryCatch(
      {
        base::print(base::paste0(sysTimePID(), " start generating DTTraits."))
        base::print(base::paste0(sysTimePID(), " before making traits table."))
        if (!is.null(session$userData$sessionVariables$traitReducedDataStructure()$clustResTraits)) {
          listTraits <- ggdendro::dendro_data(session$userData$sessionVariables$traitReducedDataStructure()$clustResTraits, type = "rectangle")$labels$label

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
        base::message("An error occurred in shiny::reactive(DTTraits):\n", e)
      },
      warning = function(w) {
        base::message("A warning occurred in shiny::reactive(DTTraits):\n", w)
      },
      finally = {
        base::print(base::paste0(sysTimePID(), " finished generating DTTraits."))
        return(result)
      }
    )
  })

  output$DTTraits <- DT::renderDataTable(as.data.frame(DTTraits()))

  output$plotDendrogramProbes <- plotly::renderPlotly(ggdendro::ggdendrogram(session$userData$sessionVariables$traitReducedDataStructure()$clustResProbes,
                                                                             rotate = TRUE, theme_dendro = FALSE))

  output$plotDendrogramTraits <- plotly::renderPlotly(ggdendro::ggdendrogram(session$userData$sessionVariables$traitReducedDataStructure()$clustResTraits,
                                                                             rotate = TRUE, theme_dendro = FALSE))

  shiny::observeEvent(input$btnPlotCombinedHM_P_Val,
    ignoreInit = TRUE,
    {
      base::tryCatch(
        {
          base::print(base::paste0(sysTimePID(), " Step 5: start plotting heatmap for P_Val. (first step in shiny::observeEvent(input$btnPlotCombinedHM_P_Val))"))
#          session <- invalidateStep5(session)
          plotCombinedHM_P_Val(input = input, output = output, session = session)
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

  shiny::observeEvent(input$btnPlotCombinedDWHM_P_Val,
    ignoreInit = TRUE,
    {
      base::tryCatch(
        {
          base::print(base::paste0(sysTimePID(), " Step 6: start plotting heatmap for distance weighted DM. (first step in shiny::observeEvent(input$btnPlotCombinedDWHM_P_Val))"))
#          session <- invalidateStep5(session)
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

  shiny::observeEvent(input$btnPlotCombinedCondDWHM_P_Val,
                      ignoreInit = TRUE,
                      {
                        base::tryCatch(
                          {
                            base::print(base::paste0(sysTimePID(), " Step 6: start plotting condensed heatmap for distance weighted DM. (first step in shiny::observeEvent(input$btnPlotCombinedCondDWHM_P_Val))"))
#                            session <- invalidateStep5(session)
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

  shiny::observeEvent(input$btnSearchCpGHM,
    ignoreInit = TRUE,
    {
      base::tryCatch(
        {
          base::print(base::paste0(sysTimePID(), " start searching CpG."))
          #find positions
          searchResult <- getSearchResultCpG(input$txtSearchCpG, session)
          length <- length(session$userData$sessionVariables$traitReducedDataStructure()$clustResProbes$labels)
          resultText <- paste0(input$txtSearchCpG, " found at position: ", searchResult, " from ", length, " CpG.")
          #write to output
          output$txtSearchResultCpG <- shiny::renderText(resultText)
          #mark in HM
          #browser()
          #plotCombinedHM_P_Val()
        },
        error = function(e) {
          base::message("An error occurred in shiny::observeEvent(input$btnSearchCpGHM):\n", e)
        },
        warning = function(w) {
          base::message("A warning occurred in shiny::observeEvent(input$btnSearchCpGHM):\n", w)
        },
        finally = {
          base::print(base::paste0(sysTimePID(), " finished searching CpG."))
        }
      )

    },
    ignoreNULL = FALSE
  )

  shiny::observeEvent(input$btnSearchTraitHM,
    ignoreInit = TRUE,
    {
      base::tryCatch(
        {
          base::print(base::paste0(sysTimePID(), " start searching trait."))
          #find positions
          searchResult <- getSearchResultTrait(input$txtSearchTrait, session)
          length <- length(session$userData$sessionVariables$traitReducedDataStructure()$clustResTraits$labels)
          resultText <- paste0(input$txtSearchTrait, " found at position: ", searchResult, " from ", length, " Traits.")
          #write to output
          output$txtSearchResultTrait <- shiny::renderText(resultText)
          #mark in HM
          #browser()
          #plotCombinedHM_P_Val()
        },
        error = function(e) {
          base:.warning("An error occurred in shiny::observeEvent(input$btnSearchTraitHM):\n", e)
        },
        warning = function(w) {
          base::message("A warning occurred in shiny::observeEvent(input$btnSearchTraitHM):\n", w)
        },
        finally = {
          base::print(base::paste0(sysTimePID(), " end search trait heatmap."))
          base::print(base::paste0(sysTimePID(), " finished searching trait."))
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
