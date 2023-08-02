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

  packageWd <<- getwd()
  session$userData$packageWd <- getwd()
  base::print(paste0(Sys.time(), " getwd: ", packageWd))
  base::print(paste0(Sys.time(), " loading configuration."))
  configFileLocation <- system.file("config.yml", package = "PatternMatchR", mustWork = TRUE)
  session$userData$config <- config::get(file = configFileLocation)
  volumes <- c(Home = paste0(getwd(), sub(".", "", session$userData$config$workDir)),
               shinyFiles::getVolumes()())
  shinyFiles::shinyFileChoose(input, "file", roots = volumes, session = session)
  shinyFiles::shinyFileSave(input, "save", roots = volumes, session = session, restrictions = system.file(package = "base"))
  #base::options(spam.force64 = TRUE)

  if (session$userData$config$debugMode == TRUE) {
    shiny::updateCheckboxInput(session, "chkDebug", value = TRUE)
  }
  else {
    shiny::updateCheckboxInput(session, "chkDebug", value = FALSE)
  }

  print(paste0(Sys.time(), " defining session variables."))
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
  result <- loadDirLists(session = session, input = input, output = output)
  dfdD1 <- result$dfdD1
  dfdD2 <- result$dfdD2
  dfdD3 <- result$dfdD3

  session$userData$sessionVariables$reactiveTestVal <- shiny::reactiveVal(value = NULL, label = "reactiveTestVal")

  session$userData$sessionVariables$resultDFListTrait1 <- shiny::reactiveVal(value = NULL, label = "resultDFListTrait1")
  session$userData$sessionVariables$resultDFListTrait2 <- shiny::reactiveVal(value = NULL, label = "resultDFListTrait2")
  session$userData$sessionVariables$resultDFListTrait3 <- shiny::reactiveVal(value = NULL, label = "resultDFListTrait3")
  session$userData$sessionVariables$combinedDFP_Val_Labels <- shiny::reactiveVal(value = NULL, label = "combinedDFP_Val_Labels")
  session$userData$sessionVariables$pReducedcombinedDFP_Val_Labels <- shiny::reactiveVal(value = NULL, label = "pReducedcombinedDFP_Val_Labels")

  session$userData$sessionVariables$matP_Val <- shiny::reactiveVal(value = NULL, label = "matP_Val")
  session$userData$sessionVariables$matP_Val.t <- shiny::reactiveVal(value = NULL, label = "matP_Val.t")
  session$userData$sessionVariables$distMatTraits <- shiny::reactiveVal(value = NULL, label = "distMatTraits")
  session$userData$sessionVariables$clustResTraits <- shiny::reactiveVal(value = NULL, label = "clustResTraits")
  session$userData$sessionVariables$traitClusters <- shiny::reactiveVal(value = NULL, label = "traitClusters")
  session$userData$sessionVariables$traitClusterMedoids <- shiny::reactiveVal(value = NULL, label = "traitClusterMedoids")
  session$userData$sessionVariables$traitDendrogram <- shiny::reactiveVal(value = NULL, label = "traitDendrogram")
  session$userData$sessionVariables$traitClustergram <- shiny::reactiveVal(value = NULL, label = "traitClustergram")
  session$userData$sessionVariables$traitReducedcombinedDFP_Val_Labels <- shiny::reactiveVal(value = NULL, label = "traitReducedcombinedDFP_Val_Labels")
  session$userData$sessionVariables$traitReducedmatP_Val <- shiny::reactiveVal(value = NULL, label = "traitReducedmatP_Val")
  session$userData$sessionVariables$traitReducedmatP_Val.t <- shiny::reactiveVal(value = NULL, label = "traitReducedmatP_Val.t")
  session$userData$sessionVariables$traitReduceddistMatTraits <- shiny::reactiveVal(value = NULL, label = "traitReduceddistMatTraits")
  session$userData$sessionVariables$traitReducedclustResTraits <- shiny::reactiveVal(value = NULL, label = "traitReducedclustResTraits")
  session$userData$sessionVariables$traitReducedDendTraits <- shiny::reactiveVal(value = NULL, label = "traitReducedDendTraits")

  session$userData$sessionVariables$distMatProbes <- shiny::reactiveVal(value = NULL, label = "distMatProbes")
  session$userData$sessionVariables$clustResProbes <- shiny::reactiveVal(value = NULL, label = "clustResProbes")
  session$userData$sessionVariables$dendProbes <- shiny::reactiveVal(value = NULL, label = "dendProbes")

  session$userData$sessionVariables$SPLOM <- FALSE

  shiny::updateSliderInput(session = session, inputId = "sldP_Val", min = 0, max = 0, value = c(0, 0))
  shiny::updateSliderInput(session = session, inputId = "sldDM", min = 0, max = 0, value = c(0, 0))
  shiny::updateSliderInput(session = session, inputId = "sldN", min = 0, max = 0, value = c(0, 0))

  base::print(paste0(Sys.time(), " starting application."))

  shinyjs::toggleClass("colRed", "red")
  shinyjs::toggleClass("colGreen", "green")
  shinyjs::toggleClass("colBlue", "blue")

################################################################################

  shiny::observe({
    shiny::invalidateLater(10000, session)
    base::print(paste0(
      Sys.time(),
      " PatternMatchR is running in idle state."
    ))
  })

  session$userData$sessionVariables$matP_Val <- shiny::reactive({
    base::print(paste0(
      Sys.time(),
      " before receiving matrix for traits."
    ))
    combinedDFP_Val_Labels <- session$userData$sessionVariables$pReducedcombinedDFP_Val_Labels()
    base::print(paste0(
      Sys.time(),
      " after receiving matrix for traits."
    ))
    return(combinedDFP_Val_Labels$dfP_Val)
  })

  session$userData$sessionVariables$matP_Val.t <- shiny::reactive({
    base::tryCatch(
      {
        base::print(paste0(
          Sys.time(),
          " before transposing matrix for traits."
        ))
        #do this only, if session$userData$sessionVariables$matP_Val() is.valid
        if (is.valid(session$userData$sessionVariables$matP_Val())) {
          dfP_Val <- session$userData$sessionVariables$matP_Val()
          dfP_Val[base::is.na(dfP_Val)] <- 1 # set missing P_VAL to 1
          result <- t(dfP_Val)
        }
        else {
          base::print(paste0(
            Sys.time(),
            " is.valid(session$userData$sessionVariables$matP_Val()) == FALSE."
          ))
          result <- FALSE
        }
      },
      error = function(e) {
        message("An error occurred in session$userData$sessionVariables$matP_Val.t():\n", e)
      },
      warning = function(w) {
        message("A warning occurred in session$userData$sessionVariables$matP_Val.t():\n", w)
      },
      finally = {
        base::print(paste0(
          Sys.time(),
          " after transposing matrix for traits."
        ))
        return(result)
      }
    )
  })

  session$userData$sessionVariables$traitReducedmatP_Val <- shiny::reactive({
    #take selected medoids as new traits for HM
    base::print(paste0(
      Sys.time(),
      " before receiving matrix for reduced traits."
    ))
    combinedDFP_Val_Labels <- session$userData$sessionVariables$traitReducedcombinedDFP_Val_Labels()
    base::print(paste0(
      Sys.time(),
      " after receiving matrix for reduced traits."
    ))
    return(combinedDFP_Val_Labels$dfP_Val)
  })


  session$userData$sessionVariables$traitReducedmatP_Val.t <- shiny::reactive({
    base::print(paste0(
      Sys.time(),
      " before transposing matrix for reduced traits."
    ))
    dfP_Val <- session$userData$sessionVariables$traitReducedmatP_Val()
    dfP_Val[base::is.na(dfP_Val)] <- 1 # set missing P_VAL to 1
    base::print(paste0(
      Sys.time(),
      " after transposing matrix for reduced traits."
    ))
    return(t(dfP_Val))
  })

  session$userData$sessionVariables$traitReduceddistMatTraits <- shiny::reactive({
    numberCores <- session$userData$numCores
    base::print(paste0(
      Sys.time(),
      " before distance matrix for reduced traits."
    ))
    result <- getDistMat(numberCores = numberCores, matrix = session$userData$sessionVariables$traitReducedmatP_Val.t())
    base::print(paste0(
      Sys.time(),
      " after distance matrix for reduced traits."
    ))
    return(result)
  })

  session$userData$sessionVariables$distMatTraits <- shiny::reactive({
    numberCores <- session$userData$numCores
    base::print(paste0(
      Sys.time(),
      " before distance matrix for traits."
    ))
    result <- getDistMat(numberCores = numberCores, matrix = session$userData$sessionVariables$matP_Val.t())
    base::print(paste0(
      Sys.time(),
      " after distance matrix for traits."
    ))
    return(result)
  })

   output$DTDebugOut1 <-
     DT::renderDataTable(as.matrix(session$userData$sessionVariables$distMatTraits()))

   output$plotDebug1 <-
     shiny::renderPlot(session$userData$sessionVariables$clustResTraits())

  session$userData$sessionVariables$clustResTraits <- shiny::reactive({
    base::print(paste0(
      Sys.time(),
      " before clustering for traits.",
      nrow(session$userData$sessionVariables$matP_Val.t())
    ))
    result <- getClustResFast(session$userData$sessionVariables$distMatTraits())
    base::print(paste0(
      Sys.time(),
      " after clustering results for traits.",
      nrow(session$userData$sessionVariables$matP_Val.t())
    ))
    return(result)
  })

  session$userData$sessionVariables$traitReducedclustResTraits <- shiny::reactive({
    base::print(paste0(
      Sys.time(),
      " before clustering for reduced traits.",
      nrow(session$userData$sessionVariables$traitReducedmatP_Val.t())
    ))
    result <- getClustResFast(session$userData$sessionVariables$traitReduceddistMatTraits())
    base::print(paste0(
      Sys.time(),
      " after clustering results for reduced traits.",
      nrow(session$userData$sessionVariables$traitReducedmatP_Val.t())
    ))
    return(result)
  })

  session$userData$sessionVariables$distMatProbes <- shiny::reactive({
    dfP_Val <- session$userData$sessionVariables$traitReducedmatP_Val()
    if (is.valid(dfP_Val)) {
    dfP_Val[dfP_Val > 0.05] <- NA # 1
    base::print(
      base::paste0(
        Sys.time(),
        " calculating distance matrix with rows= ",
        nrow(dfP_Val),
        " cols= ",
        ncol(dfP_Val)
      )
    )
    base::print(base::class(dfP_Val))
    base::print(base::class(dfP_Val))
    base::print(base::paste0(Sys.time(), " set missing p-values to 1."))
    dfP_Val[base::is.na(dfP_Val)] <- 1 # set missing P_VAL to 1
    base::print(base::class(dfP_Val))
    base::print(Cstack_info())
    if (base::nrow(dfP_Val) >= 5) {
      # clustering for rows
      base::print(
        base::paste0(
          Sys.time(),
          " before distance matrix for probes.",
          base::nrow(dfP_Val),
          " (takes some time)"
        )
      )
      #base::options(spam.force64 = TRUE)
      gc()
      numberCores <- session$userData$numCores
      result <- getDistMat(numberCores = numberCores, matrix = dfP_Val)
      base::print(
        base::paste0(
          Sys.time(),
          " after distance matrix for probes.",
          base::nrow(dfP_Val)
        )
      )
    } else {
      base::message(base::paste0(Sys.time(), " less than 5 probes remained."))
      result <- NULL
    }
    }
    else {
      result <- NULL
    }
    return(result)
  })

  session$userData$sessionVariables$clustResProbes <- shiny::reactive({
    base::print(
      base::paste0(
        Sys.time(),
        " before clustering for probes."
      )
    )
    distMat <- session$userData$sessionVariables$distMatProbes()
    if (is.valid(distMat)) {
      result <- getClustResFast(distMat)
    }
    else {
      base::print(
        base::paste0(
          Sys.time(),
          " is.valid(distMat) == FALSE."
        )
      )
      result <- NULL
    }
    base::print(
      base::paste0(
        Sys.time(),
        " after clustering for probes."
      )
    )
    return(result)
  })

  session$userData$sessionVariables$traitDendrogram <- shiny::reactive({
    base::print(
      base::paste0(
        Sys.time(),
        " before making traitDendrogram."
      )
    )
    #do this only, if a dataset is already loaded
    if (is.valid(session$userData$sessionVariables$clustResTraits())) {
      result <- getDendTraits(clustResTraits = session$userData$sessionVariables$clustResTraits(), traitClusters = input$sldNumClusters)
      base::print(
        base::paste0(
          Sys.time(),
          " after making traitDendrogram."
        )
      )
    }
    else {
      base::print(
        base::paste0(
          Sys.time(),
          " is.valid(session$userData$sessionVariables$clustResTraits() == FALSE."
        )
      )
      result <- FALSE
    }
    return(result)
  })

  session$userData$sessionVariables$traitClustergram <- shiny::reactive({
    base::print(
      base::paste0(
        Sys.time(),
        " before making traitClustergram."
      )
    )
    result <- getplotClustergramTraitsLong(matP_Val.t = session$userData$sessionVariables$matP_Val.t(),
                                           clustResTraits = session$userData$sessionVariables$clustResTraits(),
                                           traitClusters = input$sldNumClusters)
    base::print(
      base::paste0(
        Sys.time(),
        " after making traitClustergram."
      )
    )
  return(result)
  })

  session$userData$sessionVariables$traitClusterMedoids <- shiny::reactive({
    base::print(
      base::paste0(
        Sys.time(),
        " before making clusterMedoids."
      )
    )
    result <- getTraitClusterMedoids(clustResTraits = session$userData$sessionVariables$clustResTraits(),
                                     distMatTraits = session$userData$sessionVariables$distMatTraits(),
                                     numClusters = input$sldNumClusters)
    base::print(
      base::paste0(
        Sys.time(),
        " after making clusterMedoids."
      )
    )
    return(result)
  })

  session$userData$sessionVariables$traitClusters <- shiny::reactive({
    base::print(
      base::paste0(
        Sys.time(),
        " before making traitClusters."
      )
    )
    if (is.valid(session$userData$sessionVariables$clustResTraits()) && input$sldNumClusters > 1) {
      result <- cutree(session$userData$sessionVariables$clustResTraits(),
                       k = input$sldNumClusters)
    }
    else {
      result <- NULL
    }
    base::print(
      base::paste0(
        Sys.time(),
        " after making traitClusters."
      )
    )
    return(result)
  })

  output$txtLoadOut <- shiny::reactive({
    return(updateTxtLoadOut(session, session$userData$sessionVariables$resultDFListTrait1(),
                            session$userData$sessionVariables$resultDFListTrait2(),
                            session$userData$sessionVariables$resultDFListTrait3()))
  })

  output$txtMergeOut <- shiny::reactive({
    if (is.valid(session$userData$sessionVariables$combinedDFP_Val_Labels())) {
      result <- updateTxtMergeOut(session$userData$sessionVariables$combinedDFP_Val_Labels())
    }
    else {
      result <- NULL
    }
    return(result)
  })

  output$txtPReduceOut <- shiny::reactive({
    if (is.valid(session$userData$sessionVariables$pReducedcombinedDFP_Val_Labels())) {
      maxTraits <- ncol(session$userData$sessionVariables$pReducedcombinedDFP_Val_Labels()$dfP_Val)
      if (is.valid(session$userData$sessionVariables$traitReducedmatP_Val())) {
        value <- ncol(session$userData$sessionVariables$traitReducedmatP_Val())
      }
      else {
        value <- maxTraits
      }
      shiny::updateSliderInput(session = session, inputId = "sldNumClusters", max = maxTraits, min = 1, value = value, step = 1)
      result <- updateTxtpReduceOut(session$userData$sessionVariables$pReducedcombinedDFP_Val_Labels())
    }
    else {
      result <- NULL
    }
    return(result)
  })

  output$txtClusterOut <- shiny::reactive({
    minP_Val <- 5 * 10^base::as.integer(input$sldP_Val[1]) #minP_Val <- 5 * 10^-base::as.integer(input$sldP_Val[2])
    maxP_Val <- 5 * 10^base::as.integer(input$sldP_Val[2]) #maxP_Val <- 5 * 10^-base::as.integer(input$sldP_Val[1])
    minN <- base::as.integer(input$txtCases)

    sldNumClasses <- input$sldNumClusters
    return(updateTxtClusterOut(session$userData$sessionVariables$traitReducedcombinedDFP_Val_Labels(),
                               minP_Val = minP_Val, maxP_Val = maxP_Val,
                               minN = minN, sldNumClasses = sldNumClasses))
  })

  output$plotDendrogramTraitsLong <- shiny::renderPlot(getPlot(session$userData$sessionVariables$traitDendrogram()))

  output$plotClustergramTraitsLong <- shiny::renderPlot(getPlot(session$userData$sessionVariables$traitClustergram()))

  output$DTTraitsMedoids <- DT::renderDataTable(as.data.frame(getMedoidsTable(session$userData$sessionVariables$traitClusterMedoids())))

  output$DTTraitsClusters <- DT::renderDataTable(as.data.frame(getClustersTable(session$userData$sessionVariables$traitClusters(),
                                                                                session$userData$sessionVariables$traitClusterMedoids())))

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
        base::print(base::paste0(Sys.time(), " session data has been saved to ", filePath))
      }
    },
    ignoreNULL = FALSE
  )

  shiny::observeEvent(input$markingVar,
    ignoreInit = TRUE,
    {
      session$userData$sessionVariables$markingVar <- input$markingVar
      # redraw SPLOM
      height <- input$numSPLOMVSize
      width <- input$numSPLOMHSize
      fig <- getSPLOM(session$userData$sessionVariables$selectedOriginalData, session$userData$sessionVariables$markingVar, height = height, width = width)
      output$SPLOM <- plotly::renderPlotly(fig)
      session$userData$sessionVariables$SPLOM <- TRUE
    },
    ignoreNULL = FALSE
  )

  shiny::observeEvent(input$numSPLOMHSize,
                      ignoreInit = TRUE,
                      {
                        session$userData$sessionVariables$markingVar <- input$markingVar
                        # redraw SPLOM
                        height <- input$numSPLOMVSize
                        width <- input$numSPLOMHSize
                        fig <- getSPLOM(session$userData$sessionVariables$selectedOriginalData, session$userData$sessionVariables$markingVar, height = height, width = width)
                        output$SPLOM <- plotly::renderPlotly(fig)
                        session$userData$sessionVariables$SPLOM <- TRUE
                      },
                      ignoreNULL = FALSE
  )
  shiny::observeEvent(input$numSPLOMVSize,
                      ignoreInit = TRUE,
                      {
                        session$userData$sessionVariables$markingVar <- input$markingVar
                        # redraw SPLOM
                        height <- input$numSPLOMVSize
                        width <- input$numSPLOMHSize
                        fig <- getSPLOM(session$userData$sessionVariables$selectedOriginalData, session$userData$sessionVariables$markingVar, height = height, width = width)
                        output$SPLOM <- plotly::renderPlotly(fig)
                        session$userData$sessionVariables$SPLOM <- TRUE
                      },
                      ignoreNULL = FALSE
  )
  shiny::observeEvent(input$numHMHSize,
                      ignoreInit = TRUE,
                      {
                        session$userData$sessionVariables$markingVar <- input$markingVar
                        # redraw HM
                        #browser()
                        height <- input$numHMVSize
                        width <- input$numHMHSize
                        plotCombinedHM(input = input, output = output, session = session)
                      },
                      ignoreNULL = FALSE
  )
  shiny::observeEvent(input$numHMVSize,
                      ignoreInit = TRUE,
                      {
                        session$userData$sessionVariables$markingVar <- input$markingVar
                        # redraw HM
                        #browser()
                        height <- input$numHMVSize
                        width <- input$numHMHSize
                        plotCombinedHM(input = input, output = output, session = session)
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
          base::print(base::paste0(Sys.time(), " start loading trait 1 folders."))
          session$userData$sessionVariables$resultDFListTrait1(NULL)
          session$userData$sessionVariables$combinedDFP_Val_Labels(NULL)
          session$userData$sessionVariables$pReducedcombinedDFP_Val_Labels(NULL)
          base::print(base::paste0(Sys.time(), " creating empty heatmap during load process."))
          combinedHMP_VAL <- emptyHM()
          InteractiveComplexHeatmap::makeInteractiveComplexHeatmap(
            input = input,
            output = output,
            session = session,
            ht_list = combinedHMP_VAL,
            heatmap_id = "heatmap_1",
            show_layer_fun = FALSE,
            click_action = NULL,
            brush_action = NULL,
            hover_action = NULL
            )
          base::print(base::paste0(Sys.time(), " before is.numeric()."))
          if (base::is.numeric(input$trait1DirList_rows_selected)) {
            traitDirList <-
              base::as.list(dfdD1[input$trait1DirList_rows_selected, ]) #base::as.list(dfdD1[input$trait1DirList_rows_selected, ][[1]])
            base::print(base::paste0(Sys.time(), " selected folders: ", as.character(traitDirList)))
            session$userData$sessionVariables$resultDFListTrait1(loadDir(session = session, traitDirList = traitDirList))
          }
          else {
            session$userData$sessionVariables$resultDFListTrait1(NULL)
            base::message(base::paste0(Sys.time(), " no entries selected from trait1 folders."))
          }
        },
        error = function(e) {
          message("An error occurred in shiny::observeEvent(input$btnLoadDir1):\n", e)
        },
        warning = function(w) {
          message("A warning occurred in shiny::observeEvent(input$btnLoadDir1):\n", w)
        },
        finally = {
          base::print(base::paste0(Sys.time(), " end loadDir1."))
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
          base::print(base::paste0(Sys.time(), " start loading trait 2 folders."))
          session$userData$sessionVariables$resultDFListTrait2(NULL)
          session$userData$sessionVariables$combinedDFP_Val_Labels(NULL)
          session$userData$sessionVariables$pReducedcombinedDFP_Val_Labels(NULL)
          base::print(base::paste0(Sys.time(), " creating empty heatmap during load process."))
          combinedHMP_VAL <- emptyHM()
          InteractiveComplexHeatmap::makeInteractiveComplexHeatmap(
            input = input,
            output = output,
            session = session,
            ht_list = combinedHMP_VAL,
            heatmap_id = "heatmap_1",
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
            base::message(base::paste0(Sys.time(), " no entries selected from trait2 folders."))
          }
        },
        error = function(e) {
          message("An error occurred in shiny::observeEvent(input$btnLoadDir2):\n", e)
        },
        warning = function(w) {
          message("A warning occurred in shiny::observeEvent(input$btnLoadDir2):\n", w)
        },
        finally = {
          base::print(base::paste0(Sys.time(), " end loadDir2."))
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
          base::print(base::paste0(Sys.time(), " start loading trait3 folders."))
          session$userData$sessionVariables$resultDFListTrait3(NULL)
          session$userData$sessionVariables$combinedDFP_Val_Labels(NULL)
          session$userData$sessionVariables$pReducedcombinedDFP_Val_Labels(NULL)
          base::print(base::paste0(Sys.time(), " creating empty heatmap during load process."))
          combinedHMP_VAL <- emptyHM()
          InteractiveComplexHeatmap::makeInteractiveComplexHeatmap(
            input = input,
            output = output,
            session = session,
            ht_list = combinedHMP_VAL,
            heatmap_id = "heatmap_1",
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
            base::message(base::paste0(Sys.time(), " no entries selected from trait3 folders."))
          }
        },
        error = function(e) {
          message("An error occurred in shiny::observeEvent(input$btnLoadDir3):\n", e)
        },
        warning = function(w) {
          message("A warning occurred in shiny::observeEvent(input$btnLoadDir3):\n", w)
        },
        finally = {
          base::print(base::paste0(Sys.time(), " end loadDir3."))
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
          base::print(base::paste0(Sys.time(), " start loading all folders."))
          session$userData$sessionVariables$resultDFListTrait1(NULL)
          session$userData$sessionVariables$resultDFListTrait2(NULL)
          session$userData$sessionVariables$resultDFListTrait3(NULL)
          session$userData$sessionVariables$combinedDFP_Val_Labels(NULL)
          session$userData$sessionVariables$pReducedcombinedDFP_Val_Labels(NULL)
          base::print(base::paste0(Sys.time(), " creating empty heatmap."))
          combinedHMP_VAL <- emptyHM()
          InteractiveComplexHeatmap::makeInteractiveComplexHeatmap(
            input = input,
            output = output,
            session = session,
            ht_list = combinedHMP_VAL,
            heatmap_id = "heatmap_1",
            show_layer_fun = FALSE,
            click_action = NULL,
            brush_action = NULL,
            hover_action = NULL
            )
          if (base::is.numeric(input$trait1DirList_rows_selected)) {
            traitDirList <-
              base::as.list(dfdD1[input$trait1DirList_rows_selected, ])
            base::print(base::paste0(Sys.time(), " traitDirList1: ", as.character(traitDirList)))
            session$userData$sessionVariables$resultDFListTrait1(loadDir(session = session, traitDirList = traitDirList))
          } else {
            session$userData$sessionVariables$resultDFListTrait1(NULL)
            base::message(base::paste0(Sys.time(), " no entries selected from trait1 folders."))
          }
          if (base::is.numeric(input$trait2DirList_rows_selected)) {
            traitDirList <-
              base::as.list(dfdD2[input$trait2DirList_rows_selected, ])
            base::print(base::paste0(Sys.time(), " traitDirList2: ", as.character(traitDirList)))
            session$userData$sessionVariables$resultDFListTrait2(loadDir(session = session, traitDirList = traitDirList))
          } else {
            session$userData$sessionVariables$resultDFListTrait2(NULL)
            base::message(base::paste0(Sys.time(), " no entries selected from trait2 folders."))
          }
          if (base::is.numeric(input$trait3DirList_rows_selected)) {
            traitDirList <-
              base::as.list(dfdD3[input$trait3DirList_rows_selected, ])
            base::print(base::paste0(Sys.time(), " traitDirList3: ", as.character(traitDirList)))
            session$userData$sessionVariables$resultDFListTrait3(loadDir(session = session, traitDirList = traitDirList))
          } else {
            session$userData$sessionVariables$resultDFListTrait3(NULL)
            base::message(base::paste0(Sys.time(), " no entries selected from trait3 folders."))
          }
        },
        error = function(e) {
          message("An error occurred in shiny::observeEvent(input$btnLoadDirAll):\n", e)
        },
        warning = function(w) {
          message("A warning occurred in shiny::observeEvent(input$btnLoadDirAll):\n", w)
        },
        finally = {
          base::print(base::paste0(Sys.time(), " end loadDir all."))
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
          output$plotDebug1 <-
            shiny::renderPlot(session$userData$sessionVariables$clustResTraits())
        },
        error = function(e) {
          message("An error occurred in shiny::observeEvent(input$btnDebug):\n", e)
        },
        warning = function(w) {
          message("A warning occurred shiny::observeEvent(input$btnDebug):\n", w)
        },
        finally = {
          base::print(base::paste0(Sys.time(), " end debug test."))
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
          base::print(base::paste0(Sys.time(), " start merging data."))
          session$userData$sessionVariables$combinedDFP_Val_Labels(NULL)
          session$userData$sessionVariables$pReducedcombinedDFP_Val_Labels(NULL)
          base::print(base::paste0(Sys.time(), " creating empty heatmap."))
          # combinedHMP_VAL <- emptyHM()
          # InteractiveComplexHeatmap::makeInteractiveComplexHeatmap(input, output, session, combinedHMP_VAL, "heatmap_1")
          minN <- base::as.integer(input$txtCases)
          #if (session$userData$sessionVariables$LoadInitialized == TRUE) {
          if (is.valid(session$userData$sessionVariables$resultDFListTrait1()) || is.valid(session$userData$sessionVariables$resultDFListTrait2())  || is.valid(session$userData$sessionVariables$resultDFListTrait3())) {
            combinedDFP_Val_Labels <- mergeDFP_Val_Labels(session$userData$sessionVariables$resultDFListTrait1(),
                                                          session$userData$sessionVariables$resultDFListTrait2(),
                                                          session$userData$sessionVariables$resultDFListTrait3(),
                                                          minN)
          }
          else {
            base::print(base::paste0(Sys.time(), " is.valid(session$userData$sessionVariables$resultDFListTrait1()) || is.valid(session$userData$sessionVariables$resultDFListTrait2())  || is.valid(session$userData$sessionVariables$resultDFListTrait3()) == FALSE."))
            result <- NULL
          }
          session$userData$sessionVariables$combinedDFP_Val_Labels(combinedDFP_Val_Labels)
          updateSliders(session, session$userData$sessionVariables$combinedDFP_Val_Labels())
        },
        error = function(e) {
          message("An error occurred in shiny::observeEvent(input$btnMerge):\n", e)
        },
        warning = function(w) {
          message("A warning occurred in shiny::observeEvent(input$btnMerge):\n", w)
        },
        finally = {
          base::print(base::paste0(Sys.time(), " end merge data."))
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
          base::print(base::paste0(Sys.time(), " start reducing data by p-value."))
          session$userData$sessionVariables$pReducedcombinedDFP_Val_Labels(NULL)
          base::print(base::paste0(Sys.time(), " creating empty heatmap."))
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
          combinedDFP_Val_Labels <- session$userData$sessionVariables$combinedDFP_Val_Labels()
          if (is.valid(combinedDFP_Val_Labels)) {
            if (minP_Val != maxP_Val && minDM != maxDM && minN != maxN) {
              session$userData$sessionVariables$pReducedcombinedDFP_Val_Labels(getPReducedTraitData(combinedDFP_Val_Labels =
                                                                                                    combinedDFP_Val_Labels,
                                                                                                  minP_Val = minP_Val,
                                                                                                  maxP_Val = maxP_Val,
                                                                                                  minDM = minDM,
                                                                                                  maxDM = maxDM,
                                                                                                  minN = minN,
                                                                                                  maxN = maxN,
                                                                                                  debugMode = session$userData$config$debugMode))
            }
            else {
              base::print(base::paste0(Sys.time(), " minP_Val == maxP_Val && minDM == maxDM && minN == maxN."))
            }
          }
          else {
            base::print(base::paste0(Sys.time(), " is.valid(combinedDFP_Val_Labels) == FALSE."))
          }
        },
        error = function(e) {
          message("An error occurred in shiny::observeEvent(input$btnReduce):\n", e)
        },
        warning = function(w) {
          message("A warning occurred in shiny::observeEvent(input$btnReduce):\n", w)
        },
        finally = {
          base::print(base::paste0(Sys.time(), " end reduce data."))
        }
      )
    },
    ignoreNULL = FALSE
  )

  shiny::observeEvent(input$btnCluster,
    ignoreInit = TRUE,
    {
      base::tryCatch(
        {
          base::print(
            base::paste0(
              Sys.time(),
              " before making traitReducedcombinedDFP_Val_Labels."
            )
          )
          if (is.valid(session$userData$sessionVariables$pReducedcombinedDFP_Val_Labels()) &&
              is.valid(session$userData$sessionVariables$traitClusterMedoids())) {
            keys <- session$userData$config$keyAttributes
            if (is.valid(session$userData$sessionVariables$pReducedcombinedDFP_Val_Labels()) && is.valid(session$userData$sessionVariables$traitClusterMedoids())) {
              result <- getTraitReducedcombinedDFP_Val_Labels(session$userData$sessionVariables$pReducedcombinedDFP_Val_Labels(),
                                                              session$userData$sessionVariables$traitClusterMedoids(), keys)
            }
            else {
              base::print(base::paste0(Sys.time(), " (is.valid(session$userData$sessionVariables$pReducedcombinedDFP_Val_Labels()) && is.valid(session$userData$sessionVariables$traitClusterMedoids())) == FALSE."))
            }
          }
          else {
            result <- NULL
          }
        },
        error = function(e) {
          message("An error occurred in shiny::observeEvent(input$btnCluster):\n", e)
        },
        warning = function(w) {
          message("A warning occurred in shiny::observeEvent(input$btnCluster):\n", w)
        },
        finally = {
          base::print(base::paste0(Sys.time(), " after making traitReducedcombinedDFP_Val_Labels."))
          session$userData$sessionVariables$traitReducedcombinedDFP_Val_Labels(result)
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
          base::print(base::paste0(Sys.time(), " start counting probes p-value."))
          minN <- base::as.integer(input$txtCases)
          if (is.valid(session$userData$sessionVariables$combinedDFP_Val_Labels())) {
            P_VALNTable <-
              getAvailNForP_VALBorderParallel(session = session, wd = session$userData$packageWd, numCores = session$userData$numCores, DF = session$userData$sessionVariables$combinedDFP_Val_Labels()$dfP_Val)
            output$DTP_VALborder <- DT::renderDataTable(P_VALNTable)
            #insert scatterplot with table results
            #output$plotDendrogramP_VALborder <- plotly::renderPlotly(plotly::plot_ly(x = P_Val, type = "histogram", name = "histP_ValBorder"))
            plot <- plotly::plot_ly(x = P_VALNTable$P_VAL_BORDER, y = P_VALNTable$'Available CpG' , type = "scatter", mode = "lines+markers", name = "scatterP_ValBorder") %>%
              plotly::layout(xaxis = list(title = 'p-val', type = "log")) %>%
              plotly::layout(yaxis = list(title = 'n' ))
            output$plotDendrogramP_VALborder <- plotly::renderPlotly(plot)
            base::print(base::paste0(Sys.time(), " finished counting probes p-value."))
          }
          else {
            base::print(base::paste0(Sys.time(), " is.valid(session$userData$sessionVariables$combinedDFP_Val_Labels()) == FALSE."))
          }
        },
        error = function(e) {
          message("An error occurred in shiny::observeEvent(input$btnCountProbesP_ValParallel):\n", e)
        },
        warning = function(w) {
          message("A warning occurred in shiny::observeEvent(input$btnCountProbesP_ValParallel):\n", w)
        },
        finally = {
          base::print(base::paste0(Sys.time(), " count probes."))
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
          base::print(base::paste0(Sys.time(), " start counting probes delta methylation."))
          minN <- base::as.integer(input$txtCases)
          if (is.valid(session$userData$sessionVariables$combinedDFP_Val_Labels())) {
            DMNTable <-
              getAvailNForDMBorderParallel(session = session, wd = session$userData$packageWd, numCores = session$userData$numCores, DF = session$userData$sessionVariables$combinedDFP_Val_Labels()$dfDM)
            output$DTDMborder <- DT::renderDataTable(DMNTable)
            plot <- plotly::plot_ly(x = DMNTable$DM_BORDER, y = DMNTable$'Available CpG', type = "scatter", mode = "lines+markers", name = "scatterDeltaMethylationBorder") %>%
              plotly::layout(xaxis = list(title = 'DeltaMethylation', type = "linear")) %>%
              plotly::layout(yaxis = list(title = 'n' ))
            output$plotDendrogramDMborder <- plotly::renderPlotly(plot)
          }
          else {
            base::print(base::paste0(Sys.time(), " is.valid(session$userData$sessionVariables$combinedDFP_Val_Labels()) == FALSE."))
          }
        },
        error = function(e) {
          message("An error occurred in shiny::observeEvent(input$btnCountProbesDeltaMethParallel):\n", e)
        },
        warning = function(w) {
          message("A warning occurred in shiny::observeEvent(input$btnCountProbesDeltaMethParallel):\n", w)
        },
        finally = {
          base::print(base::paste0(Sys.time(), " finished counting probes delta methylation."))
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
          base::print(base::paste0(Sys.time(), " start counting probes n."))
          minN <- base::as.integer(input$txtCases)
          if (is.valid(session$userData$sessionVariables$combinedDFP_Val_Labels())) {
            NNTable <-
              getAvailNForNBorderParallel(session = session, wd = session$userData$packageWd, numCores = session$userData$numCores, DF = session$userData$sessionVariables$combinedDFP_Val_Labels()$dfN)
            output$DTNborder <- DT::renderDataTable(NNTable)
            plot <- plotly::plot_ly(x = NNTable$N_BORDER, y = NNTable$'Available CpG' , type = "scatter", mode = "lines+markers", name = "scatterNBorder") %>%
              plotly::layout(xaxis = list(title = 'n', type = "linear")) %>%
              plotly::layout(yaxis = list(title = 'n' ))
            output$plotDendrogramNborder <- plotly::renderPlotly(plot)
          }
          else {
            base::print(base::paste0(Sys.time(), " is.valid(session$userData$sessionVariables$combinedDFP_Val_Labels()) == FALSE."))
          }
        },
        error = function(e) {
          message("An error occurred in shiny::observeEvent(input$btnCountProbesNParallel):\n", e)
        },
        warning = function(w) {
          message("A warning occurred in shiny::observeEvent(input$btnCountProbesNParallel):\n", w)
        },
        finally = {
          base::print(base::paste0(Sys.time(), " finished counting probes n."))
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
          base::print(base::paste0(Sys.time(), " start counting probes."))
          minN <- base::as.integer(input$txtCases)
          #if (session$userData$sessionVariables$LoadInitialized == TRUE) {
          if (is.valid(session$userData$sessionVariables$resultDFListTrait1()) || is.valid(session$userData$sessionVariables$resultDFListTrait2())  || is.valid(session$userData$sessionVariables$resultDFListTrait3())) {
            combinedDFP_Val_Labels <- mergeDFP_Val_Labels(session$userData$sessionVariables$resultDFListTrait1(),
                                                          session$userData$sessionVariables$resultDFListTrait2(),
                                                          session$userData$sessionVariables$resultDFListTrait3(), minN)
          }
          else {
            base::print(base::paste0(Sys.time(), " is.valid(session$userData$sessionVariables$resultDFListTrait1()) || is.valid(session$userData$sessionVariables$resultDFListTrait2())  || is.valid(session$userData$sessionVariables$resultDFListTrait3()) == FALSE."))
            result <- NULL
          }
          P_VALNTable <-
            getAvailNForP_VALBorder(combinedDFP_Val_Labels$dfP_Val)
          output$DTP_VALborder <- DT::renderDataTable(P_VALNTable)
          base::print(base::paste0(Sys.time(), " finished counting probes."))
        },
        error = function(e) {
          message("An error occurred in shiny::observeEvent(input$btnCountP_ValProbes):\n", e)
        },
        warning = function(w) {
          message("A warning occurred in shiny::observeEvent(input$btnCountP_ValProbes):\n", w)
        },
        finally = {
          base::print(base::paste0(Sys.time(), " count probes."))
        }
      )
    },
    ignoreNULL = FALSE
  )

  output$histP_Val <- plotly::renderPlotly(histP_Val())

  histP_Val <- shiny::reactive({
    base::tryCatch({
      P_Val <- sort(as.numeric(unlist(session$userData$sessionVariables$traitReducedmatP_Val()))) #tbc(): here we encounter error: R character strings are limited to 2^31 bytes
      result <- plotly::plot_ly(x = P_Val, type = "histogram", name = "histP_Val")
    },
    error = function(e) {
      message("An error occurred in histP_Val <- shiny::reactive():\n", e)
      browser() #tbc() #check object.size() of session$userData$sessionVariables$traitReducedmatP_Val()
      a <- object.size(session$userData$sessionVariables$traitReducedmatP_Val())
    },
    warning = function(w) {
      message("A warning occurred in histP_Val <- shiny::reactive():\n", w)
    },
    finally = {
      base::print(base::paste0(Sys.time(), " end histP_Val <- shiny::reactive()."))
      return(result)
    })
  })

  shiny::observeEvent(input$btnExportCombinedData,
    ignoreInit = TRUE,
    {
      base::tryCatch(
        {
          fileNameCombinedHM <- base::paste0("CombinedHM.RDS")
          fileName <-
            base::paste0(session$userData$config$workDir, fileNameCombinedHM) #base::paste0(globalVariables$config$workDir, fileNameCombinedHM)
          base::print(base::paste0(Sys.time(), " start exporting session data to ", fileName, "."))
          base::saveRDS(file = fileName, session$userData)
          base::print(base::paste0(Sys.time(), " end exporting session data to ", fileName, "."))
        },
        error = function(e) {
          message("An error occurred in shiny::observeEvent(input$btnExportCombinedData):\n", e)
        },
        warning = function(w) {
          message("A warning occurred in shiny::observeEvent(input$btnExportCombinedData):\n", w)
        },
        finally = {
          base::print(base::paste0(Sys.time(), " end export combined data."))
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
          fileNameCombinedHM <- base::paste0("CombinedHM.RDS")
          fileName <-
            base::paste0(session$userData$config$workDir, fileNameCombinedHM) #base::paste0(globalVariables$config$workDir, fileNameCombinedHM)
          base::print(base::paste0(Sys.time(), " start importing data from ", fileName, "."))
          if (utils::file_test("-f", fileName) == TRUE) {
            #for whatever reason, we need to fill "session$userData$sessionVariables$resultDFListTrait1"... data structure with non-sense values (0)
            #before reading real values in order to fire reactivity, filling with NULL is not enough...
            session$userData$sessionVariables$resultDFListTrait1(0)
            session$userData$sessionVariables$resultDFListTrait2(0)
            session$userData$sessionVariables$resultDFListTrait3(0)
            session$userData$sessionVariables$combinedDFP_Val_Labels(0)
            session$userData$sessionVariables$pReducedcombinedDFP_Val_Labels(0)
            session$userData$sessionVariables$traitReducedcombinedDFP_Val_Labels(0)

            session$userData <-
              base::readRDS(file = fileName)
            base::print(base::paste0(Sys.time(), " end reading data"))
          }
          base::print(base::paste0(Sys.time(), " end importing session data from ", fileName, "."))
        },
        error = function(e) {
          message("An error occurred in shiny::observeEvent(input$btnImportCombinedData):\n", e)
        },
        warning = function(w) {
          message("A warning occurred in shiny::observeEvent(input$btnImportCombinedData):\n", w)
        },
        finally = {
          base::print(base::paste0(Sys.time(), " end import data."))
        }
      )
    },
    ignoreNULL = FALSE
  )

  output$DTP_VAL <- DT::renderDataTable(as.data.frame(session$userData$sessionVariables$traitReducedcombinedDFP_Val_Labels()$dfP_Val))

  output$DTDM <- DT::renderDataTable(as.data.frame(session$userData$sessionVariables$traitReducedcombinedDFP_Val_Labels()$dfDM))

  output$DTN <- DT::renderDataTable(as.data.frame(session$userData$sessionVariables$traitReducedcombinedDFP_Val_Labels()$dfN))

  session$userData$sessionVariables$traitReducedDendTraits <- shiny::reactive({
    base::print(base::paste0(Sys.time(), " before making dendrogram for traits"))
    #also re-read session$userData$sessionVariables$distMatProbes()
    #session$userData$sessionVariables$distMatProbes(...)
    return(stats::as.dendrogram(session$userData$sessionVariables$traitReducedclustResTraits()))
    base::print(base::paste0(Sys.time(), " after making dendrogram for traits"))
  })

  session$userData$sessionVariables$dendProbes <- shiny::reactive({
    base::print(base::paste0(Sys.time(), " before making dendrogram for probes"))
    result <- stats::as.dendrogram(session$userData$sessionVariables$clustResProbes())
    length(unlist(result))
    base::print(base::paste0(Sys.time(), " after making dendrogram for probes"))
    return(result)
  })

  DTProbes <- shiny::reactive({
    base::print(base::paste0(Sys.time(), " before making probe table."))
    if (!is.null(session$userData$sessionVariables$clustResProbes())) {
      dendProbes <- session$userData$sessionVariables$dendProbes()
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
      base::print(base::paste0(Sys.time(), " after making probe table."))
    }
    else {
      result <- NULL
    }
    return(result)
  })

  output$DTProbes <- DT::renderDataTable(as.data.frame(DTProbes()))

  DTTraits <- shiny::reactive({
    base::print(base::paste0(Sys.time(), " before making traits table."))
    if (!is.null(session$userData$sessionVariables$traitReducedclustResTraits())) {
      listTraits <- ggdendro::dendro_data(session$userData$sessionVariables$traitReducedclustResTraits(), type = "rectangle")$labels$label

      base::print(base::paste0(Sys.time(), " before rendering dendrogram tables traits"))
      DTTraits <-
        base::data.frame(row.names = seq_along(listTraits))
      DTTraits$Name <- listTraits
      DTTraits$order <- base::seq_len(base::nrow(DTTraits))
      rownames(DTTraits) <- DTTraits$Name
      DTTraits <- DTTraits[order(DTTraits$order), ]
      base::print(base::paste0(Sys.time(), " after making traits table."))
      result <- DTTraits
    }
    else {
      result <- NULL
    }
    return(result)
  })

  output$DTTraits <- DT::renderDataTable(as.data.frame(DTTraits()))

  output$plotDendrogramProbes <- plotly::renderPlotly(ggdendro::ggdendrogram(session$userData$sessionVariables$clustResProbes(),
                                                                             rotate = TRUE, theme_dendro = FALSE))

  output$plotDendrogramTraits <- plotly::renderPlotly(ggdendro::ggdendrogram(session$userData$sessionVariables$traitReducedclustResTraits(),
                                                                             rotate = TRUE, theme_dendro = FALSE))

  shiny::observeEvent(input$plotCombinedHM,
    ignoreInit = TRUE,
    {
      base::tryCatch(
        {
          plotCombinedHM(input = input, output = output, session = session)
        },
        error = function(e) {
          message("An error occurred in shiny::observeEvent(input$plotCombinedHM):\n", e)
        },
        warning = function(w) {
          message("A warning occurred in shiny::observeEvent(input$plotCombinedHM):\n", w)
        },
        finally = {
          base::print(base::paste0(Sys.time(), " end plotting heatmap."))
        }
      )
    },
    ignoreNULL = FALSE
  )

  # shiny::observeEvent(input$sldP_Val,
  #   ignoreInit = FALSE,
  #   {
  #     base::tryCatch(
  #       {
  #         if (!base::is.na(base::as.integer(input$sldP_Val[1]))) {
  #           session$userData$sessionVariables$P_ValMinBorder <-
  #             5 * 10^-base::as.integer(input$sldP_Val[2])
  #           session$userData$sessionVariables$P_ValMaxBorder <-
  #             5 * 10^-base::as.integer(input$sldP_Val[1])
  #         } else {
  #           session$userData$sessionVariables$P_ValMinBorder <- 5 * 10^-200
  #           session$userData$sessionVariables$P_ValMaxBorder <- 5 * 10^-18
  #         }
  #       },
  #       error = function(e) {
  #         message("An error occurred in shiny::observeEvent(input$sldP_Val):\n", e)
  #       },
  #       warning = function(w) {
  #         message("A warning occurred in shiny::observeEvent(input$sldP_Val):\n", w)
  #       },
  #       finally = {
  #         base::print(base::paste0(Sys.time(), " end adaptation p-values."))
  #       }
  #     )
  #   },
  #   ignoreNULL = FALSE
  # )

  shiny::observeEvent(input$chkDebug,
    ignoreInit = TRUE, #FALSE,
    {
      base::tryCatch(
        {
          if (input$chkDebug == TRUE) {
            session$userData$config$debugMode <- TRUE
            base::print(base::paste0(Sys.time(), " set debugMode = TRUE."))
          }
          else {
            session$userData$config$debugMode <- FALSE
            base::print(base::paste0(Sys.time(), " set debugMode = FALSE."))
          }
          result <- loadDirLists(session = session, input = input, output = output)
          dfdD1 <- result$dfdD1
          dfdD2 <- result$dfdD2
          dfdD3 <- result$dfdD3
        },
        error = function(e) {
          message("An error occurred in shiny::observeEvent(input$chkDebug):\n", e)
        },
        warning = function(w) {
          message("A warning occurred in shiny::observeEvent(input$chkDebug):\n", w)
        },
        finally = {
        }
      )
    },
    ignoreNULL = FALSE
  )
}
