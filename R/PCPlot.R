PCPlot_UI <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::sliderInput(ns("DMRWindow"), "", #"DMR window size",
                       1, 150, 5, 1, width = "100%"),# 1, 50, 5, 1, width = "100%"),
    shiny::selectizeInput(ns("selSelectedProbe"),
                          label = "select probe for DMR window (if no probe occurs here for selection, then there was no probe selected)",
                          choices = NULL,
                          width = "100%"
    ),
    shiny::selectizeInput(ns("selSelectedTrait"),
                          label = "select trait for color coding PC plot (if no trait occurs here for selection, then there was no trait selected)",
                          choices = NULL,
                          width = "100%"
    ),
    shiny::tabsetPanel(id = ns("tabsetPC"),
                       shiny::tabPanel(
                         "Table",
                         DT::dataTableOutput(ns("DTPCPlot"))
                       ),
                       shiny::tabPanel(
                         "Plot",
                         shiny::fluidRow(
                           shiny::column(width = 10,
                             plotly::plotlyOutput(ns("PCPlot"),
                                                  width = "100%",
                                                  height = "800px")
                           ),
                           shiny::column(width = 2,
                                         plotly::plotlyOutput(ns("ViolinPlot"),
                                                              width = "100%",
                                                              height = "800px")
                           )
                         )
                       )
      )
#    )
  )
}

PCPlot_SERVER <- function(id, session) {
  shiny::moduleServer(id, function(input, output, session) {
    output$DTPCPlot <- DT::renderDataTable(as.data.frame(DFPCplot()),
                                    options = list(pageLength = 1000, info = FALSE,
                                       lengthMenu = list(c(100, 1000, -1), c("100", "1000", "All"))))

    output$PCPlot <- plotly::renderPlotly(PCplot())

    output$ViolinPlot <- plotly::renderPlotly(Violinplot())

    getResultDataSingleTrait <- function(combinedDFP_Val_Labels, traitID) {
      result <- NULL
      rownames <- rownames(combinedDFP_Val_Labels$dfP_Val)
      if (is.valid(traitID)) {
        P_Val <- combinedDFP_Val_Labels$dfP_Val[,traitID]
        DM <- combinedDFP_Val_Labels$dfDM[,traitID]
        N <- combinedDFP_Val_Labels$dfN[,traitID]
        LogFC <- combinedDFP_Val_Labels$dfLogFC[,traitID]
        result <- cbind(P_Val, DM, N,LogFC)
        result <- as.data.frame(result)
        result$probeID <- rownames
      }
      return(result)
    }

    DMPNearRangeData <- shiny::reactive({
      shinyId <- shiny::showNotification("Creating DMP near range data...", duration = NULL, closeButton = FALSE)
      on.exit(shiny::removeNotification(shinyId), add = TRUE)
      base::tryCatch(
      {
        result <- NULL
        DMRWindow <- input$DMRWindow
        probeID <- selectedProbe() #DMP #the probeID in focus

        if(is.valid(probeID)){
          annotation <- session$userData$annotation
          if(is.valid(annotation)) {
            traitID <- selectedTrait()
            if(is.valid(traitID)) {
              probeIDs <- EpiVisR::getDMPNearRangeprobeID(annotation, probeID, DMRWindow)
              #probeIDs <- getDMPNearRangeprobeID(annotation, probeID, DMRWindow)
              probeIDs <- probeIDs[which(probeIDs %in% colnames(session$userData$Beta_tDF))] #probeIDs <- session$userData$Beta_tDF[, which(probeIDs %in% colnames(session$userData$Beta_tDF))]
              DMPNearRangeData <- as.data.frame(session$userData$Beta_tDF[, probeIDs])
              DMPNearRangeData$ID <- rownames(DMPNearRangeData)
              DMPNearRangeData <- base::merge(DMPNearRangeData, session$userData$baseData, by.x = "ID", by.y = session$userData$config$mergeAttribut)
              DMPNearRangeData$sex <- DMPNearRangeData[,session$userData$config$genderAttribut]
              DMPNearRangeData <- DMPNearRangeData[,c("ID","sex")]

              traitLocation <- session$userData$sessionVariables$probeReducedDataStructurePVal()$dfKeyShadow[as.numeric(traitID),]
              #get trait data based on traitLocation: selectedDF not found
              selectedOriginalDataTraits <- session$userData$sessionVariables$selectedOriginalDataTraits()
              #get the right variable out of selectedOriginalDataTraits
              trait <- removeAdjFromColname(traitLocation$trait)
              if (traitLocation$traitSource == "1") {
                traitSource <- "red_"
              }
              else if (traitLocation$traitSource == "2") {
                traitSource <- "green_"
              }
              else if (traitLocation$traitSource == "3") {
                traitSource <- "blue_"
              }
              else {
                browser() #should not happen
              }
              trait <- paste0(traitSource, traitLocation$trait) #trait <- paste0(trait, ".", traitLocation$traitSource)
              traitVar <- get(trait,selectedOriginalDataTraits)
              traitVar <- as.data.frame(traitVar)
              colnames(traitVar) <- trait
              traitVar$ID <- rownames(selectedOriginalDataTraits)
              rownames(traitVar) <- traitVar$ID

              DMPNearRangeData <- base::merge(DMPNearRangeData, traitVar, by.x = "ID",by.y = "ID")
              DMPNearRangeData2 <- as.data.frame(session$userData$Beta_tDF[, probeIDs])
              DMPNearRangeData2$ID <- rownames(DMPNearRangeData2)
              DMPNearRangeData <- base::merge(DMPNearRangeData, DMPNearRangeData2, by.x = "ID", by.y = "ID")
              DMPNearRangeData2 <- NULL
              result <- DMPNearRangeData
            }
          }
        }
      },
      error = function(e) {
        if(attributes(e)$class[1] != "shiny.silent.error") {
          base::message("An error occurred in DMPNearRangeData:\n", e)
          browser() #should not happen
        }
      },
      warning = function(w) {
        base::message("An error occurred in DMPNearRangeData:\n", w)
        browser() #should not happen
      },
      finally = {
        base::print(base::paste0(sysTimePID(), " finished DMPNearRangeData"))
        return(result)
      })
    })

    Violinplot <- shiny::reactive({
      shinyId <- shiny::showNotification("Creating trait reduced violin plot...", duration = NULL, closeButton = FALSE)
      on.exit(shiny::removeNotification(shinyId), add = TRUE)
      result <- NULL
      DMPNearRangeData <- DMPNearRangeData()
      result <- EpiVisR::plotlyViolinForDMP(session$userData$config$genderFemaleValue, session$userData$config$genderMaleValue, DMPNearRangeData)
      return(result)
    })

    PCplot <- shiny::reactive({
      shinyId <- shiny::showNotification("Creating trait reduced PC plot...", duration = NULL, closeButton = FALSE)
      on.exit(shiny::removeNotification(shinyId), add = TRUE)
      base::tryCatch({
        result <- NULL
#        DMRWindow <- input$DMRWindow
        probeID <- selectedProbe() #DMP #the probeID in focus
        if(is.valid(probeID)){
          annotation <- session$userData$annotation
          if(is.valid(annotation)) {
            traitID <- selectedTrait()
            if(is.valid(traitID)) {
              traitLocation <- session$userData$sessionVariables$probeReducedDataStructurePVal()$dfKeyShadow[as.numeric(traitID),]
              DMPNearRangeData <- DMPNearRangeData()
              # resultDataSingleTrait contains the regression result data for probes within the selected range for the selected trait...
              # we only extract gene.symbol, P_Val and DeltaMeth from there for legend of plot...
              resultDataSingleTrait <- getResultDataSingleTrait(session$userData$sessionVariables$combinedDataStructure()$combinedDFP_Val_Labels, traitLocation$key)
              resultDataSingleTrait$P_VAL <- resultDataSingleTrait$P_Val
              resultDataSingleTrait$P_Val <- NULL
              resultDataSingleTrait$DeltaMeth <- resultDataSingleTrait$DM
              resultDataSingleTrait$DM <- NULL
              #result <- plotlyPcPForDMP(DMPNearRange = DMPNearRangeData, probe = probeID, resultDataSingleTrait = resultDataSingleTrait, annotation = annotation, selection = session$userData$sessionVariables$selectedProbe(), shortlabel = TRUE)
              result <- EpiVisR::plotlyPcPForDMP(DMPNearRange = DMPNearRangeData, probe = probeID, resultDataSingleTrait = resultDataSingleTrait, annotation = annotation, selection = session$userData$sessionVariables$selectedProbe(), shortlabel = TRUE)
            }
          }
        }
      },
      error = function(e) {
        if (attributes(e)$class[1] != "shiny.silent.error") {
          base::message("An error occurred in PCplot:\n", e)
          browser() #should not happen
        }
      },
      warning = function(w) {
        base::message("An error occurred in PCplot:\n", w)
        browser() #should not happen
      },
      finally = {
        base::print(base::paste0(sysTimePID(), " finished PCplot"))
        return(result)
      })
    })

    DFPCplot <- shiny::reactive({
      shinyId <- shiny::showNotification("Creating trait reduced data table for PC plot...", duration = NULL, closeButton = FALSE)
      on.exit(shiny::removeNotification(shinyId), add = TRUE)
      base::tryCatch({
        result <- NULL
        result <- DMPNearRangeData()
      },
      error = function(e) {
        if (attributes(e)$class[1] != "shiny.silent.error") {
          base::message("An error occurred in DFPCplot:\n", e)
          browser() #should not happen
        }
      },
      warning = function(w) {
        base::message("An error occurred in DFPCplot:\n", w)
        browser() #should not happen
      },
      finally = {
        base::print(base::paste0(sysTimePID(), " finished  DFPCplot"))
        return(result)
      })

      return(result)
    })

    shiny::observe({
      shiny::freezeReactiveValue(input, "selSelectedTrait")
      if (is.valid(session$userData$sessionVariables$selectedTraitID())) {
        traits <- session$userData$sessionVariables$selectedTrait()
        #remove names with adj at the end, because we don't have these variables in reality...
        traits <- removeAdjFromColname(traits)
        traits <- as.data.frame(traits)
        rownames(traits) <- seq(1:nrow(traits))
        traits <- unique(traits)
        traits$ID <- rownames(traits)
        traitsIDs <- session$userData$sessionVariables$selectedTraitID()
        traitsIDs <- as.data.frame(traitsIDs)
        traitsIDs$ID <- seq(1:nrow(traitsIDs))
        traitsIDs <- traitsIDs[traitsIDs$ID %in% traits$ID, ]
        traitSource <- ifelse(traits$traitSource == "1", "red_", ifelse(traits$traitSource == "2", "green_", ifelse(traits$traitSource == "3", "blue_","")))
        traits$label <- paste0(traitSource, traits$traitName)
        traits$value <- traitsIDs$traitsIDs
      } else {
        traits <- NULL
      }
      shiny::updateSelectizeInput(
        session = session,
        inputId = "selSelectedTrait",
        choices = traits,
        server = TRUE
      )
      return(traits$label)
    })

    shiny::observe({
      shiny::freezeReactiveValue(input, "selSelectedProbe")
      if (is.valid(session$userData$sessionVariables$selectedProbe())) {
        probes <- as.list(session$userData$sessionVariables$selectedProbe())
      } else {
        probes <- NULL
      }
      shiny::updateSelectizeInput(
        session = session,
        inputId = "selSelectedProbe",
        choices = probes,
        server = TRUE
      )
      return (probes)
    })

    selectedTrait <- shiny::reactive({
      base::tryCatch({
        result <- input$selSelectedTrait
      },
      error = function(e) {
        if (attributes(e)$class[1] != "shiny.silent.error") {
          base::message("An error occurred in selectedTrait():\n", e)
          browser() #should not happen
        }
      },
      warning = function(w) {
        base::message("An error occurred in selectedTrait():\n", w)
        browser() #should not happen
      },
      finally = {
        base::print(base::paste0(sysTimePID(), " finished  selectedTrait()"))
        return(result)
      })
    })

    selectedProbe <- shiny::reactive({
      base::tryCatch({
        result <- input$selSelectedProbe
      },
      error = function(e) {
        if(attributes(e)$class[1] != "shiny.silent.error") {
          base::message("An error occurred in selectedProbe():\n", e)
          browser() #should not happen
        }
      },
      warning = function(w) {
        base::message("An error occurred in selectedProbe():\n", w)
        browser() #should not happen
      },
      finally = {
        base::print(base::paste0(sysTimePID(), " finished  selectedProbe()"))
        return(result)
      })
    })
  }) #end shiny::moduleServer
}
