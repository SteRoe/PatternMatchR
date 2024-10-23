PCPlot_UI <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    # shiny::verbatimTextOutput(ns("txtSelectedProbes"), placeholder = TRUE),
    # shiny::verbatimTextOutput(ns("txtSelectedTraits"), placeholder = TRUE),
    shiny::sliderInput(ns("DMRWindow"), "", #"DMR window size",
                       1, 50, 5, 1, width = "100%"),
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
                         DT::dataTableOutput(ns("traitReducedDTPCPlot"))
                       ),
                       shiny::tabPanel(
                         "Plot",
                         shiny::fluidRow(
                           shiny::column(width = 10,
                             plotly::plotlyOutput(ns("traitReducedPCPlot"),
                                                  width = "100%",
                                                  height = "800px")
                           ),
                           shiny::column(width = 2,
                                         plotly::plotlyOutput(ns("traitReducedViolinPlot"),
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
    # shinyjs::hide("txtSelectedProbes")
    # shinyjs::hide("txtSelectedTraits")
    #
    # output$txtSelectedProbes <- shiny::reactive({
    #   return(selectedProbes())
    # })
    #
    # output$txtSelectedTraits <- shiny::reactive({
    #   return(selectedTraits())
    # })

    output$traitReducedDTPCPlot <- DT::renderDataTable(as.data.frame(traitReducedDFPCplot()),
                                    options = list(pageLength = 1000, info = FALSE,
                                       lengthMenu = list(c(100, 1000, -1), c("100", "1000", "All"))))

    output$traitReducedPCPlot <- plotly::renderPlotly(traitReducedPCplot())

    output$traitReducedViolinPlot <- plotly::renderPlotly(traitReducedViolinplot())

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
      id <- shiny::showNotification("Creating DMP near range data...", duration = NULL, closeButton = FALSE)
      on.exit(shiny::removeNotification(id), add = TRUE)
      result <- NULL
      DMRWindow <- input$DMRWindow
      probeID <- selectedProbe() #DMP #the probeID in focus
      if(is.valid(probeID)){
        annotation <- session$userData$annotation
        if(is.valid(annotation)) {
          traitID <- selectedTrait()
          if(is.valid(traitID)) {
            probeIDs <- EpiVisR::getDMPNearRangeprobeID(annotation, probeID, DMRWindow)
            DMPNearRangeData <- as.data.frame(session$userData$Beta_tDF[, probeIDs])
            DMPNearRangeData$ID <- rownames(DMPNearRangeData)
            DMPNearRangeData <- base::merge(DMPNearRangeData, session$userData$baseData, by.x = "ID", by.y = session$userData$config$mergeAttribut)
            DMPNearRangeData$sex <- DMPNearRangeData[,session$userData$config$sexAttribut]
            DMPNearRangeData <- DMPNearRangeData[,c("ID","sex")]

            #traitLocation <- session$userData$sessionVariables$traitReducedDataStructurePVal()$dfKeyShadow[as.numeric(traitID),]
            traitLocation <- session$userData$sessionVariables$probeReducedDataStructure()$dfKeyShadow[as.numeric(traitID),]
            #get trait data based on traitLocation: selectedDF not found
            selectedOriginalDataTraits <- session$userData$sessionVariables$selectedOriginalDataTraits()
            #get the right variable out of selectedOriginalDataTraits
            trait <- removeAdjFromColname(traitLocation$trait)
            trait <- paste0(trait,".",traitLocation$traitSource)
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
      return (result)
    })

    traitReducedViolinplot <- shiny::reactive({
      id <- shiny::showNotification("Creating trait reduced violin plot...", duration = NULL, closeButton = FALSE)
      on.exit(shiny::removeNotification(id), add = TRUE)
      result <- NULL
      DMPNearRangeData <- DMPNearRangeData()
      result <- EpiVisR::plotlyViolinForDMP(session$userData$config$sexFemaleValue, session$userData$config$sexMaleValue, DMPNearRangeData)
      return(result)
    })

    traitReducedPCplot <- shiny::reactive({
      id <- shiny::showNotification("Creating trait reduced PC plot...", duration = NULL, closeButton = FALSE)
      on.exit(shiny::removeNotification(id), add = TRUE)
      result <- NULL
      DMRWindow <- input$DMRWindow
      probeID <- selectedProbe() #DMP #the probeID in focus
      if(is.valid(probeID)){
        annotation <- session$userData$annotation #session$userData$sessionVariables$selectedAnnotation()
        if(is.valid(annotation)) {
          traitID <- selectedTrait()
          if(is.valid(traitID)) {
            #traitLocation <- session$userData$sessionVariables$traitReducedDataStructurePVal()$dfKeyShadow[as.numeric(traitID),]
            traitLocation <- session$userData$sessionVariables$probeReducedDataStructure()$dfKeyShadow[as.numeric(traitID),]
            DMPNearRangeData <- DMPNearRangeData()
            # resultDataSingleTrait contains the regression result data for probes within the selected range for the selected trait...
            # we only extract gene.symbol, P_Val and DeltaMeth from there for legend of plot...
            resultDataSingleTrait <- getResultDataSingleTrait(session$userData$sessionVariables$combinedDataStructure()$combinedDFP_Val_Labels, traitLocation$key)
            resultDataSingleTrait$P_VAL <- resultDataSingleTrait$P_Val
            resultDataSingleTrait$P_Val <- NULL
            resultDataSingleTrait$DeltaMeth <- resultDataSingleTrait$DM
            resultDataSingleTrait$DM <- NULL
            result <- EpiVisR::plotlyPcPForDMP(DMPNearRange = DMPNearRangeData, probe = probeID, resultDataSingleTrait = resultDataSingleTrait, annotation = annotation)
          }
        }
      }
      return (result)
    })

    traitReducedDFPCplot <- shiny::reactive({
      id <- shiny::showNotification("Creating trait reduced data table for PC plot...", duration = NULL, closeButton = FALSE)
      on.exit(shiny::removeNotification(id), add = TRUE)
      result <- NULL
      result <- DMPNearRangeData()
      return(result)
    })

    #selectedTraits <- shiny::reactive({
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
        traits$label <- paste0(traits$traitName,".",traits$traitSource)
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

    #selectedProbes <- shiny::reactive({
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
      return(input$selSelectedTrait)
    })

    selectedProbe <- shiny::reactive({
      return(input$selSelectedProbe)
    })
  }) #end shiny::moduleServer
}
