VolcanoPlot_UI <- function(id) {
  ns <- shiny::NS(id)
  htmltools::tagList(
    shiny::tabsetPanel(id = ns("tabsetVolcano"),
      shiny::tabPanel(
        "Table",
        DT::dataTableOutput(ns("traitReducedDTVolcano"))
      ),
      shiny::tabPanel(
        "Plot",
        shiny::radioButtons(ns("RbHighlightVolcano"), "Highlight mode:",
                            choiceNames = list(
                              "Chromosome",
                              "Global Selection",
                              "DistanceNumber",
                              "DistanceMin"
                            ),
                            choiceValues = list("chr", "sel", "DistanceNumber", "DistanceMin"),
                            inline = TRUE
        ),
        plotly::plotlyOutput(ns("traitReducedPlotVolcano"), height = "150%")
      )
    )
  )
}

VolcanoPlot_SERVER <- function(id, session) {
  shiny::moduleServer(id, function(input, output, session) {

    output$traitReducedDTVolcano <- DT::renderDataTable(traitReducedDTVolcano())
    output$traitReducedPlotVolcano <- plotly::renderPlotly(traitReducedVolcano())

    shiny::observeEvent(plotly::event_data(event = "plotly_selected", source = "A"), {
      selected <- plotly::event_data(event = "plotly_selected", source = "A")
      #req(selected)
      if (is.valid(selected)) {
        session$userData$sessionVariables$selectedKey(selected$key)
        #get dfKeyShadow
        # dfKeyShadow <- session$userData$sessionVariables$traitReducedDataStructurePVal()$dfKeyShadow
        dfKeyShadow <- session$userData$sessionVariables$probeReducedDataStructure()$dfKeyShadow
        #get probe, trait and traitSource from keyShadow
        rowindex <- which(dfKeyShadow$key %in% selected$key)
        dfKeyShadow <- dfKeyShadow[rowindex,]
        probes <- dfKeyShadow$probe
        probes <- base::unique(probes)
        session$userData$sessionVariables$selectedProbe(probes)
        traits <- dfKeyShadow[,c("trait","traitSource")] #dfKeyShadow$trait
        traits <- base::unique(traits)
        traitIDs <- dfKeyShadow$traitID
        session$userData$sessionVariables$selectedTrait(traits)
        session$userData$sessionVariables$selectedTraitID(traitIDs)
        base::print(base::paste0(sysTimePID(), " last selection was from volcano plot selected event."))
        output$txtGlobalSelectOut <- shiny::renderText("last selection was from volcano plot selected event.")
      }
    })

    #leave out plotly_click event for the moment, because it fires after plotly_selected event and destroys selection made there
    # shiny::observeEvent(plotly::event_data(event = "plotly_click", source = "A"), {
    #   selected <- plotly::event_data(event = "plotly_click", source = "A")
    #   req(selected)
    #   session$userData$sessionVariables$selectedKey(selected$key)
    #   #get dfKeyShadow
    #   dfKeyShadow <- session$userData$sessionVariables$traitReducedDataStructurePVal()$dfKeyShadow
    #   #get probe, trait and traitSource from keyShadow
    #   rowindex <- which(dfKeyShadow$key %in% selected$key)
    #   dfKeyShadow <- dfKeyShadow[rowindex,]
    #   probes <- dfKeyShadow$probe
    #   probes <- base::unique(probes)
    #   session$userData$sessionVariables$selectedProbe(probes)
    #   traits <- dfKeyShadow$trait
    #   traits <- base::unique(traits)
    #   session$userData$sessionVariables$selectedTrait(traits)
    #   base::print(base::paste0(sysTimePID(), " last selection was from volcano plot click event."))
    #   output$txtGlobalSelectOut <- shiny::renderText("last selection was from volcano plot click event.")
    #   #get nearby CpG
    #
    #   #highlight them
    #   #plotly::highlight_key()
    #   #plotly::add_trace()
    #   #or highlight single point by plotly::add_trace
    # })

    traitReducedVolcano <- shiny::reactive({
      id <- shiny::showNotification("Creating trait reduced volcano plot...", duration = NULL, closeButton = FALSE)
      on.exit(shiny::removeNotification(id), add = TRUE)
      base::tryCatch({
        base::print(base::paste0(sysTimePID(), " start render plotly traitReducedVolcano()."))
        #req(session$userData$sessionVariables$traitReducedDataStructurePVal()$dfVolcano)
        highlightMode <- input$RbHighlightVolcano
        #if (is.valid(session$userData$sessionVariables$traitReducedDataStructurePVal()$dfVolcano)) {
        if (is.valid(session$userData$sessionVariables$probeReducedDataStructure()$dfVolcano)) {
          #volcano <- session$userData$sessionVariables$traitReducedDataStructurePVal()$dfVolcano
          volcano <- session$userData$sessionVariables$probeReducedDataStructure()$dfVolcano
          volcano$chr <- as.factor(volcano$chr)
          volcano$label <- paste0(volcano$probe, " ", volcano$trait, " chr=", volcano$chr, " pos=", volcano$position)
          if (highlightMode == "chr") {
            # colors <- RColorBrewer::brewer.pal(12, "Paired")
            # colors <- grDevices::colorRampPalette(colors)(22)
            colors <- viridis::inferno(22)
            #add colour
            plot <- plotly::plot_ly(volcano, x = ~LogFC, y = ~P_Val, color = ~chr, colors = colors, text = ~label, key = ~key, type = "scatter", mode = "markers", name = ~chr, source = "A") %>%
              plotly::layout(yaxis = list(title = "-log(p-val)", type = "log", autorange = "reversed", showexponent = "all", exponentformat = "e")) %>%
              plotly::layout(xaxis = list(title = "log(FC)")) %>%
              plotly::event_register("plotly_selected") %>%
              plotly::event_register("plotly_click")
            #      plotly::event_register("plotly_hover")
            #  plotly::highlight(on = "plotly_click", off = "plotly_doubleclick")
            #tbc() #add user defined hover function to show similar cg# (same trait, same cg#)
          }
          else if (highlightMode == "DistanceNumber") {
            #here we introduce the distance as color code
            volcano$DistanceNumber <- as.factor(volcano$DistanceNumber)
            colors <- viridis::inferno(length(levels(volcano$DistanceNumber)))
            plot <- plotly::plot_ly(volcano, x = ~LogFC, y = ~P_Val, color = ~DistanceNumber, colors = colors, text = ~label, key = ~key, type = "scatter", mode = "markers", name = ~DistanceNumber, source = "A") %>%
              plotly::layout(yaxis = list(title = "-log(p-val)", type = "log", autorange = "reversed", showexponent = "all", exponentformat = "e")) %>%
              plotly::layout(xaxis = list(title = "log(FC)")) %>%
              plotly::event_register("plotly_selected") %>%
              plotly::event_register("plotly_click")
          }
          else if (highlightMode == "DistanceMin") {
            volcano$DistanceMin <- as.factor(volcano$DistanceMin)
            colors <- viridis::inferno(length(levels(volcano$DistanceMin)))
            plot <- plotly::plot_ly(volcano, x = ~LogFC, y = ~P_Val, color = ~DistanceMin, colors = colors, text = ~label, key = ~key, type = "scatter", mode = "markers", name = ~DistanceMin, source = "A") %>%
              plotly::layout(yaxis = list(title = "-log(p-val)", type = "log", autorange = "reversed", showexponent = "all", exponentformat = "e")) %>%
              plotly::layout(xaxis = list(title = "log(FC)")) %>%
              plotly::event_register("plotly_selected") %>%
              plotly::event_register("plotly_click")
          }
          else if (highlightMode == "sel") {
            selectedProbe <- session$userData$sessionVariables$selectedProbe()
            positions <- base::which(volcano$probe %in% selectedProbe)
            volcano$selected <- FALSE
            if (is.valid(session$userData$sessionVariables$selectedProbe())) {
              volcano[positions,]$selected <- TRUE
            }
            volcano$selected <- base::as.factor(volcano$selected)
            colors <- viridis::inferno(length(levels(volcano$selected)))
            plot <- plotly::plot_ly(volcano, x = ~LogFC, y = ~P_Val, color = ~selected, colors = colors, text = ~label, key = ~key, type = "scatter", mode = "markers", name = ~selected, source = "A") %>%
              plotly::layout(yaxis = list(title = "-log(p-val)", type = "log", autorange = "reversed", showexponent = "all", exponentformat = "e")) %>%
              plotly::layout(xaxis = list(title = "log(FC)")) %>%
              plotly::event_register("plotly_selected") %>%
              plotly::event_register("plotly_click")
          }
          else {
            browser() #should not happen
            plot <- NULL
          }
        }
        else {
          browser() #should not happen
          plot <- NULL
        }
      },
      error = function(e) {
        base::message("An error occurred in traitReducedVolcano <- shiny::reactive():\n", e)
        browser() #should not happen
      },
      warning = function(w) {
        base::message("A warning occurred in traitReducedVolcano <- shiny::reactive():\n", w)
      },
      finally = {
        base::print(base::paste0(sysTimePID(), " finished render plotly traitReducedVolcano()."))
        plot # return(plot)
      })
    })

    traitReducedDTVolcano <- shiny::reactive({
      id <- shiny::showNotification("Creating trait reduced data table...", duration = NULL, closeButton = FALSE)
      on.exit(shiny::removeNotification(id), add = TRUE)
      selectedKey <- session$userData$sessionVariables$selectedKey()
      if(isTruthy(selectedKey)) {
        DT <- session$userData$sessionVariables$probeReducedDataStructure()$dfVolcano
        #DT <- session$userData$sessionVariables$traitReducedDataStructurePVal()$dfVolcano
        DT <- dplyr::filter(DT, key %in% selectedKey)
      }
      else {
        #DT <- as.data.frame(session$userData$sessionVariables$traitReducedDataStructurePVal()$dfVolcano)
        DT <- as.data.frame(session$userData$sessionVariables$probeReducedDataStructure()$dfVolcano)
      }
      return(DT)
    })

    probeReducedVolcano <- shiny::reactive({
      id <- shiny::showNotification("Creating probe reduced volcano plot...", duration = NULL, closeButton = FALSE)
      on.exit(shiny::removeNotification(id), add = TRUE)
      base::tryCatch({
        base::print(base::paste0(sysTimePID(), " start render plotly traitReducedVolcano()."))
        volcano <- session$userData$sessionVariables$probeReducedDataStructure()$dfVolcano
        volcano$label <- paste0(volcano$probe, " ", volcano$trait)
        plot <- plotly::plot_ly(volcano, x = ~LogFC, y = ~P_Val, text = ~label, key = ~key, type = "scatter", mode = "markers", name = "volcano", source = "B") %>%
          plotly::layout(yaxis = list(title = "-log(p-val)", type = "log", autorange = "reversed")) %>%
          plotly::layout(xaxis = list(title = "log(FC)"))
        #tbc() #add user defined hover function to show similar cg# (same trait, same cg#)
        #plotly::event_register(plot, "plotly_hover")
      },
      error = function(e) {
        base::message("An error occurred in traitReducedVolcano <- shiny::reactive():\n", e)
        browser() #should not happen
      },
      warning = function(w) {
        base::message("A warning occurred in traitReducedVolcano <- shiny::reactive():\n", w)
      },
      finally = {
        base::print(base::paste0(sysTimePID(), " finished render plotly traitReducedVolcano()."))
        return(plot)
      })
    })
  }) #end shiny::moduleServer
}
