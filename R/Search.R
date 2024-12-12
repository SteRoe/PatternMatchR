Search_Full_UI <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::fluidRow(
      shiny::verbatimTextOutput(ns("txtSearchStatus"), placeholder = TRUE),
    ),
    shiny::fluidRow(
      shiny::tabsetPanel(
        shiny::tabPanel(
              "Search for CpG in Full Heatmap (p-val)",
              shiny::textAreaInput(inputId = ns("txtSearchFullCpGPVal"), label = "Search CpG (p-val)", value = ""),
              shiny::verbatimTextOutput(outputId = ns("txtSearchResultFullCpGPVal")),
              shiny::actionButton(ns("btnSearchFullCpGHMPVal"), label = "Search CpG (p-val)"),
              "Search for Trait in Full Heatmap (p-val)",
              shiny::textAreaInput(inputId = ns("txtSearchFullTraitPVal"), label = "Search Trait (p-val)", value = ""),
              shiny::verbatimTextOutput(outputId = ns("txtSearchResultFullTraitPVal")),
              shiny::actionButton(ns("btnSearchFullTraitHMPVal"), label = "Search Trait")
        ),
        shiny::tabPanel(
          "Search for CpG in Full Heatmap (log(FC))",
          shiny::textAreaInput(inputId = ns("txtSearchFullCpGLogFC"), label = "Search CpG (log(FC))", value = ""),
          shiny::verbatimTextOutput(outputId = ns("txtSearchResultFullCpGLogFC")),
          shiny::actionButton(ns("btnSearchFullCpGHMLogFC"), label = "Search CpG (log(FC))"),
          "Search for Trait in Full Heatmap (p-val)",
          shiny::textAreaInput(inputId = ns("txtSearchFullTraitLogFC"), label = "Search Trait (log(FC))", value = ""),
          shiny::verbatimTextOutput(outputId = ns("txtSearchResultFullTraitLogFC")),
          shiny::actionButton(ns("btnSearchFullTraitHMLogFC"), label = "Search Trait (log(FC))")
        )
      )
    )
  )
}

# Search_Condensed_UI <- function(id) {
#   ns <- shiny::NS(id)
#   shinyBS::bsCollapsePanel(
#     "Search",
#     # shiny::fluidRow(
#     #   shiny::verbatimTextOutput(ns("txtSearchStatus"), placeholder = TRUE),
#     # ),
#     shiny::fluidRow(
#       shiny::tabsetPanel(
#         shiny::tabPanel(
#           "Search for CpG in Condensed Heatmap (p-val)",
#           shiny::textAreaInput(inputId = "txtSearchCondCpG", label = "Search CpG", value = ""),
#           shiny::verbatimTextOutput(outputId = "txtSearchResultCondCpG"),
#           shiny::actionButton("btnSearchCondCpGHM", label = "Search CpG"),
#           "Search for Trait in Condensed Heatmap",
#           shiny::textAreaInput(inputId = "txtSearchCondTrait", label = "Search Trait", value = ""),
#           shiny::verbatimTextOutput(outputId = "txtSearchResultCondTrait"),
#           shiny::actionButton("btnSearchCondTraitHM", label = "Search Trait")
#         ),
#         shiny::tabPanel(
#           "reserved"
#         )
#       )
#     )
#   )
# }

Search_Full_SERVER <- function(id, session) {
  shiny::moduleServer(id, function(input, output, session) {

    shiny::observeEvent(input$btnSearchFullCpGHMPVal,
                        ignoreInit = TRUE,
                        {
                          base::tryCatch(
                            {
                              base::print(base::paste0(sysTimePID(), " start searching CpG full."))
                              #find positions
                              searchResult <- getSearchResultCpG(input$txtSearchFullCpGPVal, session$userData$sessionVariables$traitReducedDataStructurePVal())
                              length <- length(session$userData$sessionVariables$traitReducedDataStructurePVal()$clustResProbes$labels)
                              resultText <- paste0(base::trimws(input$txtSearchFullCpGPVal), " found at position: ", searchResult, " from ", length, " CpG.")
                              #write to GlobalSelection:
                              if (is.valid(searchResult)) {
                                selectedProbe <- session$userData$sessionVariables$traitReducedDataStructurePVal()$clustResProbes$labels[searchResult]
                                session$userData$sessionVariables$selectedProbe(selectedProbe)
                              }
                              #write to output
                              output$txtSearchResultFullCpGPVal <- shiny::renderText(resultText)
                            },
                            error = function(e) {
                              if (attributes(e)$class[1] != "shiny.silent.error") {
                                base::message("An error occurred in shiny::observeEvent(input$btnSearchFullCpGHMPVal):\n", e)
                              }
                            },
                            warning = function(w) {
                              base::message("A warning occurred in shiny::observeEvent(input$btnSearchFullCpGHMPVal):\n", w)
                            },
                            finally = {
                              base::print(base::paste0(sysTimePID(), " finished searching CpG full."))
                            }
                          )

                        },
                        ignoreNULL = FALSE
    )

    shiny::observeEvent(input$btnSearchFullTraitHMPVal,
                        ignoreInit = TRUE,
                        {
                          base::tryCatch(
                            {
                              base::print(base::paste0(sysTimePID(), " start searching trait full."))
                              #find positions
                              searchResult <- getSearchResultTrait(input$txtSearchFullTraitPVal, session$userData$sessionVariables$traitReducedDataStructurePVal())
                              length <- length(session$userData$sessionVariables$traitReducedDataStructurePVal()$clustResTraits$labels)
                              resultText <- paste0(base::trimws(input$txtSearchFullTraitPVal), " found at position: ", searchResult, " from ", length, " traits.")
                              #write to GlobalSelection:
                              if (is.valid(searchResult)) {
                                originTrait <- session$userData$sessionVariables$traitReducedDataStructurePVal()$combinedDFP_Val_Labels$mergedOriginTrait[searchResult]
                                traitLabels <- session$userData$sessionVariables$traitReducedDataStructurePVal()$combinedDFP_Val_Labels$mergedOriginalColnames[searchResult]
                                selectedTrait <- cbind(traitLabels, originTrait)
                                colnames(selectedTrait) <- c("traitName", "traitSource")
                                session$userData$sessionVariables$selectedTrait(selectedTrait)
                              }
                              #write to output
                              output$txtSearchResultFullTraitPVal <- shiny::renderText(resultText)
                            },
                            error = function(e) {
                              if (attributes(e)$class[1] != "shiny.silent.error") {
                                base::warning("An error occurred in shiny::observeEvent(input$btnSearchFullTraitHMPVal):\n", e)
                              }
                            },
                            warning = function(w) {
                              base::message("A warning occurred in shiny::observeEvent(input$btnSearchFullTraitHMPVal):\n", w)
                            },
                            finally = {
                              base::print(base::paste0(sysTimePID(), " end search trait heatmap full."))
                              base::print(base::paste0(sysTimePID(), " finished searching trait full."))
                            }
                          )
                        },
                        ignoreNULL = FALSE
    )

    shiny::observeEvent(input$btnSearchCondCpGHM,
                        ignoreInit = TRUE,
                        {
                          base::tryCatch(
                            {
                              base::print(base::paste0(sysTimePID(), " start searching CpG condensed."))
                              #find positions
                              searchResult <- getSearchResultCpG(input$txtSearchCondCpG, session$userData$sessionVariables$probeReducedDataStructure())
                              length <- length(session$userData$sessionVariables$probeReducedDataStructure()$clustResProbes$labels)
                              resultText <- paste0(base::trimws(input$txtSearchCondCpG), " found at position: ", searchResult, " from ", length, " CpG.")
                              #write to output
                              output$txtSearchResultCondCpG <- shiny::renderText(resultText)
                            },
                            error = function(e) {
                              if (attributes(e)$class[1] != "shiny.silent.error") {
                                base::message("An error occurred in shiny::observeEvent(input$btnSearchCondCpGHM):\n", e)
                              }
                            },
                            warning = function(w) {
                              base::message("A warning occurred in shiny::observeEvent(input$btnSearchCondCpGHM):\n", w)
                            },
                            finally = {
                              base::print(base::paste0(sysTimePID(), " finished searching CpG condensed."))
                            }
                          )

                        },
                        ignoreNULL = FALSE
    )

    shiny::observeEvent(input$btnSearchCondTraitHM,
                        ignoreInit = TRUE,
                        {
                          base::tryCatch(
                            {
                              base::print(base::paste0(sysTimePID(), " start searching trait cond."))
                              #find positions
                              searchResult <- getSearchResultTrait(input$txtSearchCondTrait, session$userData$sessionVariables$probeReducedDataStructure())
                              length <- length(session$userData$sessionVariables$probeReducedDataStructure()$clustResTraits$labels)
                              resultText <- paste0(base::trimws(input$txtSearchCondTrait), " found at position: ", searchResult, " from ", length, " traits.")
                              #write to output
                              output$txtSearchResultCondTrait <- shiny::renderText(resultText)
                            },
                            error = function(e) {
                              if (attributes(e)$class[1] != "shiny.silent.error") {
                                base::warning("An error occurred in shiny::observeEvent(input$btnSearchCondTraitHM):\n", e)
                              }
                            },
                            warning = function(w) {
                              base::message("A warning occurred in shiny::observeEvent(input$btnSearchCondTraitHM):\n", w)
                            },
                            finally = {
                              base::print(base::paste0(sysTimePID(), " end search trait heatmap cond."))
                              base::print(base::paste0(sysTimePID(), " finished searching trait cond."))
                            }
                          )
                        },
                        ignoreNULL = FALSE
    )

    getSearchResultCpG <- function(txtSearchCpG, dataStructure) {
      base::tryCatch(
        {
          base::print(base::paste0(sysTimePID(), " searching CpG"))
          #find position of CpG
#          CpG <- dataStructure$clustResProbes$labels #[dataStructure$clustResProbes$order]
          CpG <- rownames(dataStructure$combinedDFP_Val_Labels$dfP_Val)
          positions <- base::which(CpG %in% unlist(base::strsplit(base::trimws(txtSearchCpG), " ")))
        },
        error = function(e) {
          if (attributes(e)$class[1] != "shiny.silent.error") {
            base::message("An error occurred in getSearchResultCpG():\n", e)
          }
        },
        warning = function(w) {
          base::message("A warning occurred in getSearchResultCpG():\n", w)
        },
        finally = {
          base::print(base::paste0(sysTimePID(), " end getSearchResultCpG()."))
          base::return(positions)
        }
      )
    }

    getSearchResultTrait <- function(txtSearchTrait, dataStructure) {
      base::tryCatch(
        {
          base::print(base::paste0(sysTimePID(), " searching Trait"))
          #find position of Trait
          Trait <- dataStructure$combinedDFP_Val_Labels$mergedColnames #clustResTraits$labels #[dataStructure$clustResTraits$order]
          positions <- base::which(Trait %in% unlist(base::strsplit(base::trimws(txtSearchTrait), " ")))
        },
        error = function(e) {
          if (attributes(e)$class[1] != "shiny.silent.error") {
            base::message("An error occurred in getSearchResultTrait():\n", e)
          }
        },
        warning = function(w) {
          base::message("A warning occurred in getSearchResultTrait():\n", w)
        },
        finally = {
          base::print(base::paste0(sysTimePID(), " end getSearchResultTrait()."))
          base::return(positions)
        }
      )
    }
  })
}
