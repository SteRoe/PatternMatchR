GlobalSelection_UI <- function(id) {
  ns <- shiny::NS(id)
  htmltools::tagList(
    #waiter::autoWaiter()
    # shinyBS::bsCollapsePanel(
    #   "Global Selection",
      shiny::fluidRow(
        shiny::verbatimTextOutput(ns("txtGlobalSelectOut"), placeholder = TRUE),
      ),
      shiny::fluidRow(
        shiny::tabsetPanel(
          shiny::tabPanel(
            "Selected Trait",
            DT::dataTableOutput(ns("DTSelectedTrait"))
          ),
          shiny::tabPanel(
            "Selected Probe",
            DT::dataTableOutput(ns("DTSelectedProbe"))
          ),
          shiny::tabPanel(
            "SPLOM of Original Data for selected area",
            shiny::selectizeInput(ns("markingVar"),
                                  label = "select variable for color marking (if no variable occurs here for selection, then there was no factor variable selected)",
                                  choices = NULL,
                                  width = "100%"
            ),
            shiny::tabsetPanel(
              shiny::tabPanel(
                "SPLOM from selected area",
                plotly::plotlyOutput(ns("SPLOM"),
                                     height = '95vh', #1200, #height = "100%",
                                     width = '95%', #100, #width = "100%",
                                     inline = TRUE)
              ),
              shiny::tabPanel(
                "SPLOM trait/trait",
                plotly::plotlyOutput(ns("SPLOMTrait"),
                                     height = '95vh', #1200, #height = "100%",
                                     width = '95%', #100, #width = "100%",
                                     inline = TRUE)
              ),
              shiny::tabPanel(
                "SPLOM probe/probe",
                shiny::tags$html("if there is no plot visible here, then probably not the full methylation data set was loaded (for debug reasons?)"),
                plotly::plotlyOutput(ns("SPLOMProbe"),
                                     height = '95vh', #1200, #height = "100%",
                                     width = '95%', #1000, #width = "100%",
                                     inline = TRUE)
              )
            )
          )
        )
      )
#    )
  )
}

GlobalSelection_SERVER <- function(id, session) {
  shiny::moduleServer(id, function(input, output, session) {

    #create DT from selectedTrait
    output$DTSelectedTrait <- DT::renderDataTable(as.data.frame(DTSelectedTraits()),
                                                  options = list(pageLength = 1000, info = FALSE,
                                                                 lengthMenu = list(c(100, 1000, -1), c("100", "1000", "All"))))

    output$DTSelectedProbe <- DT::renderDataTable(as.data.frame(DTSelectedProbes()),
                                                  options = list(pageLength = 1000, info = FALSE,
                                                                 lengthMenu = list(c(100, 1000, -1), c("100", "1000", "All"))))
    output$SPLOM <- plotly::renderPlotly({SPLOM()})
    output$SPLOMTrait <- plotly::renderPlotly({SPLOMTrait()})
    output$SPLOMProbe <- plotly::renderPlotly({SPLOMProbe()})

    session$userData$sessionVariables$selectedAnnotation <- shiny::reactive({
      result <- NULL
      if (is.valid(session$userData$sessionVariables$selectedProbe())) {
        selectedAnnotation <- session$userData$annotation[session$userData$sessionVariables$selectedProbe(),]
        nprobes <- nrow(selectedAnnotation)
        selectedAnnotation$number <- seq(1:nprobes)
        selectedAnnotation$probeID <- selectedAnnotation$name
        col_order <- c("number", "probeID", "type", "target",	"name", "chromosome",	"position", "meth.dye", "gene.symbol", "gene.accession", "gene.region", "cpg.island.name", "relation.to.island", "snp.exclude", "450k", "common", "epic", "epic2")
        selectedAnnotation <- selectedAnnotation[, col_order]
        #add links to EWAS data hub
        selectedAnnotation <- addLinkToEWASDataHubShort(selectedAnnotation, session$userData$config$baseURL_EWASDataHub, session$userData$config$probeAttribut)
        selectedAnnotation <- addLinkToMRCEWASCatalogShort(selectedAnnotation, session$userData$config$baseURL_MRCEWASCatalog, session$userData$config$probeAttribut)

        selectedAnnotation <- addLinkToEWASDataHub(selectedAnnotation, session$userData$config$baseURL_EWASDataHub, session$userData$config$probeAttribut)
        selectedAnnotation <- addLinkToMRCEWASCatalog(selectedAnnotation, session$userData$config$baseURL_MRCEWASCatalog, session$userData$config$probeAttribut)
        #selectedAnnotation$probeID <- NULL
        #session$userData$sessionVariables$selectedAnnotation(selectedAnnotation)
        result <- selectedAnnotation
      }
      return(result)
    })

    DTSelectedProbes <- shiny::reactive({
      id <- shiny::showNotification("Creating table selected probes...", duration = NULL, closeButton = FALSE)
      on.exit(shiny::removeNotification(id), add = TRUE)
      base::tryCatch(
        {
          base::print(base::paste0(sysTimePID(), " start generating DTSelectedProbes"))
          result <- NULL
          #browser()
          # if (is.valid(session$userData$sessionVariables$selectedProbe())) {
          #   selectedAnnotation <- session$userData$annotation[session$userData$sessionVariables$selectedProbe(),]
          #   nprobes <- nrow(selectedAnnotation)
          #   selectedAnnotation$number <- seq(1:nprobes)
          #   selectedAnnotation$probeID <- selectedAnnotation$name
          #   col_order <- c("number", "probeID", "type", "target",	"name", "chromosome",	"position", "meth.dye", "gene.symbol", "gene.accession", "gene.region", "cpg.island.name", "relation.to.island", "snp.exclude", "450k", "common", "epic", "epic2")
          #   selectedAnnotation <- selectedAnnotation[, col_order]
          #   #add links to EWAS data hub
          #   selectedAnnotation <- addLinkToEWASDataHubShort(selectedAnnotation, session$userData$config$baseURL_EWASDataHub, session$userData$config$probeAttribut)
          #   selectedAnnotation <- addLinkToMRCEWASCatalogShort(selectedAnnotation, session$userData$config$baseURL_MRCEWASCatalog, session$userData$config$probeAttribut)
          #
          #   selectedAnnotation <- addLinkToEWASDataHub(selectedAnnotation, session$userData$config$baseURL_EWASDataHub, session$userData$config$probeAttribut)
          #   selectedAnnotation <- addLinkToMRCEWASCatalog(selectedAnnotation, session$userData$config$baseURL_MRCEWASCatalog, session$userData$config$probeAttribut)
          #   #selectedAnnotation$probeID <- NULL
          #   session$userData$sessionVariables$selectedAnnotation(selectedAnnotation)
          #   result <- selectedAnnotation
          # }
          result <- session$userData$sessionVariables$selectedAnnotation()
        },
        error = function(e) {
          base::message("An error occurred in shiny::reactive(DTSelectedProbes):\n", e)
        },
        warning = function(w) {
          base::message("A warning occurred in shiny::reactive(DTSelectedProbes):\n", w)
        },
        finally = {
          base::print(base::paste0(sysTimePID(), " finished generating DTSelectedProbes"))
          return(result)
        }
      )
    })

    DTSelectedTraits <- shiny::reactive({
      id <- shiny::showNotification("Creating table selected traits...", duration = NULL, closeButton = FALSE)
      on.exit(shiny::removeNotification(id), add = TRUE)
      base::tryCatch(
        {
          base::print(base::paste0(sysTimePID(), " start generating DTSelectedTraits"))
          result <- NULL
          if (is.valid(session$userData$sessionVariables$selectedTrait())) {
            result <- session$userData$sessionVariables$selectedTrait()
          }
        },
        error = function(e) {
          base::message("An error occurred in shiny::reactive(DTSelectedTraits):\n", e)
        },
        warning = function(w) {
          base::message("A warning occurred in shiny::reactive(DTSelectedTraits):\n", w)
        },
        finally = {
          base::print(base::paste0(sysTimePID(), " finished generating DTSelectedTraits"))
          return(result)
        }
      )
    })

    session$userData$sessionVariables$markingVar <- shiny::reactive({
      return(input$markingVar)
    })

    session$userData$sessionVariables$selectedOriginalData <- shiny::reactive({
      id <- shiny::showNotification("Getting selected original data...", duration = NULL, closeButton = FALSE)
      on.exit(shiny::removeNotification(id), add = TRUE)
      selectedOriginalData <- getSelectedOriginalData(session$userData$sessionVariables$probeReducedDataStructure()$combinedDFP_Val_Labels, session)

      if (!is.null(selectedOriginalData)) {
        FactorialVars <- getBinaryFactorialVars(selectedOriginalData)
        if (!is.valid(FactorialVars)) {
          FactorialVars <- NULL
        }
      } else {
        FactorialVars <- NULL
      }
      shiny::updateSelectizeInput(
        session = session,
        inputId = "markingVar",
        choices = FactorialVars,
        server = TRUE
      )
      return(selectedOriginalData)
    })

    session$userData$sessionVariables$OriginalDataTraits <- shiny::reactive({
      return(getSelectedOriginalDataTraits(session$userData$sessionVariables$probeReducedDataStructure()$combinedDFP_Val_Labels, selectedOnly = FALSE, session))
    })

    session$userData$sessionVariables$selectedOriginalDataTraits <- shiny::reactive({
      return(getSelectedOriginalDataTraits(session$userData$sessionVariables$probeReducedDataStructure()$combinedDFP_Val_Labels, selectedOnly = TRUE, session))
    })

    session$userData$sessionVariables$selectedOriginalDataProbes <- shiny::reactive({
      return(getSelectedOriginalDataProbes(session$userData$sessionVariables$probeReducedDataStructure()$combinedDFP_Val_Labels, selectedOnly = TRUE, session))
    })

    SPLOM <- shiny::reactive({
      id <- shiny::showNotification("Creating SPLOM...", duration = NULL, closeButton = FALSE)
      on.exit(shiny::removeNotification(id), add = TRUE)
      if (!is.null(session$userData$sessionVariables$selectedOriginalData())) {
        base::print(base::paste0(sysTimePID(), " number of traits and probes in SPLOM (columns in selectedDF): ", ncol(session$userData$sessionVariables$selectedOriginalData()))) #thats sum of probes and traits
        base::print(base::paste0(sysTimePID(), " number of cases in SPLOM (rows in selectedDF): ", nrow(session$userData$sessionVariables$selectedOriginalData()))) #thats number of cases in data set
        base::print(base::paste0(sysTimePID(), " number of traits in SPLOM (selectedTraits): ", ncol(session$userData$sessionVariables$selectedOriginalDataTraits())))
        base::print(base::paste0(sysTimePID(), " number of probes in SPLOM (selectedProbes): ", ncol(session$userData$sessionVariables$selectedOriginalDataProbes())))
        fig <- getSPLOM(session$userData$sessionVariables$selectedOriginalData(), XVars = colnames(session$userData$sessionVariables$selectedOriginalDataTraits()), YVars = colnames(session$userData$sessionVariables$selectedOriginalDataProbes()), markingVar = session$userData$sessionVariables$markingVar())
        return(fig)
      }
    })

    SPLOMTrait <- shiny::reactive({
      id <- shiny::showNotification("Creating SPLOM trait...", duration = NULL, closeButton = FALSE)
      on.exit(shiny::removeNotification(id), add = TRUE)
      if (!is.null(session$userData$sessionVariables$selectedOriginalDataTraits())) {
        fig <- getSPLOM(session$userData$sessionVariables$selectedOriginalDataTraits(), XVars = session$userData$sessionVariables$selectedOriginalDataTraits(), YVars = session$userData$sessionVariables$selectedOriginalDataTraits(), markingVar = session$userData$sessionVariables$markingVar())
        return(fig)
      }
    })

    SPLOMProbe <- shiny::reactive({
      id <- shiny::showNotification("Creating SPLOM probe...", duration = NULL, closeButton = FALSE)
      on.exit(shiny::removeNotification(id), add = TRUE)
      if (!is.null(session$userData$sessionVariables$selectedOriginalDataProbes())) {
        fig <- getSPLOM(session$userData$sessionVariables$selectedOriginalDataProbes(), XVars = session$userData$sessionVariables$selectedOriginalDataProbes(), YVars = session$userData$sessionVariables$selectedOriginalDataProbes(), markingVar = session$userData$sessionVariables$markingVar())
        return(fig)
      }
    })

    #' getSelectedOriginalData
    #' @param combinedDFP_Val_Labels list with datastructure pointing to original data -> session$userData$sessionVariables$probeReducedDataStructure()$combinedDFP_Val_Labels
    #' @param session shiny session object
    #' @return df with merged original data
    # examples getSelectedOriginalData(combinedDFP_Val_Labels, session)
    getSelectedOriginalData <- function(combinedDFP_Val_Labels, session) {
      id <- shiny::showNotification("Getting selected original data...", duration = NULL, closeButton = FALSE)
      on.exit(shiny::removeNotification(id), add = TRUE)
      base::print(base::paste0(sysTimePID(), " start getSelectedOriginalData()"))
      base::tryCatch(
        {
          selectedColnames <- session$userData$sessionVariables$selectedTrait()[,1]
          selectedColnames <- removeAdjFromColname(selectedColnames)
          selectedTraitSources <- session$userData$sessionVariables$selectedTrait()[,2]
          selectedColnamesTrait1 <- selectedColnames[selectedTraitSources == 1]
          selectedColnamesTrait2 <- selectedColnames[selectedTraitSources == 2]
          selectedColnamesTrait3 <- selectedColnames[selectedTraitSources == 3]
          #to be sure we select only colnames, which are within PHENODF:
          selectedColnamesTrait1 <- intersect(colnames(session$userData$sessionVariables$resultDFListTrait1()$listPHENOdata[[1]]$PHENODF), selectedColnamesTrait1)
          selectedColnamesTrait2 <- intersect(colnames(session$userData$sessionVariables$resultDFListTrait2()$listPHENOdata[[1]]$PHENODF), selectedColnamesTrait2)
          selectedColnamesTrait3 <- intersect(colnames(session$userData$sessionVariables$resultDFListTrait3()$listPHENOdata[[1]]$PHENODF), selectedColnamesTrait3)
          if (!is.valid(selectedColnamesTrait1)) {
            base::message(base::paste0(sysTimePID(), " file names in trait 1 folder do not match colnames in pheno file! SPLOM will not work."))
          }
          if (!is.valid(selectedColnamesTrait2)) {
            base::message(base::paste0(sysTimePID(), " file names in trait 2 folder do not match colnames in pheno file! SPLOM will not work."))
          }
          if (!is.valid(selectedColnamesTrait3)) {
            base::message(base::paste0(sysTimePID(), " file names in trait 3 folder do not match colnames in pheno file! SPLOM will not work."))
          }
          # get selected original data from trait data
          selectedDFTrait1 <- session$userData$sessionVariables$resultDFListTrait1()$listPHENOdata[[1]]$PHENODF[selectedColnamesTrait1]
          selectedDFTrait2 <- session$userData$sessionVariables$resultDFListTrait2()$listPHENOdata[[1]]$PHENODF[selectedColnamesTrait2]
          selectedDFTrait3 <- session$userData$sessionVariables$resultDFListTrait3()$listPHENOdata[[1]]$PHENODF[selectedColnamesTrait3]
          # merge all trait data together by Kind_ID or rowname (better)
          rn <- rownames(selectedDFTrait1)
          selectedDFTrait1$Row.names <- rn
          rn <- rownames(selectedDFTrait2)
          selectedDFTrait2$Row.names <- rn
          rn <- rownames(selectedDFTrait3)
          selectedDFTrait3$Row.names <- rn
          selectedDF <- NULL
          if (!base::is.null(selectedDFTrait1) && !base::is.null(selectedDFTrait2)) {
            selectedDF <- merge(selectedDFTrait1, selectedDFTrait2, by = "Row.names", all.x = FALSE, all.y = FALSE)
          }
          else {
            if (!base::is.null(selectedDFTrait2)) {
              selectedDF <- selectedDFTrait2
            }
          }
          if (!base::is.null(selectedDF) && !base::is.null(selectedDFTrait3)) {
            selectedDF <- merge(selectedDF, selectedDFTrait3, by = "Row.names", all.x = FALSE, all.y = FALSE)
          }
          else {
            if (!base::is.null(selectedDFTrait3)) {
              selectedDF <- selectedDFTrait3
            }
          }
          rownames(selectedDF) <- selectedDF$Row.names

          # row_index <- session$userData$sessionVariables$selectedProbe()
          # get selected methylation data...
          selectedRownames <- session$userData$sessionVariables$selectedProbe()
          #subset selectedRownames to only keep those, that are loaded in $Beta_tDF
          selectedRownames <- intersect(colnames(session$userData$Beta_tDF), selectedRownames)
          selectedBeta <- as.data.frame(session$userData$Beta_tDF[, selectedRownames]) #if error
          colnames(selectedBeta) <- selectedRownames
          rownames(selectedBeta) <- rownames(session$userData$Beta_tDF)
          #"nicht definierte Spalten gewählt" occurs, this is due to debug mode,
          #where most columns in Beta_tDF are not loaded.
          #... and merge with already merged trait data
          rn <- rownames(selectedBeta)
          if (is.valid(rn)) {
            selectedBeta$Row.names <- rn
            selectedDF_Beta <- merge(selectedDF, selectedBeta, by = "Row.names", all.x = FALSE, all.y = FALSE)
            #selectedDF_Beta <- merge(selectedBeta, selectedDF, by = "Row.names", all.x = FALSE, all.y = FALSE) #beta first, then traits
            rownames(selectedDF_Beta) <- selectedDF_Beta$Row.names
            selectedDF_Beta$Row.names <- NULL
          }
          else {
            message("We miss rownames in selectedDF_Beta here... (in getSelectedOriginalData()).\n
                Reason might be, that the beta data set was not loaded in full length (debugMode == TRUE?).\n")
          }
          if (nrow(selectedDF_Beta) > 256 || ncol(selectedDF_Beta) > 256) {
            base::message(base::paste0(sysTimePID(), " Warning: nrow(selectedDF) = ",
                                       nrow(selectedDF_Beta),
                                       " || ncol(selectedDF) = ",
                                       ncol(selectedDF_Beta),
                                       " that might be too much for fast processing!"))
          }
          # result gives 3d-data structure: multiple methylation profiles!!!
          # create SPLOM... each variable a dimension
          # data structure needed by plotly is described here: https://plotly.com/r/splom/

          # get involved original data sources
          # dataSources <- combinedDFP_Val_Labels$

          # iterate over those data sources
          # append data to result data frame

          # return df with merged original data from selected area
        },
        error = function(e) {
          base::message("An error occurred in getSelectedOriginalData():\n", e)
        },
        warning = function(w) {
          base::message("A warning occurred in getSelectedOriginalData():\n", w)
        },
        finally = {
          base::print(base::paste0(sysTimePID(), " end getSelectedOriginalData()"))
          return(selectedDF_Beta)
        }
      )
    }

    #' getSelectedOriginalDataTraits
    #' @param combinedDFP_Val_Labels list with datastructure pointing to original data -> session$userData$sessionVariables$probeReducedDataStructure()$combinedDFP_Val_Labels
    #' @param selectedOnly if TRUE onloy those traits are returned, that were previously selected by the user
    #' @param session shiny session object
    #' @return df with merged original data
    # examples getSelectedOriginalDataTraits(combinedDFP_Val_Labels, selectedOnly = TRUE, session)
    getSelectedOriginalDataTraits <- function(combinedDFP_Val_Labels, selectedOnly = TRUE, session) {
      id <- shiny::showNotification("Getting selected original data traits...", duration = NULL, closeButton = FALSE)
      on.exit(shiny::removeNotification(id), add = TRUE)
      base::print(base::paste0(sysTimePID(), " start getSelectedOriginalDataTraits()"))
      base::tryCatch(
        {
          if (selectedOnly == TRUE) {
            selectedColnames <- session$userData$sessionVariables$selectedTrait()[,1]
            selectedTraitSources <- session$userData$sessionVariables$selectedTrait()[,2]
            selectedColnames <- data.frame(selectedColnames,selectedTraitSources)
            #remove adj from names
            selectedColnames$selectedColnames <- removeAdjFromColname(selectedColnames$selectedColnames)
            selectedColnames <- unique(selectedColnames)
            selectedColnamesTrait1 <- selectedColnames$selectedColnames[selectedColnames$selectedTraitSources == 1]
            selectedColnamesTrait2 <- selectedColnames$selectedColnames[selectedColnames$selectedTraitSources == 2]
            selectedColnamesTrait3 <- selectedColnames$selectedColnames[selectedColnames$selectedTraitSources == 3]
            selectedTraitID <- session$userData$sessionVariables$selectedTraitID()
            #to be sure we select only colnames, which are within PHENODF:
            selectedColnamesTrait1 <- intersect(colnames(session$userData$sessionVariables$resultDFListTrait1()$listPHENOdata[[1]]$PHENODF), selectedColnamesTrait1)
            selectedColnamesTrait2 <- intersect(colnames(session$userData$sessionVariables$resultDFListTrait2()$listPHENOdata[[1]]$PHENODF), selectedColnamesTrait2)
            selectedColnamesTrait3 <- intersect(colnames(session$userData$sessionVariables$resultDFListTrait3()$listPHENOdata[[1]]$PHENODF), selectedColnamesTrait3)
            if (!is.valid(selectedColnamesTrait1)) {
              base::message(base::paste0(sysTimePID(), " file names in trait 1 folder do not match colnames in pheno file! SPLOM will not work."))
            }
            if (!is.valid(selectedColnamesTrait2)) {
              base::message(base::paste0(sysTimePID(), " file names in trait 2 folder do not match colnames in pheno file! SPLOM will not work."))
            }
            if (!is.valid(selectedColnamesTrait3)) {
              base::message(base::paste0(sysTimePID(), " file names in trait 3 folder do not match colnames in pheno file! SPLOM will not work."))
            }
            # get selected original data from trait data
            selectedDFTrait1 <- session$userData$sessionVariables$resultDFListTrait1()$listPHENOdata[[1]]$PHENODF[selectedColnamesTrait1]
            selectedDFTrait2 <- session$userData$sessionVariables$resultDFListTrait2()$listPHENOdata[[1]]$PHENODF[selectedColnamesTrait2]
            selectedDFTrait3 <- session$userData$sessionVariables$resultDFListTrait3()$listPHENOdata[[1]]$PHENODF[selectedColnamesTrait3]
          }
          else {
            # get original data from trait data
            selectedDFTrait1 <- session$userData$sessionVariables$resultDFListTrait1()$listPHENOdata[[1]]$PHENODF
            selectedDFTrait2 <- session$userData$sessionVariables$resultDFListTrait2()$listPHENOdata[[1]]$PHENODF
            selectedDFTrait3 <- session$userData$sessionVariables$resultDFListTrait3()$listPHENOdata[[1]]$PHENODF
          }
          #place suffix after all colnames
          if (is.valid(selectedColnamesTrait1)) {
            colnames(selectedDFTrait1) <- paste0(colnames(selectedDFTrait1),".1")
          }
          if (is.valid(selectedColnamesTrait2)) {
            colnames(selectedDFTrait2) <- paste0(colnames(selectedDFTrait2),".2")
          }
          if (is.valid(selectedColnamesTrait3)) {
            colnames(selectedDFTrait3) <- paste0(colnames(selectedDFTrait3),".3")
          }
          # merge all trait data together by Kind_ID or rowname (better)
          if (is.valid(selectedColnamesTrait1)) {
            rn <- rownames(selectedDFTrait1)
            selectedDFTrait1$Row.names <- rn
          }
          if (is.valid(selectedColnamesTrait2)) {
            rn <- rownames(selectedDFTrait2)
            selectedDFTrait2$Row.names <- rn
          }
          if (is.valid(selectedColnamesTrait3)) {
            rn <- rownames(selectedDFTrait3)
            selectedDFTrait3$Row.names <- rn
          }
          selectedDF <- NULL
          if (is.valid(selectedColnamesTrait1) && is.valid(selectedColnamesTrait2)) {
            selectedDF <- merge(selectedDFTrait1, selectedDFTrait2, by = "Row.names", all.x = FALSE, all.y = FALSE)
          }
          else {
            if (is.valid(selectedColnamesTrait2)) {
              selectedDF <- selectedDFTrait2
            }
            else {
              selectedDF <- selectedDFTrait1
            }
          }
          if (is.valid(selectedDF) && is.valid(selectedColnamesTrait3)) {
            selectedDF <- merge(selectedDF, selectedDFTrait3, by = "Row.names", all.x = FALSE, all.y = FALSE)
          }
          else {
            if (is.valid(selectedColnamesTrait3)) {
              selectedDF <- selectedDFTrait3
            }
            else {
              selectedDF <- NULL
              browser() #should not happen
            }
          }
          rownames(selectedDF) <- selectedDF$Row.names

          if (length(selectedColnames) > 256) {
            base::message(base::paste0(sysTimePID(), "length(selectedColnames) = ",
                                       length(selectedColnames),
                                       " that might be too much for fast processing!"))
          }
          #remove row.names from selectedDF
          if("Row.names" %in% colnames(selectedDF)) {
            selectedDF$Row.names <- NULL
          }
          if("Row.names.2" %in% colnames(selectedDF)) {
            selectedDF$Row.names.2 <- NULL
          }
          if("Row.names.3" %in% colnames(selectedDF)) {
            selectedDF$Row.names.3 <- NULL
          }
        },
        error = function(e) {
          base::message("An error occurred in getSelectedOriginalDataTraits():\n", e)
        },
        warning = function(w) {
          base::message("A warning occurred in getSelectedOriginalDataTraits():\n", w)
        },
        finally = {
          base::print(base::paste0(sysTimePID(), " end getSelectedOriginalDataTraits()"))
          return(selectedDF)
        }
      )
    }

    #' getSelectedOriginalDataProbes
    #' @param combinedDFP_Val_Labels list with datastructure pointing to original data -> session$userData$sessionVariables$probeReducedDataStructure()$combinedDFP_Val_Labels
    #' @param session shiny session object
    #' @return df with merged original data
    # examples getSelectedOriginalDataProbes(combinedDFP_Val_Labels, session)
    getSelectedOriginalDataProbes <- function(combinedDFP_Val_Labels, selectedOnly, session) {
      id <- shiny::showNotification("Getting selected original data probes...", duration = NULL, closeButton = FALSE)
      on.exit(shiny::removeNotification(id), add = TRUE)
      base::print(base::paste0(sysTimePID(), " start getSelectedOriginalDataProbes()"))
      base::tryCatch(
        {
          if (selectedOnly == TRUE) {
            row_index <- session$userData$sessionVariables$selectedProbe()
            # get selected methylation data...
            rowInd <- which(rownames(combinedDFP_Val_Labels$dfP_Val) %in% row_index)
            selectedRownames <- rownames(combinedDFP_Val_Labels$dfP_Val)[rowInd]
            #subset selectedRownames to only keep those, that are loaded in $Beta_tDF
            selectedRownames <- intersect(colnames(session$userData$Beta_tDF), selectedRownames)
            selectedBeta <- session$userData$Beta_tDF[, selectedRownames] #if error
          }
          else {
            selectedBeta <- session$userData$Beta_tDF
          }
          #"nicht definierte Spalten gewählt" occurs, this is due to debug mode,
          #where most columns in Beta_tDF are not loaded.
          #... and merge with trait data from markingVar
          #select only ID# and markingVar
          # traits$id <- rownames(traits)
          # Vars <- c("id", markingVar)
          # if (all(Vars %in% colnames(traits))) {
          #   traits <- traits[, Vars]
          #   traits$markingVar <- traits[, markingVar]
          #   traits[, markingVar] <- NULL
          #   #merge
          #   selectedBeta$id <- rownames(selectedBeta)
          #   selectedBeta <- merge(selectedBeta, traits, by.x = "id", by.y = "id")
          #   rownames(selectedBeta) <- selectedBeta$id
          #   selectedBeta$id <- NULL
          # }
          rn <- rownames(selectedBeta)
          if (!is.valid(rn)) {
            message("We miss rownames in selectedDF here... (in getSelectedOriginalDataProbes()).\n
                Reason might be, that the beta data set was not loaded in full length (debugMode == TRUE?).\n")
          }
          if (length(selectedRownames) > 256) {
            base::message(base::paste0(sysTimePID(), "length(selectedRownames) = ",
                                       length(selectedRownames),
                                       " that might be too much for fast processing!"))
          }
          #remove row.names from selectedDF
          if ("Row.names" %in% colnames(selectedBeta)) {
            selectedBeta$Row.names <- NULL
          }
        },
        error = function(e) {
          base::message("An error occurred in getSelectedOriginalDataProbes():\n", e)
        },
        warning = function(w) {
          base::message("A warning occurred in getSelectedOriginalDataProbes():\n", w)
        },
        finally = {
          base::print(base::paste0(sysTimePID(), " end getSelectedOriginalDataProbes()"))
          return(selectedBeta)
        }
      )
    }
  }) #end shiny::moduleServer
}
