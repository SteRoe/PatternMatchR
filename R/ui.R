ui <- shiny::shinyUI(
  shiny::fluidPage(
    # Some custom CSS for a smaller font for preformatted text
    shiny::tags$head(
      shiny::tags$style(shiny::HTML("      pre, table.table {        font-size: smaller;      }    "))
    ),
    shiny::fluidRow(
      shiny::textInput("inpWorkDir", "Work Directory", "", placeholder = TRUE),
      shiny::verbatimTextOutput("txtDebugOut", placeholder=TRUE),
    ),
    shiny::tabsetPanel(
      shiny::tabPanel("HM",
        shinyBS::bsCollapse(id = "collapse", open = c("files", "#probes", "borders", "HM"), multiple = TRUE,
          shinyBS::bsCollapsePanel("# cases",
          shiny::fluidRow(
             shiny::textInput("txtCases", "minimum # cases for each trait", 100, placeholder = TRUE),
            )
          ),
          shinyBS::bsCollapsePanel("folders", #label = "P_VAL border vs. # of probes",
            shiny::fluidRow(
              shiny::column(width = 4,
                DT::dataTableOutput("trait1DirList"),
                shiny::actionButton("btnLoadDir1", label = "Load trait 1 dir"),
              ),
              shiny::column(width = 4,
                 DT::dataTableOutput("trait2DirList"),
                 shiny::actionButton("btnLoadDir2", label = "Load trait 2 dir"),
              ),
              shiny::column(width = 4,
#                shiny::textInput("txtInDir3", "trait 3 directory", dataDir3, placeholder = TRUE),
                DT::dataTableOutput("trait3DirList"),
                shiny::actionButton("btnLoadDir3", label = "Load trait 3 dir"),
              )
            )
          ),
          shinyBS::bsCollapsePanel("#probes", #label = "P_VAL border vs. # of probes",
            shiny::fluidRow(
              shiny::actionButton("btnCountProbes", label = "Count Probes"),
              shiny::actionButton("btnBrowser", label = "Break to Browser()")
            ),
            shiny::fluidRow(
              shiny::column(width = 5,
                DT::dataTableOutput("DTP_VALborder")
              )
            )
          ),
          shinyBS::bsCollapsePanel("borders", #label = "Border Settings",
            shiny::fluidRow(
              shiny::column(width = 2,
#                shiny::textInput("txtP_VALborder", "P_VAL Border", 0.000000000000000005, placeholder = TRUE),
#                shiny::textInput("txtP_VALborder", "P_VAL Border 5*10^-", 18, placeholder = TRUE),
                shiny::textInput("txtMaxProbes", "maximum Probes", 32000, placeholder = TRUE), #sqrt(2^31)=46340
              ),
              shiny::column(width = 2,
                shiny::textInput("txtMaxClassesProbes", "maximum colored Classes Probes (in dendrogram)", 10, placeholder = TRUE),
                shiny::textInput("txtMaxClassesTraits", "maximum colored Classes Traits (in dendrogram)", 10, placeholder = TRUE)
              ),
              shiny::column(width = 4,
                shiny::tags$html("n"),
                shiny::verbatimTextOutput("txtResultingN", placeholder = TRUE),
                shiny::tags$html("minimum p-value"),
                shiny::verbatimTextOutput("txtMinP_Val", placeholder = TRUE)
              )
            )
          )
        ),
        shinyBS::bsCollapsePanel("HM", #label = "Border Settings",
          shiny::fluidRow(
            shiny::column(width = 4,
              shiny::actionButton("plotCombinedHM", label = "Make HM/ Dendrograms/ Tables")
            )
          ),
          shiny::fluidRow(
            shiny::column(width = 8,
              shiny::verbatimTextOutput("txtHMDescription", placeholder = TRUE)
            )
          ),
          shiny::fluidRow(
            shiny::tabsetPanel(
            # tabPanel("Full",
            #   tabsetPanel(
                shiny::tabPanel("HM",
                  shiny::fluidRow(
                    shiny::column(width = 6,
                      shiny::sliderInput("sldMinP_Val", "minimum p-val for colour scale, 5e-x", 0, 200, step = -1, 200, 1)
                    ),
                    shiny::column(width = 6,
                      shiny::sliderInput("sldMaxP_Val", "maximum p-val for colour scale, 5e-x", 0, 200, step = -1, 2, 1)
                    )
                  ),
                  shiny::fluidRow(
                    shiny::column(width = 12,
                      InteractiveComplexHeatmap::InteractiveComplexHeatmapOutput("heatmap_1",
                                                       height1 = 2000, width1 = 1500, height2 = 2000, width2 = 1500, inline = F),
                    )
                  )
                ),
                shiny::tabPanel("Table P_VAL",
                         #DT tableout
                         DT::dataTableOutput("DTP_VAL")
                ),
                shiny::tabPanel("Table Delta Methylation",
                         #DT tableout
                         DT::dataTableOutput("DTDM")
                ),
                shiny::tabPanel("Table N",
                         #DT tableout
                         DT::dataTableOutput("DTN")
                ),
                # tabPanel("Normalized and Filtered P_VAL",
                #          #DT tableout
                #          DT::dataTableOutput("DTNormP_VAL")
                # ),
                shiny::tabPanel("Dendrogram Probes",
    #              column(width = 12, plotly::plotlyOutput("plotDendrogramProbes"))
                  plotly::plotlyOutput("plotDendrogramProbes", height = "800px")
                ),
                shiny::tabPanel("Annotated Table Probes",
                         DT::dataTableOutput("DTProbes")
                ),
                shiny::tabPanel("Dendrogram Traits",
    #              column(width = 12, plotly::plotlyOutput("plotDendrogramTraits"))
                   plotly::plotlyOutput("plotDendrogramTraits", height = "80%")
                ),
                shiny::tabPanel("Table Traits",
                         DT::dataTableOutput("DTTraits")
                )
              )
          )
        ),
      ),
      shiny::tabPanel("Settings",
        shiny::fluidRow(
#           shiny::textInput("inpTraitsFile", "File with Traits", "FileName", placeholder = TRUE),
           shiny::textInput("inpDNAmFile", "File with DNAm", "FileName", placeholder = TRUE),
        )
      )
    )
  )
)
