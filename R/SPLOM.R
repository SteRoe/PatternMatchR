#' structure of sessionVariables for SPLOM:
#' detailed original data for SPLOM analysis
#' session$userData$sessionVariables$SPLOM - data structure for scatter plot matrix
#' session$userData$sessionVariables$selectedOriginalData - merged original data for creation of SPLOM
#' session$userData$sessionVariables$markingVar - variable, which should be marked (stratified for) in SPLOM

#' getDimensionsForPlotlySPLOM
#' @param df data frame
#' @param omitVar Variable to leave out of nested list
#' @return nested list for plotly with variable labels and formulae
#' examples getDimensionsForPlotlySPLOM(df, omitVar)
getDimensionsForPlotlySPLOM <- function(df, omitVar) {
  dimensionList <- list()
  for (i in seq_along(df)) {
    if (is.valid(omitVar)) {
      if (colnames(df[i]) != omitVar) { # omit grouping variabe for plotly
        l <- list(label = colnames(df[i]), values = as.formula(paste("~", colnames(df)[i]))) # ~ as formula operator for plotly
        dimensionList[[i]] <- l
      }
    } else {
      l <- list(label = colnames(df[i]), values = as.formula(paste("~", colnames(df)[i]))) # ~ as formula operator for plotly
      dimensionList[[i]] <- l
    }
  }
  return(dimensionList)
}

#' getBinaryFactorialVars
#' @param df data frame
#' @return list with variables that can be treated as factor in df
# examples getBinaryFactorialVars(df)
getBinaryFactorialVars <- function(df) {
  shinyId <- shiny::showNotification("Getting binary factor variables...", duration = NULL, closeButton = FALSE)
  on.exit(shiny::removeNotification(shinyId), add = TRUE)
  base::print(base::paste0(sysTimePID(), " start getBinaryFactorialVars()."))
  result <- list()
  for (i in base::seq_along(df)) {
    if (is.integer(df[[i]])) {
      if (length(levels(as.factor(df[[i]]))) == 2) {
        result <- c(result, colnames(df)[[i]])
      }
    }
  }
  if (length(result) == 0) result <- NULL
  base::print(base::paste0(sysTimePID(), " found ", length(result), " binary factorial vars."))
  base::print(base::paste0(sysTimePID(), " end getBinaryFactorialVars()."))
  return(result)
}

#' getSPLOM
#' @param df data frame, from which the SPLOM will be drawn
#' @param markingVar variable inside the data frame, which will be not drawn in SPLOM, instead it will be used for visual stratification
#' @return a plotly figure, which can be drawn using renderPlotly
# examples getSPLOM(df, markingVar)
getSPLOM <- function(df, XVars, YVars, markingVar) {
  shinyId <- shiny::showNotification("Getting SPLOM...", duration = NULL, closeButton = FALSE)
  on.exit(shiny::removeNotification(shinyId), add = TRUE)
  if (!is.null(df)) {
    dimensions <- getDimensionsForPlotlySPLOM(df, markingVar)
    pl_colorscale <- list(
      c(0.0, "#119dff"),
      c(0.5, "#119dff"),
      c(0.5, "#ef553b"),
      c(1, "#ef553b")
    )
    #25900x3200 is possible for dimensions
    # if (is.numeric(height)) {
    #   h <- base::paste0(height, "px")
    # }
    # if (is.numeric(width)) {
    #   w <- base::paste0(width, "px")
    # }
    if ("markingVar" %in% colnames(df)) {
      df[, "markingVar"] <- as.factor(df[, "markingVar"])
    }
    fig <- df %>%
      #plotly::plot_ly(height = height, width = width)
      plotly::plot_ly()
    if ("markingVar" %in% colnames(df)) {
      fig <- fig %>%
        plotly::add_trace(
          type = "splom",
          dimensions = dimensions,
          text = as.formula("~factor(markingVar, labels=c(\"control\",\"case\"))"),
          diagonal = list(visible = FALSE),
          marker = list(
            color = as.formula(paste0("~markingVar")),
            colorscale = pl_colorscale,
            size = 5,
            line = list(
              width = 1,
              color = "rgb(230,230,230)"
            )
          ),
          xaxes = XVars
          #yaxes = YVars
          # xaxis = list(title = "trait", tickangle = 90),
          # yaxis = list(title = "probe", tickangle = 90)
        )
    }
    else if (is.valid(markingVar) && any(markingVar %in% colnames(df))) {
      fig <- fig %>%
        plotly::add_trace(
        type = "splom",
        dimensions = dimensions,
        text = as.formula(paste0("~factor(", markingVar, ", labels=c(\"control\",\"case\"))")),
        diagonal = list(visible = FALSE),
        marker = list(
          color = as.formula(paste0("~", markingVar)),
          colorscale = pl_colorscale,
          size = 5,
          line = list(
            width = 1,
            color = "rgb(230,230,230)"
          )
        ),
        xaxes = XVars
        #yaxes = YVars
        # xaxis = list(title = "trait", tickangle = 90),
        # yaxis = list(title = "trait", tickangle = 90)
      )
    }
    else {
      fig <- fig %>%
        plotly::add_trace(
        type = "splom",
        dimensions = dimensions,
        diagonal = list(visible = FALSE),
        marker = list(
          colorscale = pl_colorscale,
          size = 5,
          line = list(
            width = 1,
            color = "rgb(230,230,230)"
          )
        ),
        xaxes = XVars
        #yaxes = YVars
        # xaxis = list(title = "probe", tickangle = 90),
        # yaxis = list(title = "probe", tickangle = 90)
      )
    }
  } else {
    fig <- NULL
  }
  return(fig)
}
