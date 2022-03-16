server <- function(input, output, session) {
###because <InteractiveComplexHeatmap> does not support modules, we have to do this here outside the heatmap module plotHM.R
###nevertheless, we refer to the functions defined in plotHM.R
# define globally used sessionVariables here
  globalVariables <- list()
  print(paste0(Sys.time(), " loading configuration."))
  globalVariables$config <- config::get(file = "config.yml")
  globalVariables = loadObjects(globalVariables)
  print(paste0(Sys.time(), " defining session variables."))
#  sessionVariables <- shiny::reactiveValues(trait = list(), probe = list(), df_plot = data.frame(), resultDataSingleTrait = data.frame(), resultDataTrait1 = data.frame(), resultDataTrait2 = data.frame(), resultDataTrait3 = data.frame(), gP_ValMaxBorder = double(), gP_ValMinBorder = double(), gMaxProbes = integer())
  sessionVariables <- shiny::reactiveValues(gP_ValMaxBorder = double(), gP_ValMinBorder = double(), gMaxProbes = integer(), numberVariables = integer())
  print(paste0(Sys.time(), " starting application."))
#read config
  output$inpWorkDir <- shiny::renderText(globalVariables$config$workDir)
  dfdD1 = data.table::as.data.table(unlist(globalVariables$config$dataDir1))
  dfdD2 = data.table::as.data.table(unlist(globalVariables$config$dataDir2))
  dfdD3 = data.table::as.data.table(unlist(globalVariables$config$dataDir3))
  #output$trait1DirList = DT::renderDataTable(dfdD, editable = 'cell')
  output$trait1DirList = DT::renderDataTable(dfdD1)
  output$trait2DirList = DT::renderDataTable(dfdD2)
  output$trait3DirList = DT::renderDataTable(dfdD3)
  shiny::observeEvent(input$btnBrowser, ignoreInit = TRUE, {
    browser()
  }, ignoreNULL = FALSE)
  shiny::observeEvent(input$btnLoadDir1, ignoreInit = TRUE, {
    print(paste0(Sys.time(), " start loading trait 1 folders."))
#browser()
    if (is.numeric(input$trait1DirList_rows_selected)) {
      trait1DirList = as.list(dfdD1[input$trait1DirList_rows_selected,][[1]])
      print(paste0(Sys.time(), " before merge folders1()"))
      traitDFs = lapply(trait1DirList, loadResultDF)

      resultDFList = loadtraitDFs(traitDFs)
      resultDFP_Val = resultDFList$resultDFP_Val
      resultDFDM = resultDFList$resultDFDM
      resultDFN = resultDFList$resultDFN
#browser()
      sessionVariables$resultP_ValTrait1 = resultDFP_Val #traitDF
      sessionVariables$resultDMTrait1 = resultDFDM
      sessionVariables$resultNTrait1 = resultDFN
      sessionVariables$numberVariables = length(sessionVariables$resultP_ValTrait1) + length(sessionVariables$resultP_ValTrait2) + length(sessionVariables$resultP_ValTrait3)
      print(paste0(Sys.time(), " '", trait1DirList, "' successfully loaded."))
    }
  }, ignoreNULL = FALSE)

  shiny::observeEvent(input$btnLoadDir2, ignoreInit = TRUE, {
    print(paste0(Sys.time(), " start loading trait 2 folders."))
    if (is.numeric(input$trait2DirList_rows_selected)) {
      trait2DirList = as.list(dfdD2[input$trait2DirList_rows_selected,][[1]])
      print(paste0(Sys.time(), " before merge folders2()"))
      traitDFs = lapply(trait2DirList, loadResultDF)

      resultDFList = loadtraitDFs(traitDFs)
      resultDFP_Val = resultDFList$resultDFP_Val
      resultDFDM = resultDFList$resultDFDM
      resultDFN = resultDFList$resultDFN

      sessionVariables$resultP_ValTrait2 = resultDFP_Val #traitDF
      sessionVariables$resultDMTrait2 = resultDFDM
      sessionVariables$resultNTrait2 = resultDFN
      sessionVariables$numberVariables = length(sessionVariables$resultP_ValTrait1) + length(sessionVariables$resultP_ValTrait2) + length(sessionVariables$resultP_ValTrait3)
      print(paste0(Sys.time(), " '", trait2DirList, "' successfully loaded."))
    }
  }, ignoreNULL = FALSE)

  shiny::observeEvent(input$btnLoadDir3, ignoreInit = TRUE, {
    print(paste0(Sys.time(), " start loading trait3 folders."))
    if (is.numeric(input$trait3DirList_rows_selected)) {
      trait3DirList = as.list(dfdD3[input$trait3DirList_rows_selected,][[1]])
      print(paste0(Sys.time(), " before merge folders3()"))
      traitDFs = lapply(trait3DirList, loadResultDF)

      resultDFList = loadtraitDFs(traitDFs)
      resultDFP_Val = resultDFList$resultDFP_Val
      resultDFDM = resultDFList$resultDFDM
      resultDFN = resultDFList$resultDFN

      sessionVariables$resultP_ValTrait3 = resultDFP_Val #traitDF
      sessionVariables$resultDMTrait3 = resultDFDM
      sessionVariables$resultNTrait3 = resultDFN
      sessionVariables$numberVariables = length(sessionVariables$resultP_ValTrait1) + length(sessionVariables$resultP_ValTrait2) + length(sessionVariables$resultP_ValTrait3)
      print(paste0(Sys.time(), " '", trait3DirList, "' successfully loaded."))
    }
    else {
      message(paste0(Sys.time(), " no entries selected from trait3 folders."))
    }
  }, ignoreNULL = FALSE)

  shiny::observeEvent(input$btnCountProbes, ignoreInit = TRUE, {
    print(paste0(Sys.time(), " start counting probes."))
    minN = as.integer(input$txtCases)
    combinedDFP_Val_Labels = getCombinedDFP_Val_Labels(sessionVariables, minN)
    P_VALNTable = getAvailNForP_VALBorder(combinedDFP_Val_Labels$dfP_Val)
    output$DTP_VALborder = DT::renderDataTable(P_VALNTable)
    print(paste0(Sys.time(), " finished counting probes."))
  }, ignoreNULL = FALSE)

  output$txtMinP_Val = shiny::reactive({
    a = sessionVariables$numberVariables
    combinedDFP_Val_Labels = getCombinedDFP_Val_Labels(sessionVariables, minN)
    dfP_Val = combinedDFP_Val_Labels$dfP_Val
    minP_Val = min(dfP_Val[which(dfP_Val>0, arr.ind = TRUE)])
    exponent = extractMantissaExponent(minP_Val)$exponent
    if (exponent<0) {exponent = exponent * -1}
    if (is.finite(exponent)) {
      if (exponent == 0) {exponent = 200}
      shiny::updateSliderInput(session, "sldMaxP_Val", max = exponent, value = 2)
      shiny::updateSliderInput(session, "sldMinP_Val", max = exponent, value = exponent)
    }
    return(minP_Val)
  })

  shiny::observeEvent(input$plotCombinedHM, ignoreInit = TRUE, {

    sessionVariables$gMaxProbes = input$txtMaxProbes
    print(paste0(Sys.time(), " start plotting heatmap."))
    output$txtHMDescription = shiny::renderText(paste0("calculating heatmap..., current plot is not valid"))
    while (!is.null(grDevices::dev.list())) grDevices::dev.off()
    print(paste0(Sys.time(), " creating empty heatmap."))
    combinedHMP_VAL = emptyHM()
    InteractiveComplexHeatmap::makeInteractiveComplexHeatmap(input, output, session, combinedHMP_VAL, "heatmap_1")
    minN = as.integer(input$txtCases)
    fileNameCombinedHM = paste0("CombinedHM.RDS")
    if (utils::file_test("-f", paste0(globalVariables$config$workDir,fileNameCombinedHM)) == TRUE) {
      result = readRDS (file=paste0(globalVariables$config$workDir,fileNameCombinedHM))
      combinedDFP_Val_Labels = result$combinedDFP_Val_Labels
      distMatProbes = result$distMatProbes
      clustResProbes = result$clustResProbes
      distMatTraits = result$distMatTraits
      clustResTraits = result$clustResTraits
      dfP_Val = combinedDFP_Val_Labels$dfP_Val
      dfDM = combinedDFP_Val_Labels$dfDM
      dfN = combinedDFP_Val_Labels$dfN
    }
    else {
      print(paste0(Sys.time(), " calculating combined heatmap."))
      combinedDFP_Val_Labels = getCombinedDFP_Val_Labels(sessionVariables, minN)
      dfP_Val = combinedDFP_Val_Labels$dfP_Val
      dfDM = combinedDFP_Val_Labels$dfDM
      dfN = combinedDFP_Val_Labels$dfN
print(class(dfP_Val))
      print(paste0(Sys.time()," nrow before reduce ",nrow(dfP_Val)))
#check numRows...
      if (sessionVariables$gMaxProbes>32000) {
        #in theory, ComplexHeatmap is limited to sqrt(2^31) = 46340 by 46340 items;
        #in reality, limitation is around 30000 items (sometimes more, unfortunately no constant value; checked this in August 2021); with matrices larger than that,
        #Error in Cairo: Failed to create Cairo backend!
        message(paste0("probelist for heatmap has ", sessionVariables$gMaxProbes, " probes. This is more than 32000 probes, which is due to the sqr(2^31) limit for ComplexHeatmap."))
        sessionVariables$gMaxProbes=32000
      }
      dfP_Val = getReducedP_Valdf(dfP_Val, numRows = as.integer(sessionVariables$gMaxProbes), upperP_VALborder = as.numeric(sessionVariables$gP_ValMaxBorder), lowerP_VALborder = as.numeric(sessionVariables$gP_ValMinBorder))
print(class(dfP_Val))
      nRowReduced = nrow(dfP_Val)
      print(paste0(Sys.time()," nrow after reduce ",nrow(dfP_Val)))
#      dfP_Val[which(is.na(dfP_Val))] = 1 #set missing P_VAL to 1
      print(paste0(Sys.time(), " set missing p-values to 1."))
      dfP_Val[is.na(dfP_Val)] = 1 #set missing P_VAL to 1
print(class(dfP_Val))
      output$txtResultingN = shiny::renderText(paste0("distance matrix will be computed with n probes: ", nrow(dfP_Val)))
      print(Cstack_info())
      if (nrow(dfP_Val)>=5) {
        if (globalVariables$config$debugMode == TRUE && nrow(dfP_Val) > 1000) {
          print(paste0(Sys.time(), " debug mode n probes=1000."))
          dfP_Val = dfP_Val[1:1000,]
        }

        print(paste0(Sys.time(), " set infinite p-values to .Machine$integer.max."))
        indexdfP_Val_infinite = is.infinite(as.matrix(dfP_Val))
        dfP_Val[indexdfP_Val_infinite] = .Machine$integer.max
        if (!is.data.frame(dfDM)) {
          dfDM = as.data.frame(dfDM)
        }
        print(paste0(Sys.time(), " shortening dfDM"))
        dfDM = dfDM[rownames(dfP_Val),colnames(dfP_Val)]
        if (!is.data.frame(dfN)) {
          dfN = as.data.frame(dfN)
        }
        print(paste0(Sys.time(), " shortening dfN"))
        dfN = dfN[rownames(dfP_Val),colnames(dfP_Val)]
        combinedDFP_Val_Labels$dfP_Val = dfP_Val
        combinedDFP_Val_Labels$dfDM = dfDM
        combinedDFP_Val_Labels$dfN = dfN
        #clustering for rows
        distMatProbes = stats::as.dist(matrix(nrow = nrow(dfP_Val), ncol = nrow(dfP_Val)))
        print(paste0(Sys.time()," before distance matrix for probes ",nrow(dfP_Val)))
        if (!is.null(globalVariables$numCores)) {
          nC = globalVariables$numCores - 10
        }
        else {
          nC = 1
        }
        distMatProbes <- parallelDist::parallelDist(as.matrix(dfP_Val), method = "euclidean", labels = TRUE, threads = nC)
        print(paste0(Sys.time()," before hierarchical clustering for probes ",nrow(dfP_Val)))
        print(paste0(Sys.time(), " mode distMatProbes: ", mode(distMatProbes), "."))
        print(paste0(Sys.time(), " NA in distMatProbes: ", which(is.na(distMatProbes), arr.ind = TRUE),".")) #identify NA if given

        #clustResProbes <- fastcluster::hclust(distMatProbes, method = "complete")
        clustResProbes <- fastcluster::hclust(distMatProbes, method = "ward.D2")

        #clustering for cols
        matP_Val.t = base::t(as.matrix(dfP_Val))
        distMatTraits = stats::as.dist(matrix(nrow = nrow(matP_Val.t), ncol = nrow(matP_Val.t)))
        #    dist.mat <- philentropy::distance(mat.t, method = "euclidean")
        print(paste0(Sys.time()," before distance matrix for traits ",nrow(matP_Val.t)))
        distMatTraits <- parallelDist::parallelDist(matP_Val.t, method = "euclidean", labels = TRUE, threads = nC)
        print(paste0(Sys.time()," before hierarchical clustering for traits ",nrow(matP_Val.t)))
  #      clustResTraits <- fastcluster::hclust(distMatTraits, method = "complete")
        clustResTraits <- fastcluster::hclust(distMatTraits, method = "ward.D2")

        result = list()
        result$combinedDFP_Val_Labels = combinedDFP_Val_Labels
        result$distMatProbes = distMatProbes
        result$clustResProbes = clustResProbes
        result$distMatTraits = distMatTraits
        result$clustResTraits = clustResTraits
  #      saveRDS(result, file = paste0(workDir,fileNameCombinedHM))
      }
      else {
        message(paste0(Sys.time()," less than 5 probes remained."))
      }
    }
    if (nrow(dfP_Val)>5) {
      start_time <- Sys.time()
      output$txtResultingN = shiny::renderText(paste0("number of resulting probes: ",nrow(dfP_Val)))
      print(paste0(Sys.time()," gc()"))
      gc()
      print(paste0(Sys.time()," before making dendrogram for probes"))

#check clustResProbes > 8
length(clustResProbes)
  #even with:
  #RStudio\bin\rstudio.exe /high --max-ppsize 500000
  #Error in : protect(): protection stack overflow
      options(expression = 500000)
      dendProbes = stats::as.dendrogram(clustResProbes)
      dendProbes = dendextend::color_branches(dendProbes, input$txtMaxClassesProbes)
      listProbes = labels(dendProbes)
      print(paste0(Sys.time()," before making dendrogram for traits"))
      dendTraits = stats::as.dendrogram(clustResTraits)
      dendTraits = dendextend::color_branches(dendTraits, input$txtMaxClassesTraits)
      listTraits = labels(dendTraits)
      print(paste0(Sys.time()," before reordering matrices"))
      # reorder matrices
      dfP_Val = dfP_Val[listProbes,listTraits]
      dfDM = dfDM[listProbes,listTraits]
      dfN = dfN[listProbes,listTraits]

      print(paste0(Sys.time()," before plotting trait dendrogram"))
      plotTraits <- ggplot2::ggplot(dendTraits, horiz = TRUE, offset_labels = -1)
      print(paste0(Sys.time()," before rendering trait dendrogram"))
      output$plotDendrogramTraits <- plotly::renderPlotly(plotTraits)
      #build annotated tables for probes in order of dendrogram
      print(paste0(Sys.time()," before rendering dendrogram tables probes"))
      DTProbes = data.frame(row.names = 1:length(listProbes))
      DTProbes$probeID = listProbes
      DTProbes$order = seq(nrow(DTProbes))
      rownames(DTProbes) = DTProbes$probeID
      #add annotation
      DTProbes = base::merge(x = DTProbes, y = globalVariables$annotation, by.x = "probeID", by.y = "name", all.x = TRUE, all.y = FALSE)
      #sort
      DTProbes = DTProbes[order(DTProbes$order),]
      rownames(DTProbes) = DTProbes$probeID
      output$DTProbes <- DT::renderDataTable(DTProbes)
      print(paste0(Sys.time()," before rendering dendrogram tables traits"))
      DTTraits = data.frame(row.names = 1:length(listTraits))
      DTTraits$Name = listTraits
      DTTraits$order = seq(nrow(DTTraits))
      rownames(DTTraits) = DTTraits$Name
      DTTraits = DTTraits[order(DTTraits$order),]
      output$DTTraits <- DT::renderDataTable(DTTraits)
      print(paste0(Sys.time()," before rendering table P_Val"))
      output$DTP_VAL = DT::renderDataTable(dfP_Val)

      print(paste0(Sys.time()," before rendering table DeltaMeth"))
      output$DTDM = DT::renderDataTable(dfDM)

      print(paste0(Sys.time()," before rendering table N"))
      output$DTN = DT::renderDataTable(dfN)

      print(paste0(Sys.time()," before calculating heatmap"))
      #check dimensions of combinedDFP_Val_Labels, dendProbes, dendTraits
      l = combinedDFInteractiveHeatMapP_Val(combinedDFP_Val_Labels, dendProbes, dendTraits)
      print(paste0(Sys.time()," before plotting heatmap"))

      while (!is.null(grDevices::dev.list())) grDevices::dev.off()
      combinedHMP_VAL = l$combinedHMP_VAL
      subHM_layer_fun = l$layer_fun
      end_time <- Sys.time()
      elapsed_time <- end_time - start_time
#      InteractiveComplexHeatmap::makeInteractiveComplexHeatmap(input = input, output = output, session = session,ht_list = combinedHMP_VAL, heatmap_id = "heatmap_1")
      InteractiveComplexHeatmap::makeInteractiveComplexHeatmap(input = input, output = output, session = session,ht_list = combinedHMP_VAL, heatmap_id = "heatmap_1", sub_heatmap_layer_fun = subHM_layer_fun)
      output$txtHMDescription = shiny::renderText(paste0(Sys.time(), " done calculating heatmap..., current plot is valid. n(probe) = ", nrow(as.matrix(combinedDFP_Val_Labels[[1]])), "n(reduced): ", nRowReduced, "; n(trait) = ", ncol(as.matrix(combinedDFP_Val_Labels[[1]])), "; elapsed time: ", elapsed_time))
      print(paste0(Sys.time()," finished plotting heatmap with n(probes)=", nrow(dfP_Val), " n(traits)=", ncol(dfP_Val)))

      orderedDistMatTraits <- as.matrix(result$distMatTraits)
      orderedDistMatTraits[upper.tri(orderedDistMatTraits)] <- NA
      diag(orderedDistMatTraits) <- NA
      orderedDistMatTraits <- reshape2::melt(orderedDistMatTraits)
      #sort by distance, closest first
      orderedDistMatTraits <- orderedDistMatTraits[order(orderedDistMatTraits$value),]
      colnames(orderedDistMatTraits) <- c("Var1", "Var2", "distance")
      orderedDistMatTraits<-stats::na.omit(orderedDistMatTraits)
      #apply names/colour of variable source
      if (ncol(sessionVariables$resultP_ValTrait1)>0) {
        trait1names <- colnames(sessionVariables$resultP_ValTrait1)
      }
      if (ncol(sessionVariables$resultP_ValTrait2)>0) {
        trait2names <- colnames(sessionVariables$resultP_ValTrait2)
      }
      if (ncol(sessionVariables$resultP_ValTrait3)>0) {
        trait3names <- colnames(sessionVariables$resultP_ValTrait3)
      }
      #annotation = annotation[which(annotation$name %in% rownames(beta)),]
      orderedDistMatTraits$Var1Source<-0
      orderedDistMatTraits$Var2Source<-0

      orderedDistMatTraits[which(orderedDistMatTraits$Var1 %in% trait1names),]$Var1Source <- 1
      orderedDistMatTraits[which(orderedDistMatTraits$Var1 %in% trait2names),]$Var1Source <- 2
      orderedDistMatTraits[which(orderedDistMatTraits$Var1 %in% trait3names),]$Var1Source <- 3

      orderedDistMatTraits[which(orderedDistMatTraits$Var2 %in% trait1names),]$Var2Source <- 1
      orderedDistMatTraits[which(orderedDistMatTraits$Var2 %in% trait2names),]$Var2Source <- 2
      orderedDistMatTraits[which(orderedDistMatTraits$Var2 %in% trait3names),]$Var2Source <- 3

      output$DTPairwiseDistTraits <- DT::renderDataTable(orderedDistMatTraits)

      DifferentOrderedDistMatTraits<- orderedDistMatTraits[which((orderedDistMatTraits$Var1Source != orderedDistMatTraits$Var2Source)),]

      output$DTDifferentOrderedDistMatTraits <- DT::renderDataTable(DifferentOrderedDistMatTraits)

      orderedDistMatProbes <- as.matrix(result$distMatProbes)
      orderedDistMatProbes[upper.tri(orderedDistMatProbes)] <- NA
      diag(orderedDistMatProbes) <- NA
      orderedDistMatProbes <- reshape2::melt(orderedDistMatProbes)
      #sort by distance, closest first
      orderedDistMatProbes <- orderedDistMatProbes[order(orderedDistMatProbes$value),]
      colnames(orderedDistMatProbes) <- c("probeID1", "probeID2", "distance")
      orderedDistMatProbes<-stats::na.omit(orderedDistMatProbes)
      output$DTPairwiseDistProbes<-DT::renderDataTable(orderedDistMatProbes)
    }
  }, ignoreNULL = FALSE)

  shiny::observeEvent(input$sldMinP_Val, ignoreInit = FALSE, {
    if (!is.na(as.integer(input$sldMinP_Val))) {
      sessionVariables$gP_ValMinBorder = 5*10^-as.integer(input$sldMinP_Val)
#      shiny::updateSliderInput(session, "sldMinP_Val", value = as.integer(input$sldMinP_Val))
#      output$txtDebugOut = shiny::renderText(sessionVariables$gP_ValMinBorder);
    }
    else{
      sessionVariables$gP_ValMinBorder = 5*10^-200
    }

  }, ignoreNULL = FALSE)

  shiny::observeEvent(input$sldMaxP_Val, ignoreInit = FALSE, {
    if (!is.na(as.integer(input$sldMaxP_Val))) {
      sessionVariables$gP_ValMaxBorder = 5*10^-as.integer(input$sldMaxP_Val)
#      shiny::updateSliderInput(session, "sldMaxP_Val", value = as.integer(input$sldMaxP_Val))
#      output$txtDebugOut = shiny::renderText(sessionVariables$gP_ValMaxBorder);
    }
    else{
      sessionVariables$gP_ValMaxBorder = 5*10^-18
    }

  }, ignoreNULL = FALSE)

  shiny::observeEvent(input$txtMaxProbes, ignoreInit = TRUE, {
    sessionVariables$gMaxProbes = input$txtMaxProbes
  }, ignoreNULL = FALSE)
}

