#' calculate models
#' @description function to calculate models defined in the config.yml of the working directory outside the shiny interactive mode using multip
#' @param workDir working directory to look for config.yml in
#' @export
CalculateModels <- function(workDir) {
  #read config.yml
  configFileLocation <- paste0(workDir, "/config.yml")
  base::print(base::paste0(sysTimePID(), " read config from: ", configFileLocation))
  config <- config::get(file = configFileLocation)
  loadObjectsForCalculation(config)
  numCores <- getNumCores()
  if (numCores < config$no_cores) {
    config$no_cores <- numCores
  }
  originalDir <- getwd()
  base::print(base::paste0(sysTimePID(), " starting calculations."))
  baseDir <-  paste0(workDir,"/all")
  if(!base::dir.exists(baseDir)) {
    base::dir.create(baseDir)
  }
  setwd(baseDir)
  ExposureVariableList<-PHENO[,config$firstPHENOVar:config$lastPHENOVar]
  res<-mclapply_preserve_names(ExposureVariableList, calculateEpigeneticEffectsAdjusted, config = config)

  base::print(base::paste0(sysTimePID(), " starting calculations for females."))
  baseDir <-  paste0(workDir,"/female")
  if(!base::dir.exists(baseDir)) {
    base::dir.create(baseDir)
  }
  setwd(baseDir)
  P <- PHENO[PHENO$sex=="w",]
  ExposureVariableList<-P[,config$firstPHENOVar:config$lastPHENOVar]
  res<-mclapply_preserve_names(ExposureVariableList, calculateEpigeneticEffectsAdjusted, config = config)

  base::print(base::paste0(sysTimePID(), " starting calculations for males."))
  baseDir <-  paste0(workDir,"/male")
  if(!base::dir.exists(baseDir)) {
    base::dir.create(baseDir)
  }
  setwd(baseDir)
  P <- PHENO[PHENO$sex=="m",]
  ExposureVariableList<-P[,config$firstPHENOVar:config$lastPHENOVar]
  res<-mclapply_preserve_names(ExposureVariableList, calculateEpigeneticEffectsAdjusted, config = config)

  setwd(originalDir)
  base::print(base::paste0(sysTimePID(), " finished calculations in ", workDir, "."))
}

#' load required data for model calculation
#' @description function to load all data to calculate regression model into memory
loadObjectsForCalculation <- function(config) {
  library(foreach)
  library(psych)
  annotation <- meffil::meffil.get.features("450k")
  assign("annotation",annotation,envir=globalenv())

  beta <- data.table::fread(file = config$betaFileName, stringsAsFactors=FALSE, header=TRUE, sep="\t")
  #  beta<-data.frame(column_to_rownames(beta, var = "PROBEID"))
  #rownames(beta) <- beta$probeID
  #beta$probeID <- NULL
  beta<-as.data.frame(beta)
  rownames(beta) <- beta$probeID
  beta$probeID <- NULL
  beta_wo_outliers <- removeOutliers(as.matrix(beta))
  beta <- as.data.frame(beta_wo_outliers[[1]])
  rm(beta_wo_outliers)
  #head(beta)
  #typeof(beta)
  #typeof(beta[5,5])
  assign("beta", beta, envir=globalenv())

  PHENO<-data.table::fread(file=config$PHENOFileName, sep="\t", dec=".", header=TRUE)
  #PHENO[,3] <- NULL
  rownames(PHENO) <- PHENO[[config$mergeAttribut]]
  PHENO <- as.data.frame(PHENO)
  rownames(PHENO) <- PHENO[[config$mergeAttribut]]
  #typeof(PHENO[5,5])
  adjustDF <- data.table::fread(file= config$adjustFileName, sep="\t", dec=".")
  adjustDF <- base::subset(adjustDF, select=c("ID_Kind", "PC1_cp", "PC2_cp", "PC3_cp", "PC4_cp", "PC5_cp", "PC6_cp", "sex", "age_mother_yrs", "gestational_age", "mat_smoking", "ME", "CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran", "nRBC", "batch"))
  PHENO <- base::merge(PHENO, adjustDF, by.x = config$mergeAttribut, by.y = config$mergeAttribut)
  #winsorize PHENO
  PHE<-PHENO
  foreach(i = config$firstPHENOVar:config$lastPHENOVar, .combine = cbind, .verbose=FALSE) %do% {
    PHE[,i] <- psych::winsor(PHENO[, i])
  }
  PHENO<-PHE
  rm(PHE)
  assign("PHENO", PHENO, envir=globalenv())
  #typeof(PHENO[15,15])

  #load list of multimodal CpG
  MultiModCpG <- data.table::fread(file = config$MultiModProbesFileName, sep="\t", dec=".")
  assign("MultiModCpG",MultiModCpG,envir=globalenv())
  beta <- removeMultiModelCpGFromBeta(beta, MultiModCpG, config)
  assign("beta", beta, envir=globalenv())
}

removeOutliers <- function(probes) {
  require(matrixStats)
  if(nrow(probes) < ncol(probes)) {
    base::message(base::paste0(sysTimePID(), " expect probes in rows (long dataset)"))
  }
  rowIQR <- matrixStats::rowIQRs(probes, na.rm = T)
  row2575 <- matrixStats::rowQuantiles(probes, probs = c(0.25, 0.75), na.rm = T)
  maskL <- probes < row2575[,1] - 3 * rowIQR
  maskU <- probes > row2575[,2] + 3 * rowIQR
  initial_NAs <- base::rowSums(is.na(probes))
  probes[maskL] <- NA
  removed_lower <- base::rowSums(base::is.na(probes))-initial_NAs
  probes[maskU] <- NA
  removed_upper <- base::rowSums(base::is.na(probes)) - removed_lower - initial_NAs
  N_for_probe <- base::rowSums(!is.na(probes))
  Log <- base::data.frame(initial_NAs, removed_lower, removed_upper, N_for_probe)
  return(list(probes, Log))
}

removeMultiModelCpGFromBeta <- function(probes, MultiModCpG, config) {
  if(nrow(probes) < ncol(probes)) {
    base::message(base::paste0(sysTimePID(), " expect probes in rows (long dataset)"))
  }
  probes[, config$probeAttribut] <- rownames(probes)
  result <- base::merge(probes, MultiModCpG, by.x = config$probeAttribut, by.y = "CpG")
  row.names(result) <- result[, config$probeAttribut]
  result[, config$probeAttribut] <- NULL
  #select only CpG with NumModes=1
  result <- result[result$NumModes<2, ]
  result$NumModes <- NULL
  result$NormalP <- NULL
  return (result)
}

mclapply_preserve_names <- function(list, fun, ...){
  parallel::mclapply(seq_along(list), function(i) {
    obj = list[i]
    names(obj) = names(list)[i]
    fun(obj, ...)
  })
}

#' calculate multiple models in parallel
#' @description function to calculate multiple models in parallel for given ExposureVariables
calculateEpigeneticEffectsAdjusted <- function(ExposureVariable, config) {
  ExV <- ExposureVariable
  Name <- (names(ExposureVariable)[1])
  base::message(base::paste0(sysTimePID(), " start calculating regression models for ", Name))
  analysisname <- paste0(as.character(Name),"adj")
  fileName <- paste0(analysisname,".csv")
  if (is.list(ExposureVariable)) {
    ExposureVariableValues<-unlist(ExposureVariable)
  }
  else {
    ExposureVariableValues<-as.numeric(ExposureVariable)
  }
  if (file_test("-f", fileName) != TRUE || file.size(fileName) == 0) {
    if ((is.numeric(ExposureVariableValues)) && (min(ExposureVariableValues, na.rm=TRUE)!=Inf && max(ExposureVariableValues, na.rm=TRUE)!=-Inf) && (min(ExposureVariableValues, na.rm=TRUE) != max(ExposureVariableValues, na.rm=TRUE))) {
      print(paste0("start ",Name))
      beta <- globalenv()$beta
      PHENO <- globalenv()$PHENO
      cl <- parallel::makeCluster(config$no_cores)
      doParallel::registerDoParallel(cl)
      parallel::clusterExport(cl, varlist = c("beta", "PHENO", "robregCTAdjustedWithMerge2"))
      result <- foreach::foreach(i = 1:nrow(beta), .combine = rbind, .packages = c("tidyverse", "MASS", "sandwich", "lmtest"), .verbose=FALSE) %dopar% {
      #result <- foreach::foreach(i = 1:10, .combine = rbind, .packages = c("tidyverse", "MASS", "sandwich", "lmtest"), .verbose=FALSE) %do% {
        result <- robregCTAdjustedWithMerge2(meth_matrix = beta,
                                            methcol = i,
                                            pheno_df = PHENO,
                                            exposure_variable = ExposureVariable,
                                            mergeIDPheno = config$mergeAttribut,
                                            mergeIDBeta = config$probeAttribut)
        return(result)
      }
      parallel::stopCluster(cl)

      #set rownames
      rownames(result) <- result[, "probeID"]
      result <- data.table::as.data.table(result)
      result[,1 ] <- base::as.numeric(base::unlist(result[, 1]))
      result[,2 ] <- base::as.numeric(base::unlist(result[, 2]))
      result[,3 ] <- base::as.numeric(base::unlist(result[, 3]))
      result[,4 ] <- base::as.numeric(base::unlist(result[, 4]))
      result[,6 ] <- base::as.numeric(base::unlist(result[, 6]))
      result[,7 ] <- base::as.numeric(base::unlist(result[, 7]))
      data.table::setattr(result, 'class', 'data.frame')
      data.table::setattr(result, "row.names", c(NA_integer_,4))
      data.table::setattr(result, "names", make.names(names(result), unique=TRUE))
      result <- data.table::data.table(result)

      # rename columns
      data.table::setnames(result, c("BETA","SE", "P_VAL", "N", "probeID", "DeltaMeth", "logFC"))
      data.table::setcolorder(result, c("probeID","BETA","SE", "P_VAL", "DeltaMeth", "logFC", "N"))

      # Adjustment for multiple testing
      result$FDR <- p.adjust(result$P_VAL, method="fdr")
      data.table::setcolorder(result, c("probeID","BETA","SE", "P_VAL", "FDR", "DeltaMeth", "logFC", "N"))

      result <- result[which(!is.na(result$BETA)), ]
      result <- result[order(result$P_VAL), ]
      result.write <- result
      data.table::fwrite(result.write, file=fileName, sep="\t", dec=".", col.names=TRUE)

      # Extract P_VAL columns
      P_VALS <- result$P_VAL
      lambda <- lambda_fun(P_VAL)
      N<-dim(beta)[1]
      CpGFails <- sum(is.na(result[, 2]))
      Tablestats <- t(c(N, CpGFails, lambda))
      colnames(Tablestats) <- c("Sample_Count", "CpGFails", "lambda")
      #  Write file
      fileName <- paste0(analysisname,"SampleSizeLambda.csv")
      data.table::fwrite(Tablestats, file=fileName, sep="\t", dec=".", col.names=TRUE)
    }
  }
  base::message(base::paste0(sysTimePID(), " finished calculating regression models for ", Name))
}

robregCTAdjustedWithMerge2 <- function(meth_matrix, methcol, pheno_df, exposure_variable, mergeIDPheno, mergeIDBeta) {
  name <- names(exposure_variable)
  Probe <- meth_matrix[methcol,]
  probeID <- rownames(Probe)
  #transpose
  Probe <- as.data.frame(t(Probe))
  Probe <- as.data.frame(rownames_to_column(Probe, var = mergeIDBeta))
  if (class(exposure_variable) == "character") {
    pheno <- subset(pheno_df, select=c(mergeIDPheno, exposure_variable, "CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran", "nRBC"))
  }
  else {
    pheno <- subset(pheno_df, select=c(mergeIDPheno, names(exposure_variable), "CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran", "nRBC"))
  }
  #merge CpG from beta with PHENO
  Pr <- base::merge(Probe, pheno, by.x = mergeIDBeta, by.y = mergeIDPheno)
  cf <- tryCatch({
    mod <- rlm(as.double(Pr[,2])~Pr[,name]+Pr$CD8T+Pr$CD4T+Pr$NK+Pr$Bcell+Pr$Mono+Pr$Gran+Pr$nRBC,maxit=200)
    cf <- coeftest(mod, vcov=vcovHC(mod, type="HC0"))
    deltaMeth <- max(mod$fitted.values)-min(mod$fitted.values)
    if(cf[2,1]<0) deltaMeth=deltaMeth*-1
    if (class(exposure_variable) == "character") {
      exposure_variable_name <- exposure_variable
    }
    else {
      exposure_variable_name <- names(exposure_variable)
    }
    #if binary trait var: take both values and calculate log2(mean(class1)/mean(class2))
    #if numeric trait var: class1 = lower median, class2 = upper median
    if (length(levels(as.factor(Pr[,exposure_variable_name]))) == 2) {
      level1 <- levels(as.factor(Pr[,exposure_variable_name]))[1]
      level2 <- levels(as.factor(Pr[,exposure_variable_name]))[2]
      mean1 <- base::mean(Pr[(which(Pr[,exposure_variable_name] == level1)), 2])
      mean2 <- base::mean(Pr[(which(Pr[,exposure_variable_name] == level2)), 2])
    }
    else {
      #calculate median:
      median <- stats::median(Pr[,exposure_variable_name], na.rm = TRUE)
      mean1 <- base::mean(Pr[(which(Pr[,exposure_variable_name] < median)), 2])
      mean2 <- base::mean(Pr[(which(Pr[,exposure_variable_name] > median)), 2])
    }
    #logConst <- 10e256
    #logFC <- log2((mean1 + logConst) / (mean2 + logConst))
    logFC <- log2((mean1) / (mean2))
    if (is.infinite(logFC)) {
      logFC <- "-Inf"
    }
    N <- length(mod$fitted.values)
    #    print(N)
    cf <- cf[2, c("Estimate", "Std. Error", "Pr(>|z|)")]
    cf <- append(cf, N) #; names(cf)[4] <- "N"
    cf <- append(cf, probeID) #; names(cf)[5] <- "probeID"
    cf <- append(cf, deltaMeth)
    cf <- append(cf, logFC)
    cf
  }, error=function(err){
    cf <- c(NaN, NaN, NaN, NaN, probeID, NaN, NaN)
  });
  names(cf) <- c("Estimate", "Std. Error", "Pr(>|z|)", "N", "probeID", "DeltaMeth", "logFC")
  return(cf)
}

# calculate lambdas
lambda_fun<-function(x){
  qchisq(median(x,na.rm=T), df =1, lower.tail=F)/qchisq(0.5,1)
}
