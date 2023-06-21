#' structure of sessionVariables:
#'
#' result data to be analyzed loaded from regression results
#' session$userData$sessionVariables$resultDFListTrait1()$resultDFP_Val - list of resulting p-values
#' session$userData$sessionVariables$resultDFListTrait1()$resultDFDM - list of resulting delta methylations
#' session$userData$sessionVariables$resultDFListTrait1()$resultDFN - list of resulting n
#' session$userData$sessionVariables$resultDFListTrait1()$listPHENOdata[[x]]$PHENODF - original data, the regression is based on
#' session$userData$sessionVariables$resultDFListTrait1()$listPHENOdata[[x]]$PHENOFileName - file name of the original data
#' PHENOdata <- list(PHENODF = getPHENODF(PHENOFileName),PHENOFileName = PHENOFileName)

#' loadResultDF
#' loads data.frame from folder
#' @param folder folder containing all files to read to list
#' @param session session object
#' @param loadRDS loads RDS file or not (if only RDS files should be written)
#' @return result list()
#' data.frame P_Val,
#' data.frame DM (DeltaMethylation),
#' data.frame N,
#' data.frame PHENOdata,
#' chr PHENOFileName
#'
# examples loadResultDF(session, folder, FALSE)
loadResultDF <- function(session, folder, loadRDS = FALSE) {
  tryCatch(
    {
      base::print(base::paste0(Sys.time(), " start loadResultDF() for ", folder))
      fileNameLPV <- "listResultP_Val_DeltaMeth_N.RDS"
      fileNameLPV <- base::paste0(folder, fileNameLPV)
      doLoadFolderDFList <- TRUE
      if (utils::file_test("-f", fileNameLPV) == TRUE && getYoungestFile(folder) == fileNameLPV) {
        if (loadRDS != FALSE) {
          base::print(base::paste0(Sys.time(), " read loadResultDF() from ", fileNameLPV))
          listResultP_Val_DeltaMeth_N <- base::readRDS(file = fileNameLPV)
          doLoadFolderDFList <- FALSE
        } else {
          base::print(base::paste0(Sys.time(), " can't read RDS in loadResultDF() from ", fileNameLPV))
        }
      }
      else {
        doLoadFolderDFList <- TRUE
      }
      # check for original data in .RDS
      if (doLoadFolderDFList == FALSE) {
        if (checklistResultP_Val_DeltaMeth_NValidity(listResultP_Val_DeltaMeth_N) == FALSE) {
          doLoadFolderDFList <- TRUE
        }
      }
#      else {
      if (doLoadFolderDFList == TRUE) {
        # load or generate list of DFs with results
        base::print(base::paste0(Sys.time(), " read from loadFolderDFList()."))
        listOfResultDF <- loadFolderDFList(session = session, folder = folder) #listOfResultDF <- loadFolderDFList(folder, globalVariables)
        if (base::length(listOfResultDF) != 0) {
          # read original data for each member of listofResultDF
          # look in data folder into all config.yml files for PHENOFileName:
          PHENOFileNameLine <-
            findInFiles(
              "PHENOFileName:",
              folder,
              in_files = "\\config.yml$",
              recursive = FALSE
            ) #look into data folder
          if (base::length(PHENOFileNameLine) == 0 || PHENOFileNameLine == FALSE) {
          PHENOFileNameLine <-
            findInFiles(
              "PHENOFileName:",
              paste0(folder, ".."),
              in_files = "\\config.yml$",
              recursive = FALSE
            ) #look one folder above
          }
          if (base::length(PHENOFileNameLine) == 0 || PHENOFileNameLine == FALSE) {
            base::message(base::paste0(Sys.time(), " config.yml in data folder or
                                       PHENOFileNameLine not found. Consider placing
                                       a config.yml file with PHENOFileName:
                                       <location of PHENO file> in each data folder."))
          }
          else {
            PHENOFileName <-
              stringr::str_match(PHENOFileNameLine, "\"(.*?)\"")[2]
            PHENOFileName <- paste0(folder, PHENOFileName)
            if (base::length(PHENOFileName) == 0) {
              base::message(base::paste0(Sys.time(), " PHENOFileName in PHENOFileNameLine not found."))
            }
            else {
              listPrimaryKeys <- as.list(session$userData$config$keyAttributes)
              PHENOdata <- list(PHENODF = getPHENODF(PHENOFileName, listPrimaryKeys), PHENOFileName = PHENOFileName)
              listResultP_Val <- getResultDfP_D_N(listOfResultDF, "P")
              rn <- rownames(listResultP_Val)
              listResultP_Val <- base::as.data.frame(listResultP_Val)
              rownames(listResultP_Val) <- rn
              listResultDeltaMeth <- getResultDfP_D_N(listOfResultDF, "D")
              rn <- rownames(listResultDeltaMeth)
              listResultDeltaMeth <- base::as.data.frame(listResultDeltaMeth)
              rownames(listResultDeltaMeth) <- rn
              listResultN <- getResultDfP_D_N(listOfResultDF, "N")
              rn <- rownames(listResultN)
              listResultN <- base::as.data.frame(listResultN)
              rownames(listResultN) <- rn
              listResultP_Val_DeltaMeth_N <-
                list(
                  P_Val = listResultP_Val,
                  DM = listResultDeltaMeth,
                  N = listResultN,
                  PHENOdata = PHENOdata
                )
            }
          }
          if (base::exists("listResultP_Val_DeltaMeth_N")) {
            rn <- rownames(listResultP_Val_DeltaMeth_N$P_Val)
            if (!base::is.data.frame(listResultP_Val_DeltaMeth_N$P_Val)) {
              listResultP_Val_DeltaMeth_N$P_Val <-
                base::as.data.frame(listResultP_Val_DeltaMeth_N$P_Val)
            }
            if (!base::is.data.frame(listResultP_Val_DeltaMeth_N$DM)) {
              listResultP_Val_DeltaMeth_N$DM <-
                base::as.data.frame(listResultP_Val_DeltaMeth_N$DM)
            }
            if (!base::is.data.frame(listResultP_Val_DeltaMeth_N$N)) {
              listResultP_Val_DeltaMeth_N$N <-
                base::as.data.frame(listResultP_Val_DeltaMeth_N$N)
            }
            #        rownames(listResultP_Val_DeltaMeth_N$P_Val) <- rn
            rownames(listResultP_Val_DeltaMeth_N$DM) <- rn
            rownames(listResultP_Val_DeltaMeth_N$N) <- rn
          }
          base::print(base::paste0(Sys.time(), " saveRDS loadFolderDFList()")) # , fileNameLPV, "."))
          base::saveRDS(listResultP_Val_DeltaMeth_N, file = fileNameLPV)
        } else {
          base::message(base::paste0(Sys.time(), " length(listOfResultDF) == 0."))
        }
      }
      if (loadRDS != FALSE) {
        if (session$userData$config$debugMode == TRUE) { #if (globalVariables$config$debugMode == TRUE) {
          # it is mandatory to reassign rownames and also colnames... Thanks R!
          rn <- base::rownames(listResultP_Val_DeltaMeth_N$P_Val)[1:1000]
          cn <- base::colnames(listResultP_Val_DeltaMeth_N$P_Val)
          listResultP_Val_DeltaMeth_N$P_Val <-
            base::as.data.frame(listResultP_Val_DeltaMeth_N$P_Val[1:1000, ])
          listResultP_Val_DeltaMeth_N$DM <-
            base::as.data.frame(listResultP_Val_DeltaMeth_N$DM[1:1000, ])
          listResultP_Val_DeltaMeth_N$N <-
            base::as.data.frame(listResultP_Val_DeltaMeth_N$N[1:1000, ])
          rownames(listResultP_Val_DeltaMeth_N$P_Val) <- rn
          rownames(listResultP_Val_DeltaMeth_N$DM) <- rn
          rownames(listResultP_Val_DeltaMeth_N$N) <- rn
          colnames(listResultP_Val_DeltaMeth_N$P_Val) <- cn
          colnames(listResultP_Val_DeltaMeth_N$DM) <- cn
          colnames(listResultP_Val_DeltaMeth_N$N) <- cn
        }
        result <- listResultP_Val_DeltaMeth_N
        base::print(base::paste0(Sys.time(), " finished read from loadFolderDFList()", "."))
      } else {
        result <- FALSE
      }
    },
    error = function(e) {
      message("An error occurred in loadResultDF():\n", e)
    },
    warning = function(w) {
      message("A warning occurred in loadResultDF():\n", w)
    },
    finally = {
      base::print(base::paste0(Sys.time(), " end loadResultDF()."))
      return(result)
    }
  )
}

#' loadtraitDFs
#' loads and concatenates corresponding data.frames for p_val (DFP_Val), delta methylation (DFDM) and n (DFN)
#' @param traitDFs data.frames
#' @return result list()
#' data.frame P_Val,
#' data.frame DM (DeltaMethylation),
#' data.frame N,
#' data.frame PHENOdata,
#' chr PHENOFileName,
#' list resultOriginDF
#' list resultColnames
# examples loadtraitDFs(traitDFs)
loadtraitDFs <- function(traitDFs) {
  tryCatch(
    {
      # browser() #everything seems fine here (all DFs become loaded)
      #listPHENOdata <- base::list(1:base::length(traitDFs))
      listPHENOdata <- base::list(seq_along(traitDFs))
      i <- NULL
      #foreach(i = 1:base::length(traitDFs)) %do% {
      foreach::foreach(i = seq_along(traitDFs)) %do% {
        rownames(traitDFs[[i]]$P_Val)
        colnames(traitDFs[[i]]$P_Val)
        if (i == 1) {
          resultDFP_Val <- traitDFs[[i]]$P_Val # [[1]]
          resultDFDM <- traitDFs[[i]]$DM # [[2]]
          resultDFN <- traitDFs[[i]]$N # [[3]]
          resultDFP_Val$Row.names <- rownames(traitDFs[[i]]$P_Val)
          resultDFDM$Row.names <- rownames(traitDFs[[i]]$DM)
          resultDFN$Row.names <- rownames(traitDFs[[i]]$N)
          resultColnames <- base::colnames(traitDFs[[i]]$P_Val)
          resultOriginDF <- base::rep(i, length(resultColnames))
        } else {
          base::print(base::paste0(Sys.time(), " merge trait ", i, "."))
          traitDFs[[i]]$P_Val$Row.names <-
            base::rownames(traitDFs[[i]]$P_Val)
          traitDFs[[i]]$DM$Row.names <- base::rownames(traitDFs[[i]]$DM)
          traitDFs[[i]]$N$Row.names <- base::rownames(traitDFs[[i]]$N)
          OriginColnames <- base::colnames(traitDFs[[i]]$P_Val)
          OriginDF <- base::rep(i, length(OriginColnames))

          # Cannot merge with rownames due to:
          # A non-empty vector of column names is required for `by.x` and `by.y`.
          # resultDFP_Val = base::merge(x = resultDFP_Val, y = traitDFs[[i]]$P_Val, by.x = 0, by.y = 0, all.x = TRUE, all.y = TRUE)
          resultDFP_Val <-
            base::merge(
              x = resultDFP_Val,
              y = traitDFs[[i]]$P_Val,
              by.x = "Row.names",
              by.y = "Row.names",
              all.x = TRUE,
              all.y = TRUE
            )
          resultDFDM <-
            base::merge(
              x = resultDFDM,
              y = traitDFs[[i]]$DM,
              by.x = "Row.names",
              by.y = "Row.names",
              all.x = TRUE,
              all.y = TRUE
            )
          resultDFN <-
            base::merge(
              x = resultDFN,
              y = traitDFs[[i]]$N,
              by.x = "Row.names",
              by.y = "Row.names",
              all.x = TRUE,
              all.y = TRUE
            )
          resultOriginDF <- base::c(resultOriginDF, OriginDF)
          resultColnames <- base::c(resultColnames, OriginColnames)
        }
        listPHENOdata[[i]] <- traitDFs[[i]]$PHENOdata
      }
      # omit column Row.names
      if ("Row.names" %in% colnames(resultDFP_Val)) {
        rownames(resultDFP_Val) <- resultDFP_Val$Row.names
        resultDFP_Val$Row.names <- NULL
      }
      if ("Row.names" %in% colnames(resultDFDM)) {
        rownames(resultDFDM) <- resultDFDM$Row.names
        resultDFDM$Row.names <- NULL
      }
      if ("Row.names" %in% colnames(resultDFN)) {
        rownames(resultDFN) <- resultDFN$Row.names
        resultDFN$Row.names <- NULL
      }
      result <- base::list()
      result$resultDFP_Val <- resultDFP_Val
      result$resultDFDM <- resultDFDM
      result$resultDFN <- resultDFN
      result$rownames <- base::rownames(resultDFP_Val)
      result$colnames <- base::colnames(resultDFP_Val)
      result$listPHENOdata <- listPHENOdata
      result$resultOriginDF <- resultOriginDF
      result$resultColnames <- resultColnames
    },
    error = function(e) {
      message("An error occurred in loadtraitDFs():\n", e)
    },
    warning = function(w) {
      message("A warning occurred in loadtraitDFs():\n", w)
    },
    finally = {
      base::print(base::paste0(Sys.time(), " end loadtraitDFs()."))
      return(result)
    }
  )
}

#' getlistOfResultsDF
#' gets list of result data.frame from folder
#' @param folder folder containing all files to read to list
#' @param session delivers session variables to this function
#' @return list of data.frame from folder
# examples getlistOfResultsDF(session, folder)
#getlistOfResultsDF <- function(session, folder, globalVariables) {
getlistOfResultsDF <- function(session, folder) {
  tryCatch(
    {
      base::print(base::paste0(Sys.time(), " before list.files()"))
      if (base::dir.exists(folder)) {
        temp <- base::list.files(path = folder, pattern = "*.csv")
        result <- list()
        #for (i in 1:base::length(temp)) {
        for (i in seq_along(temp)) {
          trait <-
            stringr::str_sub(temp[i], 1, stringr::str_length(temp[i]) - 4)
          fileName <-
            base::paste0(folder, base::as.character(trait), ".csv")
          if (base::file.exists(fileName)) {
            firstlines <-
              utils::read.table(
                file = fileName,
                sep = "\t",
                dec = ".",
                header = TRUE,
                nrows = 5
              )
            if (colnames(firstlines)[1] == "probeID") {
              if (nrow(firstlines) >= 5) {
                #          if (grepl("adj", temp[i], fixed = TRUE) == TRUE) {
                base::print(base::paste0(Sys.time(), " file ", i, " of ", base::length(temp)))
                base::print(base::paste0(Sys.time(), " load result file ", folder, trait))
                # read results into DF
                resultDF <-
                  loadResultFile(session = session, folder = folder, fileName = trait) #loadResultFile(folder, trait, globalVariables)
                # omit unneccesary variables
                resultDF <-
                  resultDF[, c("probeID", "P_VAL", "DeltaMeth", "N")]
                tryCatch({
                  if (base::min(resultDF$P_VAL, na.rm = TRUE) <
                      base::as.numeric(session$userData$config$P_VALWarningThreshold)) {
                    base::message(
                      base::paste0(
                        Sys.time(),
                        " p-value below warning threshold found in ",
                        folder
                      )
                    )
                  }
                  resultDF <- list(resultDF)
                  names(resultDF) <- trait
                  if (length(result) > 0) {
                    result <- base::append(result, resultDF)
                  } else {
                    result <- resultDF
                  }
                },
                error = function(e) {
                  base::message("An error occurred in getlistOfResultsDF(), inner tryCatch:\n", e)
                },
                warning = function(w) {
                  base::message(
                    base::paste0(
                      Sys.time(),
                      " no other p-values than NA found in ",
                      trait
                    )
                  )
                },
                finally = {

                }
                )
              }
            }
            else {
              result <- FALSE
              base::message(
                base::paste0(
                  Sys.time(),
                  " only 1 column found in file ", fileName, ". Are you using \t as column separator?"
                )
              )
            }
          }
        }
      } else {
        result <- FALSE
        base::message(base::paste0("folder ", folder, " does not exist."))
      }
    },
    error = function(e) {
      message("An error occurred in getlistOfResultsDF():\n", e)
    },
    warning = function(w) {
      message("A warning occurred in getlistOfResultsDF():\n", w)
    },
    finally = {
      base::print(base::paste0(Sys.time(), " end getlistOfResultsDF()."))
      return(result)
    }
  )
}

#' loadFolderDFList
#' loads list of data.frame from folder
#' @param folder folder containing all files to read to list
#' @param session session object
#' @return list of data.frame from folder
# examples loadFolderDFList(session, folder)
loadFolderDFList <- function(session, folder) {
  tryCatch(
    {
      fileNameLR <- "listOfResultsDF.RDS"
      fileNameLR <- base::paste0(folder, fileNameLR)
      if (utils::file_test("-f", fileNameLR) == TRUE && getYoungestFile(folder) == fileNameLR) {
        listOfResultDF <- readRDS(file = fileNameLR)
        getList <- FALSE
        if (base::length(listOfResultDF) == 0) {
          getList <- TRUE
        } else {
          getList <- FALSE
        }
      } else {
        getList <- TRUE
      }
      #check for original data in RDS
#browser()
      if (getList == TRUE) {
        #        listOfResultDF = getlistOfResultsDF(dataDir())
        base::print(base::paste0(Sys.time(), " before getlistOfResultsDF()"))
        listOfResultDF <- getlistOfResultsDF(session = session, folder = folder) #listOfResultDF <- getlistOfResultsDF(folder, globalVariables)
        if (listOfResultDF != FALSE) {
          base::print(base::paste0(Sys.time(), " saveRDS getlistOfResultsDF()", fileNameLR))
          base::saveRDS(listOfResultDF, file = fileNameLR)
        }
      }
    },
    error = function(e) {
      message("An error occurred in loadFolderDFList():\n", e)
    },
    warning = function(w) {
      message("A warning occurred in loadFolderDFList():\n", w)
    },
    finally = {
      base::print(base::paste0(Sys.time(), " end loadFolderDFList()."))
      return(listOfResultDF)
    }
  )
}

#' getPHENODF
#' loads a pheno DF for a certain PHENOFileName
#' @param PHENOFileName filename for pheno DF to load
#' @param listPrimaryKeys list with possible primary keys in original data
#' @return pheno DF for a certain PHENOFileName
# examples getPHENODF(PHENOFileName)
getPHENODF <- function(PHENOFileName, listPrimaryKeys) {
  tryCatch(
    {
      if (utils::file_test("-f", PHENOFileName) == TRUE) {
        #  if (length(base::readLines(PHENOFileName))>0) {
        PHENO <- data.table::fread(
          file = PHENOFileName,
          sep = "\t",
          dec = ".",
          header = TRUE,
          quote = ""
        )
#try different primary keys for rownames^and take the first
        #foreach::foreach(key = seq_along(listPrimaryKeys)) %do% {
        #foreach::foreach(key in listPrimaryKeys) %do% {
        found <- FALSE
        for (key in listPrimaryKeys) {
          if (key %in% colnames(PHENO)) {
            found <- TRUE
            keyvar <- key
            break
          }
        }
        if (found == TRUE) {
          # PHENO <-
          #   data.frame(tibble::column_to_rownames(PHENO, var = "ID_Kind2"))
          PHENO <-
            data.frame(tibble::column_to_rownames(PHENO, var = keyvar))
        }
        else {
          base::message(base::paste0(Sys.time(), " no key found in data. (from list:", as.character(listPrimaryKeys), "."))
        }
        PHENO <- as.data.frame(PHENO)
        result <- PHENO
      } else {
        base::message(base::paste0(Sys.time(), " file not found: ", PHENOFileName))
        result <- FALSE
      }
    },
    error = function(e) {
      message("An error occurred in getPHENODF():\n", e)
    },
    warning = function(w) {
      message("A warning occurred in getPHENODF():\n", w)
    },
    finally = {
      base::print(base::paste0(Sys.time(), " end getPHENODF()."))
      return(result)
    }
  )
}

#' getAvailNForP_VALBorder
#' counts traits with at least 2 elements > P_VAL_BORDER
#' @param DF data.frame with P_Val
#' @return data.frame with reduced data set
# examples getAvailNForP_VALBorder(data.frame)
getAvailNForP_VALBorder <- function(DF) {
  tryCatch(
    {
      numRows <- 300
      result <- base::matrix(nrow = numRows, ncol = 2)
      for (i in 1:numRows) {
        mat <- DF
        P_VAL_BORDER <- 5 * 10^-i
        mat[mat > P_VAL_BORDER] <- NA
        mat <-
          delete.na(mat, ncol(mat) - 1) # -1, because we need at least 2 traits to associate
        n <- base::nrow(mat)
        if (!base::is.numeric(n)) {
          break()
        }
        if (n <= 1) {
          break()
        }
        base::print(
          base::paste0(
            Sys.time(),
            " counting remaining probes at p = ",
            P_VAL_BORDER,
            " remaining n = ",
            n
          )
        )
        result[i, 1] <- P_VAL_BORDER
        result[i, 2] <- n
      }
      colnames(result) <- base::c("P_VAL_BORDER", "Available n")
      result <- result[1:i - 1, ]
      result <- base::as.data.frame((result))
      result <- result[base::order(result[1]), ]
    },
    error = function(e) {
      message("An error occurred in getAvailNForP_VALBorder():\n", e)
    },
    warning = function(w) {
      message("A warning occurred in getAvailNForP_VALBorder():\n", w)
    },
    finally = {
      base::print(base::paste0(Sys.time(), " end getAvailNForP_VALBorder()."))
      return(result)
    }
  )
}

#' delete.na
#' deletes rows with all NA from data.frame DF
#' @param df data.frame
#' @param n minimum number of NA
#' @return data.frame
# examples delete.na(df, 0)
delete.na <- function(df, n = 0) {
  # source: https://stackoverflow.com/questions/4862178/remove-rows-with-all-or-some-nas-missing-values-in-data-frame
  df[base::rowSums(base::is.na(df)) <= n, ]
}

#' getNForP_ValBorder
#' counts number of remaining models below a defined n
#' @param mat matrix for which to calculate results
#' @param n minimum n for which to calculate results
#' @return data.frame with results
# examples getNForP_ValBorder(mat, n)
getNForP_ValBorder <- function(mat, n) {
  tryCatch({
    result <- base::matrix(nrow = 1, ncol = 2)
    P_VAL_BORDER <- 5 * 10^-n
    mat[mat > P_VAL_BORDER] <- NA
    mat <-
      delete.na(mat, ncol(mat) - 1) # -1, because we need at least 2 traits to associate
    nrow <- base::nrow(mat)
    # base::print(
    #   base::paste0(
    #     Sys.time(),
    #     " counting remaining probes at p = ",
    #     P_VAL_BORDER,
    #     " remaining n = ",
    #     nrow
    #   )
    # )
    if (base::is.numeric(nrow)) {
      result[1, 1] <- P_VAL_BORDER
      result[1, 2] <- nrow
    }
    else {
      result[1, 1] <- P_VAL_BORDER
      result[1, 2] <- 0
    }
  },
  error = function(e) {
    message("An error occurred in getNForP_ValBorder():\n", e)
  },
  warning = function(w) {
    message("A warning occurred in getNForP_ValBorder():\n", w)
  },
  finally = {
    base::print(base::paste0(Sys.time(), " end getNForP_ValBorder()."))
    return(result)
  }
  )
}

#' getAvailNForP_VALBorderParallel
#' counts number of remaining models below a defined n; uses parallel processing for faster results
#' @param numCores number of CPU cores to use
#' @param DF data frame for which to calculate results
#' @return data.frame with results
# examples getAvailNForP_VALBorderParallel(numCores, DF)
getAvailNForP_VALBorderParallel <- function(numCores, DF) {
  tryCatch(
    {
      base::print(base::paste0(Sys.time(), " start getAvailNForP_VALBorderParallel()."))
      i <- NULL
      numRows <- 300
      result <- base::matrix(nrow = numRows, ncol = 2)
      cl <- parallel::makeCluster(numCores) #makeCluster(no_cores)
      doParallel::registerDoParallel(cl)
      parallel::clusterExport(cl, varlist = c("DF", "result", "numRows"), envir = environment())
      #browser() #check, if result is NA, then data or exported function is not known in the compute nodes
      result <- foreach::foreach(i = 1:numRows, .combine = rbind, .packages = c("base"),
                                 .export = c("getNForP_ValBorder", "delete.na"),
                                 .verbose = TRUE) %dopar% {
        base::source(paste0(getwd(), "/R/ResultData.R")) #this is necessary for foreach %dopar% to run properly
        result <- getNForP_ValBorder(mat = DF, n = i)
        return(result)
      }
      parallel::stopCluster(cl)
      colnames(result) <- base::c("P_VAL_BORDER", "Available n")
      result <- base::as.data.frame((result))
      result <- result[base::order(result[[1]], decreasing = TRUE), ]
    },
    error = function(e) {
      message("An error occurred in getAvailNForP_VALBorderParallel():\n", e)
    },
    warning = function(w) {
      message("A warning occurred in getAvailNForP_VALBorderParallel():\n", w)
    },
    finally = {
      base::print(base::paste0(Sys.time(), " end getAvailNForP_VALBorderParallel()."))
      return(result)
    }
  )
}

#' getReducedP_ValMatrix
#' counts and gets back traits with at least 2 elements > P_VAL_BORDER
#' @param df data.frame with P_Val, rows are probes, cols are traits
#' @param numRows scalar with maximum allowed row number
#' @param upperP_VALborder scalar with maximum allowed P_Val
#' @param lowerP_VALborder optionalscalar with minimum allowed P_Val
#' @return matrix with reduced data set
# examples getReducedP_ValMatrix(matrix)
getReducedP_Valdf <-
  function(df,
           numRows,
           upperP_VALborder,
           lowerP_VALborder) {
    tryCatch(
      {
        if (!base::is.data.frame(df)) {
          df <- base::as.data.frame(df)
        }
        if (base::missing(upperP_VALborder)) {
          upperP_VALborder <- 0.05
        }
        orig_df <- df
        # remove rowsProbes with all p-val > 0.05
        df[df > upperP_VALborder] <- NA
        #  if (!missing(lowerP_VALborder) && length(lowerP_VALborder) == 0) {
        if (!base::missing(lowerP_VALborder)) {
          df[df < lowerP_VALborder] <- NA
        }
        df <-
          delete.na(df, base::ncol(df) - 1) # -1, because we need at least 2 traits to associate
        df[base::is.na(df)] <- 1
        result_df <- orig_df[base::rownames(df), base::colnames(df)]
        base::print(base::paste0("upperP_VALborder: ", upperP_VALborder))
        base::print(base::paste0("lowerP_VALborder: ", lowerP_VALborder))
        base::print(base::paste0("nrow result matrix: ", nrow(df)))
        # browser()
        if (!base::missing(numRows)) {
          if (base::is.numeric(numRows)) {
            #      if (numRows < nrow(mat)) {
            if (numRows < base::nrow(result_df)) {
              #        mat = mat[1:numRows,]
              result_df <- result_df[1:numRows, ]
            }
          }
        }
        #  return(mat)
        base::print(dim(result_df))
      },
      error = function(e) {
        message("An error occurred  in getReducedP_Valdf():\n", e)
      },
      warning = function(w) {
        message("A warning occurred in getReducedP_Valdf():\n", w)
      },
      finally = {
        base::print(base::paste0(Sys.time(), " end getReducedP_Valdf()."))
        return(result_df)
      }
    )
  }

#' removeTraitsMinN
#' removes CpG with casecount < minN from corresponding data.frames for P_Val, for DeltaMethylation and for N
#' @param dfList list containing corresponding data.frames for P_Val, for DeltaMethylation and for N
#' @param minN minimum value for n
#' @return named list of data.frames, one df for P_Val, one for DeltaMethylation, one for N as well as labels
#' @return result$dfP_Val for p-values
#' @return result$dfDM for delta methylation values
#' @return result$dfN for n
# examples removeTraitsMinN(dfList, minN)
removeTraitsMinN <- function(dfList, minN) {
  # check for minimum n in each trait
  # also remove from <resultOriginDF> and <resultColnames>
  # check for valid minN
  tryCatch(
    {
      if (exists("minN")) {
        dfN <- dfList$dfN
        rn <- base::rownames(dfN)
        cn <- base::colnames(dfN)
        dfP_Val <- dfList$dfP_Val
        dfDM <- dfList$dfDM
        resultOriginDF <- dfList$resultOriginDF
        resultColnames <- dfList$resultColnames
        listPHENOdata <- dfList$listPHENOdata
        dfN <- base::as.data.frame(dfN)
        dfP_Val <- base::as.data.frame(dfP_Val)
        rownames(dfP_Val) <- rn
        dfDM <- base::as.data.frame(dfDM)
        rownames(dfDM) <- rn
        traitNames <-
          stats::na.omit(colnames(dfN)[matrixStats::colMins(as.matrix(dfN),
                                                            na.rm = TRUE) > base::as.integer(minN)])
        positions <- which(colnames(dfN) %in% traitNames)
        dfN <-
          base::as.data.frame(dfN[, traitNames])
        rownames(dfN) <- rn
        colnames(dfN) <- cn

        # select the same content than in dfN
        dfP_Val <- as.data.frame(dfP_Val[base::rownames(dfN), base::colnames(dfN)]) # does not work with data.table
        rownames(dfP_Val) <- rn
        colnames(dfP_Val) <- cn
        dfDM <- as.data.frame(dfDM[base::rownames(dfN), base::colnames(dfN)])
        rownames(dfDM) <- rn
        colnames(dfDM) <- cn
        resultOriginDF <- resultOriginDF[positions]
        resultColnames <- resultColnames[positions]
        dfList <- base::list(
          dfP_Val = NULL,
          dfDM = NULL,
          dfN = NULL,
          resultOriginDF = NULL,
          resultColnames = NULL,
          listPHENOdata = NULL
        )
        dfList$dfP_Val <- dfP_Val
        dfList$dfDM <- dfDM
        dfList$dfN <- dfN
        dfList$resultOriginDF <- resultOriginDF
        dfList$resultColnames <- resultColnames
        dfList$listPHENOdata <- listPHENOdata
      }
      else {
        base::print(base::paste0(Sys.time(), " minN does not exist."))
      }
    },
    error = function(e) {
      message("An error occurred in removeTraitsMinN():\n", e)
    },
    warning = function(w) {
      message("A warning occurred in removeTraitsMinN():\n", w)
    },
    finally = {
      base::print(base::paste0(Sys.time(), " end removeTraitsMinN()."))
      #browser() #check for rownames and colnames at the end
      return(dfList)
    }
  )
}

#' checklistResultP_Val_DeltaMeth_NValidity
#' checks for validity of structure
#' @param listResultP_Val_DeltaMeth_N structure to check
#' @return TRUE or FALSE
# examples checklistResultP_Val_DeltaMeth_NValidity(listResultP_Val_DeltaMeth_N)
checklistResultP_Val_DeltaMeth_NValidity <- function(listResultP_Val_DeltaMeth_N) {
  result <- TRUE
  if (!is.valid(listResultP_Val_DeltaMeth_N$P_Val)) {
    result <- FALSE
  }
  if (!is.valid(listResultP_Val_DeltaMeth_N$DM)) {
    result <- FALSE
  }
  if (!is.valid(listResultP_Val_DeltaMeth_N$N)) {
    result <- FALSE
  }
  # if (!is.valid(listResultP_Val_DeltaMeth_N$resultOriginDF)) {
  #   result <- FALSE
  # }
  # if (!is.valid(listResultP_Val_DeltaMeth_N$resultColnames)) {
  #   result <- FALSE
  # }
  if (!is.valid(listResultP_Val_DeltaMeth_N$PHENOdata$PHENODF)) {
    result <- FALSE
  }
  if (!is.valid(listResultP_Val_DeltaMeth_N$PHENOdata$PHENOFileName)) {
    result <- FALSE
  }
  return(result)

}

#' checkResultP_Val_cg
#' checks whether rownames of results start with 'cg'
#' @param listResultP_Val structure to check
#' @return TRUE or FALSE
# examples checkResultP_Val_cg(listResultP_Val)
checkResultP_Val_cg <- function(listResultP_Val) {
  if (base::exists("listResultP_Val")) {
    if (length(rownames(listResultP_Val)[1]) != 0) {
      if (stringr::str_starts(
        rownames(
          listResultP_Val
        )[1],
        "cg"
      ) != TRUE) {
        base::message(
          base::paste0(
            Sys.time(),
            "warning: rownames(listResultP_Val)[1] does not start with 'cg'"
          )
        )
        result <- FALSE
      }
      else {
        result <- TRUE
      }
    }
    else {
      base::message(
        base::paste0(
          Sys.time(),
          "warning: length(listResultP_Val)[1] == 0"
        )
      )
      result <- FALSE
    }
  }
  return(result)
}
