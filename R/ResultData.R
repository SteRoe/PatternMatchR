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
  base::tryCatch(
    {
      fileNameLPV <- "listResultP_Val_DeltaMeth_N.RDS"
      fileNameLPV <- base::paste0(folder, fileNameLPV)
#browser()
      doLoadFolderDFList <- TRUE
      if (utils::file_test("-f", fileNameLPV) == TRUE && getYoungestFile(folder) == fileNameLPV) {
        if (loadRDS != FALSE) {
#          base::print(base::paste0(sysTimePID(), " read loadResultDF() from ", fileNameLPV))
          listResultP_Val_DeltaMeth_N <- base::readRDS(file = fileNameLPV)
          doLoadFolderDFList <- FALSE
        } else {
          base::print(base::paste0(sysTimePID(), " can't read RDS in loadResultDF() from ", fileNameLPV))
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
        base::print(base::paste0(sysTimePID(), " read from loadFolderDFList()."))
        listOfResultDF <- loadFolderDFList(session = session, folder = folder) #listOfResultDF <- loadFolderDFList(folder, globalVariables)
        if (base::length(listOfResultDF) != 0) {
          # read traitListName
          traitListName <-
            findInFiles(
              "traitListName:",
              folder,
              in_files = "\\config.yml$",
              recursive = FALSE
            ) #look into data folder
          if (base::length(traitListName) == 0 || traitListName == FALSE) {
            traitListName <-
              findInFiles(
                "traitListName:",
                paste0(folder, ".."),
                in_files = "\\config.yml$",
                recursive = FALSE
              ) #look one folder up
          }
          if (base::length(traitListName) == 0 || traitListName == FALSE) {
            base::message(base::paste0(sysTimePID(), " config.yml in data folder or
                                       traitListName not found. Consider placing
                                       a config.yml file with traitListName:
                                       <name of trait list> in each data folder."))
          }
          else {
            traitListName <-
              stringr::str_match(traitListName, "\"(.*?)\"")[2]
            if (base::length(traitListName) == 0) {
              base::message(base::paste0(sysTimePID(), " traitListName in config.yml not found."))
              traitListName <- folder
            }
          }
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
            base::message(base::paste0(sysTimePID(), " config.yml in data folder or
                                       PHENOFileNameLine not found. Consider placing
                                       a config.yml file with PHENOFileName:
                                       <location of PHENO file> in each data folder."))
          }
          else {
            PHENOFileName <-
              stringr::str_match(PHENOFileNameLine, "\"(.*?)\"")[2]
            if (base::length(PHENOFileName) == 0) {
              base::message(base::paste0(sysTimePID(), " PHENOFileName in PHENOFileNameLine not found."))
            }
            else {
              if (!file.exists(PHENOFileName)) {
                PHENOFileName <- paste0(folder, PHENOFileName)
              }
              listPrimaryKeys <- as.list(session$userData$config$keyAttributes)
              PHENOdata <- list(PHENODF = getPHENODF(PHENOFileName, listPrimaryKeys), PHENOFileName = PHENOFileName)
              listResultP_Val <- getResultDfP_D_N(listOfResultDF, "P")
              rn <- base::rownames(listResultP_Val)
              #add traitListName to colnames
              cn <- base::colnames(listResultP_Val)
              OriginalColnames <- cn
              cn <- paste0(traitListName, "_", cn)
              listResultP_Val <- base::as.data.frame(listResultP_Val)
              rownames(listResultP_Val) <- rn
              colnames(listResultP_Val) <- cn
              listResultDeltaMeth <- getResultDfP_D_N(listOfResultDF, "D")
              rn <- rownames(listResultDeltaMeth)
              listResultDeltaMeth <- base::as.data.frame(listResultDeltaMeth)
              rownames(listResultDeltaMeth) <- rn
              colnames(listResultDeltaMeth) <- cn
              listResultN <- getResultDfP_D_N(listOfResultDF, "N")
              rn <- rownames(listResultN)
              listResultN <- base::as.data.frame(listResultN)
              rownames(listResultN) <- rn
              colnames(listResultN) <- cn

              listResultP_Val_DeltaMeth_N <-
                list(
                  P_Val = listResultP_Val,
                  DM = listResultDeltaMeth,
                  N = listResultN,
                  OriginalColnames = OriginalColnames,
                  PHENOdata = PHENOdata,
                  traitListName = traitListName,
                  folder = folder
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
            listResultP_Val_DeltaMeth_N$folder <- folder
          }
          base::print(base::paste0(sysTimePID(), " saveRDS loadFolderDFList()")) # , fileNameLPV, "."))
          base::saveRDS(listResultP_Val_DeltaMeth_N, file = fileNameLPV)
        } else {
          base::message(base::paste0(sysTimePID(), " length(listOfResultDF) == 0."))
        }
      }
      if (loadRDS != FALSE) {
        if (session$userData$config$debugMode == TRUE) { #if (globalVariables$config$debugMode == TRUE) {
          # it is mandatory to reassign rownames and also colnames... Thanks R!
          #rn <- base::rownames(listResultP_Val_DeltaMeth_N$P_Val)[1:session$userData$sessionVariables$debugNumber]
          rn <- head(base::rownames(listResultP_Val_DeltaMeth_N$P_Val),session$userData$sessionVariables$debugNumber)
          cn <- base::colnames(listResultP_Val_DeltaMeth_N$P_Val)
          listResultP_Val_DeltaMeth_N$P_Val <-
            base::as.data.frame(head(listResultP_Val_DeltaMeth_N$P_Val,session$userData$sessionVariables$debugNumber))
          listResultP_Val_DeltaMeth_N$DM <-
            base::as.data.frame(head(listResultP_Val_DeltaMeth_N$DM,session$userData$sessionVariables$debugNumber))
          listResultP_Val_DeltaMeth_N$N <-
            base::as.data.frame(head(listResultP_Val_DeltaMeth_N$N,session$userData$sessionVariables$debugNumber))
          rownames(listResultP_Val_DeltaMeth_N$P_Val) <- rn
          rownames(listResultP_Val_DeltaMeth_N$DM) <- rn
          rownames(listResultP_Val_DeltaMeth_N$N) <- rn
          colnames(listResultP_Val_DeltaMeth_N$P_Val) <- cn
          colnames(listResultP_Val_DeltaMeth_N$DM) <- cn
          colnames(listResultP_Val_DeltaMeth_N$N) <- cn
          listResultP_Val_DeltaMeth_N$folder <- folder
        }
        result <- listResultP_Val_DeltaMeth_N
        base::print(base::paste0(sysTimePID(), " finished read from loadFolderDFList()", "."))
      } else {
        result <- FALSE
      }
    },
    error = function(e) {
      base::message("An error occurred in loadResultDF():\n", e)
    },
    warning = function(w) {
      base::message("A warning occurred in loadResultDF():\n", w)
    },
    finally = {
      base::print(base::paste0(sysTimePID(), " end loadResultDF()."))
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
  base::tryCatch(
    {
      # browser() #everything seems fine until here (all DFs become loaded)
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
          resultColnames <- base::colnames(traitDFs[[i]]$P_Val)
          resultOriginalColnames <- traitDFs[[i]]$OriginalColnames
          resultOriginDF <- base::rep(i, length(resultColnames))
          resultDFP_Val$Row.names <- rownames(traitDFs[[i]]$P_Val)
          resultDFDM$Row.names <- rownames(traitDFs[[i]]$DM)
          resultDFN$Row.names <- rownames(traitDFs[[i]]$N)
          resultFolder <- traitDFs[[i]]$folder
        } else {
          base::print(base::paste0(sysTimePID(), " merge trait ", i, "."))
          cn <- base::colnames(traitDFs[[i]]$P_Val)
          OriginalColnames <- traitDFs[[i]]$OriginalColnames
          OriginDF <- base::rep(i, length(cn))
          traitDFs[[i]]$P_Val$Row.names <-
            base::rownames(traitDFs[[i]]$P_Val)
          traitDFs[[i]]$DM$Row.names <- base::rownames(traitDFs[[i]]$DM)
          traitDFs[[i]]$N$Row.names <- base::rownames(traitDFs[[i]]$N)

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
          resultColnames <- base::c(resultColnames, cn)
          resultOriginalColnames <- base::c(resultOriginalColnames, OriginalColnames)
          resultFolder <- traitDFs[[i]]$folder
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
      result$resultOriginalColnames <- resultOriginalColnames
      result$folder <- resultFolder
    },
    error = function(e) {
      base::message("An error occurred in loadtraitDFs():\n", e)
    },
    warning = function(w) {
      base::message("A warning occurred in loadtraitDFs():\n", w)
    },
    finally = {
      base::print(base::paste0(sysTimePID(), " end loadtraitDFs()."))
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
  base::tryCatch(
    {
      base::print(base::paste0(sysTimePID(), " before list.files()"))
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
                base::print(base::paste0(sysTimePID(), " file ", i, " of ", base::length(temp)))
                base::print(base::paste0(sysTimePID(), " load result file ", folder, trait))
                # read results into DF
                resultDF <-
                  loadResultFile(session = session, folder = folder, fileName = trait) #loadResultFile(folder, trait, globalVariables)
                # omit unneccesary variables
                resultDF <-
                  resultDF[, c("probeID", "P_VAL", "DeltaMeth", "N")]
                base::tryCatch({
                  if (base::min(resultDF$P_VAL, na.rm = TRUE) <
                      base::as.numeric(session$userData$config$P_VALWarningThreshold)) {
                    base::message(
                      base::paste0(
                        sysTimePID(),
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
                  base::message("An error occurred in getlistOfResultsDF(), inner base::tryCatch:\n", e)
                },
                warning = function(w) {
                  base::message(
                    base::paste0(
                      sysTimePID(),
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
              #result <- FALSE
              base::message(
                base::paste0(
                  sysTimePID(),
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
      base::message("An error occurred in getlistOfResultsDF():\n", e)
    },
    warning = function(w) {
      base::message("A warning occurred in getlistOfResultsDF():\n", w)
    },
    finally = {
      base::print(base::paste0(sysTimePID(), " end getlistOfResultsDF()."))
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
  base::tryCatch(
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
      if (getList == TRUE) {
        #        listOfResultDF = getlistOfResultsDF(dataDir())
        base::print(base::paste0(sysTimePID(), " before getlistOfResultsDF()"))
        listOfResultDF <- getlistOfResultsDF(session = session, folder = folder) #listOfResultDF <- getlistOfResultsDF(folder, globalVariables)
        if (is.list(listOfResultDF)) { #if (listOfResultDF != FALSE) {
          base::print(base::paste0(sysTimePID(), " saveRDS getlistOfResultsDF()", fileNameLR))
          base::saveRDS(listOfResultDF, file = fileNameLR)
        }
      }
    },
    error = function(e) {
      base::message("An error occurred in loadFolderDFList():\n", e)
    },
    warning = function(w) {
      base::message("A warning occurred in loadFolderDFList():\n", w)
    },
    finally = {
      base::print(base::paste0(sysTimePID(), " end loadFolderDFList()."))
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
  base::tryCatch(
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
          base::message(base::paste0(sysTimePID(), " no key found in data. (from list:", as.character(listPrimaryKeys), "."))
        }
        PHENO <- as.data.frame(PHENO)
        result <- PHENO
      } else {
        base::message(base::paste0(sysTimePID(), " file not found: ", PHENOFileName))
        result <- FALSE
      }
    },
    error = function(e) {
      base::message("An error occurred in getPHENODF():\n", e)
    },
    warning = function(w) {
      base::message("A warning occurred in getPHENODF():\n", w)
    },
    finally = {
      base::print(base::paste0(sysTimePID(), " end getPHENODF()."))
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
  base::tryCatch(
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
            sysTimePID(),
            " counting remaining probes at p = ",
            P_VAL_BORDER,
            " remaining n = ",
            n
          )
        )
        result[i, 1] <- P_VAL_BORDER
        result[i, 2] <- n
      }
      colnames(result) <- base::c("maximum P_VAL_BORDER", "Available n")
      result <- result[1:i - 1, ]
      result <- base::as.data.frame((result))
      result <- result[base::order(result[1]), ]
    },
    error = function(e) {
      base::message("An error occurred in getAvailNForP_VALBorder():\n", e)
    },
    warning = function(w) {
      base::message("A warning occurred in getAvailNForP_VALBorder():\n", w)
    },
    finally = {
      base::print(base::paste0(sysTimePID(), " end getAvailNForP_VALBorder()."))
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
#' @description counts number of remaining models below a defined n based on p-value
#' @param mat matrix for which to calculate results
#' @param n minimum n for which to calculate results
#' @return data.frame with results
# examples getNForP_ValBorder(mat, n)
getNForP_ValBorder <- function(mat, n) {
  base::tryCatch({
    result <- base::matrix(nrow = 1, ncol = 2)
    P_VAL_BORDER <- 5 * 10^-n
    mat[mat > P_VAL_BORDER] <- NA
    mat <-
      delete.na(mat, ncol(mat) - 1) # -1, because we need at least 2 traits to associate
    nrow <- base::nrow(mat)
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
    base::message("An error occurred in getNForP_ValBorder():\n", e)
  },
  warning = function(w) {
    base::message("A warning occurred in getNForP_ValBorder():\n", w)
  },
  finally = {
    base::print(base::paste0(sysTimePID(), " end getNForP_ValBorder()."))
    return(result)
  }
  )
}

#' getNForDMBorder
#' @description counts number of remaining models below a defined n based on delta methylation
#' @param mat matrix for which to calculate results
#' @param n minimum n for which to calculate results
#' @return data.frame with results
# examples getNForDMBorder(mat, n)
getNForDMBorder <- function(mat, DMBorder) {
  base::tryCatch({
    result <- base::matrix(nrow = 1, ncol = 2)
    if (DMBorder > 0) {
      mat[mat < DMBorder] <- NA #mat[mat > DMBorder] <- NA
    }
    else if (DMBorder < 0) {
      mat[mat > DMBorder] <- NA
    }
    mat <-
      delete.na(mat, ncol(mat) - 1) # -1, because we need at least 2 traits to associate
    nrow <- base::nrow(mat)
    if (base::is.numeric(nrow)) {
      result[1, 1] <- DMBorder
      result[1, 2] <- nrow
    }
    else {
      result[1, 1] <- DMBorder
      result[1, 2] <- 0
    }
  },
  error = function(e) {
    base::message("An error occurred in getNForDMBorder():\n", e)
  },
  warning = function(w) {
    base::message("A warning occurred in getNForDMBorder():\n", w)
  },
  finally = {
    base::print(base::paste0(sysTimePID(), " end getNForDMBorder()."))
    return(result)
  }
  )
}

#' getNForNBorder
#' @description counts number of remaining models below a defined n based on n
#' @param mat matrix for which to calculate results
#' @param n minimum n for which to calculate results
#' @return data.frame with results
# examples getNForNBorder(mat, n)
getNForNBorder <- function(mat, NBorder) {
  base::tryCatch({
    result <- base::matrix(nrow = 1, ncol = 2)
    mat[mat > NBorder] <- NA
    mat <-
      delete.na(mat, ncol(mat) - 1) # -1, because we need at least 2 traits to associate
    nrow <- base::nrow(mat)
    if (base::is.numeric(nrow)) {
      result[1, 1] <- NBorder
      result[1, 2] <- nrow
    }
    else {
      result[1, 1] <- NBorder
      result[1, 2] <- 0
    }
  },
  error = function(e) {
    base::message("An error occurred in getNForNBorder():\n", e)
  },
  warning = function(w) {
    base::message("A warning occurred in getNForNBorder():\n", w)
  },
  finally = {
    base::print(base::paste0(sysTimePID(), " end getNForNBorder()."))
    return(result)
  }
  )
}

#' getAvailNForP_VALBorderParallel
#' counts number of remaining models below a defined n; uses parallel processing for faster results
#' @param session session
#' @param wd working directory
#' @param numCores number of CPU cores to use
#' @param DF data frame for which to calculate results
#' @return data.frame with results
# examples getAvailNForP_VALBorderParallel(session, wd, numCores, DF)
getAvailNForP_VALBorderParallel <- function(session, wd, numCores, DF) {
  base::tryCatch({
    base::print(base::paste0(sysTimePID(), " start getAvailNForP_VALBorderParallel()."))
    i <- NULL
    DF <- as.matrix(DF)
    #minP <- base::min(DF, na.rm = TRUE)
    minP <- base::apply(DF, 2, FUN = function(x) {base::min(x[x > 0], na.rm = TRUE)})
    minP <- base::min(minP)
    minP <- extractMantissaExponent(minP)$exponent
    if (minP > -1) {
      base::print(base::paste0(sysTimePID(), "Warning: minP > -1. Please check your data.")) #that should not be the case, please check data!
      browser()
    }
    maxP <- base::apply(DF, 2, FUN = function(x) {base::max(x[x > 0], na.rm = TRUE)})
    maxP <- base::max(maxP)
    maxP <- extractMantissaExponent(maxP)$exponent
    # if (maxP < 1) {
    #   base::print(base::paste0(sysTimePID(), "Warning: maxP < 1. Please check your data.")) #that should not be the case, please check data!
    #   browser()
    # }
    numRows <- maxP - minP
    shiny::updateSliderInput(session = session, inputId = "sldP_Val", min = minP, max = maxP, value = c(minP, maxP))
    result <- base::matrix(nrow = numRows, ncol = 2)
    #check size of exported global DF
    lengthDF <- length(DF)
    #limit <- 3000*1024^2 # for 3 GB
    #limit <- lengthDF*1024^2 # for real DF
    limit <- lengthDF
    limit <- limit * 1.1
    options(future.globals.maxSize = limit * 1024^2)
    #check, whether limit * numores exceeds memory limit of compute nodes...
    maxMemory <- limit * numCores
    #memorySize <- 512 * (1024)^2 # 512MB, hard coded, because memory.size() and memory.limit() are no longer supported from R 4.2 on...
    memorySize <- 5120 * (1024)^2 # 5120MB, hard coded
    multiple <- base::as.integer(memorySize/limit)
    if (multiple >= 1) {
      numCoresMemSize <- multiple
    }
    else {
      base::message(base::paste0(sysTimePID(), " size of DF is too big for computers memory: ", memorySize, "MB."))
      browser()
    }
    numCores <- base::min(numCores, numCoresMemSize)
    nWorkers <- parallelly::availableCores(constraints = "connections")
    numCores <- base::min(numCores,nWorkers)
    base::tryCatch(
      {
        future::plan(strategy = future::multisession, workers = numCores)
      },
      error = function(e) {
        warning("An error occurred in future::plan():\n", e)
        ##extract # of available connections
        e <- "Cannot create 112 parallel PSOCK nodes. Each node needs one connection, but there are only 75 connections left out of the maximum 128 available on this R installation"
        pattern <- "(Cannot create) ([[:digit:]]+) (parallel PSOCK nodes. Each node needs one connection, but there are only) ([[:digit:]]+) (connections left out of the maximum) ([[:digit:]]+) (available on this R installation)"
        # empty dataframe with columns to match after split
        proto <- data.frame(a = character(), numWanted = integer(), c = character(),
                            numLeft = integer(), e = character(), numMax = integer(), g = character())
        # extract
        numConnections <- utils::strcapture(pattern, e, proto)
        numCores <- numConnections$numLeft / 2
        if (numCores < 1) {
          numCores <- 1
        }
        session$userData$numCores <- numCores
        future::plan(strategy = future::multisession, workers = numCores)
      },
      warning = function(w) {
        warning("A warning occurred future::plan():\n", w)
      },
      finally = {
        base::print(base::paste0(sysTimePID(), " end error handler future::plan()."))
      }
    )
    library(future) #we have this already in DESCRIPTION file, but without "library(future)" here, it won't work. Strange.
    library(doFuture)
    result <- foreach::foreach(i = 1:numRows, .combine = rbind, #.packages = c("base"),
                               #.export = c("getNForP_ValBorder", "delete.na", "is.numeric", "nrow"),
                               .verbose = TRUE) %do% { #.verbose = TRUE) %dofuture% {
      #for some reason neither %dopar% nor %dofuture% find getNForDMBorder(), so the solution with source() is not very elegant, but works
      base::source(paste0(wd, "/R/ResultData.R")) #this is necessary for foreach %dopar% to run properly
      result <- getNForP_ValBorder(mat = DF, n = i)
      return(result)
    }
    colnames(result) <- base::c("maximum P_VAL_BORDER", "Available CpG")
    result <- base::as.data.frame((result))
    result <- result[base::order(result[[1]], decreasing = TRUE), ]
    },
    error = function(e) {
      warning("An error occurred in getAvailNForP_VALBorderParallel():\n", e)
    },
    warning = function(w) {
      warning("A warning occurred in getAvailNForP_VALBorderParallel():\n", w)
    },
    finally = {
      base::print(base::paste0(sysTimePID(), " end getAvailNForP_VALBorderParallel()."))
      return(result)
    }
  )
}

getAvailNForDMBorderParallel <- function(session, wd, numCores, DF) {
  base::tryCatch({
    base::print(base::paste0(sysTimePID(), " start getAvailNForP_VALBorderParallel()."))
    result <- NULL
    i <- NULL
    #check min DM and max DM
    #minDM <- base::round(base::min(DF, na.rm = TRUE), 5)
    DF <- as.matrix(DF)
    minDM <- base::apply(DF, 2, FUN = function(x) {base::min(x[x > 0], na.rm = TRUE)})
    minDM <- base::min(minDM)
    if (minDM < 0) {
      base::print(base::paste0(sysTimePID(), "Warning: minDM < 0. Please check your data.")) #that should not be the case, please check data!
      browser() #that should not be the case, please check data!
    }
    #maxDM <- base::round(base::max(DF, na.rm = TRUE), 5)
    maxDM <- base::apply(DF, 2, FUN = function(x) {base::max(x[x > 0], na.rm = TRUE)})
    maxDM <- base::max(maxDM)
    if (maxDM > 1) {
      base::print(base::paste0(sysTimePID(), "Warning: maxDM > 1. Please check your data.")) #that should not be the case, please check data!
      browser() #that should not be the case, please check data!
    }
    #shiny::updateSliderInput(session = session, inputId = "sldDM", min = minDM, max = maxDM, value = c(minDM, maxDM), step = NULL)
    shiny::updateSliderInput(session = session, inputId = "sldDM", min = minDM, max = maxDM, value = c(minDM, maxDM), step = 0.001)
    minDM <- as.integer(minDM * 100) #0
    maxDM <- as.integer(maxDM * 100)
    numRows <- maxDM - minDM
    result <- base::matrix(nrow = numRows, ncol = 2)
    #check size of exported global DF
    lengthDF <- length(DF)
    #limit <- 3000*1024^2 # for 3 GB
    #limit <- lengthDF*1024^2 # for real DF
    limit <- lengthDF
    limit <- limit * 1.1
    options(future.globals.maxSize = limit * 1024^2)
    #check, whether limit * numores exceeds memory limit of compute nodes...
    maxMemory <- limit * numCores
    #memorySize <- 512 * (1024)^2 # 512MB, hard coded, because memory.size() and memory.limit() are no longer supported from R 4.2 on...
    memorySize <- 5120 * (1024)^2 # 5120MB, hard coded
    multiple <- base::as.integer(memorySize/limit)
    if (multiple >= 1) {
      numCoresMemSize <- multiple
    }
    else {
      base::message(base::paste0(sysTimePID(), " size of DF is too big for computers memory: ", memorySize, "MB."))
      browser()
    }
    numCores <- base::min(numCores, numCoresMemSize)
    nWorkers <- parallelly::availableCores(constraints = "connections")
    numCores <- base::min(numCores,nWorkers)
    base::tryCatch(
      {
        future::plan(strategy = future::multisession, workers = numCores)
      },
      error = function(e) {
        browser()
        base::message("An error occurred in future::plan():\n", e)
        ##extract # of available connections
        e <- "Cannot create 112 parallel PSOCK nodes. Each node needs one connection, but there are only 75 connections left out of the maximum 128 available on this R installation"
        pattern <- "(Cannot create) ([[:digit:]]+) (parallel PSOCK nodes. Each node needs one connection, but there are only) ([[:digit:]]+) (connections left out of the maximum) ([[:digit:]]+) (available on this R installation)"
        # empty dataframe with columns to match after split
        proto <- data.frame(a = character(), numWanted = integer(), c = character(),
                            numLeft = integer(), e = character(), numMax = integer(), g = character())
        # extract
        numConnections <- utils::strcapture(pattern, e, proto)
        numCores <- numConnections$numLeft / 2
        if (numCores < 1) {
          numCores <- 1
        }
        session$userData$numCores <- numCores
        future::plan(strategy = future::multisession, workers = numCores)
      },
      warning = function(w) {
        base::message("A warning occurred future::plan():\n", w)
      },
      finally = {
        base::print(base::paste0(sysTimePID(), " end error handler future::plan()."))
      }
    )
    library(future) #we have this already in DESCRIPTION file, but without "library(future)" here, it won't work. Strange.
    library(doFuture)
    result <- foreach::foreach(i = minDM:maxDM, .combine = rbind,
                               .verbose = TRUE) %do% { # .verbose = TRUE) %dofuture% {
      #for some reason neither %dopar% nor %dofuture% find getNForDMBorder(), so the solution with source() is not very elegant, but works
      base::source(paste0(wd, "/R/ResultData.R")) #this is necessary for foreach %dopar% to run properly
      result <- getNForDMBorder(mat = DF, DMBorder = i/100)
      return(result)
    }
    colnames(result) <- base::c("DM_BORDER", "Available CpG")
    result <- base::as.data.frame((result))
    result <- result[base::order(result[[1]], decreasing = TRUE), ]
  },
  error = function(e) {
    base::message("An error occurred in getAvailNForDMBorderParallel():\n", e)
  },
  warning = function(w) {
    base::message("A warning occurred in getAvailNForDMBorderParallel():\n", w)
  },
  finally = {
    base::print(base::paste0(sysTimePID(), " end getAvailNForDMBorderParallel()."))
    return(result)
  })
}

getAvailNForNBorderParallel <- function(session, wd, numCores, DF) {
  base::tryCatch({
    base::print(base::paste0(sysTimePID(), " start getAvailNForP_VALBorderParallel()."))
    result <- NULL
    i <- NULL
    DF <- as.matrix(DF)
    minN <- base::apply(DF, 2, FUN = function(x) {base::min(as.integer(x[x > 0]), na.rm = TRUE)})
    minN <- base::min(minN)
    if (minN < 1) {minN <- 1}
    if (minN != as.integer(minN)) {
      base::print(base::paste0(sysTimePID(), "Warning: minN != as.integer(minN). Please check your data.")) #that should not be the case, please check data!
      browser()
    }
    maxN <- base::apply(DF, 2, FUN = function(x) {base::max(as.integer(x[x > 0]), na.rm = TRUE)})
    maxN <- base::max(maxN)
    numRows <- maxN - minN
    if (maxN != as.integer(maxN)) {
      base::print(base::paste0(sysTimePID(), "Warning: maxN != as.integer(maxN). Please check your data.")) #that should not be the case, please check data!
      browser()
    }
    if (maxN < 1) {
      base::print(base::paste0(sysTimePID(), "Warning: maxN < 1. Please check your data.")) #that should not be the case, please check data!
      browser()
    }
    shiny::updateSliderInput(session = session, inputId = "sldN", min = minN, max = maxN, value = c(minN, maxN))
    result <- base::matrix(nrow = numRows, ncol = 2)
    #check size of exported global DF
    lengthDF <- length(DF)
    #limit <- 3000*1024^2 # for 3 GB
    #limit <- lengthDF*1024^2 # for real DF
    limit <- lengthDF
    limit <- limit * 1.1
    options(future.globals.maxSize = limit * 1024^2)
    #check, whether limit * numores exceeds memory limit of compute nodes...
    maxMemory <- limit * numCores
    #memorySize <- 512 * (1024)^2 # 512MB, hard coded, because memory.size() and memory.limit() are no longer supported from R 4.2 on...
    memorySize <- 5120 * (1024)^2 # 5120MB, hard coded
    multiple <- base::as.integer(memorySize/limit)
    if (multiple >= 1) {
      numCoresMemSize <- multiple
    }
    else {
      base::print(base::paste0(sysTimePID(), " size of DF is too big for computers memory: ", memorySize, "MB."))
      browser()
    }
    numCores <- base::min(numCores, numCoresMemSize)
    nWorkers <- parallelly::availableCores(constraints = "connections")
    numCores <- base::min(numCores,nWorkers)
    base::tryCatch(
      {
        future::plan(strategy = future::multisession, workers = numCores)
      },
      error = function(e) {
        browser()
        base::message("An error occurred in future::plan():\n", e)
        ##extract # of available connections
        e <- "Cannot create 112 parallel PSOCK nodes. Each node needs one connection, but there are only 75 connections left out of the maximum 128 available on this R installation"
        pattern <- "(Cannot create) ([[:digit:]]+) (parallel PSOCK nodes. Each node needs one connection, but there are only) ([[:digit:]]+) (connections left out of the maximum) ([[:digit:]]+) (available on this R installation)"
        # empty dataframe with columns to match after split
        proto <- data.frame(a = character(), numWanted = integer(), c = character(),
                            numLeft = integer(), e = character(), numMax = integer(), g = character())
        # extract
        numConnections <- utils::strcapture(pattern, e, proto)
        numCores <- numConnections$numLeft / 2
        if (numCores < 1) {
          numCores <- 1
        }
        session$userData$numCores <- numCores
        future::plan(strategy = future::multisession, workers = numCores)
      },
      warning = function(w) {
        base::message("A warning occurred future::plan():\n", w)
      },
      finally = {
        base::print(base::paste0(sysTimePID(), " end error handler future::plan()."))
      }
    )
    library(future) #we have this already in DESCRIPTION file, but without "library(future)" here, it won't work. Strange.
    library(doFuture)
    result <- foreach::foreach(i = minN:maxN, .combine = rbind, #.packages = c("base"),
                               #.export = c("getNForNBorder", "delete.na", "is.numeric", "nrow"),
                               .verbose = TRUE) %do% { #.verbose = TRUE) %dofuture% {
      #for some reason neither %dopar% nor %dofuture% find getNForDMBorder(), so the solution with source() is not very elegant, but works
      base::source(paste0(wd, "/R/ResultData.R")) #this is necessary for foreach %dopar% to run properly
      result <- getNForNBorder(mat = DF, NBorder = i)
      return(result)
    }
    colnames(result) <- base::c("N_BORDER", "Available CpG")
    result <- base::as.data.frame((result))
    result <- result[base::order(result[[1]], decreasing = TRUE), ]
    },
  error = function(e) {
    base::message("An error occurred in getAvailNForNBorderParallel():\n", e)
  },
  warning = function(w) {
    base::message("A warning occurred in getAvailNForNBorderParallel():\n", w)
  },
  finally = {
    base::print(base::paste0(sysTimePID(), " end getAvailNForNBorderParallel()."))
    return(result)
  })
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
    base::tryCatch(
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
        base::message("An error occurred  in getReducedP_Valdf():\n", e)
      },
      warning = function(w) {
        base::message("A warning occurred in getReducedP_Valdf():\n", w)
      },
      finally = {
        base::print(base::paste0(sysTimePID(), " end getReducedP_Valdf()."))
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
  base::tryCatch(
    {
      if (exists("minN")) {
        dfN <- dfList$dfN
        rn <- base::rownames(dfN)
        cn <- base::colnames(dfN)
        dfP_Val <- dfList$dfP_Val
        dfDM <- dfList$dfDM
        resultOriginDF <- dfList$resultOriginDF
        resultColnames <- dfList$resultColnames
        resultOriginalColnames <- dfList$resultOriginalColnames
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
        cn <- colnames(dfN)

        # select the same content than in dfN
        dfP_Val <- as.data.frame(dfP_Val[base::rownames(dfN), base::colnames(dfN)]) # does not work with data.table
        rownames(dfP_Val) <- rn
        colnames(dfP_Val) <- cn
        dfDM <- as.data.frame(dfDM[base::rownames(dfN), base::colnames(dfN)])
        rownames(dfDM) <- rn
        colnames(dfDM) <- cn
        resultOriginDF <- resultOriginDF[positions]
        resultColnames <- resultColnames[positions]
        resultOriginalColnames <- resultOriginalColnames[positions]
        dfList <- base::list(
          dfP_Val = NULL,
          dfDM = NULL,
          dfN = NULL,
          resultOriginDF = NULL,
          resultColnames = NULL,
          resultOriginalColnames = NULL,
          listPHENOdata = NULL
        )
        dfList$dfP_Val <- dfP_Val
        dfList$dfDM <- dfDM
        dfList$dfN <- dfN
        dfList$resultOriginDF <- resultOriginDF
        dfList$resultColnames <- resultColnames
        dfList$resultOriginalColnames <- resultOriginalColnames
        dfList$listPHENOdata <- listPHENOdata
      }
      else {
        base::message(base::paste0(sysTimePID(), " minN does not exist."))
      }
    },
    error = function(e) {
      base::message("An error occurred in removeTraitsMinN():\n", e)
    },
    warning = function(w) {
      base::message("A warning occurred in removeTraitsMinN():\n", w)
    },
    finally = {
      base::print(base::paste0(sysTimePID(), " end removeTraitsMinN()."))
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
            sysTimePID(),
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
        base::message(
          sysTimePID(),
          "warning: length(listResultP_Val)[1] == 0"
        )
      )
      result <- FALSE
    }
  }
  return(result)
}
