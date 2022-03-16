#rm(list = ls(all.names = TRUE))
#options(expression = 500000)
#start "Y:\Home\roeder\RStudio\bin\rstudio.exe" --max-ppsize=500000
#libDir = "Y:/Home/roeder/R-4.1.0.3/library"
#appDir = "Y:/Home/roeder/daten/Forschungsplanung/EWAS_HM/"

#dataDir = "F:/roeder/EWAS_DM_development/all/"
#"F:/roeder/EWAS_DrugMetabolites/all/",
#"Z:/EWAS_MVOC/all/",
#"Z:/EWAS_VOC/all/",
# dataDir1 = c("Z:/EWAS_Biocrates_1LJ/all/",
#                "Z:/EWAS_MercapturicAcid/all/",
#                "Z:/EWAS_Parabene/all/",
#                "Z:/EWAS_Phthalate/all/",
#                "Z:/EWAS_Stress/",
#                "Y:/Home/roeder/daten/Forschungsplanung/EWAS_Zytokine/",
#                "Z:/EWAS_VitD/",
#                "Y:/Home/roeder/daten/Forschungsplanung/EWAS_IgE/",
#                "Z:/EWAS_Test/all/",
#                "Z:/EWAS_Test2/all/",
#                 "Z:/EWAS_Test3/all/")

# dataDir2 = c("Z:/EWAS_Biocrates_1LJ/all/",
#              "Z:/EWAS_MercapturicAcid/all/",
#              "Z:/EWAS_Parabene/all/",
#              "Z:/EWAS_Phthalate/all/",
#              "Z:/EWAS_Stress/",
#              "Y:/Home/roeder/daten/Forschungsplanung/EWAS_Zytokine/",
#              "Z:/EWAS_VitD/",
#              "Y:/Home/roeder/daten/Forschungsplanung/EWAS_IgE/",
#              "Z:/EWAS_Test/all/",
#              "Z:/EWAS_Test2/all/",
#              "Z:/EWAS_Test3/all/")
# dataDir3 = c("Z:/EWAS_Outcomes/all/",
#              "F:/roeder/EWAS_OutcomesLong/all/")
#dataDir3 = "Z:/EWAS_Outcomes/all/"
#dataDir3 = "F:/roeder/EWAS_OutcomesLong/all/"
#workDir = "Y:/Home/roeder/daten/Forschungsplanung/EWAS_DrugMetabolites_OutcomeShort/"
#library(compiler)
#setwd(appDir)
#setwd(dataDir)
#source(paste0(appDir,"util.R"))

#installLibraries()
#loadLibraries()
#loadObjects()

#setwd(appDir)
#setwd(workDir)
#setwd(dataDir)

#debugMode=FALSE #TRUE
#debugMode=TRUE

#source(paste0(appDir,"util.R"))
#source(paste0(appDir,"ui.R"))
#source(paste0(appDir,"server.R"))
#shiny::shinyApp(ui,server)
