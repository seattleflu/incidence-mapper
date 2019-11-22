library(testthat)


library('dbViewR')
library('incidenceMapR')
library('modelVisualizeR')
library('modelServR')


context("local packages content")

ls("package:dbViewR")
test_that("package:dbViewR content", {
  expect_setequal(ls("package:dbViewR"), c("addCensusData", "expandDB", "loadBingData", "masterSpatialDB", "selectFromDB"))
})

ls("package:incidenceMapR")
test_that("package:incidenceMapR content", {
  expect_setequal(ls("package:incidenceMapR"), c("appendCatchmentModel", "constructAdjacencyNetwork", "effectsModel", "fluVaxEfficacyModel", "forecast_ILI", "latentFieldModel", "modelTrainR", "plotforecasts" , "smoothModel" ))
})

ls("package:modelServR")
test_that("pakcage:modelServR content", {
  expect_setequal(ls("package:modelServR"), c("getHumanReadableModelIdFromModel", "getHumanReadableModelIdFromQuery", "getModelIdFromModel", "getModelIdFromQuery", "getModelQueryObjectFromModel", "getModelQueryObjectFromQuery", "loadModelFileById", "queryLoadedModel", "queryModelById", "returnModel", "saveModel" ))
})

ls("package:modelVisualizeR")
test_that("pakcage:modelVisualizeR content", {
  expect_setequal(ls("package:modelVisualizeR"), c("ggplotDiagnosticPlots", "ggplotFixedEffects", "ggplotLatentMap", "ggplotSmoothEffects", "ggplotSmoothMap", "ggplotSmoothSequential"))
})

