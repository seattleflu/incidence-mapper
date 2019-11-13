library(testthat)
library('dbViewR')

context("test dbViewR packages")

shp <- masterSpatialDB(shape_level = 'census_tract', source = 'simulated_data', rm_files = TRUE)

test_that("masterSpatialDB output vars", {
  expect_setequal(colnames(shp), c("STATEFP", "COUNTYF", "TRACTCE", "AFFGEOI", "GEOID", "NAME", "LSAD", "ALAND", "AWATER", "rowID", "CRA_NAM", "NEIGHBO", "PUMA5CE", "geometry", "residence_census_tract", "work_census_tract", "residence_cra_name", "work_cra_name", "residence_neighborhood_district_name", "work_neighborhood_district_name", "residence_puma", "work_puma", "domain"))
})

test_that("masterSpatialDB output rows", {
  expect_equal(nrow(shp), 397)
})

files <- list.files(getwd(), pattern="2016_CensusTracts_KingCountyWa")

test_that("masterSpatialDB remove_local_files", {
  expect_equal(length(files), 0)
})

queryIn <- list(
    SELECT   =list(COLUMN=c('site_type','residence_census_tract')),
    WHERE    =list(COLUMN='site_type', IN = c('kiosk')),
    GROUP_BY =list(COLUMN=c('site_type','residence_census_tract')),
    SUMMARIZE=list(COLUMN='site_type', IN= c('kiosk'))
)
db <- expandDB( selectFromDB(  queryIn, source = 'simulated_data' ) )

test_that("expandDB output vars", {
  expect_setequal(colnames(db$observedData), c("residence_census_tract", "site_type", "pathogen", "n", "positive" ))
})

test_that("expandDB output rows", {
  expect_equal(nrow(db$observedData), 397)
})

