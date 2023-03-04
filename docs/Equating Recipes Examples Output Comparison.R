library(jsonlite)
library(stringr)

originalExamplesPath = "~/Developer/EquatingRecipes/EqRCpp/tests/examples/resources/json/original_results/"
cppExamplesPath = "~/Developer/EquatingRecipes/EqRCpp/build/tests/examples/"

importCh2 = function(filename) {}

importCh3 = function(filename) {}

compareCh4 = function() {
  originalResults = fromJSON(paste0(originalExamplesPath, "chapter4.json"))
  
  cppResults = fromJSON(paste0(cppExamplesPath, "chapter4.json"))
  
  originalResults$summary$msx
  
  cppResults$EquatedRawScoreResults$equatedRawScoreMoments
}



