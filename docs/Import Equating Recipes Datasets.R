library(jsonlite)
library(stringr)

equatingRecipesExamplesFoldername = "~/Developer/EquatingRecipes/ER for distribution 9-10-09/Examples/"
equatingRecipesExamplesJsonFoldername = "~/Developer/EquatingRecipes/ER for distribution 9-10-09/Examples/JSON/"

importFixedColumnsTextFile = function(filename) {
  linesRead = trimws(readLines(filename))
  
  values = str_split(linesRead[1], pattern = " +|\t", simplify = TRUE)
  
  if (ncol(values) == 1) {
    numberOfRows = as.numeric(linesRead[1])
    
    linesRead = linesRead[2:(numberOfRows + 1)] 
  }
  
  values = str_split(linesRead, pattern = " +|\t", simplify = TRUE)
  
  columnsToKeep = colSums(values == "") < nrow(values)
  rowsToKeep = rowSums(values == "") < ncol(values)
  
  values = values[rowsToKeep, columnsToKeep]
  
  dataset = NULL
  for (columnIndex in 1:ncol(values)) {
    if (sum(is.na(as.numeric(values[, columnIndex]))) == 0) {
      if (columnIndex == 1) {
        dataset = data.frame(V1 = as.numeric(values[, columnIndex]))
      } else {
        dataset = cbind(dataset,
                        data.frame(V1 = as.numeric(values[, columnIndex])))  
      }
      
    } else {
      if (columnIndex == 1) {
        dataset = data.frame(V1 = values[, columnIndex])
      } else {
        dataset = cbind(dataset,
                   data.frame(V1 = values[, columnIndex]))
      }
    }
    
    colnames(dataset)[columnIndex] = paste0("V", columnIndex)
  }
  
  return(dataset)
}

datasetFilenames = c("actmathfreq.dat",
                     "Conversion-reordered.DAT",
                     "ConvRawSS2.in",
                     "dummyV1items",
                     "dummyv2items",
                     "dummyV3items",
                     "dummyXdist.ST",
                     "dummyXdist",
                     "dummyXitems.ST",
                     "dummyXitems",
                     "dummyYdist",
                     "dummyYitems",
                     "fdexample.dat",
                     "KB6Vitems",
                     "KB6Xdist.ST",
                     "KB6Xdist",
                     "KB6Xitems.ST",
                     "KB6Xitems",
                     "KB6Ydist",
                     "KB6Yitems",
                     "yctmath.TXT")

for (index in 1:length(datasetFilenames)) {
  dataset = importFixedColumnsTextFile(paste0(equatingRecipesExamplesFoldername, datasetFilenames[index]))
 
  json = toJSON(as.data.frame(dataset), dataframe = "values",matrix = "rowmajor", pretty = TRUE) 
  
  jsonFilename = paste0(equatingRecipesExamplesJsonFoldername, datasetFilenames[index], ".json")
  
  write(json, jsonFilename, append = FALSE)
}

dataset = read.fortran(paste0(equatingRecipesExamplesFoldername, "mondatx.dat"), 
                       format = c("36I1", "2I2"))
json = toJSON(as.data.frame(dataset), dataframe = "values",matrix = "rowmajor", pretty = TRUE) 
jsonFilename = paste0(equatingRecipesExamplesJsonFoldername, "mondatx.dat", ".json")
write(json, jsonFilename, append = FALSE)

dataset = read.fortran(paste0(equatingRecipesExamplesFoldername, "mondaty.dat"), 
                       format = c("36I1", "2I2"))
json = toJSON(as.data.frame(dataset), dataframe = "values",matrix = "rowmajor", pretty = TRUE) 
jsonFilename = paste0(equatingRecipesExamplesJsonFoldername, "mondaty.dat", ".json")
write(json, jsonFilename, append = FALSE)

