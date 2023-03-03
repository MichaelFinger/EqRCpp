equatingRecipesExamplesFoldername = "~/Developer/EquatingRecipes/ER for distribution 9-10-09/Examples/"

datasetFilenames = c("ConvRawSS2.in",
                     "Conversion-reordered.DAT",
                     "KB6Vitems",
                     "KB6Xdist",
                     "KB6Xdist.ST",
                     "KB6Xitems",
                     "KB6Xitems.ST",
                     "KB6Ydist",
                     "KB6Yitems",
                     "actmathfreq.dat",
                     "dummyV1items",
                     "dummyV3items",
                     "dummyXdist",
                     "dummyXdist.ST",
                     "dummyXitems",
                     "dummyXitems.ST",
                     "dummyYdist",
                     "dummyYitems",
                     "dummyv2items",
                     "fdexample.dat",
                     "mondatx.dat",
                     "mondaty.dat",
                     "yctmath.TXT")

index = 1
readLines(paste0(equatingRecipesExamplesFoldername, datasetFilenames[index]))
