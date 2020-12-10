library(Chicago)
library(PCHiCdata)
vignette("Chicago")

testDesignDir <- file.path("~/Project/Data/Des")
dir(testDesignDir)
testDataPath <- file.path("~/Project/Data/Final")
dir(testDataPath)
files <- file.path(testDataPath, "5.chinput")
settingsFile <- file.path(system.file("extdata", package="PCHiCdata"),
                          "sGM12878Settings", "sGM12878.settingsFile")
cd <- setExperiment(designDir = testDesignDir, settingsFile = settingsFile)
cd <- readAndMerge(files=files, cd=cd)
cd <- chicagoPipeline(cd)