library(bsseq)
smooth <- HDF5Array::loadHDF5SummarizedExperiment()

pData <- pData(smooth)
pData$col <- c(rep("red", 9), rep("blue", 9), rep("green", 8))
pData(smooth) <- pData

dmr <- <dmr of interest> # output of callDMR function (DSS-Package)

pdf(file = "dmrs_top50.pdf", width = 10, height = 5)
plotManyRegions(smooth, dmr[1:50,], extend = 2000, 
                addRegions = dmrs)
dev.off()