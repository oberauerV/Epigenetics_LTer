library(DSS)
require(bsseq)
#loading all samples
#control 9
c_1 = read.table(file.path("samples/c12-1.tsv"), header = TRUE)
c_10 = read.table(file.path("samples/c12-10.tsv"), header = TRUE)
c_11 = read.table(file.path("samples/c12-11.tsv"), header = TRUE)
c_2 = read.table(file.path("samples/c12-2.tsv"), header = TRUE)
c_3 = read.table(file.path("samples/c12-3.tsv"), header = TRUE)
c_5 = read.table(file.path("samples/c12-5.tsv"), header = TRUE)
c_6 = read.table(file.path("samples/c12-6.tsv"), header = TRUE)
c_7 = read.table(file.path("samples/c12-7.tsv"), header = TRUE)
c_9 = read.table(file.path("samples/c12-9.tsv"), header = TRUE)
#cd10 9
l_1 = read.table(file.path("samples/10cd12-1.tsv"), header = TRUE)
l_10 = read.table(file.path("samples/10cd12-10.tsv"), header = TRUE)
l_2 = read.table(file.path("samples/10cd12-2.tsv"), header = TRUE)
l_3 = read.table(file.path("samples/10cd12-3.tsv"), header = TRUE)
l_4 = read.table(file.path("samples/10cd12-4.tsv"), header = TRUE)
l_5 = read.table(file.path("samples/10cd12-5.tsv"), header = TRUE)
l_6 = read.table(file.path("samples/10cd12-6.tsv"), header = TRUE)
l_7 = read.table(file.path("samples/10cd12-7.tsv"), header = TRUE)
l_9 = read.table(file.path("samples/10cd12-9.tsv"), header = TRUE)
#cd25 8
h_1 = read.table(file.path("samples/25cd12-1.tsv"), header = TRUE)
h_11 = read.table(file.path("samples/25cd12-11.tsv"), header = TRUE)
h_2 = read.table(file.path("samples/25cd12-2.tsv"), header = TRUE)
h_3 = read.table(file.path("samples/25cd12-3.tsv"), header = TRUE)
h_5 = read.table(file.path("samples/25cd12-5.tsv"), header = TRUE)
h_6 = read.table(file.path("samples/25cd12-6.tsv"), header = TRUE)
h_7 = read.table(file.path("samples/25cd12-7.tsv"), header = TRUE)
h_9 = read.table(file.path("samples/25cd12-9.tsv"), header = TRUE)


BS_c_10 = makeBSseqData(list(c_1,c_2,c_3,c_5,c_6,c_7,c_9,c_10,c_11,l_1,l_2,l_3,l_4,l_5,l_6,l_7,l_9,l_10), c("C1","C2","C3","C4","C5","C6","C7","C8","C9","L1","L2","L3","L4","L5","L6","L7","L8","L9"))

BS_c_25 = makeBSseqData(list(c_1,c_2,c_3,c_5,c_6,c_7,c_9,c_10,c_11,h_1,h_11,h_2,h_3,h_5,h_6,h_7,h_9), c("C1","C2","C3","C4","C5","C6","C7","C8","C9","H1","H2","H3","H4","H5","H6","H7","H8"))

#DML test and DMR
dmltest_c_10 = DMLtest(BS_c_10, group1=c("C1","C2","C3","C4","C5","C6","C7","C8","C9"), group2=c("L1","L2","L3","L4","L5","L6","L7","L8","L9"), smoothing = TRUE, smoothing.span = 300)
dmls_c_10 = callDML(dmltest_c_10)
dmrs_c_10 <- callDMR(dmltest_c_10, p.threshold = 0.0001)

dmltest_c_25 = DMLtest(BS_c_25, group1=c("C1","C2","C3","C4","C5","C6","C7","C8","C9"),  group2=c("H1","H2","H3","H4","H5","H6","H7","H8"), smoothing = TRUE, smoothing.span = 300)
dmls_c_25 = callDML(dmltest_c_25)
dmrs_c_25 <- callDMR(dmltest_c_25, p.threshold = 0.0001)

#only for Methylation-plots and PCA
BS_all = makeBSseqData(list(c_1,c_2,c_3,c_5,c_6,c_7,c_9,c_10,c_11,l_1,l_2,l_3,l_4,l_5,l_6,l_7,l_9,l_10,h_1,h_11,h_2,h_3,h_5,h_6,h_7,h_9), c("C1","C2","C3","C4","C5","C6","C7","C8","C9","L1","L2","L3","L4","L5","L6","L7","L8","L9","H1","H2","H3","H4","H5","H6","H7","H8"))

BSHDF <- realize(BSobj, "HDF5Array")
smooth <- BSmooth(BSHDF,
                    BPPARAM = MulticoreParam(workers = 4))
saveRDS(smooth, file = "smooth_all.RDS")
HDF5Array::saveHDF5SummarizedExperiment(smooth)