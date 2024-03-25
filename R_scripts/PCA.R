#genomic PCA (needs data from Angsd/PCAngsd)
library(ggplot2)
LTer <- as.matrix(read.table("LumTer.cov"))
individuals <- read.table("indlist")

ind <- paste0(individuals[,1])
pop <- read.table("pop")
LTerEigen <- eigen(LTer)
sum_eigenvalues <- sum(LTerEigen$values)
explained_variance <- (LTerEigen$values / sum_eigenvalues) * 100
print(explained_variance)

theme_set(theme_bw())

ggplot(LTer, aes(x = V1, y = V2, colour = samples)) + 
  geom_point() +
  labs(title = "Angsd-based PCA, min 50% ind, maf 0.05", 
       x = "PC1 (14,7%)", 
       y = "PC2 (4,7%)") +
  theme(legend.position = "none") +
  theme(plot.title = element_text(size = 16))

#epigenetic PCA (loads smoothed data usding HDF5Array)
library(DSS)
library(ggplot2)

smooth <- HDF5Array::loadHDF5SummarizedExperiment()
chr <- c("OX457036.1", "OX457037.1", "OX457038.1", "OX457039.1", "OX457040.1",
        "OX457041.1", "OX457042.1", "OX457043.1", "OX457044.1", "OX457045.1",
        "OX457046.1", "OX457047.1", "OX457048.1", "OX457049.1", "OX457050.1",
        "OX457051.1", "OX457052.1", "OX457053.1", "OX457054.1")


Z <- getMeth(chrSelectBSseq(smooth, seqnames = chr), type = "smooth", what = "perBase")
t<-as.data.frame(Z)
cova <- cov(t)
individuals <- read.table("indlist_epi")
ind <- paste0(individuals[,1])
pop <- read.table("popfull_epi")
LTerEigen <- eigen(cova)
LTerEigen$values[1:5]
sum_eigenvalues <- sum(LTerEigen$values)
explained_variance <- (LTerEigen$values / sum_eigenvalues) * 100
print(explained_variance)
Conditions <- paste0(pop[,1])

theme_set(theme_bw())

ggplot(df, aes(x = V1, y = V2, colour = Conditions)) + 
  geom_point(size = 3) +
  labs(title = "Smoothed epigenetic PCA", 
       x = "PC1 (80.5%)", 
       y = "PC2 (1.2%)" ) +
#  theme(legend.position = "none") +
  theme(legend.text = element_text(size = 14)) +  
  theme(plot.title = element_text(size = 16)) +
  theme(legend.position=c(0.82,0.22))