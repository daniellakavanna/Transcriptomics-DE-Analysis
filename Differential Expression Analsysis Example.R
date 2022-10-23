# Differential Expression Analaysis Example 
# Install required packages 
#To install Bioconductor:
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()

#  install the GEOquery package:
BiocManager::install(c("GEOquery"))

library(limma)
library(dplyr)
library(pheatmap)
library(ggplot2)
library(ggrepel)

rm(list = ls()) # clear the R environment
library(GEOquery) # load GEOquery

# load the aforementioned DENV expression dataset with accession number GDS5093 
# using the getGEO() function:
gse <- getGEO("GDS5093", GSEMatrix = TRUE)  
# Make a dataframe to get the expression objects using Table()
X <- Table(gse)

# Convert GSE to an ExpressionSet object
eset <- GDS2eSet(gse, do.log2=TRUE)

# The expression data come from microarray experiments

# Get the phenotypic state of each sample (labels)
pDat <- pData(eset)
fdat <- fData(eset)

geneNames <- as.character(X$IDENTIFIER) 
X <- exprs(eset) 
rownames(X) <- geneNames

X <- avereps(X)

dim(X)

pDat <- select(pDat, group = disease.state, patient = individual)
design <- model.matrix(~0+pDat$group)

design

colnames(design) <- c("Convalescent", "DengueHemorrhagicFever", "DengueFever","healthycontrol")


#calculate median expression level
cutoff <- median(X)

is_expressed <- X > cutoff


model <- lmFit(X, design)

head(model$coefficients)

contrasts.matrix <- makeContrasts(Convalescent - DengueHemorrhagicFever,
                                  DengueHemorrhagicFever - DengueFever,
                                  DengueFever - healthycontrol,
                                  Convalescent - DengueFever,
                                  Convalescent - healthycontrol,
                                  DengueHemorrhagicFever - healthycontrol,
                                  levels=design)

# Obtain differential expression statistics and p-values.
fit2 <- contrasts.fit(model, contrasts.matrix)
fit2 <- eBayes(fit2)
topTable(fit2)

# The outcome of each hypothesis test
outcome <- decideTests(fit2)


topTable(fit2, coef=6)

results <- topTable(fit2, number = Inf, coef = 6,  adjust="BH")

# Volcano Plot Visualisations for DE analysis 
BiocManager::install('EnhancedVolcano')
library(EnhancedVolcano)
png(paste(format(Sys.time(),"%Y-%m-%d"), "volcano.png", sep = ". "), width = 10, height = 10, units = 'in', res = 300)
setwd("~/Downloads")
EnhancedVolcano(results,
                lab = rownames(results),
                x = 'logFC',
                y = 'P.Value',
                pCutoff = 0.05,
                FCcutoff = 0.6,
                pointSize = 2.0,
                xlim = c(-1,1.3),
                ylim = c(0,25),
                labSize = 3.0)
dev.off()
