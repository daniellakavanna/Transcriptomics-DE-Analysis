# Transcriptomics-DE-Analysis

## Exploratory Data Analysis for Transcriptomic Data

R-studio was used to analyse a gene expression data set collected using DNA Affymetrix microarrays.
The dataset was sourced from the Gene Expression Omnibus; a genomics data repository.
The work was performed using R version 4.1.3 and R-studio. The dataset was imported into R using the Bioconductor’s GEOquery package. For background, data was collected using DNA microarrays from the analysis of blood from patients with acute Dengue virus (DENV) infection and during convalescence. The GEOquery package contains the *getGEO* function which allows seamless download and parsing of a GEO SOFT format file into an R data structure. The transcriptomics dataset was in the form of single channel microarray Affymetrix data.

Differential expression (DE) analysis was carried out on the full Dengue gene expression dataset to show relationships between the gene expression profiles and distinguish the genes that are differentially regulated and show the most significant difference in expression between the clinical groups. Differentially expressed genes (DEG) were revealed by specifying a linear model, creating a design matrix and implementing the Empirical Bayes statistic. The Limma package was used to implement the Empirical Bayes statistic to share information across all the genes.

A linear model was derived using the function *model.matrix* to create a design matrix for every gene for each of the population groups (Phipson et al., 2016). The coefficients of the model were fitted by the function *lmFit* (Phipson et al., 2016). All pair-wise comparisons between the four population groups were made (Phipson et al., 2016). The *topTable* function provides the summary statistics of the linear model for the gene expression profiles (Phipson et al., 2016). The summary statistics include the (log2) fold changes which infers the quantitative value of the contrast. A volcano plot, a type of scatter plot was used to visualise and interpret the the most significant difference of individual genes.

## Volcano Plot for DE-analysis 
![image](https://user-images.githubusercontent.com/93345220/197392084-f3cb9677-d731-42e0-aae4-59d10958fe84.png)


Volcano plot where the (log2) fold change cut off is 0.6 and p-value cut off is 0.05 to visualise DEGs. The position of each gene is determined by its p-value and (log2) fold change. The x-axis is the (log2) fold change and y-axis is the (log10) transformed p-value. The upper right and left hand corners show the most statistical significant genes, marked as red circles. Black dots indicate the remaining genes that were not statistically significant.
