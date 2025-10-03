if (!requireNamespace("BiocManager", quietly = TRUE)) 
  install.packages("BiocManager")

BiocManager::install(c("GEOquery","affy","arrayQualityMetrics"))
install.packages("dplyr")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GEOquery")
# load GEOquery
library(GEOquery)

# load the dataset (returns a list if multiple platforms exist)
gse_list <- getGEO("GSE16558", GSEMatrix = TRUE)

names(gse_list)
# "GSE16558-GPL6244_series_matrix.txt.gz" "GSE16558-GPL8695_series_matrix.txt.gz"

# access each platform separately
gse6244 <- gse_list[[1]]   # or gse_list[["GPL6244"]] if names match
gse8695 <- gse_list[[2]]   # or gse_list[["GPL8695"]]

# get expression data
expr6244 <- exprs(gse6244)
expr8695 <- exprs(gse8695)

# get feature data (genes/probes)
feature6244 <- fData(gse6244)
feature8695 <- fData(gse8695)

# get phenotype/metadata (samples, disease vs control, etc.)
pheno6244 <- pData(gse6244)
pheno8695 <- pData(gse8695)

# check missing values in sample annotation
sum(is.na(pheno6244$source_name_ch1))
sum(is.na(pheno8695$source_name_ch1))

# raw data required full preprocessing (e.g., RMA normalization, QC)

# CEL files are large, and downloads may fail even with a good connection. 
#It's recommended to download raw data directly from NCBI GEO

# skip this step if you already downloaded data from NCBI

# fetch GEO supplementry files
dir.create("raw_data/GSE16558", recursive = TRUE, showWarnings = FALSE)

options(timeout = 1200)  # increase timeout (seconds)

url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE16nnn/GSE16558/suppl/GSE16558_RAW.tar"
dest <- "raw_data/GSE16558/GSE16558_RAW.tar"

# resume if partially downloaded (libcurl supports -C -)
download.file(url,
              destfile = dest,
              method   = "libcurl",
              extra    = "-C -",
              mode     = "wb")

# extract the archive
untar(dest, exdir = "raw_data/GSE16558/CEL")

# read CEL files into R as an AffyBatch object
BiocManager::install("affy")
library(affy)

# verify the folder path
getwd()   # shows your current working directory
list.files("raw_data/GSE16558/CEL")   # shows files in that folder

raw_data <- ReadAffy(celfile.path = "raw_data/GSE16558/CEL")

# inspect raw data
raw_data

# install arrayQualityMetrics package
BiocManager::install("arrayQualityMetrics")
library(arrayQualityMetrics)

# QC identifies outlier arrays, hybridization problems, or technical biases.
# arrayQualityMetrics: # this package generates automated QC reports for microarray data.
arrayQualityMetrics(expressionset = raw_data,
                    outdir = "results/QC_Raw_Data",
                    force = TRUE,
                    do.logtransform = TRUE)

# RMA (robust multi-array average) normalization is a popular method for normalizing Affymetrix microarray data
normalized_data <- rma(raw_data)

# QC after data normalization 
arrayQualityMetrics(expressionset = normalized_data,
                    outdir = "Results/QC_Normalized_Data",
                    force = TRUE)

# extract normalized expression values into a data frame
processed_data <- as.data.frame(exprs(normalized_data))

dim(processed_data)   # dimensions: number of probes × number of samples

# filter low-variance transcripts (“soft” intensity based filtering)
# filtering removes probes with low or uninformative expression signals

# calculate median intensity per probe across samples
row_median <- rowMedians(as.matrix(processed_data))

# visualize distribution of probe median intensities
par(mfrow = c(1,1), mar = c(5,4,2,1), oma = c(0,0,0,0))
hist(row_median,
     breaks = 100,
     freq = FALSE,
     main = "median intensity distribution")

# set a threshold to remove low variance probes (dataset-specific, adjust accordingly)
threshold <- 3
abline(v = threshold, col = "black", lwd = 2) 

# select probes above threshold
indx <- row_median > threshold 
filtered_data_1 <- processed_data[indx, ]
filtered_data_2 <- processed_data[indx, ]

# rename filtered expression data with sample metadata
colnames(filtered_data_1) <- rownames(pheno6244)
colnames(filtered_data_2) <- rownames(pheno8695)

# overwrite processed data with filtered dataset
processed_data <- filtered_data_1

# phenotype data preparation
class(pheno6244$source_name_ch1)
class(pheno8695$source_name_ch1)

colnames(pheno6244)

head(pheno6244$characteristics_ch1, 10)
head(pheno6244$characteristics_ch1.1, 11)
head(pheno6244$characteristics_ch1.2, 12)

table(pheno6244$source_name_ch1)

# define experimental groups (normal vs cancer)
vals <- trimws(as.character(pheno6244$source_name_ch1))
groups_1 <- factor(vals,
                   levels = c("Bone marrow, healthy", "Bone marrow, MM"),
                   labels = c("control", "myeloma"))
table(groups_1, useNA = "ifany")

class(groups_1)
levels(groups_1)

table(groups_1)