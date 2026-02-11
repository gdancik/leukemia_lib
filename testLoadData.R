# Example to load expression for a gene
source("functions/leukemia_lib/getData.R")

# get expression values selected gene across all
# datasets
GENE <- "ALDH1A1"
datasets <- get_all_data(GENE)

# names of the datasets
names(datasets)

# View expression of gene in TCGA
View(datasets$TCGA$X)

# View risk information for TCGA
View(datasets$TCGA$Y)


