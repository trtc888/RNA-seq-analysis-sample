# TCGA breast cancer RNA seq data analysis: an example procedure 
Date: 2025-03-10
## Purpuse
This document is to show a sample procedure to analyze TCGA RNA-seq results.
Metaplastic breast cancer is a rare and aggresive cancer type with poor prognosis. Understanding its gene expression pattern is helpful to deciper the behavior of this cancer type. This procedure intends to compare the transcriptomes of metaplastic breast cancer with the infiltrating duct carcinoma from surgical samples. Results can be further aligned with experimenental findings and modelings.
## Step 1: download RNA-seq results from TCGA
1. goto TCGA [database](https://portal.gdc.cancer.gov/).
2. Open Projects, and check breast in primary sites.
3. Find Project TCGA-BRCA, which is the collection of breast cancer, then save the cohort as Breast.
4. Close the Project view, and open the Repository view.
5. Check RNA-Seq and STAR-Counts in the left filter panel.
6. Add all files to cart, and download cart by selecting Manifest. Note, TCGA does not allow downloading files directly larger than 5 GB, therefore, we must use the manifest file along with the [DGDC Data Transfer Tool](https://gdc.cancer.gov/access-data/gdc-data-transfer-tool) to download. In the code below, -m indicates using the manifest file followed by the file name in the same folder.
```
$ ./dgc-client download -m gdc_manifest.[date].[session].txt
```
7. Download associated clinical data, biospecimen data, and metadata.
8. Downloaded count files are saved as .tsv in separate subfolders. Use the Python code below to consolidate them into a single folder.

```Python
# Python code
import os
import shutil

def gather_tsv_files(source_dir, target_dir):
    
    # detect if the directories are exist
        if not os.path.exists(target_dir):
        os.makedirs(target_dir)

    # Parse all subfolders
    for root, dirs, files in os.walk(source_dir):
        for file in files:
            if file.endswith(".tsv"):
                
                # obtain full directories
                file_path = os.path.join(root, file)
                
                # Make full paths
                target_file_path = os.path.join(target_dir, file)
                
                # Change extension to avoid file replacement
                if os.path.exists(target_file_path):
                    base, extension = os.path.splitext(file)
                    counter = 1
                    while os.path.exists(target_file_path):
                        target_file_path = os.path.join(target_dir, f"{base}_{counter}{extension}")
                        counter += 1
                        
                # Copy files
                shutil.copy2(file_path, target_file_path)
                print(f"Copied {file_path} to {target_file_path}")

if __name__ == "__main__":
    
    # define source and target directories
    source_directory = "/Path/To/TCGA_breast/gdc_data"
    target_directory = "/Path/To/TCGA_breast/TCGA_counts"

    # Run the code
    gather_tsv_files(source_directory, target_directory)


```
## Step 2: Organize data
In this procedure, I am interested in the transcriptomic differences between metaplastic breast cancer and infiltrating ductal breast cancer. Therefore, only these two types of samples will be selectd and organized for further analysis.
Only surgical samples were analyzed to avoid potential interference of therapeutics.
To save time, only 18 samples from each group were picked up for the analysis. Related codes can be removed for full analysis.
```r
library("dplyr")
library("stringr")
library("jsonlite")

#Function to extract case_id
extract_file <- function(x) {
  file_id <- x  %>% as.data.frame %>% select(case_id) %>% as.character
  return(file_id)
}

# Prepare samples
setwd("/media/rt/code_data/Codes/rnaseq/TCGA_breast")
sample_info <- read.csv("./TCGA_info/clinical.cart.2025-03-10/clinical.tsv",sep="\t",header=TRUE)

## Select specific cancer types for the analysis
sample_select <- sample_info[sample_info$primary_diagnosis %in% c("Infiltrating duct carcinoma, NOS","Lobular carcinoma, NOS","Metaplastic carcinoma, NOS"),]

## Generate sample information for DESeqDataSetFromMatrix() function
sample_groups <- sample_select %>% select(case_id,primary_diagnosis, treatment_type) %>% filter(treatment_type == "Surgery, NOS")

## Obtain data file names stored in the folder after Python combining
#file_names <- list.files(path = "./TCGA_counts",full.names = TRUE)
file_names <- list.files(path = "./TCGA_counts")

## Read json metadata
json_data <- fromJSON("./TCGA_info/metadata.cart.2025-03-10.json")

#Extract case_id from associated entities to entry
json_data$case_id <- sapply(json_data$associated_entities, extract_file)
json_data_use<- json_data %>% filter(endsWith(file_name, ".tsv")) %>% select(case_id,file_name)

#Make colData for DESeqDataSetFromMatrix function
sample_groups <-merge(sample_groups, json_data_use, by="case_id")

sample_groups_final <- sample_groups %>% 
  mutate(is_file_exist = ifelse(file_name %in% file_names, sample_groups$file_name, NA)) %>%
  filter(!is.na(is_file_exist)) %>%
  select(-is_file_exist) %>%
  group_by(primary_diagnosis) %>%     #Randomly select 18 samples from each diagnosis group
  slice_sample(n=18)             #Used 18 because there are 18 samples in metaplastic group

```
## Step 3: Analysis
Analysis of the results is straightforward. Here, the DESeq2 package in Bioconductor is used.
For first time use, install the package following the [tutorial](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("DESeq2")
```
Since the required data have been prepared, simply call the *DESeqDataSetFromMatrix()* and *DESeq()* functions will do the work. Note that *DESeqDataSetFromMatrix()* function requires direct count results, which is the "unstranded" column in the downloaded files.

```r
library("DESeq2")
dds1 <- DESeqDataSetFromMatrix(countData = merged_table,
                              colData = sample_groups_final,
                              design= ~ primary_diagnosis)
dds2 <- DESeq(dds1)
res <- results(dds2, contrast = c("primary_diagnosis", "Metaplastic carcinoma, NOS", "Infiltrating duct carcinoma, NOS"))

## Re-order the results and select rows with adjusted p value less than 0.05
res_significant <- as.data.frame(subset(res[order(res$padj),], padj < 0.05)) 
```
In the *results()* function, arg *contrast* is the way how the comparison should be analyzed. 

# Step 4: Finalize the data comparison and save
In this step, package $biomaRt$ will be used to retrieve the gene names and descriptions. The tutorial of this pacckage is [here](https://bioconductor.org/packages/release/bioc/html/biomaRt.html).
After gene info are retrieved, they can be merged with the result table. Below is an example. After that, genes expressions with adjusted p value less than 0.05 are extracted, annotated, sorted and saved for further analysis and plotting
.
```r
# Search gene information if needed.
library(biomaRt)
mart <- useMart('ENSEMBL_MART_ENSEMBL')
mart <- useDataset('hsapiens_gene_ensembl', mart)

## Note that the gene_id was stored as "ENSG00000053108.17" in res file. 
## Only the first part before "." will be used to search the gene information.
values = str_split_i(row.names(res), "\\.", 1) 
lookup <- getBM(
  mart = mart,
  attributes = c('entrezgene_id', 'ensembl_gene_id',
                 'gene_biotype','hgnc_symbol',  
                 "definition_1006"), # Adding definition may result in multiplication of result rows
  filter = 'ensembl_gene_id',
  values = values,
  uniqueRows = TRUE)

# listAttributes(mart)  # Show attributes of mart

res_significant$ensembl_gene_id <- str_split_i(row.names(res_significant), "\\.", 1)
res_significant <- merge(lookup, res_significant, by = "ensembl_gene_id", all.x = FALSE, all.y = TRUE)
res_significant <- res_significant[order(res_significant$log2FoldChange),]

write.csv(res_significant, "significant_results.csv")
```

## Distribution and Citation
This document and related codes can be freely distributed and reused. Cite the following reference if you want:

Rui Tang, Aixiang Ding, Marvin Rivera, Eben Alsberg. Modeling breast cancer tissue in vitro using extracted native collagen fibers. [Cancer Res., 15 February 2021; 81 (4_Supplement): PS17–53.](https://doi.org/10.1158/1538-7445.SABCS20-PS17-53)

