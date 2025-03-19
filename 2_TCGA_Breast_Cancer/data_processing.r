# Prepare samples

library("dplyr")
library("stringr")
library("jsonlite")

#Function to extract case_id
extract_file <- function(x) {
  file_id <- x  %>% as.data.frame %>% dplyr::select(case_id) %>% as.character
  return(file_id)
}

setwd("/media/rt/code_data/Codes/rnaseq/TCGA_breast")
sample_info <- read.csv("./TCGA_info/clinical.cart.2025-03-10/clinical.tsv",sep="\t",header=TRUE)

## Select specific cancer types for the analysis
sample_select <- sample_info[sample_info$primary_diagnosis %in% c("Infiltrating duct carcinoma, NOS","Metaplastic carcinoma, NOS"),]

## Generate sample information for DESeqDataSetFromMatrix() function
sample_groups <- sample_select %>% dplyr::select(case_id,primary_diagnosis, treatment_type) %>% filter(treatment_type == "Surgery, NOS")

## Obtain data file names stored in the folder after Python combining
#file_names <- list.files(path = "./TCGA_counts",full.names = TRUE)
file_names <- list.files(path = "./TCGA_counts",pattern="*.tsv")

## Read json metadata
json_data <- fromJSON("./TCGA_info/metadata.cart.2025-03-10.json")

#Extract case_id from associated entities to entry
json_data$case_id <- sapply(json_data$associated_entities, extract_file)
json_data_use<- json_data %>% filter(endsWith(file_name, ".tsv")) %>% dplyr::select(case_id,file_name)

#Make colData for DESeqDataSetFromMatrix function
sample_groups <-merge(sample_groups, json_data_use, by="case_id")

sample_groups_final <- sample_groups %>% 
  mutate(is_file_exist = ifelse(file_name %in% file_names, sample_groups$file_name, NA)) %>%
  filter(!is.na(is_file_exist)) %>%
  dplyr::select(-is_file_exist) %>%
  group_by(file_name) %>%   #Remove duplicatte file names
  filter(row_number()==1) %>%
  group_by(primary_diagnosis) %>%      #Randomly select 18 samples from each diagnosis group
  slice_sample(n=18)     #Used 18 because there are 18 samples in metaplastic group

## Create count data for DESeqDataSetFromMatrix
merged_table <- NULL
for (file in sample_groups_final$file_name) {
  full_path = paste("./TCGA_counts/", file, sep="")
  print(file)
  file_content <- read.csv(full_path, header=TRUE, skip = 1, sep="\t")
  new_name <- sample_groups_final$case_id[sample_groups_final$file_name==file]
  merge_content <-file_content %>%
    dplyr::select(gene_id,unstranded) %>%
    #    rename(unstranded = sample_groups_final$case_id[sample_groups_final$file_name==file])
    rename_with(~new_name, "unstranded")
  if (is.null(merged_table)){
    merged_table <- merge_content
  }
  else{
    merged_table <- merge(merged_table, merge_content, by.x="gene_id",by.y="gene_id", all=TRUE)
  }
}
merged_table<- merged_table %>%
  filter(startsWith(gene_id,"EN"))
row.names(merged_table) <- merged_table$gene_id
merged_table <- merged_table[,c(2:ncol(merged_table))]

# Analysis

library("DESeq2")

dds1 <- DESeqDataSetFromMatrix(countData = merged_table,
                               colData = sample_groups_final,
                               design= ~ primary_diagnosis)
dds2 <- DESeq(dds1)
res <- results(dds2, contrast = c("primary_diagnosis", "Metaplastic carcinoma, NOS", "Infiltrating duct carcinoma, NOS"))

## Re-order the results and select rows with adjusted p value less than 0.05
res_significant <- as.data.frame(subset(res[order(res$padj),], padj < 0.05)) 

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
                 'gene_biotype','hgnc_symbol'),  
  #                "definition_1006"), # Adding definition may result in multiplication of result rows
  filter = 'ensembl_gene_id',
  values = values,
  uniqueRows = TRUE)

# listAttributes(mart)  # Show attributes of mart

res_significant$ensembl_gene_id <- str_split_i(row.names(res_significant), "\\.", 1)
res_significant <- merge(lookup, res_significant, by = "ensembl_gene_id", all.x = FALSE, all.y = TRUE)
res_significant <- res_significant[order(res_significant$log2FoldChange),]

write.csv(res_significant, "significant_results.csv")

# Draw volcano plot
library("EnhancedVolcano")
res_merge <- as.data.frame(res)
res_merge$ensembl_gene_id <- str_split_i(row.names(res_merge), "\\.", 1)
res_merge <- merge(lookup, res_merge, by = "ensembl_gene_id", all.x = FALSE, all.y = TRUE)
jpeg(filename = "volcano.jpg", width = 1500, height = 1000, pointsize = 20)
EnhancedVolcano(res_merge,
                lab = res_merge$hgnc_symbol,
                x = 'log2FoldChange',
                y = 'padj')
dev.off()


#Draw heatmap
library("pheatmap")
res_significant$ensembl_gene_id <-row.names(res_significant)
merged_table$ensembl_gene_id <- row.names(merged_table)
res_significant <- merge(res_significant, merged_table, by = "ensembl_gene_id", all.x = TRUE, all.y = FALSE)
row.names(res_significant) <- res_significant$ensembl_gene_id
res_heatmap <- res_significant[,c(8:ncol(res_significant))]
jpeg(filename = "heatmap.jpg", width = 1500, height = 1000)
pheatmap(res_heatmap, scale="row",color=colorRampPalette(c("blue","white","red"))(100))
dev.off()