library("dplyr")
library("enrichplot")
library("clusterProfiler")
library("pathview")
library("ReactomePA")

setwd("/media/rt/code_data/Codes/rnaseq/TCGA_breast")
gene_full <- read.csv("significant_results.csv",header=TRUE)
gene <- gene_full %>%
  filter(log2FoldChange >0) %>%  # select upregulated genes in metaplastic breast cancer
  dplyr::select(entrezgene_id)

#Gene ontology analysis

GO<-enrichGO( gene$entrezgene_id,
              OrgDb = 'org.Hs.eg.db',
              keyType = "ENTREZID",
              ont = "MF",  #one or more in  c("BP", "CC", "MF"), or "ALL"
              pvalueCutoff = 1, # this p value indicates the difference compared to background.
              qvalueCutoff = 1,# This q value is also in comparison with background
              readable = TRUE) # True means converting ids to gene names

# Enrichment plot
png(filename="GO_enrichment.png", width = 1500, height = 1000, pointsize = 30)
cnetplot(GO,circular=FALSE,color.params = list(edge = TRUE))
dev.off()


#KEGG
KEGG_database <- 'hsa'
KEGG<-enrichKEGG(gene$entrezgene_id,
                 organism = KEGG_database,
                 pvalueCutoff = 1,
                 qvalueCutoff = 1
)
png(filename="KEGG_enrichment.png",width = 1500, height = 1000, pointsize = 30)
cnetplot(KEGG, showCategory = 10, circular=FALSE,color.params = list(edge = TRUE))
dev.off()
browseKEGG(KEGG, "hsa04110")

# Reactome
PA <- enrichPathway(gene=gene$entrezgene_id, pvalueCutoff = 1, readable=TRUE)
png(filename="Reactome_enrichment.png",width = 1500, height = 1000, pointsize = 30)
cnetplot(PA, showCategory = 5, circular=FALSE,colorEdge = TRUE)
dev.off()