library(pathview)

setwd("/Users/willtrim/Documents/projs/bedrest/outputs/kegg_pathviews/DE")

#### data load
deseq2_results <- read.csv("/Users/willtrim/Documents/projs/bedrest/outputs/star_rsem_tmp_post_pre.csv", row.names = 1)
deseq2_results$gene_ids <- sapply(strsplit(rownames(deseq2_results), "\\."), head, 1)
rownames(deseq2_results) <- deseq2_results$gene_ids
deseq2_stats <- deseq2_results[, "stat"]
names(deseq2_stats) <- deseq2_results$gene_ids



### KEGG pathways
library(clusterProfiler)

kegg_pathways_df <- download_KEGG(species = "hsa", keggType = "KEGG", keyType = "kegg")
View(kegg_pathways_df[["KEGGPATHID2NAME"]])
library(data.table)
setwd("/Users/willtrim/Documents/projs/bedrest/outputs/kegg_pathviews")

fwrite(kegg_pathways_df[["KEGGPATHID2NAME"]], "kegg_pathways.csv",
       row.names = FALSE)

####
library(pathview)
library(tidyverse)

out_suffix <- "star_rsem_tmp_post_pre"

limit <- floor(max(abs(deseq2_stats)))

draw_pathway <- function(pathway_name) {
  
  pathway_id <- str_replace_all(kegg_pathways_df[["KEGGPATHID2NAME"]][kegg_pathways_df[["KEGGPATHID2NAME"]]$to == pathway_name, "from"], pattern="hsa", replacement="")
  out_suffix <- str_replace_all(paste(str_replace_all(pathway_name, pattern = " ", replacement = "_"), "star_rsem_tmp_post_pre", collapse="_", sep = "."), pattern="/", replacement="_")
  pv.out <- pathview(gene.data = deseq2_stats, pathway.id = pathway_id, gene.idtype = "ENSEMBL", limit = list(gene = limit, cpd = 1),
                     species = "hsa", out.suffix = out_suffix,  low = list(gene = "blue", cpd = "green"), mid =
                       list(gene = "gray", cpd = "gray"), high = list(gene = "red", cpd =
                                                                        "yellow"))
}
  
####
## Insulin signaling pathway, glucose transport, AKT, PTPN1

pathways <- c(
  "Insulin signaling pathway",
  "Insulin resistance",
  "Type II diabetes mellitus",
  "AGE-RAGE signaling pathway in diabetic complications",
  "Glucagon signaling pathway",

####
## TAG storage/Lipolysis, FOXO, PPARg, G3P synthesis
  "Regulation of lipolysis in adipocytes",
  "Lipoic acid metabolism",
  "Cholesterol metabolism",
  "Lipid and atherosclerosis",
  "Non-alcoholic fatty liver disease",
  "Glycerolipid metabolism",

  "Fatty acid biosynthesis",
  "Fatty acid elongation",
  "Fatty acid degradation",
#   "Fatty acid metabolism",

  "ABC transporters",
  "Thermogenesis",

####
## PPP, NADPH, ROS
  "Pentose phosphate pathway",
  "Glutathione metabolism",

  "Ferroptosis", # (hsa04216): Iron-dependent oxidative cell death.
  "Apoptosis", #(hsa04210): Programmed cell death pathways triggered by stress.
  "Mitophagy - animal",
  "Autophagy - animal", # (hsa04137, hsa04140): Removal of damaged mitochondria and cellular components.
  "Cellular senescence", # (hsa04114): Age-related stress response.
  "HIF-1 signaling pathway", # (hsa04066): Responds to hypoxic stress.
  "Chemical carcinogenesis - reactive oxygen species",

####
## TCA, glycolysis
  "Glycolysis / Gluconeogenesis",
  "Citrate cycle (TCA cycle)",
  "Pyruvate metabolism",
#   "2-Oxocarboxylic acid metabolism",
####
## Mitochondria and b-oxidation
  "Oxidative phosphorylation",
####
## adipokines
  "Adipocytokine signaling pathway",

  "Apelin signaling pathway",
####
## peroxisomes and fatty acid oxidation
  "Glyoxylate and dicarboxylate metabolism",
  "Peroxisome",
####
## ribosomes and mTOR?
  "mTOR signaling pathway",
  "Ribosome biogenesis in eukaryotes",
  "Ribosome",
  

  "PPAR signaling pathway",
  "FoxO signaling pathway",



  "Spliceosome",
  "Proteasome",
  "Protein export",
  "Cell cycle",
  "Ubiquitin mediated proteolysis",
  "Protein processing in endoplasmic reticulum",
  "Lysosome",
  "Endocytosis",
  "Phagosome",



  "MAPK signaling pathway",
  "ErbB signaling pathway",
  "Ras signaling pathway",
  "Rap1 signaling pathway",
  "cGMP-PKG signaling pathway",
  "cAMP signaling pathway",
  "NF-kappa B signaling pathway",
  "AMPK signaling pathway",
  "Wnt signaling pathway",
  "Notch signaling pathway",
  "Hedgehog signaling pathway",
  "TGF-beta signaling pathway",
  "Hippo signaling pathway",
  "JAK-STAT signaling pathway",
  "TNF signaling pathway",
  "Toll-like receptor signaling pathway",
  "B cell receptor signaling pathway",
  "Fc gamma R-mediated phagocytosis",
  "Antigen processing and presentation",
  "Cytokine-cytokine receptor interaction",
  "Leukocyte transendothelial migration",
  "Natural killer cell mediated cytotoxicity",
  "Cell adhesion molecule (CAM) interaction",
  "T cell receptor signaling pathway"
)

for (pathway in pathways) {
  draw_pathway(pathway)
}
draw_pathway("Chemical carcinogenesis - reactive oxygen species")
