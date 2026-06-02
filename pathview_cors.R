if(!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("pathview")
BiocManager::install("tidyverse")


library(pathview)

setwd("/Users/willtrim/Documents/projs/bedrest/outputs/kegg_pathviews")

##### example
data(gse16873.d)
pv.out <- pathview(gene.data = gse16873.d[, 1:2], pathway.id = "04110",
                   species = "hsa", out.suffix = "gse16873")

#### data load
rs <- read.csv("/Users/willtrim/Documents/projs/bedrest/outputs/cors2/rs.csv", row.names = 1)

### KEGG pathways
library(clusterProfiler)

kegg_pathways_df <- download_KEGG(species = "hsa", keggType = "KEGG", keyType = "kegg")
View(kegg_pathways_df[["KEGGPATHID2NAME"]])

# 
# kegg_pathways_list <- kegg_pathways_df["KEGGPATHID2EXTID"][[1]] %>%
#   group_by(from) %>%
#   group_split()
# 
# library(purrr)
# kegg_pathways <- map(kegg_pathways_list, ~pull(.x, to))
# kegg_pathway_ids <-  map(kegg_pathways_list, ~ .x[[1, 1]])
# kegg_pathways_names_map <- kegg_pathways_df[["KEGGPATHID2NAME"]]$to
# names(kegg_pathways_names_map) <- kegg_pathways_df[["KEGGPATHID2NAME"]]$from
# kegg_pathways_names <- lapply(kegg_pathway_ids, function(x) {kegg_pathways_names_map[[x]]})
# names(kegg_pathways) <- kegg_pathways_names

library(pathview)
library(tidyverse)

setwd("/Users/willtrim/Documents/projs/bedrest/outputs/kegg_pathviews/cors")


draw_pathway_cors <- function(pathway_name, phens) {
  
  limit <- 1

  pathway_id <- str_replace_all(kegg_pathways_df[["KEGGPATHID2NAME"]][kegg_pathways_df[["KEGGPATHID2NAME"]]$to == pathway_name, "from"], pattern="hsa", replacement="")
  # print(pathway_name)
  # print(pathway_id)
  out_suffix <- str_replace_all(paste(str_replace_all(pathway_name, pattern = " ", replacement = "_"), paste(phens, collapse="_"), collapse="_", sep = "."), pattern="/", replacement="_")
  
  pv.out <- pathview(
    gene.data = rs[, phens], pathway.id = pathway_id, gene.idtype = "ENSEMBL", limit = list(gene = limit, cpd = 1),
    species = "hsa", out.suffix = out_suffix,  low = list(gene = "blue", cpd = "green"), mid =
    list(gene = "gray", cpd = "gray"), high = list(gene = "red", cpd = "yellow"))
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
  "Fat digestion and absorption",
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
  
  
  ####
  ## TCA, glycolysis
  "Glycolysis / Gluconeogenesis",
  "Citrate cycle (TCA cycle)",
  "Pyruvate metabolism",
  #   "2-Oxocarboxylic acid metabolism",
  "Carbohydrate digestion and absorption",
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
  "Protein digestion and absorption",
  
  
  
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

phen_names <- c("insulin_art", "insulin_ven", "FFA", "TAG", "fat_ox", "cho_ox", "glucose_disposal", "glucose_art", "glucose_ven")

for (pathway in pathways) {
  # print(pathway)
  for (phen_name in phen_names) {
    pre_post_phen_names <- c(paste(phen_name, "pre", sep = "_"), paste(phen_name, "post", sep="_"))
    # print(pre_post_phen_names)
    draw_pathway_cors(pathway, pre_post_phen_names)
  }
}

for (phen_name in phen_names) {
  pre_post_phen_names <- c(paste(phen_name, "pre", sep = "_"), paste(phen_name, "post", sep="_"))
  # print(pre_post_phen_names)
  draw_pathway_cors("", pre_post_phen_names)
}
