#All required packages 
library(GOSemSim)
library(AnnotationHub)
library(GO.db)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggarchery)
library(enrichplot)
library(stats)
library(DOSE)
library(ReactomePA)
library(tidyr)
library(dplyr)
library(tibble)
library(pheatmap)
library(ComplexHeatmap)
library(ggtree)
library(stringr)
library(ggplot2)
library(ggrepel)
library(scales)
library(cowplot)
library(patchwork)
library(igraph)
library(ggraph)
library(tidygraph)
library(corrplot)
library(RColorBrewer)

#setup the working directory that has all your input files and would have your outputs stored
setwd("E:/AAT/OneDrive/1 - ONGOING PUBLICATION/EN PREPARATION/4- Espourteille et al-2023/WGCNA/PSP_PICK_WGCNA")

########### PREPROCESSING ###################################################

#load the datafile contaning the proteomic data with sample details
data_mito_PSP = read.csv("PiD_fold_change_GSEA.csv")

#load background list (usually all the proteins names in your data)
#background = read.csv("Background_list.csv")


#Converting the gene symbol ids to entrez ids
entrez_mapping <- bitr(data_mito_PSP$Proteins,
                       fromType = "SYMBOL",
                       toType = "ENTREZID",
                       OrgDb = org.Hs.eg.db)

# Merge the Entrez IDs back into your dataframe
df_merged <- merge(data_mito_PSP, entrez_mapping, by.x = "Proteins", by.y = "SYMBOL", all.x = TRUE)

## assume 1st column is ID
## 2nd column is FC (Fold change)

## feature 1: numeric vector: FC
geneList = df_merged$Fold_change

## feature 2: named vector: entrez ID
names(geneList) = as.charcater(df_merged$ENTREZID)

## feature 3: decreasing order
geneList = sort(geneList, decreasing = TRUE)

#remove NA's 
geneList <- geneList[!is.na(names(geneList))]


########### GO ENRICHMENT ###################################################
#GO classification first; initiate with BP  
#level 3 is better than level 2 because we are not getting the specific terms

#Biological Pathways
ego_bp <- gseGO(geneList   = geneList,
                OrgDb        = org.Hs.eg.db,
                ont          = "BP",
                minGSSize    = 2,
                maxGSSize    = 100,
                pvalueCutoff = 0.05,
                verbose      = FALSE)

#write.csv(ego_bp@result, "GSEA_mito_BP_all.csv")

#Molecular Function
ego_mf <- gseGO(geneList   = geneList,
                OrgDb        = org.Hs.eg.db,
                ont          = "MF",
                minGSSize    = 2,
                maxGSSize    = 100,
                pvalueCutoff = 0.05,
                verbose      = FALSE)
#write.csv(ego_mf@result, "GSEA_mito_MF_all.csv")

#Cellular Components
ego_cc <- gseGO(geneList     = geneList,
                OrgDb        = org.Hs.eg.db,
                ont          = "CC",
                minGSSize    = 2,
                maxGSSize    = 100,
                pvalueCutoff = 0.05,
                verbose      = FALSE)
#write.csv(ego_cc@result, "GSEA_mito_CC_all.csv")

#KEGG
ego_kk <- gseKEGG(geneList     = geneList,
                  organism     = 'hsa',
                  minGSSize    = 2,
                  pvalueCutoff = 0.05,
                  verbose      = FALSE)
#write.csv(ego_kk@result, "GSEA_mito_KEGG_all.csv")

#Wiki Pathways
ego_wp = gseWP(geneList, organism = "Homo sapiens")
#write.csv(ego_wp@result, "GSEA_mito_Wiki_all.csv")

#Reactome
ego_react = gsePathway(geneList, pvalueCutoff = 0.2, pAdjustMethod = "BH", verbose=FALSE)

########### GROUP ANNOTATION ###################################################

# Pick top N from each
extract_top_gsea <- function(gsea_obj, db_name, top_n) {
  gsea_obj@result %>%
    dplyr::arrange(p.adjust) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::select(ID, Description, NES, p.adjust, core_enrichment) %>%
    dplyr::mutate(Database = db_name)
}


df_go_bp     <- extract_top_gsea(ego_bp, "GOBP", 50)
df_go_mf     <- extract_top_gsea(ego_mf, "GOMF", 50)
df_go_cc     <- extract_top_gsea(ego_cc, "GOCC", 50)
df_kegg   <- extract_top_gsea(ego_kk, "KEGG",  50)
df_react  <- extract_top_gsea(ego_react, "Reactome", 50)
df_wp     <- extract_top_gsea(ego_wp, "WP", 50)

combined_df <- bind_rows(df_go_bp,df_go_mf,df_go_cc, df_kegg, df_react, df_wp)
write.csv(combined_df,"Enrich_50_PiD_ombine_all.csv")

# Create a lowercase version of the term description
df <- combined_df %>%
  mutate(Description_lc = tolower(Description))

#I tried to assign the name of module manually here 
# Define a function to assign modules
assign_module <- function(desc) {
  if (str_detect(desc, "synapse|synaptic|neuro|glutamate|gaba|presynaptic|exocytosis|recycling")) {
    return("Synaptic regulation")
  } else if (str_detect(desc, "immune|b cell|t cell|immunoglobulin|inflammatory|cytokine|neutrophil|complement")) {
    return("Immune Response and Regulation")
  } else if (str_detect(desc, "peptidase|hydrolase|endopeptidase|enzyme inhibitor|protease")) {
    return("Peptidase activity")
  } else if (str_detect(desc, "extracellular matrix|ecm|integrin|collagen|adhesion|basement membrane")) {
    return("Extracellular matrix and adhesion")
  } else if (str_detect(desc, "mitochondria|oxidoreductase|tryptophan|metabolism|vitamin d")) {
    return("Mitochondria and metabolism")
  } else if (str_detect(desc, "endoplasmic reticulum|er lumen|golgi|secretory|retrograde")) {
    return("Endoplasmic reticulum and secretory pathway")
  } else {
    return("Unclassified")
  }
}

# Apply the function to assign modules
df$Functional_Module <- sapply(df$Description_lc, assign_module)

write.csv(df, "Assign_PiD_50.csv")



###########################################################################
#GAP43 enriched terms 

# Convert symbol to ENTREZ (assuming human)

gap43_entrez <- AnnotationDbi::select(org.Hs.eg.db, keys = "GAP43", keytype = "SYMBOL", columns = "ENTREZID")$ENTREZID


#Helper function to extract terms 
get_gap43_terms <- function(gsea_obj, gap43_id) {
  df <- gsea_obj@result
  df$includes_GAP43 <- sapply(df$core_enrichment, function(gene_str) {
    genes <- unlist(strsplit(gene_str, "/"))
    return(gap43_id %in% genes)
  })
  return(df[df$includes_GAP43, ])
}

#apply the function
gap43_terms_go_bp     <- get_gap43_terms(ego_bp, gap43_entrez)
gap43_terms_go_mf     <- get_gap43_terms(ego_mf, gap43_entrez)
gap43_terms_go_cc     <- get_gap43_terms(ego_cc, gap43_entrez)
gap43_terms_kegg   <- get_gap43_terms(ego_kk, gap43_entrez)
gap43_terms_react  <- get_gap43_terms(ego_react, gap43_entrez)
gap43_terms_wp     <- get_gap43_terms(ego_wp, gap43_entrez)


combined_gap43 <- bind_rows(
  gap43_terms_go_bp %>% mutate(Source = "GOBP"),
  gap43_terms_go_mf %>% mutate(Source = "GOMF"),
  gap43_terms_go_cc %>% mutate(Source = "GOCC"),
  gap43_terms_kegg %>% mutate(Source = "KEGG"),
  gap43_terms_react %>% mutate(Source = "Reactome"),
  gap43_terms_wp %>% mutate(Source = "WP")
)

# Plot NES of terms involving GAP43

ggplot(combined_gap43, aes(x = reorder(Description, NES), y = NES, fill = Source)) +
  geom_col() +
  coord_flip() +
  labs(
    title = "Enriched Pathways Involving GAP43",
    x = NULL, y = "Normalized Enrichment Score (NES)"
  ) +
  theme_minimal(base_size = 13)

gap43_in_bp <- ego_bp@result %>%
  filter(str_detect(core_enrichment, gap43_entrez))


gap43_cc <- setReadable(ego_cc, 'org.Hs.eg.db', 'ENTREZID')


# Let's say you have enrichResult object `ego`
# You want to keep only specific GO terms
terms_to_keep <- c("GO:0098984", "GO:0032279", "GO:0099572", "GO:0014069", "GO:0150034")

# Subset the object manually
filtered_res <- ego_cc
filtered_res@result <- ego_cc@result %>%
  dplyr::filter(ID %in% terms_to_keep)

filtered_res <- setReadable(filtered_res, 'org.Hs.eg.db', "ENTREZID")


foldChange <- setNames(filtered_res$log2FC, filtered_res$ENTREZID)

# Plot
h1 <- heatplot(filtered_res)
heatmap <- h1$data

gene_specific <- unique(heatmap$Gene)

write.csv(gene_specific, "Genes.csv")


#############################################################

###############################################################
#Leading Edge Analysis 

# Extract core enrichment data from each gseaResult
get_gsea_core_df <- function(gsea_obj, source_label) {
  df <- gsea_obj@result %>%
    dplyr::select(ID, Description,core_enrichment) %>%
    dplyr::mutate(Database = source_label)
}
core_go_bp     <- get_gsea_core_df(ego_bp, "GOBP")
core_go_mf     <- get_gsea_core_df(ego_mf, "GOMF")
core_go_cc     <- get_gsea_core_df(ego_cc, "GOCC")
core_kegg   <- get_gsea_core_df(ego_kk, "KEGG")
core_react  <- get_gsea_core_df(ego_react, "Reactome")
core_wp     <- get_gsea_core_df(ego_wp, "WP")

# Combine into one long dataframe
core_all <- bind_rows(core_go_bp ,core_go_mf ,core_go_cc , core_kegg, core_react, core_wp)

#Expand gene list into long format
core_long <- core_all %>%
  tidyr::separate_rows(core_enrichment, sep = "/") %>%
  dplyr::rename(ENTREZID = core_enrichment)

#MErge the data with Assigned_data 
# Merge assigned_df (your file with Functional_Module column)
merged <- Assigned_data %>%
  dplyr::select(ID, Functional_Module, Description, Database, NES, p.adjust) %>%
  dplyr::inner_join(core_long, by = c("ID", "Database"))

#Annotate the gene symbol with Entrez id
merged_annot <- merged %>%
  left_join(
    AnnotationDbi::select(org.Hs.eg.db, keys = .$ENTREZID, keytype = "ENTREZID", columns = "SYMBOL"),
    by = "ENTREZID" 
  )

write.csv(merged_annot, "Module_genes_leading.csv")

leading_edge_by_module <- merged_annot %>%
  group_by(Functional_Module, Description.x,SYMBOL, Database, NES, p.adjust) %>%
  summarise(
    n_terms = n(),              # How many terms the gene appears in
    gene_ids = paste(unique(ID), collapse = ", "),
    .groups = "drop"
  ) %>%
  arrange(Functional_Module, desc(n_terms))

leading_edge_by_module <- leading_edge_by_module %>%
  filter(Functional_Module != "Unclassified")

write.csv(leading_edge_by_module, "gene_specific_term.csv")

top_leading_genes <- leading_edge_by_module %>%
  group_by(Functional_Module) %>%
  slice_max(order_by = n_terms, n = 5)

###############################################################################

Corr <- read.csv("Sig_correlated.csv", stringsAsFactors = FALSE)

colnames(Corr) = c("SYMBOL", "R-coef", "P-value")

merge_lead <- leading_edge_by_module %>% left_join(Corr, by = "SYMBOL")

write.csv(merge_lead, "Correlated_leads.csv")

###############################################################################

########################################################################
#gene-Concept Network around GAP43

#filter correlated data 
network_data <- merge_lead %>%
  filter(abs(`R-coef`) >= 0.6)

#create a egde list
edge_list <- network_data %>%
  select(SYMBOL, gene_ids, Description.x, Functional_Module, NES, p.adjust)

edge_list <- edge_list %>%
  mutate(TermLabel = paste(Description.x, "\n[", Functional_Module, "]"))

#Add GAP43 as the central hub (connect it to correlated genes)

gap43_edges <- data.frame(
  from = "GAP43",
  to = unique(edge_list$SYMBOL),
  type = "correlated"
)

gene_term_edges <- network_data %>%
  mutate(TermLabel = paste0(Description.x, "\n[", Functional_Module, "]")) %>%
  select(from = SYMBOL, to = TermLabel, NES) %>%
  distinct()

gene_term_edges <- gene_term_edges %>%
  mutate(NES_direction = case_when(
    NES > 0 ~ "Positive",
    NES < 0 ~ "Negative",
    TRUE ~ "Neutral"
  ))

gap43_edges <- data.frame(
  from = "GAP43",
  to = unique(gene_term_edges$from),
  NES = NA,
  NES_direction = "GAP43"
)

all_edges <- bind_rows(gap43_edges, gene_term_edges)

all_edges <- all_edges %>%
  mutate(
    edge_type = case_when(
      from == "GAP43" ~ "GAP43_edge",
      TRUE ~ NES_direction
    ),
    edge_color = case_when(
      edge_type == "Positive" ~ "blue",
      edge_type == "Negative" ~ "red",
      edge_type == "GAP43_edge" ~ "black",
      TRUE ~ "grey70"
    )
  )


# Get unique node names
node_df <- data.frame(name = unique(c(all_edges$from, all_edges$to)))

node_df <- node_df %>%
  left_join(
    network_data %>% select(SYMBOL, `R-coef`) %>% distinct(),
    by = c("name" = "SYMBOL")
  ) %>%
  mutate(
    color_group = case_when(
      name == "GAP43" ~ "GAP43",
      `R-coef` > 0 ~ "Pos",
      `R-coef` < 0 ~ "Neg",
      TRUE ~ "Other"
    ),
    display_name = name
  )

#Create the graph with ggraph

graph <- tbl_graph(nodes = node_df, edges = all_edges, directed = FALSE)

ggraph(graph, layout = "fr") +
  geom_edge_link(aes(color = edge_color), alpha = 0.7, width = 0.9) +
  geom_node_point(aes(color = color_group), size = 4) +
  geom_node_text(aes(label = display_name), repel = TRUE, size = 3) +
  scale_color_manual(values = c(
    "Pos" = "blue", "Neg" = "red", "GAP43" = "black", "Other" = "grey60"
  )) +
  scale_edge_color_identity() +  # use manual colors from edge_color
  theme_void() +
  ggtitle("GAP43-Centered Gene-Term Network with NES and Correlation Edges")

# Convert to igraph (if not already)
cyto_graph <- graph_from_data_frame(
  all_edges,     # edges with columns: from, to
  vertices = node_df, # nodes with metadata (e.g. correlation, type)
  directed = FALSE
)

write_graph(cyto_graph, file = "gap43_network_new.graphml", format = "graphml")


#Correlation with GAP43

gap = read.csv("synapse_correlation.csv")

gap = gap[,-1]

gap_corr = cor(gap)

gap_p = cor.mtest(gap, conf.level = 0.95)


file2 = corrplot(gap_corr, p.mat = gap_p$p, 
                 method = "square", type = 'lower', 
                 tl.col =  "black",
                 tl.srt = 45,
                 sig.level = 0.01,
                 insig = 'blank',
                 order = 'hclust',
                 number.cex = 0.7,
                 number.digits = 1,
                 addCoef.col = 'grey20',
                 diag = FALSE)

gap_inf = read.csv("Inflammation_correlation.csv")

gap_inf = gap_inf[,-1]

gap_inf_corr = cor(gap_inf)

gap_p = cor.mtest(gap_inf, conf.level = 0.95)

file2 = corrplot(gap_inf_corr, p.mat = gap_p$p, 
                 method = "square", type = 'lower', 
                 tl.col =  "black",
                 tl.srt = 45,
                 sig.level = 0.05,
                 insig = 'blank',
                 order = 'hclust',
                 number.cex = 0.7,
                 number.digits = 1,
                 addCoef.col = 'grey20',
                 diag = FALSE)



