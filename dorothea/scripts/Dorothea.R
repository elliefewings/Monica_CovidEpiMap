# Time stamp
Sys.time()

# Clean up
rm(list=ls())

# Load libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(magrittr)
library(Seurat)
library(stringr)
library(dorothea)
library(tibble)
library(pheatmap)

# Folders
setwd("/net/data.isilon/ag-saez/bq_efewings/data/dorothea/")

#Load data
data <- readRDS("integrated.RNA.Tcells.annotated.rds")

#Separate metadata for exploration
metadata <- data@meta.data

rundor <- function(x){
  dorothea_regulon_human <- get(data("dorothea_hs", package = "dorothea"))
  
  ## We obtain the regulons based on interactions with confidence level A, B and C
  regulon <- dorothea_regulon_human %>%
    dplyr::filter(confidence %in% c("A","B","C"))
  
  
  ## We compute Viper Scores 
  x <- run_viper(x, regulon,
                 options = list(method = "scale", minsize = 4, 
                                eset.filter = FALSE, cores = 1, 
                                verbose = FALSE))
  
  #Extract viper scores
  viper_scores_df <- GetAssayData(x, assay = "dorothea") %>%
    data.frame(check.names = F) %>%
    t()
  
  ## We create a data frame containing the cells and their clusters
  CellsClusters <- data.frame(cell = names(Idents(x)), 
                              cell_type = as.character(Idents(x)),
                              check.names = F)
  
  ## We create a data frame with the Viper score per cell and its clusters
  viper_scores_clusters <- viper_scores_df  %>%
    data.frame() %>% 
    rownames_to_column("cell") %>%
    gather(tf, activity, -cell) %>%
    inner_join(CellsClusters)
  
  ## We summarize the Viper scores by cellpopulation
  summarized_viper_scores <- viper_scores_clusters %>% 
    group_by(tf, cell_type) %>%
    summarise(avg = mean(activity),
              std = sd(activity))
  
  ## For visualization purposes, we select the 20 most variable TFs across clusters according to our scores.
  
  #Count number of populations
  n <- length(unique(CellsClusters$cell_type))
  
  ## We select the 20 most variable TFs. (20*9 populations = 180)
  highly_variable_tfs <- summarized_viper_scores %>%
    group_by(tf) %>%
    mutate(var = var(avg))  %>%
    ungroup() %>%
    top_n((20*n), var) %>%
    distinct(tf)
  
  ## We prepare the data for the plot
  summarized_viper_scores_df <- summarized_viper_scores %>%
    semi_join(highly_variable_tfs, by = "tf") %>%
    dplyr::select(-std) %>%   
    spread(tf, avg) %>%
    data.frame(row.names = 1, check.names = FALSE) 
  palette_length = 100
  my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)
  
  my_breaks <- c(seq(min(summarized_viper_scores_df), 0, 
                     length.out=ceiling(palette_length/2) + 1),
                 seq(max(summarized_viper_scores_df)/palette_length, 
                     max(summarized_viper_scores_df), 
                     length.out=floor(palette_length/2)))
  
  viper_hmap <- pheatmap(t(summarized_viper_scores_df),fontsize=14, 
                         fontsize_row = 10, 
                         color=my_color, breaks = my_breaks, 
                         main = "DoRothEA (ABC)", angle_col = 45,
                         treeheight_col = 0,  border_color = NA) 
  
  out <- list(scores=summarized_viper_scores, heat=viper_hmap)
  
  return(out)
}

#alldata <- rundor(data)

#Create condition per cell type annotation and rereun
data$CTI <- paste(data@meta.data$condition, data@meta.data$integrated_annotation, sep=".")

data <- SetIdent(data, value="CTI")

cti <- rundor(data)

#Count number of populations
n <- length(unique(cti$scores$cell_type))

## We select the 20 most variable TFs. (20*9 populations = 180)
highly_variable_tfs <- cti$scores %>%
  group_by(tf) %>%
  mutate(var = var(avg))  %>%
  ungroup() %>%
  top_n((20*n), var) %>%
  distinct(tf)

summarized_viper_scores_df <- cti$scores %>%
  semi_join(highly_variable_tfs, by = "tf") %>%
  dplyr::select(-std) %>%   
  spread(tf, avg) %>%
  data.frame(row.names = 1, check.names = FALSE) 
palette_length = 100
my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)

annotation <- data.frame(condition=sapply(strsplit(row.names(summarized_viper_scores_df),"\\."), `[`, 1),
                         row.names = row.names(summarized_viper_scores_df))

celltype <- sapply(strsplit(row.names(summarized_viper_scores_df),"\\."), `[`, 2)

my_breaks <- c(seq(min(summarized_viper_scores_df), 0, 
                   length.out=ceiling(palette_length/2) + 1),
               seq(max(summarized_viper_scores_df)/palette_length, 
                   max(summarized_viper_scores_df), 
                   length.out=floor(palette_length/2)))

CTI_hmap <- pheatmap(t(summarized_viper_scores_df),fontsize=14, 
                     fontsize_row = 10, 
                     color=my_color, breaks = my_breaks, 
                     main = "DoRothEA (ABC)", angle_col = 45,
                     treeheight_col = 0,  border_color = NA,
                     annotation_col = annotation, labels_col=celltype,
                     cluster_cols = FALSE)


#Run per condition
data <- SetIdent(data, value="condition")

#Subset data
am <- subset(data, idents="active_mild")
am <- SetIdent(am, value="integrated_annotations")
am.p <- rundor(am)

#Subset data
as <- subset(data, idents="active_severe")
as <- SetIdent(as, value="integrated_annotations")
as.p <- rundor(as)

#Subset data
h <- subset(data, idents="healthy")
h <- SetIdent(h, value="integrated_annotations")
h.p <- rundor(h)

#Subset data
rm <- subset(data, idents="recovered_mild")
rm <- SetIdent(rm, value="integrated_annotations")
rm.p <- rundor(rm)

#Subset data
rs <- subset(data, idents="recovered_severe")
rs <- SetIdent(rs, value="integrated_annotations")
rs.p <- rundor(rs)


#Write scores
write.csv(am.p$scores, "dorothea_am_scores.csv")
write.csv(as.p$scores, "dorothea_as_scores.csv")
write.csv(h.p$scores, "dorothea_h_scores.csv")
write.csv(rm.p$scores, "dorothea_rm_scores.csv")
write.csv(rs.p$scores, "dorothea_rs_scores.csv")
write.csv(cti$scores, "dorothea_compareconditons_scores.csv")

#Create heatmap names
cti.h <- "heatmap_compareconditons.pdf"
am.h <- "heatmap_am.pdf"
as.h <- "heatmap_as.pdf"
h.h <- "heatmap_h.pdf"
rm.h <- "heatmap_rm.pdf"
rs.h <- "heatmap_rs.pdf"

#Write heatmap
pdf(cti.h, width=18, height=9)
print(CTI_hmap)
dev.off()

pdf(am.h, width=18, height=9)
print(am.p$heat)
dev.off()

pdf(as.h, width=18, height=9)
print(as.p$heat)
dev.off()

pdf(h.h, width=18, height=9)
print(h.p$heat)
dev.off()

pdf(rm.h, width=18, height=9)
print(rm.p$heat)
dev.off()

pdf(rs.h, width=18, height=9)
print(rs.p$heat)
dev.off()

#Save environment
save.image(file = "dorothea.RData")
