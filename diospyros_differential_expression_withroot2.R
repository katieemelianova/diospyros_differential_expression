library(dplyr)
library(magrittr)
library(readr)
library(stringr)
library(tidyverse)
library(DESeq2)
library(ape)
library(topGO)
library(UpSetR)
library(ggnewscale)
library("wesanderson")
library(egg)
library(eulerr)
library(pheatmap)
library(ggtree)
library(reshape2)
library(cowplot)
library(phytools)
source("differential_expression_functions.R")
source("/Users/katieemelianova/Desktop/Diospyros/diospyros_R_functions/diospyros_soil_leaf_species_element_dataset.R")
source("/Users/katieemelianova/Desktop/Diospyros/diospyros_R_functions/diospyros_orthogroups_long.R")
source("/Users/katieemelianova/Desktop/Diospyros/diospyros_R_functions/diospyros_genome_annotation.R")



##############################################################
#        load in the GO ID mappings and GO annotations       #
##############################################################

# for the GO term enrichment tests
mp<-readMappings("/Users/katieemelianova/Desktop/Diospyros/vieillardii_functionalAnnotation_transcript.topGO.txt")
mp_go <- read_delim("/Users/katieemelianova/Desktop/Diospyros/vieillardii_functionalAnnotation_transcript.txt")


##############################################
#         load in braker annotations         #
##############################################


braker_annotation <- get_braker_annotation()
braker_vieillardii <- braker_annotation$vieillardii_braker
braker_revolutissima <- braker_annotation$revolutissima_braker
braker_impolita <- braker_annotation$impolita_braker


##################################################
#       function for reading in featurecounts    #
##################################################

read_in_featurecounts<-function(input_path, strings_to_remove){
  # read in each featurecounts file and join them together
  df <- list.files(path=input_path, full.names = TRUE) %>% 
    lapply(read_tsv, skip=1) %>% 
    purrr::reduce(left_join, by = c("Geneid", "Chr", "Start", "End", "Strand", "Length"))
  # get and store the first standard featurecounts column names
  fc_cols<-colnames(df)[1:6]
  #get and store the egbne lengths to return for DE later
  gene_lengths<-df$Length
  # take the 7th until the last colnames of the data frame
  sample_names<-colnames(df)[7:length(colnames(df))]
  # remove any prefixes that need to be removed
  for (strings in strings_to_remove){
    sample_names<-str_replace(sample_names, strings, "")
  }
  # now set the df colnames to the shortened sample names for ease of reading
  colnames(df)<-c(fc_cols, sample_names)
  # get the relevant columns (geneid and sample counts) and set the geneid column to be rownames
  df_counts<-df %>% 
    dplyr::select(c(1,7:length(colnames(df)))) %>% 
    column_to_rownames("Geneid")
  return(list(counts=df_counts, lengths=gene_lengths))
}

#################################################################################
#    diospyros notation specific function for annotating based on sample name   #
#################################################################################

annotate_diospyros_samples<-function(sample_tibble){
  # expected format of samples:
  # #(species)_(mother)(plant_id)(pseudorep)(techrep)
  # I remove the "r" for technical replicate, this will be evident from repeated lines in resulting tibble
  # where there is no pseudoreplicate, list "none"
  to_return<-sample_tibble %>%
    mutate(to_rownames=sample_id) %>%
    separate(sample_id, 
             into = c("species", "data"), 
             sep = "_") %>%
    separate(data, 
             into = c("data", "to_remove"), 
             sep = "r") %>%
    dplyr::select(-to_remove) %>%
    separate(data, 
             into = c("mother", "pseudoreplicate"), 
             sep = "(?=[0-9])(?<=[a-z])") %>%
    replace_na(list(pseudoreplicate = "none")) %>%
    separate(mother, 
             into = c("mother", "offspring"), 
             sep = "(?=[a-z])(?<=[0-9])") %>%
    column_to_rownames("to_rownames")
  return(to_return)
}

#####################################################################
#     read in feauturecounts and annotate samples based on names    #
#####################################################################

strings_to_remove<-c("hequetiae/", "impolita/", "calciphila/", "sp_pic_nga/", "revolutissima/", "labillardierei/",  "Aligned.sortedByCoord.out.bam")
all_featurecounts<-read_in_featurecounts('leaf_root_counts', strings_to_remove)

# remove outlier
all_featurecounts$counts %<>% dplyr::select(-"l11a1_MKDL250006837-1A_22VTLFLT4_L7")





all_samples<-all_featurecounts$counts %>% 
  colnames() %>% 
  tibble() %>% 
  set_colnames("sample_id") %>% 
  mutate(Tissue = ifelse(grepl("MKDL", sample_id ), "Root", "Leaf")) %>% 
  mutate(sample_id = case_when(Tissue == "Root" & substr(sample_id,1,1) == "c" ~ str_replace(sample_id, "c", "cal_"),
                               Tissue == "Root" & substr(sample_id,1,1) == "h" ~ str_replace(sample_id, "h", "heq_"),
                               Tissue == "Root" & substr(sample_id,1,1) == "i" ~ str_replace(sample_id, "i", "imp_"),
                               Tissue == "Root" & substr(sample_id,1,1) == "l" ~ str_replace(sample_id, "l", "lab_"),
                               Tissue == "Root" & substr(sample_id,1,1) == "r" ~ str_replace(sample_id, "r", "rev_"),
                               Tissue == "Root" & substr(sample_id,1,1) == "p" ~ str_replace(sample_id, "p", "spn_"),
                               Tissue == "Leaf" ~ sample_id)) %>% 
  annotate_diospyros_samples()

# set soil preferences
all_samples %<>% mutate(soil=case_when(species == "cal" ~ "nonultramafic",
                                       species == "heq" ~ "ultramafic",
                                       species == "imp" ~ "nonultramafic",
                                       species == "lab" ~ "nonultramafic",
                                       species == "rev" ~ "ultramafic",
                                       species == "spn" ~ "ultramafic",))



all_counts<-all_featurecounts$counts
all_lengths<-all_featurecounts$lengths

# set colnames of counts to updated ones in samples
colnames(all_counts) <- rownames(all_samples)

all_comparison<-specify_comparison(all_samples, all_counts, "species %in% c('cal', 'spn', 'heq', 'imp', 'lab', 'rev')")
all_dds <- DESeqDataSetFromMatrix(countData = all_comparison[["counts"]],
                                  colData = all_comparison[["samples"]],
                                  design = ~ species) %>%
  varianceStabilizingTransformation()
pca_all <-plotPCA(all_dds, intgroup=c("species", "pseudoreplicate", "Tissue"), ntop = 5000, returnData = TRUE)
species_colours <- c("#68C7AA", "#F1A1FD", "#AD1640", "#FF007E", "#22660D", "#FE841C")

pca_all %<>% mutate(species=case_when(species == "cal" ~ "calciphila",
                                      species == "heq" ~ "hequetiae",
                                      species == "imp" ~ "impolita",
                                      species == "rev" ~ "revolutissima",
                                      species == "spn" ~ "sp. Pic N'ga",
                                      species == "lab" ~ "labillardierei"))



ggplot(pca_all, aes(PC1, PC2)) +          
  geom_point(size = 1, stroke = 5, aes(fill = species, colour=species, shape=Tissue)) + 
  scale_color_manual(values=species_colours) +
  theme(axis.text = element_text(size=15),
        axis.title = element_text(size=18),
        legend.title = element_blank(),
        legend.text = element_text(size = 15))






########################
#       plot PCA       #
########################

all_comparison<-specify_comparison(all_samples, all_counts, "species %in% c('cal', 'spn', 'heq', 'imp', 'lab', 'rev')")
all_dds <- DESeqDataSetFromMatrix(countData = all_comparison[["counts"]],
                                  colData = all_comparison[["samples"]],
                                  design = ~ species) %>%
  varianceStabilizingTransformation()
pca_all <-plotPCA(all_dds, intgroup=c("species", "pseudoreplicate"), ntop = 5000, returnData = TRUE)
species_colours <- c("#68C7AA", "#F1A1FD", "#AD1640", "#FF007E", "#22660D", "#FE841C")

pca_all %<>% mutate(species=case_when(species == "cal" ~ "calciphila",
                                      species == "heq" ~ "hequetiae",
                                      species == "imp" ~ "impolita",
                                      species == "rev" ~ "revolutissima",
                                      species == "spn" ~ "sp. Pic N'ga",
                                      species == "lab" ~ "labillardierei"))



##########################################
#          plot PCA with phylo tree      #
##########################################

pic1<-ggplot(pca_all, aes(PC1, PC2)) +          
  geom_point(size = 1, stroke = 5, aes(fill = species, colour=species)) + 
  scale_color_manual(values=species_colours) +
  theme(axis.text = element_text(size=15),
        axis.title = element_text(size=18),
        legend.title = element_blank(),
        legend.text = element_text(size = 15))

#species_tree<-get_species_tree()
species_tree<-ape::read.tree("/Users/katieemelianova/Desktop/Diospyros/TeernaDiospyrosTrees/species_tree_small_caledit.txt")
species_tree$tip.label<-case_when(startsWith(species_tree$tip.label, "vie") ~"vieillardii",
                                  startsWith(species_tree$tip.label, "lab") ~"labillardierei",
                                  startsWith(species_tree$tip.label, "rev") ~"revolutissima",
                                  startsWith(species_tree$tip.label, "heq") ~"hequetiae",
                                  startsWith(species_tree$tip.label, "imp") ~"impolita",
                                  startsWith(species_tree$tip.label, "cal") ~"calciphila",
                                  startsWith(species_tree$tip.label, "spn") ~"sp. Pic N'ga")

species_tree <- root(species_tree, which(species_tree$tip.label == "vieillardii"))

# make a list where item name is species name and objects within are the tip labels belonging to that species
groupInfo<-split(species_tree$tip.label, species_tree$tip.label)

# use groupOTU to group the tips by species and plot
species_tree<-groupOTU(species_tree, groupInfo, group_name = "species")

#outgroup<-c("sandwicensis")
outgroup<-c("vieillardii")

ultramafic<-c("hequetiae", "revolutissima", "sp. Pic N'ga")
volcanic<-c("calciphila", "impolita", "labillardierei")

colours_tips <- case_when(species_tree$tip.label %in% ultramafic ~"Ultramafic",
                          species_tree$tip.label %in% volcanic ~"Volcanic",
                          !(species_tree$tip.label %in% c(outgroup, ultramafic, volcanic)) ~ "No data")

colours_labels <- case_when(species_tree$tip.label == "calciphila" ~ "#68C7AA",
                            species_tree$tip.label == "impolita" ~ "#AD1640",
                            species_tree$tip.label == "labillardierei" ~ "#FF007E",
                            species_tree$tip.label == "hequetiae" ~ "#F1A1FD",
                            species_tree$tip.label == "revolutissima" ~ "#22660D",
                            species_tree$tip.label == "sp. Pic N'ga" ~ "#FE841C",
                            !(species_tree$tip.label %in% c("calciphila", "impolita", "labillardierei", "hequetiae", "revolutissima", "sp. Pic N'ga")) ~ "grey77")



colours_labels_pairs <- case_when(species_tree$tip.label == "calciphila" ~ "dodgerblue2",
                                  species_tree$tip.label == "impolita" ~ "hotpink",
                                  species_tree$tip.label == "labillardierei" ~ "darkolivegreen3",
                                  species_tree$tip.label == "hequetiae" ~ "darkolivegreen3",
                                  species_tree$tip.label == "revolutissima" ~ "hotpink",
                                  species_tree$tip.label == "sp. Pic N'ga" ~ "dodgerblue2",
                                  !(species_tree$tip.label %in% c("calciphila", "impolita", "labillardierei", "hequetiae", "revolutissima", "sp. Pic N'ga")) ~ "grey77")

dd <- data.frame(taxa=species_tree$tip.label, tipcols=colours_tips, labelcols=colours_labels)
p<-ggtree(species_tree, size=1)
p <- p %<+% dd
p2<-p + new_scale_color() + geom_tiplab(size=5.5, aes(color=species), offset=0.005, show.legend=FALSE) + 
  scale_color_manual(values=colours_labels, limits=species_tree$tip.label) + 
  theme(legend.title = element_blank(),
        legend.text = element_text(size=12)) +
  xlim_tree(5)



elements<-get_element_dataset()
# match tiplab name with element name
elements$species[elements$species == "sp PicN'ga"] <- "sp. Pic N'ga"
test<-elements[elements$species %in% c("calciphila", "impolita", "labillardierei", "hequetiae", "revolutissima", "sp. Pic N'ga"),]  %>% dplyr::select(species, Cr_soil, Ni_soil, Co_soil) %>% data.frame()
pic2<-facet_plot(p2, panel="Relative Soil Nickel Concentration", data=test, geom_boxplot, mapping = aes(x=Ni_soil, group = label, fill=species)) + 
  theme_bw() +
  theme(legend.position="none",
        strip.text.x = element_text(size = 15),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()
  ) + 
  scale_fill_manual(values=c("#68C7AA", "#F1A1FD", "#AD1640", "#FF007E", "#22660D", "#FE841C"))

grid.arrange(pic1, pic2, nrow = 1)

all_samples %<>% mutate(soil=case_when(species == "cal" ~ "nonultramafic",
                                       species == "heq" ~ "ultramafic",
                                       species == "imp" ~ "nonultramafic",
                                       species == "lab" ~ "nonultramafic",
                                       species == "rev" ~ "ultramafic",
                                       species == "spn" ~ "ultramafic",))


######################################################################
#        DEGs for each ultramafic-nonultramafic pair per tissue      #
######################################################################

# per tissue and species pair DEG tests
cal_spn_root <-specify_comparison(all_samples, all_counts, "species %in% c('cal', 'spn') & Tissue == 'Root'") %>% run_diffexp("species", "spn", "cal", all_lengths, cpm_threshold=5, min_samples=3)
heq_lab_root <-specify_comparison(all_samples, all_counts, "species %in% c('heq', 'lab') & Tissue == 'Root'") %>% run_diffexp("species", "heq", "lab", all_lengths, cpm_threshold=5, min_samples=3)
rev_imp_root <-specify_comparison(all_samples, all_counts, "species %in% c('rev', 'imp') & Tissue == 'Root'") %>% run_diffexp("species", "rev", "imp", all_lengths, cpm_threshold=5, min_samples=3)
cal_spn_leaf <-specify_comparison(all_samples, all_counts, "species %in% c('cal', 'spn') & Tissue == 'Leaf'") %>% run_diffexp("species", "spn", "cal", all_lengths, cpm_threshold=5, min_samples=3)
heq_lab_leaf <-specify_comparison(all_samples, all_counts, "species %in% c('heq', 'lab') & Tissue == 'Leaf'") %>% run_diffexp("species", "heq", "lab", all_lengths, cpm_threshold=5, min_samples=3)
rev_imp_leaf <-specify_comparison(all_samples, all_counts, "species %in% c('rev', 'imp') & Tissue == 'Leaf'") %>% run_diffexp("species", "rev", "imp", all_lengths, cpm_threshold=5, min_samples=3)

# get significant genes
cal_spn_root_sig <- cal_spn_root$results %>% data.frame() %>% filter(padj < 0.05 & abs(log2FoldChange) > 2)
heq_lab_root_sig <- heq_lab_root$results %>% data.frame() %>% filter(padj < 0.05 & abs(log2FoldChange) > 2)
rev_imp_root_sig <- rev_imp_root$results %>% data.frame() %>% filter(padj < 0.05 & abs(log2FoldChange) > 2)
cal_spn_leaf_sig <- cal_spn_leaf$results %>% data.frame() %>% filter(padj < 0.05 & abs(log2FoldChange) > 2)
heq_lab_leaf_sig <- heq_lab_leaf$results %>% data.frame() %>% filter(padj < 0.05 & abs(log2FoldChange) > 2)
rev_imp_leaf_sig <- rev_imp_leaf$results %>% data.frame() %>% filter(padj < 0.05 & abs(log2FoldChange) > 2)

# GO enrichment 
hq_enrich_root <- get_enriched_terms(rownames(cal_spn_root_sig), mp, return_sample_GOData=TRUE)
ri_enrich_root <- get_enriched_terms(rownames(heq_lab_root_sig), mp, return_sample_GOData=TRUE)
cs_enrich_root <- get_enriched_terms(rownames(rev_imp_root_sig), mp, return_sample_GOData=TRUE)
hq_enrich_leaf <- get_enriched_terms(rownames(cal_spn_leaf_sig), mp, return_sample_GOData=TRUE)
ri_enrich_leaf <- get_enriched_terms(rownames(heq_lab_leaf_sig), mp, return_sample_GOData=TRUE)
cs_enrich_leaf <- get_enriched_terms(rownames(rev_imp_leaf_sig), mp, return_sample_GOData=TRUE)

# make results object contain significant only genes
hq_enrich_root$result %<>% filter(as.numeric(classicFisher) < 0.05)
ri_enrich_root$result %<>% filter(as.numeric(classicFisher) < 0.05)
cs_enrich_root$result %<>% filter(as.numeric(classicFisher) < 0.05)
hq_enrich_leaf$result %<>% filter(as.numeric(classicFisher) < 0.05)
ri_enrich_leaf$result %<>% filter(as.numeric(classicFisher) < 0.05)
cs_enrich_leaf$result %<>% filter(as.numeric(classicFisher) < 0.05)





c(hq_enrich_root$GO.ID, ri_enrich_root$GO.ID, cs_enrich_root$GO.ID) %>% table() %>% data.frame() %>% filter(Freq > 2) %>% pull(".")
c(hq_enrich_leaf$GO.ID, ri_enrich_leaf$GO.ID, cs_enrich_leaf$GO.ID) %>% table() %>% data.frame() %>% filter(Freq > 2) %>% pull(".")



get_de_genes_in_term<-function(degs, go_term, go_data){
  genes_in_term<-genesInTerm(go_data, go_term)[[1]]
  degs_in_term<-intersect(genes_in_term, degs)
  return(degs_in_term)
}


a <- get_de_genes_in_term(rownames(hq_leaf), c("GO:0006468"), hq_enrich_leaf$goData)

hq_leaf[a, ]


c("GO:0006468", "GO:0006952")


