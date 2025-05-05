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
all<-read_in_featurecounts('all_counts', strings_to_remove)
all_samples<-all$counts %>% colnames() %>% tibble() %>% set_colnames("sample_id") %>% annotate_diospyros_samples()
all_counts<-all$counts
all_lengths<-all$lengths


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

species_tree<-get_species_tree()
outgroup<-c("sandwicensis")
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
  expand_limits(x = 0.15)

##############################
#    load in elements data   #
##############################

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

#pdf("PCA_with_tree.pdf", width = 15, height = 7)
#grid.arrange(pic1, pic2, nrow = 1)
#dev.off()


############################################################################
#     add in soil preference to sample info for differential expression   #
############################################################################

all_samples %<>% mutate(soil=case_when(species == "cal" ~ "nonultramafic",
                                       species == "heq" ~ "ultramafic",
                                       species == "imp" ~ "nonultramafic",
                                       species == "lab" ~ "nonultramafic",
                                       species == "rev" ~ "ultramafic",
                                       species == "spn" ~ "ultramafic",))


##########################
#     Run DE analysis    #
##########################

cal_spn <-specify_comparison(all_samples, all_counts, "species %in% c('cal', 'spn')") %>% run_diffexp("species", "spn", "cal", all_lengths, cpm_threshold=5, min_samples=3)
heq_lab <-specify_comparison(all_samples, all_counts, "species %in% c('heq', 'lab')") %>% run_diffexp("species", "heq", "lab", all_lengths, cpm_threshold=5, min_samples=3)
rev_imp<-specify_comparison(all_samples, all_counts, "species %in% c('rev', 'imp')") %>% run_diffexp("species", "rev", "imp", all_lengths, cpm_threshold=5, min_samples=3)

# write results to file
#cal_spn$results %>% data.frame() %>% arrange(padj) %>% rownames_to_column() %>% write_tsv("spPicNga_calciphila.DEG.tsv")
#heq_lab$results %>% data.frame() %>% arrange(padj) %>% rownames_to_column() %>% write_tsv("hequetiae_labillardierei.DEG.tsv")
#rev_imp$results %>% data.frame() %>% arrange(padj) %>% rownames_to_column() %>% write_tsv("revolutissimae_impolita.DEG.tsv")


# write sample info to file
all_samples %>% rownames_to_column() %>% set_colnames(c("sample_id", "species", "mother", "offspring", "pseudoreplicate", "soil_type")) %>% write_tsv("diospyros_common_garden_samples.tsv")



############################################################
#     Use pairwise FSTs to decide how to pair up species   #
############################################################

#read.table("avg_fst_lib12345.txt", header=TRUE) %>% filter(sp1 %in% c("impolita") & sp2 %in% c("hequetiae", "revolutissima", "impolita"))
#read.table("avg_fst_lib12345.txt", header=TRUE) %>% filter(sp1 %in% c("labillardierei") & sp2 %in% c("hequetiae", "revolutissima", "impolita"))
#read.table("avg_fst_lib12345.txt", header=TRUE) %>% filter(sp1 %in% c("revolutissima", "spPicNga", "hequetiae") & sp2 %in% c("impolita", "calciphila", "labillardierei"))


############################################################
#     Upset plot intersection between species pair DEGs    #
############################################################

# create the intersect list
listInput<-list(`calciphila vs sp. Pic N'Ga` = cal_spn$results %>% data.frame() %>% filter(padj < 0.05 & abs(log2FoldChange) > 2) %>% rownames(),
                `hequetiae vs labillardierei` = heq_lab$results %>% data.frame() %>% filter(padj < 0.05 & abs(log2FoldChange) > 2) %>% rownames(),
                `revolutissima vs impolita` = rev_imp$results %>% data.frame() %>% filter(padj < 0.05 & abs(log2FoldChange) > 2) %>% rownames())

# plot upset plot
upset(fromList(listInput), order.by = "freq", nsets = 5)

# get genes DE between each ultramafic-nonultramafic pair
ultramaf_nonultramaf_intersect <- intersect(intersect(listInput$`calciphila vs sp. Pic N'Ga`, listInput$`hequetiae vs labillardierei`), listInput$`revolutissima vs impolita`)


#################################################################################################
#     get the annotations for the DEGs found to intersect between ultramaf-nonultramaf pairs    #
#################################################################################################

mp_go %>% filter(SeqName %in% ultramaf_nonultramaf_intersect) %>%
  data.frame() %>% dplyr::select(SeqName, Description, GO.Names)


############################################################
#    Eulerr plot intersection between species pair DEGs    #
############################################################

fit2<-euler(listInput) 

euler_pairs<-plot(fit2,
     fills = c("dodgerblue2", "darkolivegreen3", "hotpink2"),
     edges = TRUE,
     fontsize = 50,
     labels = NULL,
     quantities = list(fontsize = 13),
     legend = list(labels = c("calciphila vs sp. Pic N'Ga", "hequetiae vs labillardierei", "revolutissima vs impolita"), fontsize=20, side = "bottom"),
     alpha=0.8)


tree_pairs<-p2 + scale_color_manual(values=colours_labels_pairs, limits=species_tree$tip.label)


pdf("euler_with_tree.pdf", width = 8, height = 7)
grid.arrange(tree_pairs, euler_pairs, nrow = 1)
dev.off()


############################################################
#               plot average expression of DEGs            #
############################################################


cal_av<-fpm(cal_spn$dds) %>% data.frame() %>% dplyr::select(cal_9a, cal_9b, cal_9c1, cal_9c2, cal_9d1, cal_9d2) %>% rowMeans()
pic_av<-fpm(cal_spn$dds) %>% data.frame() %>% dplyr::select(spn_18b, spn_18c, spn_19a, spn_19b) %>% rowMeans()

heq_av<-fpm(heq_lab$dds) %>% data.frame() %>% dplyr::select(heq_13a, heq_13b, heq_13br, heq_13c, heq_13d, heq_14a, heq_14b, heq_14c, heq_14cr) %>% rowMeans()
lab_av<-fpm(heq_lab$dds) %>% data.frame() %>% dplyr::select(lab_11a1, lab_11a2, lab_12a, lab_12b, lab_12c1, lab_12c2) %>% rowMeans()

imp_av<-fpm(rev_imp$dds) %>% data.frame() %>% dplyr::select(imp_3ba1, imp_3ba2, imp_3bb, imp_3bc, imp_3bd1, imp_3bd2, imp_3bd2r) %>% rowMeans()
rev_av<-fpm(rev_imp$dds) %>% data.frame() %>% dplyr::select(rev_24a1, rev_24a1r, rev_24a2, rev_24b1, rev_24b2, rev_25a1, rev_25a2, rev_25a3, rev_25aq3r) %>% rowMeans()




average_expression_plot<-inner_join(inner_join(
  data.frame(cal_av, pic_av) %>% rownames_to_column(), 
  data.frame(heq_av, lab_av) %>% rownames_to_column()), 
  data.frame(imp_av, rev_av) %>% rownames_to_column()) %>% 
  filter(rowname %in% ultramaf_nonultramaf_intersect) %>% 
  melt() %>%
  mutate(rowname=str_replace(rowname, ".t1", "")) %>%
  mutate(rowname=str_replace(rowname, ".t2", "")) %>%
  mutate(rowname=str_replace(rowname, ".t3", "")) %>%
  mutate(pair_name=case_when(variable %in% c("cal_av", "pic_av") ~ "Sp. Pic N'Ga - Calciphila",
                             variable %in% c("heq_av", "lab_av") ~ "Hequetiae - Labillardierei",
                             variable %in% c("rev_av", "imp_av") ~ "Revolutissima - Impolita"),
         soiltype=case_when(variable %in% c("rev_av", "heq_av", "pic_av") ~ "Ultramafic",
                            variable %in% c("cal_av", "lab_av", "imp_av") ~ "Non-Ultramafic")) %>%
  ggplot(aes(x=rowname, y=sqrt(value))) + 
  geom_point(size= 6, alpha=0.55, aes(shape=pair_name, color=soiltype)) +
  geom_point(size= 6, alpha=0.55, aes(shape=pair_name, color=soiltype)) + 
  scale_shape_manual(values=c(17, 19, 15)) +
  scale_colour_manual(values=c("darkolivegreen3", "hotpink2")) +
  ylab("sqrt(Counts per Million)") +
  xlab("Gene") + 
  theme(axis.text = element_text(size=10), 
        axis.title = element_text(size=15),
        legend.text = element_text(size=13),
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.3),
        panel.background = element_rect(fill = 'white', colour = 'grey72'),
        legend.position=c(.85, 0.8))

#pdf("average_expressionDEGs.pdf", width=12, height=6)
average_expression_plot
#dev.off()

#######################################################################################################################
#      How many genes have an average expression greater than 100 counts per million (cpm) in at least one species?    #                      #
#######################################################################################################################

inner_join(inner_join(
  data.frame(cal_av, pic_av) %>% rownames_to_column(), 
  data.frame(heq_av, lab_av) %>% rownames_to_column()), 
  data.frame(imp_av, rev_av) %>% rownames_to_column()) %>% 
  filter(rowname %in% ultramaf_nonultramaf_intersect) %>% 
  melt() %>%
  mutate(rowname=str_replace(rowname, ".t1", "")) %>%
  mutate(rowname=str_replace(rowname, ".t2", "")) %>%
  mutate(rowname=str_replace(rowname, ".t3", "")) %>%
  mutate(pair_name=case_when(variable %in% c("cal_av", "pic_av") ~ "Sp. Pic N'Ga - Calciphila",
                             variable %in% c("heq_av", "lab_av") ~ "Hequetiae - Labillardierei",
                             variable %in% c("rev_av", "imp_av") ~ "Revolutissima - Impolita"),
         soiltype=case_when(variable %in% c("rev_av", "heq_av", "pic_av") ~ "Ultramafic",
                            variable %in% c("cal_av", "lab_av", "imp_av") ~ "Non-Ultramafic")) %>%
  filter(value > 100) %>%
  pull(rowname) %>%
  unique()

######################################################################################################
#                  get subsets of DEGs - all DEGs and those of large effect                          #
######################################################################################################

# get the deseq results and get genes to focus on using baseMean as filter

# DEGs
ultramafic_nonultramafic_DEGs_large_effect <- results(rev_imp$dds)[ultramaf_nonultramaf_intersect,] %>% data.frame() %>% filter(baseMean > 1000 & abs(log2FoldChange) > 2) %>% rownames()
ultramafic_nonultramafic_DEGs <- results(rev_imp$dds)[ultramaf_nonultramaf_intersect,] %>% data.frame() %>% rownames()

# DEG counts
rev_imp_ultramafic_nonultramafic_counts_large_effect <- fpm(rev_imp$dds) %>% data.frame() %>% rownames_to_column() %>% filter(rowname %in% ultramafic_nonultramafic_DEGs_large_effect)
rev_imp_ultramafic_nonultramafic_counts <- fpm(rev_imp$dds) %>% data.frame() %>% rownames_to_column() %>% filter(rowname %in% ultramafic_nonultramafic_DEGs)



#####################################################################################
#    get all counts no filter and draw heatmap of different gene copy expression    #
#####################################################################################

# get unfiltered deseq objects
cal_spn_nothresh <- specify_comparison(all_samples, all_counts, "species %in% c('cal', 'spn')") %>% run_diffexp("species", "spn", "cal", all_lengths, cpm_threshold=0, min_samples=0)
heq_lab_nothresh <- specify_comparison(all_samples, all_counts, "species %in% c('heq', 'lab')") %>% run_diffexp("species", "heq", "lab", all_lengths, cpm_threshold=0, min_samples=0)
rev_imp_nothresh <- specify_comparison(all_samples, all_counts, "species %in% c('rev', 'imp')") %>% run_diffexp("species", "rev", "imp", all_lengths, cpm_threshold=0, min_samples=0)

# creeate annotation dataframe
annotation_values <-rbind(counts(cal_spn_nothresh$dds, normalized=TRUE) %>% t(),
         counts(heq_lab_nothresh$dds, normalized=TRUE) %>% t(),
         counts(rev_imp_nothresh$dds, normalized=TRUE) %>% t()) %>% 
  rownames()
annotation_species<-c(rep("calciphila", 6), rep("Sp. Pic N'Ga", 4), rep("Hequetiae", 9), 
                      rep("Labillardiei", 6), rep("impolita", 7), rep("revolutissima", 9)) %>% data.frame()
rownames(annotation_species) <- annotation_values

# draw heatmap
rbind(counts(cal_spn_nothresh$dds, normalized=TRUE)[c("g10152.t1", "g10154.t1", "g10156.t1", "g10157.t1"),] %>% t(),
      counts(heq_lab_nothresh$dds, normalized=TRUE)[c("g10152.t1", "g10154.t1", "g10156.t1", "g10157.t1"),] %>% t(),
      counts(rev_imp_nothresh$dds, normalized=TRUE)[c("g10152.t1", "g10154.t1", "g10156.t1", "g10157.t1"),] %>% t()) %>% 
  pheatmap(scale="row", cluster_rows = FALSE, annotation_row=annotation_species)

# draw heatmap feature scaled
rbind(t(apply(fpm(cal_spn_nothresh$dds)[ultramafic_nonultramafic_DEGs_large_effect,] %>% t(), 1, function(x) x/sum(x))),
      t(apply(fpm(heq_lab_nothresh$dds)[ultramafic_nonultramafic_DEGs_large_effect,] %>% t(), 1, function(x) x/sum(x))),
      t(apply(fpm(rev_imp_nothresh$dds)[ultramafic_nonultramafic_DEGs_large_effect,] %>% t(), 1, function(x) x/sum(x)))) %>% 
  pheatmap(scale="row", cluster_rows = FALSE, annotation_row=annotation_species)

a <- apply(fpm(cal_spn_nothresh$dds)[ultramafic_nonultramafic_DEGs_large_effect,], 1, function(x) x/sum(x))
b <- apply(fpm(heq_lab_nothresh$dds)[ultramafic_nonultramafic_DEGs_large_effect,], 1, function(x) x/sum(x))
c <- apply(fpm(rev_imp_nothresh$dds)[ultramafic_nonultramafic_DEGs_large_effect,], 1, function(x) x/sum(x))


#################################################################################
#       get calciphila-spn genes to overlay with those found by Florian         # 
#################################################################################

cal_spn_degs <- cal_spn_nothresh$results %>% 
  data.frame() %>% 
  filter(abs(log2FoldChange) > 2 & padj < 0.05) %>% 
  rownames()

braker_vieillardii %>% 
  filter(annotation %in% cal_spn_degs) %>% 
  dplyr::select(seqname, start, end, annotation) %>% 
  data.frame() %>%
  write.table("cal_spn_degs.bed", quote=FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")


# found the following genes to overap between florians set and mine:
florian_mine_genes <- c("g931.t1", "g1801.t1", "g1852.t1", "g1855.t1", "g1856.t1", "g1861.t1", "g1877.t2", "g1878.t1", "g1893.t1", "g11936.t1", "g11962.t1", "g12622.t1", "g12622.t1", "g14809.t1", "g14824.t1", "g14835.t1", "g16635.t1", "g20731.t1", "g20783.t1", "g20791.t1", "g20832.t1", "g21722.t1", "g23287.t1", "g23932.t2")

florian_mine_genes %>% get_enriched_terms(mp, return_sample_GOData=TRUE)

mp[florian_mine_genes] %>% unlist() %>% unname()


#############################################################################################################################
#                       get genes DE between all ultramafic vs nonultramafic species                                        #
#     Then need to find their homologs in each genome rev and imp using orthofinder output                                  #
#############################################################################################################################

# get orthogroup info
orthogroups_long<-get_long_orthogroups() %>% drop_na()

# get DEG orthogroups and the vieillardii gene that was DE in each
ultramafic_OGs<-orthogroups_long %>% filter(gene %in% (ultramaf_nonultramaf_intersect %>% str_replace(".t1", "")) & species == "D. vieillardii")

# read in orthogroup gene trees for DEGs and large effect DEGs
ultramafic_OG_trees <- list.files(path="/Users/katieemelianova/Desktop/Diospyros/diospyros_gene_family_analysis/diospyros_gene_family_analysis/fastas/OrthoFinder/Results_Feb17/Gene_Trees", 
           full.names = TRUE, 
           pattern=paste(ultramafic_OGs$Orthogroup, collapse="|"))
# somehow not all the orthogroups have a tree file, so can't guarantee that I can purrr::set_names using the OG names used to fish them out
# instead getting the names of the orthogroups In was ablw to retrieve and using that to set names
ultramafic_OG_trees_names <- str_split_i(ultramafic_OG_trees, "/", 12) %>% str_split_i("_", 1)
ultramafic_OG_trees %<>% 
  lapply(ggtree::read.tree) %>%
  purrr::set_names(ultramafic_OG_trees_names)


##########################################################################################################
#           draw trees of orthogroups containing DEGs, hghlighting monophyletic homologs                 #
##########################################################################################################


draw_highlighted_genetree<-function(orthogroup, focal_gene){
  focal_orthogroup <- ultramafic_OG_trees[[orthogroup]]
  
  # get the vieillardii gene which is the DEG
  # removed this line because sometomes 2 DEGs are in a single orthogroup
  # now instead passing focal tip from the ultramafic_OGs data frame explicitly
  focal_tip <<- focal_orthogroup$tip.label %>% str_subset(focal_gene)

  # get the parent node of the vieillardii DEG in tree
  focal_parent_node <- focal_orthogroup %>% as_tibble() %>% filter(label == focal_tip) %>% pull(parent)
  focal_grandparent_node <- focal_orthogroup %>% as_tibble() %>% filter(node == focal_parent_node) %>% pull(parent)
  
  
  # get descendant nodes and tip names of DE viellardii gene
  descendant<-getDescendants(focal_orthogroup, focal_grandparent_node)
  
  # use these to get the labels (i.e. tip labels) which have revolutissima or impolita in the name (i.e. homologs of the viellardii DEGs)
  focal_homolog_tips_impolita <<- focal_orthogroup %>% 
    as_tibble() %>% 
    filter(node %in% descendant) %>% 
    filter(grepl(c("impolita"), label)) %>%
    pull(label)
  
  focal_homolog_tips_revolutissima <<- focal_orthogroup %>% 
    as_tibble() %>% 
    filter(node %in% descendant) %>% 
    filter(grepl(c("revolutissima"), label)) %>%
    pull(label)
  
  # plot gene tree with highlighted tips
  output_tree <- ggtree(focal_orthogroup) + 
    xlim(0, 1) +
    geom_tiplab(size=3, aes(subset = isTip & !(label %in% c(focal_tip, focal_homolog_tips_impolita, focal_homolog_tips_revolutissima)))) +
    theme(legend.position="none") +
    geom_tiplab(aes(subset = label == focal_tip), size=3, colour = 'firebrick', fontface="bold") +
    geom_tiplab(aes(subset = (label %in% focal_homolog_tips_impolita)), size=3, colour = 'dodgerblue2', fontface="bold") +
    geom_tiplab(aes(subset = (label %in% focal_homolog_tips_revolutissima)), size=3, colour = 'forestgreen', fontface="bold")
  
  return(list(output_tree=output_tree, focal_tip=focal_tip, focal_homolog_tips_impolita=focal_homolog_tips_impolita, focal_homolog_tips_revolutissima=focal_homolog_tips_revolutissima))
}


draw_highlighted_genetree("OG0000363", "g19147")
draw_highlighted_genetree("OG0000066", "g4374")


# this loop gets the orthogroup and the vieillardii gene which is found to be DE between impolita and revolutissima
# It then records the impolita and vieillardii orthologs associated with the vieillardii gene (anything within the clade stemming from the grandparental node of vieillardii gene), if there are any
##### here, should I exclude any orthogroups whcih do no have members of both revolutissima and impoita there? I think I should #####

homologs <- data.frame()
for (i in 1:nrow(ultramafic_OGs)) {
  current_row <- ultramafic_OGs[i,]
  tree <- draw_highlighted_genetree(current_row$Orthogroup, current_row$gene)
  homologs_impolita <- data.frame(homolog=tree$focal_homolog_tips_impolita)
  homologs_revolutissima <- data.frame(homolog=tree$focal_homolog_tips_revolutissima)
  
  if(nrow(homologs_impolita) > 0){
    homologs_impolita %<>% mutate(orthogroup=current_row$Orthogroup, vieillardii_homolog=current_row$gene, species="impolita")
    homologs <- rbind(homologs, homologs_impolita)
  }
  
  if(nrow(homologs_revolutissima) > 0){
    homologs_revolutissima %<>% mutate(orthogroup=current_row$Orthogroup, vieillardii_homolog=current_row$gene, species="revolutissima")
    homologs <- rbind(homologs, homologs_revolutissima)
  }
}

# this bit converts the format of revolutissima and impolita homologs to be able to filter GRanges objects with them
homologs$homolog %<>% str_split_i("_", 4)

# get impolita and revolutissima homologs of DEGs common to ultramafic-nonultramafic pairs
imp_hom <- homologs %>% filter(species == "impolita")
rev_hom <- homologs %>% filter(species == "revolutissima")

# get the GRanges in the genomes of impolita and revolutissima of DEG homologs
imp_hom_granges <- gene_granges$impolita %>% plyranges::filter(annotation %in% imp_hom$homolog)
rev_hom_granges <- gene_granges$revolutissima %>% plyranges::filter(annotation %in% rev_hom$homolog)

# using the granges objects, ask which DEG homologs in impolita and revolutissima have a TE within 1000bp
# get the IDs of the impolita and revolutissima genes tha fit this croteria
imp_te_proximal <- imp_hom_granges[findOverlaps(imp_hom_granges, te_intact_granges$impolita,  maxgap = 1000) %>% data.frame() %>% pull(queryHits)]
rev_te_proximal <- rev_hom_granges[findOverlaps(rev_hom_granges, te_intact_granges$revolutissima,  maxgap = 1000) %>% data.frame() %>% pull(queryHits)]


############################################################################################################################
#           now we can plot the same plot as before, highlighting only genes which have a TE next to them                  #
#             yellow is te within 1kb of impolita, blue is te within 1kb of revolutissima, pink is te within 1kb of both   #
############################################################################################################################

average_expression_plot_teprox<-inner_join(inner_join(
  data.frame(cal_av, pic_av) %>% rownames_to_column(), 
  data.frame(heq_av, lab_av) %>% rownames_to_column()), 
  data.frame(imp_av, rev_av) %>% rownames_to_column()) %>% 
  filter(rowname %in% ultramaf_nonultramaf_intersect) %>% 
  melt() %>%
  mutate(rowname=str_replace(rowname, ".t1", "")) %>%
  mutate(rowname=str_replace(rowname, ".t2", "")) %>%
  mutate(rowname=str_replace(rowname, ".t3", "")) %>%
  mutate(pair_name=case_when(variable %in% c("cal_av", "pic_av") ~ "Sp. Pic N'Ga - Calciphila",
                             variable %in% c("heq_av", "lab_av") ~ "Hequetiae - Labillardierei",
                             variable %in% c("rev_av", "imp_av") ~ "Revolutissima - Impolita"),
         soiltype=case_when(variable %in% c("rev_av", "heq_av", "pic_av") ~ "Ultramafic",
                            variable %in% c("cal_av", "lab_av", "imp_av") ~ "Non-Ultramafic")) %>%
  ggplot(aes(x=rowname, y=sqrt(value))) + 
  
  geom_rect(xmin = 2.5, xmax = 3.5, ymin = 0, ymax = 800, fill = 'darkseagreen1', alpha = 0.02) +
  geom_rect(xmin = 16.5, xmax = 18.5, ymin = 0, ymax = 800, fill = 'darkseagreen1', alpha = 0.02) +
  geom_rect(xmin = 43.5, xmax = 44.5, ymin = 0, ymax = 800, fill = 'darkseagreen1', alpha = 0.02) +
  geom_rect(xmin = 28.5, xmax = 29.5, ymin = 0, ymax = 800, fill = 'mistyrose', alpha = 0.1) +
  geom_rect(xmin = 30.5, xmax = 31.5, ymin = 0, ymax = 800, fill = 'mistyrose', alpha = 0.1) +
  geom_rect(xmin = 7.5, xmax = 8.5, ymin = 0, ymax = 800, fill = 'cornsilk2', alpha = 0.1) +
  geom_rect(xmin = 19.5, xmax = 20.5, ymin = 0, ymax = 800, fill = 'cornsilk2', alpha = 0.1) +
  geom_rect(xmin = 33.5, xmax = 34.5,, ymin = 0, ymax = 800, fill = 'cornsilk2', alpha = 0.1) +
  geom_point(size= 6, alpha=0.55, aes(shape=pair_name, color=soiltype)) +
  geom_point(size= 6, alpha=0.55, aes(shape=pair_name, color=soiltype)) + 
  scale_shape_manual(values=c(17, 19, 15)) +
  scale_colour_manual(values=c("darkolivegreen3", "hotpink2")) +
  ylab("sqrt(Counts per Million)") +
  xlab("Gene") + 
  theme(axis.text = element_text(size=10), 
        axis.title = element_text(size=15),
        legend.text = element_text(size=10),
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.3),
        panel.background = element_rect(fill = 'white', colour = 'grey72'),
        legend.position=c(.85, 0.8))

#pdf("average_expressionDEGs_te_proximal.pdf", width=12, height=6)
average_expression_plot_teprox
#dev.off()



###########################################################################################################
#     Now that we have genes of interest, we can plot the replicate expression values (not averages)
#        this can be plotted fin a boxplot for imp and rev, alongside 
###########################################################################################################



# get vieillardii homologs of impolita and revolutissima which have a TE proximal
# **NOTE** that I am combining veillardii homologs of imp and rev together
# therefore some items in the list will not have a proximal TE in one species (but will in the other)
te_prox_vieillardii_homologs <- c(homologs %>% filter(species == "impolita" & homolog %in% imp_te_proximal$annotation) %>% pull(vieillardii_homolog), 
                                  homologs %>% filter(species == "revolutissima" & homolog %in% rev_te_proximal$annotation) %>% pull(vieillardii_homolog)) %>%
                                  unique()


te_prox_expression <- rbind(fpm(rev_imp$dds) %>% data.frame() %>% 
             dplyr::select(imp_3ba1, imp_3ba2, imp_3bb, imp_3bc, imp_3bd1, imp_3bd2, imp_3bd2r) %>% 
             rownames_to_column(var="gene_id") %>% 
             mutate(gene_id = str_split_i(gene_id, "\\.", 1)) %>%
             filter(gene_id %in% te_prox_vieillardii_homologs) %>%
               melt() %>%
               mutate(species = "revolutissima"),
           fpm(rev_imp$dds) %>% data.frame() %>% 
             dplyr::select(rev_24a1, rev_24a1r, rev_24a2, rev_24b1, rev_24b2, rev_25a1, rev_25a2, rev_25a3, rev_25aq3r) %>% 
             rownames_to_column(var="gene_id") %>% 
             mutate(gene_id = str_split_i(gene_id, "\\.", 1)) %>%
             filter(gene_id %in% te_prox_vieillardii_homologs) %>%
             melt() %>%
             mutate(species = "impolita")) %>%
  set_colnames(c("gene_id", "individual", "cpm", "species"))


plot_te_prox_expression <- function(gene_id){
  expression_boxplot <- te_prox_expression %>% 
    filter(gene_id == !!gene_id) %>% 
    ggplot(aes(x=species, y=log(cpm), fill=species)) + 
    geom_boxplot() +
    facet_wrap(~gene_id) +
    xlab("") + 
    ylab("") + 
    theme(strip.background = element_rect(fill = "grey95"),
      #strip.text.x = element_blank(),
      strip.text.x = element_text(size=20),
      axis.title = element_text(size=20),
      axis.text = element_text(size=18),
      legend.position = "none",
      panel.background = element_rect(fill = "white"),
      panel.grid = element_line(color = "grey85")) +
    scale_fill_manual(values=c("darkolivegreen3", "hotpink2"))
  
  return(expression_boxplot)
}


g11472_expression <- plot_te_prox_expression("g11472")
g13660_expression <- plot_te_prox_expression("g13660")
g19147_expression <- plot_te_prox_expression("g19147")
g19148_expression <- plot_te_prox_expression("g19148")
g21937_expression <- plot_te_prox_expression("g21937")
g594_expression <- plot_te_prox_expression("g594")
g7857_expression <- plot_te_prox_expression("g7857")
g9428_expression <- plot_te_prox_expression("g9428")














###################################################
#                  gggenes plots                  #
###################################################


# need to find the correct homologs which cant be done automatically
# use these as input along with vieillardii homolog
# outout can be used to plot gggenes

setup_gggenes <- function(impolita_homolog, revolutissima_homolog, vieillardii_homolog){
  
  # get the granges of the imp and rev homologs
  revolutissima_gene <- gene_granges$revolutissima %>% plyranges::filter(annotation == revolutissima_homolog)
  impolita_gene <- gene_granges$impolita %>% plyranges::filter(annotation == impolita_homolog)
  
  # overlap these with their respective TE granges objects
  revolutissima_te <- findOverlaps(revolutissima_gene, te_intact_granges$revolutissima,  maxgap = 1000)
  impolita_te <- findOverlaps(impolita_gene, te_intact_granges$impolita,  maxgap = 1000)
  
  # extract the grange of any TEs within 1000bp using the subject of the overlaps
  revolutissima_te <- te_intact_granges$revolutissima[revolutissima_te@to]
  impolita_te <- te_intact_granges$impolita[impolita_te@to]
  
  # combine the TE granges of revolutissima and impolita, setting gene slot to "TE"
  te <- rbind(revolutissima_te %>% data.frame() %>% mutate(gene="TE"), 
              impolita_te %>% data.frame() %>% mutate(gene="TE")) %>% 
    dplyr::select(c(seqnames, gene, start, end, ph2))
  
  # combine gene granges of revolutissima and impolita, setting the gene slot to the name of the gene
  gene <- rbind(revolutissima_gene %>% data.frame() %>% mutate(gene=vieillardii_homolog),
                impolita_gene %>% data.frame() %>% mutate(gene=vieillardii_homolog)) %>% 
    dplyr::select(c(seqnames, gene, start, end, ph2))
  
  ## now combine both TE and gene granges, using the gene granges seqnames to infer which species is which
  toplot <- rbind(te, gene) %>% 
    mutate(species = case_when(seqnames == as.character(revolutissima_gene@seqnames) ~ "revolutissima", 
                               seqnames == as.character(impolita_gene@seqnames) ~ "impolita"))
  
  # make a dummy so that the gggenes plot is lined up
  dummies <- make_alignment_dummies(
    toplot,
    aes(xmin = start, xmax = end, y = species, id = gene),
    on = unique(vieillardii_homolog))  
  #  
  to_return=list(dummy=dummies, toplot=toplot)
  return(to_return)
}

###################################
#       gggenes plot function     #
###################################



plot_gggenes <- function(input_list){
  ggplot2::ggplot(input_list$toplot, ggplot2::aes(xmin = start, xmax = end, y = species, fill = interaction(gene, species, sep=':'), label = gene)) +
  geom_gene_arrow(arrow_body_height = grid::unit(10, "mm"),
                    arrowhead_height = grid::unit(12, "mm")) +
    #geom_gene_label(height = grid::unit(6, "mm"), grow = TRUE) +
    ggplot2::facet_wrap(~ species, ncol = 1, scales = "free") +
    theme_genes() +
    theme(legend.text = element_text(size=18),
          axis.text.x = element_blank(),
          line = element_blank(),
          #axis.text.x = element_text(size=11),
          axis.text.y = element_text(size=20),
          #axis.title = element_text(size=20),
          legend.title = element_blank(),
          legend.position = "none",
          
          axis.title.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.line.x = element_blank(),
          plot.margin = unit(c(2,0,2,0), "cm")) +
    geom_blank(data = input_list$dummy) +
    ylab("")
}

#################################
#           g11472              #
#################################

vieillardii_homolog <- "g11472"
impolita_homolog <- "g14194.t1"
revolutissima_homolog <- "g18897.t1"
g11472_toplot <- setup_gggenes(impolita_homolog, revolutissima_homolog, vieillardii_homolog)
g11472_gggenes <- plot_gggenes(g11472_toplot) + scale_fill_manual(values=c("darkolivegreen3", "grey76", "hotpink2", "grey76"))


#################################
#           g13660              #
#################################


#ultramafic_OGs %>% filter(gene == vieillardii_homolog)
#test <- draw_highlighted_genetree("OG0000292", "g13660")

vieillardii_homolog <- "g13660"
impolita_homolog <- "g13031.t1"
revolutissima_homolog <- "g1268.t1"
g13660_toplot <- setup_gggenes(impolita_homolog, revolutissima_homolog, vieillardii_homolog)
g13660_gggenes <- plot_gggenes(g13660_toplot) + scale_fill_manual(values=c("darkolivegreen3", "hotpink2", "grey76"))



#################################
#           g19147              #
#################################

#ultramafic_OGs %>% filter(gene == vieillardii_homolog)
#draw_highlighted_genetree("OG0000363", "g19147")

vieillardii_homolog <- "g19147"
impolita_homolog <- "g5404.t1"
revolutissima_homolog <- "g18409.t1"
g19147_toplot <- setup_gggenes(impolita_homolog, revolutissima_homolog, vieillardii_homolog)
g19147_gggenes <- plot_gggenes(g19147_toplot) + scale_fill_manual(values=c("darkolivegreen3", "grey76", "hotpink2"))

#################################
#           g19148              #
#################################

# same homologs as previous

#ultramafic_OGs %>% filter(gene == vieillardii_homolog)
#draw_highlighted_genetree("OG0000363", "g19148")

vieillardii_homolog <- "g19148"
impolita_homolog <- "g5404.t1"
revolutissima_homolog <- "g18409.t1"
g19148_toplot <- setup_gggenes(impolita_homolog, revolutissima_homolog, vieillardii_homolog)
g19148_gggenes <- plot_gggenes(g19148_toplot) + scale_fill_manual(values=c("darkolivegreen3", "grey76", "hotpink2"))


#################################
#           g21937              #
#################################

#ultramafic_OGs %>% filter(gene == vieillardii_homolog)
#draw_highlighted_genetree("OG0011205", "g21937")

vieillardii_homolog <- "g21937"
impolita_homolog <- "g17749.t1"
revolutissima_homolog <- "g6984.t1"
g21937_toplot <- setup_gggenes(impolita_homolog, revolutissima_homolog, vieillardii_homolog)
g21937_gggenes <- plot_gggenes(g21937_toplot) + scale_fill_manual(values=c("darkolivegreen3", "grey76", "hotpink2", "grey76"))

#################################
#           g4374               # unresolveable
#################################

#ultramafic_OGs %>% filter(gene == vieillardii_homolog)
#draw_highlighted_genetree("OG0000066", "g4374")
#
#vieillardii_homolog <- "g4374"
#impolita_homolog <- ""
#revolutissima_homolog <- ""

#################################
#           g594                #
#################################

#ultramafic_OGs %>% filter(gene == vieillardii_homolog)
#draw_highlighted_genetree("OG0006346", "g594")

vieillardii_homolog <- "g594"
impolita_homolog <- "g3611.t1"
revolutissima_homolog <- "g14832.t1"
g594_toplot <- setup_gggenes(impolita_homolog, revolutissima_homolog, vieillardii_homolog)
g594_gggenes <- plot_gggenes(g594_toplot) + scale_fill_manual(values=c("darkolivegreen3", "hotpink2", "grey76"))


#################################
#           g7857               #
#################################

#ultramafic_OGs %>% filter(gene == vieillardii_homolog)
#draw_highlighted_genetree("OG0000219", "g7857")

vieillardii_homolog <- "g7857"
impolita_homolog <- "g25993.t1"
#impolita_homolog <- "g26603.t1"
revolutissima_homolog <- "g24808.t1"
g7857_toplot <- setup_gggenes(impolita_homolog, revolutissima_homolog, vieillardii_homolog)
# not clear which homolog has been flagged as being TE proximal
homologs %>% filter(vieillardii_homolog == "g7857" & homolog %in% imp_te_proximal$annotation)
g7857_gggenes <- plot_gggenes(g7857_toplot) + scale_fill_manual(values=c("darkolivegreen3", "grey76", "hotpink2", "grey76"))

#################################
#           g9428               #
#################################

#ultramafic_OGs %>% filter(gene == vieillardii_homolog)
#draw_highlighted_genetree("OG0000291", "g9428")

vieillardii_homolog <- "g9428"
impolita_homolog <- "g7821.t1"
revolutissima_homolog <- "g25092.t1"
g9428_toplot <- setup_gggenes(impolita_homolog, revolutissima_homolog, vieillardii_homolog)
g9428_gggenes <- plot_gggenes(g9428_toplot) + scale_fill_manual(values=c("darkolivegreen3", "grey76", "hotpink2"))






g11472_grid <- plot_grid(g11472_expression, g11472_gggenes, rel_widths = c(1.2,2))
g13660_grid <- plot_grid(g13660_expression, g13660_gggenes, rel_widths = c(1.2,2))
g19147_grid <- plot_grid(g19147_expression, g19147_gggenes, rel_widths = c(1.2,2))
g19148_grid <- plot_grid(g19148_expression, g19148_gggenes, rel_widths = c(1.2,2))
g21937_grid <- plot_grid(g21937_expression, g21937_gggenes, rel_widths = c(1.2,2))
g594_grid <- plot_grid(g594_expression, g594_gggenes, rel_widths = c(1.2,2))
g7857_grid <- plot_grid(g7857_expression, g7857_gggenes, rel_widths = c(1.2,2))
g9428_grid <- plot_grid(g9428_expression, g9428_gggenes, rel_widths = c(1.2,2))



all_grid <- plot_grid(g11472_grid,
          g13660_grid,
          g19147_grid,
          g19148_grid,
          g21937_grid,
          g594_grid,
          g7857_grid,
          g9428_grid, ncol = 1)

pdf("test_panel.pdf", height=18, width=10)
all_grid
dev.off()



































revolutissima_g11314_intact_TE <- findOverlaps(revolutissima_g11314_gene, te_intact_granges$revolutissima,  maxgap = 35000)
vieillardii_g1861_intact_TE <- findOverlaps(vieillardii_g1861_gene, te_intact_granges$vieillardii,  maxgap = 35000)
impolita_g4621_intact_TE <- findOverlaps(impolita_g4621_gene, te_intact_granges$impolita,  maxgap = 35000)

revolutissima_g11314_intact_TE <- te_intact_granges$revolutissima[revolutissima_g11314_intact_TE@to]
vieillardii_g1861_intact_TE <- te_intact_granges$vieillardii[vieillardii_g1861_intact_TE@to]
impolita_g4621_intact_TE <- te_intact_granges$impolita[impolita_g4621_intact_TE@to]

test_te <- rbind(revolutissima_g11314_intact_TE %>% data.frame() %>% mutate(gene="TE"), 
                 vieillardii_g1861_intact_TE %>% data.frame() %>% mutate(gene="TE"), 
                 impolita_g4621_intact_TE %>% data.frame()%>% mutate(gene="TE")) %>% 
  dplyr::select(c(seqnames, gene, start, end, ph2))


test_gene <- rbind(revolutissima_g11314_gene %>% data.frame() %>% mutate(gene="Annexin"),
                   vieillardii_g1861_gene %>% data.frame() %>% mutate(gene="Annexin"),
                   impolita_g4621_gene %>% data.frame() %>% mutate(gene="Annexin")) %>% 
  dplyr::select(c(seqnames, gene, start, end, ph2))

test <- rbind(test_te, test_gene) %>% mutate(species = case_when(seqnames == "ptg000017l" ~ "revolutissima",
                                                                 seqnames == "ptg000002l" ~ "vieillardii",
                                                                 seqnames == "Scaffolds_1151" ~ "impolita"))

dummies <- make_alignment_dummies(
  test,
  aes(xmin = start, xmax = end, y = species, id = gene),
  on = "Annexin"
)

#pdf("annexin_OG0000336_geneplot.pdf", height=5.5, width=9)
ggplot2::ggplot(test, ggplot2::aes(xmin = start, xmax = end,
                                   y = species, fill = gene, label = gene)) +
  geom_gene_arrow(arrow_body_height = grid::unit(10, "mm"),
                  arrowhead_height = grid::unit(12, "mm")) +
  #geom_gene_label(height = grid::unit(6, "mm"), grow = TRUE) +
  ggplot2::facet_wrap(~ species, ncol = 1, scales = "free") +
  theme_genes() +
  theme(legend.text = element_text(size=20),
        axis.text.x = element_text(size=11),
        axis.text.y = element_text(size=14),
        axis.title = element_text(size=20),
        legend.title = element_blank()) +
  geom_blank(data = dummies) +
  ylab("Species")
#dev.off()





















## to return to later when looking at diostance to TE for homoogs of DE genes
distanceToNearest(imp_hom_granges, te_intact_granges$impolita, ignore.strand=FALSE) 
distanceToNearest(rev_hom_granges, te_intact_granges$revolutissima, ignore.strand=FALSE)
rbind(distanceToNearest(imp_hom_granges, te_intact_granges$impolita, ignore.strand=FALSE) %>% data.frame() %>% mutate(species="impolita"),  
      distanceToNearest(rev_hom_granges, te_intact_granges$revolutissima, ignore.strand=FALSE) %>% data.frame() %>% mutate(species="revolutissima")) %>%
  ggplot(aes(x=species, y=distance)) + 
  geom_boxplot()
rbind(distanceToNearest(imp_hom_granges, te_intact_granges$impolita, ignore.strand=FALSE) %>% data.frame() %>% mutate(species="impolita"),  
      distanceToNearest(rev_hom_granges, te_intact_granges$revolutissima, ignore.strand=FALSE) %>% data.frame() %>% mutate(species="revolutissima")) %>%
  group_by(species) %>%
  summarise(meandist=median(distance))
te_intact_granges$impolita[distanceToNearest(imp_hom_granges, te_intact_granges$impolita, ignore.strand=FALSE)$subjectHits]




















braker_annotation <- get_braker_annotation()



brakr_vieillardii <- braker_annotation$vieillardii_braker
brakr_revolutissima <- braker_annotation$revolutissima_braker
brakr_impolita <- braker_annotation$impolita_braker



brakr_vieillardii %>% filter(annotation == "g1861")
brakr_revolutissima %>% filter(annotation == "g11314")
brakr_impolita %>% filter(annotation == "g4621")




library(gggenes)

ggplot2::ggplot(example_genes, ggplot2::aes(xmin = start, xmax = end,
                                            y = molecule, fill = gene, label = gene)) +
  geom_gene_arrow() +
  geom_gene_label() +
  ggplot2::facet_wrap(~ molecule, ncol = 1, scales = "free") +
  theme_genes()


example_genes %>% filter(gene == "genA")












impolita_DE_orthogroup_genes_High[impolita_DE_orthogroup_genes_High %in% c("g4621", "g22665", "g3611")]
revolutissima_DE_orthogroup_genes_High[revolutissima_DE_orthogroup_genes_High %in% c("g11314", "g4181", "g482", "g14832")]









df<-data.frame(start=c(594198,596540,598457,600085,983488,984345),
               stop=c(596450,598423,600070,601182,984336,986495),
               species=rep("Ferriphaselus amnicola",6),
               gene=c("gene1","gene2","gene3","gene4","gene5","gene6"))


ggplot(df, aes(xmin = start, xmax = stop, y = species, fill = gene)) +
  geom_gene_arrow() +
  facet_wrap(~ species, scales = "free", ncol = 1) +
  scale_fill_brewer(palette = "Set3") +
  theme_genes()





konecna<-c("g10708.t1", "g10727.t1", "g11122.t1", "g11665.t1", "g11665.t2", 
           "g12060.t1", "g1210.t1", "g12994.t1", "g13037.t1", "g13175.t1", 
           "g13379.t1", "g13737.t1", "g13737.t2", "g13958.t1", "g14806.t1", 
           "g15288.t1", "g16171.t1", "g16602.t1", "g17371.t1", "g18223.t1", 
           "g21039.t1", "g21512.t1", "g22372.t1", "g22372.t2", "g22540.t1", 
           "g23219.t1", "g23867.t1", "g24061.t1", "g24061.t2", "g24344.t1", 
           "g24433.t1", "g24434.t1", "g24634.t1", "g24804.t1", "g24971.t1", 
           "g3080.t1", "g3736.t1", "g3827.t1", "g3943.t1", "g4026.t1", 
           "g4106.t1", "g4179.t1", "g4724.t1", "g4959.t1", "g5140.t1", 
           "g5140.t2", "g5712.t1", "g6209.t1", "g6518.t1", "g6519.t1", 
           "g6664.t1", "g6687.t1", "g6880.t1", "g739.t1", "g9014.t1", "g9247.t1")

turner<-c("g10554.t1", "g11922.t1", "g1383.t1", "g16134.t1", "g21710.t1", "g23837.t1",
          "g23837.t2", "g24433.t1", "g24434.t1", "g264.t1", "g3633.t1", "g4959.t1", 
          "g6362.t1", "g7194.t1", "g7406.t1", "g9069.t1")







# none of these genes come up as DE - I think because theyre all zero?
intersect((assay(heq_imp$dds)), te_gene_mine)

# lets check - not all of them are zeo!
all_counts[te_gene_mine,] %>% 
  pheatmap::pheatmap(cluster_rows=T, cluster_cols=F)

pheatmap::pheatmap(assay(cal_spn$dds)[te_gene_mine,], 
                   scale = "row", 
                   cluster_cols = FALSE,
                   cluster_rows = TRUE,
                   show_rownames = FALSE,
                   treeheight_row = 0, treeheight_col = 0)

heq_imp_k<-data.frame(results(heq_imp$dds, contrast=c("species", "heq", "imp")))[konecna_de,] %>% dplyr::select(log2FoldChange) %>% mutate(pair="heq_imp")
cal_spn_k<-data.frame(results(cal_spn$dds, contrast=c("species", "spn", "cal")))[konecna_de,] %>% dplyr::select(log2FoldChange) %>% mutate(pair="cal_spn")

imp_cal_k<-data.frame(results(imp_cal$dds, contrast=c("species", "cal", "imp")))[konecna_de,] %>% dplyr::select(log2FoldChange) %>% mutate(pair="imp_cal")
heq_spn_k<-data.frame(results(heq_spn$dds, contrast=c("species", "spn", "heq")))[konecna_de,] %>% dplyr::select(log2FoldChange) %>% mutate(pair="heq_spn")


# ask how many genes are DE between each species pair, starting with ultramafic-volcanic pairs
# then do volcanic-volcanic pairs
# ask how many of these genes are in the konecna gene set
data.frame(results(heq_imp$dds, contrast=c("species", "heq", "imp"))) %>% filter(abs(log2FoldChange) > 1 & as.numeric(padj) < 0.05) %>% rownames_to_column(var="gene_id") %>% filter(gene_id %in% konecna) 
data.frame(results(cal_spn$dds, contrast=c("species", "spn", "cal"))) %>% filter(abs(log2FoldChange) > 1 & as.numeric(padj) < 0.05) %>% rownames_to_column(var="gene_id") %>% filter(gene_id %in% konecna) 
data.frame(results(imp_cal$dds, contrast=c("species", "cal", "imp"))) %>% filter(abs(log2FoldChange) > 1 & as.numeric(padj) < 0.05) %>% rownames_to_column(var="gene_id") %>% filter(gene_id %in% konecna) 
data.frame(results(heq_spn$dds, contrast=c("species", "spn", "heq"))) %>% filter(abs(log2FoldChange) > 1 & as.numeric(padj) < 0.05) %>% rownames_to_column(var="gene_id") %>% filter(gene_id %in% konecna)

# ask how many genes are DE between each species pair, starting with ultramafic-volcanic pairs
# then do volcanic-volcanic pairs
# ask how many of these genes are in the TURNER gene set
data.frame(results(heq_imp$dds, contrast=c("species", "heq", "imp"))) %>% filter(abs(log2FoldChange) > 1 & as.numeric(padj) < 0.05) %>% rownames_to_column(var="gene_id") %>% filter(gene_id %in% turner) 
data.frame(results(cal_spn$dds, contrast=c("species", "spn", "cal"))) %>% filter(abs(log2FoldChange) > 1 & as.numeric(padj) < 0.05) %>% rownames_to_column(var="gene_id") %>% filter(gene_id %in% turner) 
data.frame(results(imp_cal$dds, contrast=c("species", "cal", "imp"))) %>% filter(abs(log2FoldChange) > 1 & as.numeric(padj) < 0.05) %>% rownames_to_column(var="gene_id") %>% filter(gene_id %in% turner) 
data.frame(results(heq_spn$dds, contrast=c("species", "spn", "heq"))) %>% filter(abs(log2FoldChange) > 1 & as.numeric(padj) < 0.05) %>% rownames_to_column(var="gene_id") %>% filter(gene_id %in% turner)



# check expression levcels of the one DEG bertween calciphila and spn
assay(cal_spn$dds)[c("g9247.t1"),]

# get the peptide sequence from D. vieillardii

#>g9247.t1
MAGMHAAAAVRMMAAGPSAGDEIEFGTAWWFMYAGISCFLIVFAGVMSGLTLGLMSLGRV
ELEILQSSGTSAEMRQAAAILPVVQKQYQLLVTLLLCNAGAMEALPIYLDKLCNQYVAIL
LSVTFVLAFGEVIPQAICTRHGLAVGANFVWLVKILMTICYPIAYPIGKILDCVLGHDEA
LFRRAQLKALVSIHSQEAGKGGELTHDETTIISGALDLSEKTAEAAMTPIESTFSLDVNS
KLDWEVMENILARGHSRVPVYSGNPKNIIGLLLVKCLLTVLPEKETPVSAVSIRRIPRVP
ADMPLYDILNEFQKGSSHMAAVVKAKGRKKNPPSSVTEGYTEEVRIYSANTEVNIPLLSE
QDEKSECFVVDIDKASKASTISKKTLSQYNDAARNGLLHMSEDIEDSEIIGIITLEDVFE
ELLQEEILDETDEFVDVHNRTSVAAAAAASLAVGCSPAVWRLSCYKGTGCQSKQGQAQAV
KKSAEEDPFTPKLTGEVGEPLLENKR*

# put the sequence into domain search, setting pval threshold to 0.0005
#https://www.genome.jp/tools/motif/
  
# this returns one hit: "PF01595, Cyclin M transmembrane N-terminal domain"

# google this, it comes up with iron ion homeostasis:
# https://www.ebi.ac.uk/interpro/entry/pfam/PF01595/

# normalise by numer of total DEGs
# Get BAMs from Teerna to get illumina mapped to see population wide TE insertions





# npw do the same for species pairs, asking how many genes are found in the turner gene set
data.frame(results(heq_imp$dds, contrast=c("species", "heq", "imp"))) %>% filter(abs(log2FoldChange) > 1 & as.numeric(padj) < 0.05) %>% rownames_to_column(var="gene_id") %>% filter(gene_id %in% turner) 
data.frame(results(cal_spn$dds, contrast=c("species", "spn", "cal"))) %>% filter(abs(log2FoldChange) > 1 & as.numeric(padj) < 0.05) %>% rownames_to_column(var="gene_id") %>% filter(gene_id %in% turner) 
data.frame(results(imp_cal$dds, contrast=c("species", "cal", "imp"))) %>% filter(abs(log2FoldChange) > 1 & as.numeric(padj) < 0.05) %>% rownames_to_column(var="gene_id") %>% filter(gene_id %in% turner) 
data.frame(results(heq_spn$dds, contrast=c("species", "spn", "heq"))) %>% filter(abs(log2FoldChange) > 1 & as.numeric(padj) < 0.05) %>% rownames_to_column(var="gene_id") %>% filter(gene_id %in% turner)

# get the konecna genes and ask what go terms are enriched in these genes
heq_imp_k_go<-data.frame(results(heq_imp$dds, contrast=c("species", "heq", "imp"))) %>% filter(abs(log2FoldChange) > 1.5 & as.numeric(padj) < 0.05) %>% rownames_to_column(var="gene_id") %>% filter(gene_id %in% konecna) %>% pull(gene_id) %>% get_enriched_terms(mp, return_sample_GOData=TRUE)
cal_spn_k_go<-data.frame(results(cal_spn$dds, contrast=c("species", "spn", "cal"))) %>% filter(abs(log2FoldChange) > 1.5 & as.numeric(padj) < 0.05) %>% rownames_to_column(var="gene_id") %>% filter(gene_id %in% konecna) %>% pull(gene_id) %>% get_enriched_terms(mp, return_sample_GOData=TRUE)
imp_cal_k_go<-data.frame(results(imp_cal$dds, contrast=c("species", "cal", "imp"))) %>% filter(abs(log2FoldChange) > 1.5 & as.numeric(padj) < 0.05) %>% rownames_to_column(var="gene_id") %>% filter(gene_id %in% konecna) %>% pull(gene_id) %>% get_enriched_terms(mp, return_sample_GOData=TRUE)
heq_spn_kj_go<-data.frame(results(heq_spn$dds, contrast=c("species", "spn", "heq"))) %>% filter(abs(log2FoldChange) > 1.5 & as.numeric(padj) < 0.05) %>% rownames_to_column(var="gene_id") %>% filter(gene_id %in% konecna) %>% pull(gene_id) %>% get_enriched_terms(mp, return_sample_GOData=TRUE)



heq_imp_k_go$result
# get mean log2fold change between species pairs
rbind(heq_imp_k,
      cal_spn_k,
      imp_cal_k,
      heq_spn_k) %>%
  drop_na() %>%
  group_by(pair) %>%
  summarise(Mean=mean(log2FoldChange))


# plot log2gold change between süecies pairs
ggplot(all_k, aes(x=log2FoldChange, fill = pair)) + geom_density(aes(y = after_stat(count)), alpha = 0.25) 

ggplot(all_k, aes(x=log2FoldChange, fill = pair)) + 
  geom_density(aes(y = after_stat(count)), alpha = 0.25) + facet_wrap(~pair)


# we need to plug in divergence into this to get an understanding of whether this is 
# more DE gegenes than would be expected given divergence time
get_significant_genes(cal_spn) %>% length()
get_significant_genes(heq_imp) %>% length()
get_significant_genes(imp_cal) %>% length()
get_significant_genes(heq_spn) %>% length()




###################################
#     Perform GO term enrichment  #
###################################

get_de_genes_in_term<-function(degs, go_term, go_data){
  genes_in_term<-genesInTerm(go_data, go_term)[[1]]
  degs_in_term<-intersect(genes_in_term, degs)
  return(degs_in_term)
}



cal_spn_enrichedGO<-cal_spn$results %>% data.frame() %>% filter(log2FoldChange > 2 & padj < 0.05) %>% rownames()  %>% get_enriched_terms(mp, return_sample_GOData=TRUE)
heq_imp_enrichedGO<-heq_imp$results %>% data.frame() %>% filter(log2FoldChange > 2 & padj < 0.05) %>% rownames()  %>% get_enriched_terms(mp, return_sample_GOData=TRUE)

imp_cal_enrichedGO<-imp_cal$results %>% data.frame() %>% filter(log2FoldChange > 2 & padj < 0.05) %>% rownames() %>% get_enriched_terms(mp, return_sample_GOData=TRUE)
imp_spn_enrichedGO<-imp_spn$results %>% data.frame() %>% filter(log2FoldChange > 2 & padj < 0.05) %>% rownames() %>% get_enriched_terms(mp, return_sample_GOData=TRUE)

heq_cal_enrichedGO<-heq_cal$results %>% data.frame() %>% filter(log2FoldChange > 2 & padj < 0.05) %>% rownames() %>% get_enriched_terms(mp, return_sample_GOData=TRUE)
heq_spn_enrichedGO<-heq_spn$results %>% data.frame() %>% filter(log2FoldChange > 2 & padj < 0.05) %>% rownames() %>% get_enriched_terms(mp, return_sample_GOData=TRUE)



head(cal_spn_enrichedGO$result)
head(heq_imp_enrichedGO$result)

head(imp_cal_enrichedGO$result)
head(heq_spn_enrichedGO$result)




intersect(intersect(heq_imp_degs$Term, imp_cal_degs$Term), cal_spn_degs$Term)



a<-get_de_genes_in_term((cal_spn_de$results %>% data.frame() %>% filter(log2FoldChange > 2 & padj < 0.05) %>% rownames()), "GO:0002215", cal_spn_degs$goData)
b<-get_de_genes_in_term((imp_cal_de$results %>% data.frame() %>% filter(log2FoldChange > 2 & padj < 0.05) %>% rownames()), "GO:0002215", imp_cal_degs$goData)


intersect(a, b)

test<-heq_imp$results %>% rownames() %>% sample(500)
get_enriched_terms(test, mp) %>% filter(as.numeric(classicFisher) < 0.05)

heq_imp_de$results %>% rownames() %>% sample(500)


assay(heq_imp$dds) %>% rownames() %>% tibble() %>% filter(. %in% konecna)






imp_test<-tibble(sample_id=c("imp_3ba1",
                             "imp_3ba2",
                             "imp_3bb",
                             "imp_3bc",
                             "imp_3bd1",
                             "imp_3bd2",
                             "imp_3bd2r"))

imp_test %>%
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
           into = c("mother", "seed"), 
           sep = "(?=[a-z])(?<=[0-9])") %>%
  column_to_rownames("to_rownames")






