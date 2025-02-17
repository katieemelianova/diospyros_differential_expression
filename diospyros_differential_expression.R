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
source("differential_expression_functions.R")
source("/Users/katieemelianova/Desktop/Diospyros/diospyros_R_functions/diospyros_soil_leaf_species_element_dataset.R")
source("/Users/katieemelianova/Desktop/Diospyros/diospyros_R_functions/diospyros_orthogroups_long.R")




############################################
#        load in the GO ID mappings        #
############################################

# for the GO term enrichment tests
mp<-readMappings("vieillardii_topGO_annotation.txt")
mp<-readMappings("/Users/katieemelianova/Desktop/Diospyros/vieillardii_functionalAnnotation_transcript.topGO.txt")

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

pdf("PCA_with_tree.pdf", width = 15, height = 7)
grid.arrange(pic1, pic2, nrow = 1)
dev.off()


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
listInput<-list(`calciphila vs sp. Pic N'Ga` = cal_spn$results %>% data.frame() %>% filter(padj < 0.05 & log2FoldChange > 2) %>% rownames(),
                `hequetiae vs labillardierei` = heq_lab$results %>% data.frame() %>% filter(padj < 0.05 & log2FoldChange > 2) %>% rownames(),
                `revolutissima vs impolita` = rev_imp$results %>% data.frame() %>% filter(padj < 0.05 & log2FoldChange > 2) %>% rownames())


# plot upset plot
upset(fromList(listInput), order.by = "freq", nsets = 5)

# get genes DE between each ultramafic-nonultramafic pair
ultramaf_nonultramaf_intersect <- intersect(intersect(listInput$`calciphila vs sp. Pic N'Ga`, listInput$`hequetiae vs labillardierei`), listInput$`revolutissima vs impolita`)


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


#pdf("euler_with_tree.pdf", width = 8, height = 7)
#grid.arrange(tree_pairs, euler_pairs, nrow = 1)
#dev.off()


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
  mutate(pair_name=case_when(variable %in% c("cal_av", "pic_av") ~ "Sp. Pic N'Ga - Calciphila",
                             variable %in% c("heq_av", "lab_av") ~ "Hequetiae - Labillardierei",
                             variable %in% c("rev_av", "imp_av") ~ "Revolutissima - Impolita"),
         soiltype=case_when(variable %in% c("rev_av", "heq_av", "pic_av") ~ "Ultramafic",
                            variable %in% c("cal_av", "lab_av", "imp_av") ~ "Non-Ultramafic")) %>%
  ggplot(aes(x=rowname, y=(value), fill=soiltype)) + 
  geom_point(size= 6, alpha=0.55, aes(shape=soiltype, color=pair_name)) +
  scale_shape_manual(values=c(17, 19)) +
  scale_colour_manual(values=c("darkolivegreen3", "hotpink2", "dodgerblue2")) +
  ylab("Counts per Million") +
  xlab("Gene") + 
  theme(axis.text = element_text(size=10), 
        axis.title = element_text(size=15),
        legend.text = element_text(size=13),
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 25, vjust = 0.5, hjust=0.3),
        panel.background = element_rect(fill = 'white', colour = 'grey72'),
        legend.position=c(.75, 0.8))

#pdf("average_expressionDEGs.pdf", width=6, height=6)
#average_expression_plot 
#dev.off()

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




######################################################################################################
#                  get subsets of DEGs - all DEGs and those of large effect                          #
######################################################################################################

# get the deseq results and get genes to focus on using baseMean as filter

# DEGs
ultramafic_nonultramafic_DEGs_large_effect <- results(rev_imp$dds)[ultramaf_nonultramaf_intersect,] %>% data.frame() %>% filter(baseMean > 1000 & log2FoldChange > 2) %>% rownames()
ultramafic_nonultramafic_DEGs <- results(rev_imp$dds)[ultramaf_nonultramaf_intersect,] %>% data.frame() %>% rownames()

# DEG counts
rev_imp_ultramafic_nonultramafic_counts_large_effect <- fpm(rev_imp$dds) %>% data.frame() %>% rownames_to_column() %>% filter(rowname %in% ultramafic_nonultramafic_DEGs_large_effect)
rev_imp_ultramafic_nonultramafic_counts <- fpm(rev_imp$dds) %>% data.frame() %>% rownames_to_column() %>% filter(rowname %in% ultramafic_nonultramafic_DEGs)



#############################################################################################################################
#                       get genes DE between all ultramafic vs nonultramafic species                                        #
#     Then need to find their homologs in each genome rev and imp using orthofinder output                                  #
#############################################################################################################################

# get orthogroup info
orthogroups_long<-get_long_orthogroups()

# get genes DE between all ultramafic vs non-ultramafic pairs
ultramafic_DEGs <- ultramaf_nonultramaf_intersect %>%
  str_replace(".t1", "")



# plot copy number of each ortogrop per species
ultramafic_OGs<-orthogroups_long %>% filter(gene %in% ultramafic_DEGs & species == "D. vieillardii") %>% pull(Orthogroup)
ultramafic_OG_copy_number<-orthogroups_long %>% filter(Orthogroup %in% ultramafic_OGs) %>% 
  group_by(Orthogroup, species) %>%
  summarise(gene_count = length(gene)) %>% 
  data.frame()
ultramafic_OG_copy_number$species <-factor(ultramafic_OG_copy_number$species, levels=c("D. sandwicensis", "D. vieillardii", "D. pancheri", "D. yahouensis", "D. impolita", "D. revolutissima"))
ggplot(ultramafic_OG_copy_number, aes(x=species, y=gene_count, group=Orthogroup, colour=Orthogroup)) +
  geom_line(size=1) +
  geom_point(size=6) +
  scale_colour_manual(values=c("red", "orange", "pink", "green", "blue", "grey", "black", "purple", "yellow", "brown", "dodgerblue2", "gold", "aquamarine"))





# get genes DE in all ulramafic-nonultramafic species apirs and their orthogrouops
impolita_ultramafic_orthogroup_genes<- orthogroups_long %>% filter(Orthogroup %in% ultramafic_OGs & species == "D. impolita" & gene != "NA") %>% pull(gene)
revolutissima_ultramafic_orthogroup_genes<- orthogroups_long %>% filter(Orthogroup %in% ultramafic_OGs & species == "D. revolutissima" & gene != "NA") %>% pull(gene)





# get DEG homologs in impolita and revolutissima
imp_rev_OGs<-orthogroups_long %>% filter(gene %in% ultramafic_DEGs & species == "D. vieillardii") %>% pull(Orthogroup)
impolita_DE_orthogroup_genes<-orthogroups_long %>% filter(Orthogroup %in% imp_rev_OGs & species == "D. impolita" & gene != "NA") %>% pull(gene)
revolutissima_DE_orthogroup_genes<-orthogroups_long %>% filter(Orthogroup %in% imp_rev_OGs & species == "D. revolutissima" & gene != "NA") %>% pull(gene)










impolita.gene_te %>% filter(gene %in% impolita_DE_orthogroup_genes) %>% 
  set_colnames(c("annotation", "gene", "gene_dist", "te_length", "TE")) %>% 
  left_join(impolita_te_annotation, by="TE", relationship = "many-to-many") %>% 
  left_join(orthogroups_long %>% filter(Orthogroup %in% imp_rev_OGs & species == "D. impolita" & gene != "NA"), by="gene")





revolutissima.gene_te %>% filter(gene %in% revolutissima_DE_orthogroup_genes) %>%
  set_colnames(c("annotation", "gene", "gene_dist", "te_length", "TE")) %>% 
  left_join(revolutissima_te_annotation, by="TE", relationship = "many-to-many") %>% 
  left_join(orthogroups_long %>% filter(Orthogroup %in% imp_rev_OGs & species == "D. revolutissima" & gene != "NA"), by="gene") %>%
  filter(gene == "g11314")



te_annotation


# next do a fishers test for this to put a p-value on it

ultramafic_DEGs_High<-ultramafic_DEGs[ultramafic_DEGs %in% c("g10152", "g1861", "g22083", "g594", "g8963")]

# get DEG homologs in impolita and revolutissima
imp_rev_OGs_High<-orthogroups_long %>% filter(gene %in% ultramafic_DEGs_High & species == "D. vieillardii") %>% pull(Orthogroup)
impolita_DE_orthogroup_genes_High<-impolita_DE_orthogroup_genes<-orthogroups_long %>% filter(Orthogroup %in% imp_rev_OGs & species == "D. impolita" & gene != "NA") %>% pull(gene)
revolutissima_DE_orthogroup_genes_High<-revolutissima_DE_orthogroup_genes<-orthogroups_long %>% filter(Orthogroup %in% imp_rev_OGs & species == "D. revolutissima" & gene != "NA") %>% pull(gene)

vie_g1861 impolita_g4621  rev_g11314
vie_8963  impolita_absent rev_g4181
vie_10152 impolita_22665  rev_482
vie_594   impolita_3611   rev_14832



impolita_DE_orthogroup_genes_High[impolita_DE_orthogroup_genes_High %in% c("g4621", "g22665", "g3611")]
revolutissima_DE_orthogroup_genes_High[revolutissima_DE_orthogroup_genes_High %in% c("g11314", "g4181", "g482", "g14832")]















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






