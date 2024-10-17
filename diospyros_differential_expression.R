library(dplyr)
library(magrittr)
library(readr)
library(stringr)
library(tidyverse)
library(DESeq2)
library(topGO)
library("wesanderson")
library(egg)
source("differential_expression_functions.R")

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


species_colours <- c("yellowgreen", "royalblue3", "plum1", "steelblue1", "salmon3", "orange2")



pic1<-ggplot(pca_all, aes(PC1, PC2)) +          
  geom_point(size = 1, stroke = 5, aes(fill = species, colour=species)) + 
  scale_color_manual(values=species_colours) 

species_tree<-get_species_tree()
outgroup<-c("sandwicensis")
ultramafic<-c("hequetiae", "revolutissima", "sp. Pic N'ga")
volcanic<-c("calciphila", "impolita", "labillardierei")

colours_tips <- case_when(species_tree$tip.label %in% ultramafic ~"Ultramafic",
                          species_tree$tip.label %in% volcanic ~"Volcanic",
                          !(species_tree$tip.label %in% c(outgroup, ultramafic, volcanic)) ~ "No data")

colours_labels <- case_when(species_tree$tip.label == "calciphila" ~ "yellowgreen",
                            species_tree$tip.label == "impolita" ~ "plum1",
                            species_tree$tip.label == "labillardierei" ~ "steelblue1",
                            species_tree$tip.label == "hequetiae" ~ "royalblue3",
                            species_tree$tip.label == "revolutissima" ~ "salmon3",
                            species_tree$tip.label == "sp. Pic N'ga" ~ "orange2",
                            !(species_tree$tip.label %in% c("calciphila", "impolita", "labillardierei", "hequetiae", "revolutissima", "sp. Pic N'ga")) ~ "grey77")


dd <- data.frame(taxa=species_tree$tip.label, tipcols=colours_tips, labelcols=colours_labels)
p<-ggtree(species_tree, size=1)
p <- p %<+% dd + geom_tippoint(aes(color=tipcols), size=6, show.legend=FALSE) + scale_color_manual(values=c("magenta", "chartreuse2"), limits = c("Ultramafic", "Volcanic"), na.value = "grey77") 
#p <- p + geom_tippoint(size=6, show.legend=FALSE, colour="blue")
p2<-p + new_scale_color() + geom_tiplab(size=8, aes(color=species), offset=0.01, show.legend=FALSE) + 
  scale_color_manual(values=colours_labels, limits=species_tree$tip.label) + 
  theme(legend.title = element_blank(),
        legend.text = element_text(size=12)) +
  expand_limits(x = 0.15)


pic2<-p2

pdf("testplot.pdf", width = 12, height = 7)
grid.arrange(pic1, pic2, nrow = 1)
dev.off()













ggplot(pca_all, aes(PC1, PC2)) +          
  geom_point(size = 1, stroke = 5, aes(fill = name, colour=name)) +
  scale_fill_manual(values=c(test)) +
  scale_colour_manual(values=c(test))



library(RColorBrewer)
library(ggh4x)
library(Polychrome)
P36 = createPalette(38,  c("red", "white", "green")) %>% as.character()



test<-c("#9DC022", "#35B600", "#0DFDA5", "#C3EEAC", "#66C3A7", "#ABF700",
"#89660D", "#AA8A8A", "#893D47", "#F8A193", "#A3004D", "#42494F", "#ED7149", "red", "tomato3",
"#FED5FB", "mistyrose1", "#FE35F5", "#FD7FBB", "#D000FE", "#FE16C5", "#831677",
"#00EAFC", "#BBE6FA", "#2A47FF", "#1C94FB", "#C0BBFE", "#2E85BB", "#5D51A5", "blue")

# get as many colours as you need, minus one to put grey at the end for Other
col_vector<-c(P36[1:length(levels(all$B2)) -1], "gray63")
##########################
#     Run DE analysis    #
##########################

cal_spn <-specify_comparison(all_samples, all_counts, "species %in% c('cal', 'spn')") %>% run_diffexp("species", all_lengths, cpm_threshold=1, min_samples=3)
heq_imp <-specify_comparison(all_samples, all_counts, "species %in% c('heq', 'imp')") %>% run_diffexp("species", all_lengths, cpm_threshold=1, min_samples=3)

imp_cal <-specify_comparison(all_samples, all_counts, "species %in% c('imp', 'cal')") %>% run_diffexp("species", all_lengths, cpm_threshold=1, min_samples=3)
imp_spn <-specify_comparison(all_samples, all_counts, "species %in% c('imp', 'spn')") %>% run_diffexp("species", all_lengths, cpm_threshold=1, min_samples=3)

heq_cal <-specify_comparison(all_samples, all_counts, "species %in% c('heq', 'cal')") %>% run_diffexp("species", all_lengths, cpm_threshold=1, min_samples=3)
heq_spn <-specify_comparison(all_samples, all_counts, "species %in% c('heq', 'spn')") %>% run_diffexp("species", all_lengths, cpm_threshold=1, min_samples=3)



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

te_gene_mine<-dna_met_vie




pheatmap::pheatmap(assay(heq_imp$dds)[konecna_de,], 
                   scale = "row", 
                   cluster_cols = FALSE,
                   cluster_rows = TRUE,
                   show_rownames = FALSE,
                   treeheight_row = 0, treeheight_col = 0)




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






