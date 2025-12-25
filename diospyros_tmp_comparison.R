library(dplyr)
library(magrittr)
library(readr)
library(stringr)
library(tidyverse)
library(DESeq2)

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


test <- read_in_featurecounts('revolutissima_counts', strings_to_remove)


calculate_tpm <- function(counts, lengths){
  scaling_factor <- colSums(counts/lengths)/1000000
  tpm <- (counts/lengths)/scaling_factor
  return(tpm)
}

calculate_tpm(test$counts, test$lengths)




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
                                       species == "spn" ~ "ultramafic"))









