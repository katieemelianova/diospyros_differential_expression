library(dplyr)
library(magrittr)
library(readr)
library(stringr)
library(tidyverse)


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



strings_to_remove<-c("calciphila/", "Aligned.sortedByCoord.out.bam")
test<-read_in_featurecounts('calciphila_counts', strings_to_remove)

test_samples <- test$counts %>% colnames() %>% tibble() %>% set_colnames("sample_id")

annotate_diospyros_samples(test_samples) 



data.frame(mycol="cal_19c2r") %>%
  separate(mycol, 
         into = c("text", "num"), 
         sep = "(?<=[A-Za-z])(?=[0-9])")

data.frame(mycol="cal_19c2r") %>%
  separate(mycol, 
           into = c("species", "data"), 
           sep = "_")

#species mother plant_id pseudorep techrep

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



