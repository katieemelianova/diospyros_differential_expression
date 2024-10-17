


# make a function which takes the original full data and gives back the counts and sample info in a named list based on selection statement
specify_comparison<-function(samples, counts, selection_statement){
  samples_return<-samples %>% filter(rlang::eval_tidy(rlang::parse_expr(selection_statement)))
  counts_return <- counts %>% dplyr::select(rownames(samples_return))
  return(list(samples=samples_return, counts=counts_return))
}


run_diffexp<-function(comparison, design_term, gene_lengths, cpm_threshold=1, min_samples=3){
  formula_parsed<-paste("~", design_term)
  dds <- DESeqDataSetFromMatrix(countData = comparison[["counts"]],
                                colData = comparison[["samples"]],
                                design = as.formula(formula_parsed))
  mcols(dds)$basepairs<-gene_lengths
  dds <- estimateSizeFactors(dds)
  print(sizeFactors(dds))
  idx <- rowSums( counts(dds, normalized=TRUE) >= cpm_threshold ) >= min_samples
  dds <- dds[idx,]
  dds <- DESeq(dds)
  res <- results(dds)
  return(list(dds=dds, results=res))
}


get_dds_object<-function(comparison, design_term, gene_lengths, cpm, min_samples, vst=FALSE){
  formula_parsed<-paste("~", design_term)
  dds <- DESeqDataSetFromMatrix(countData = comparison[["counts"]],
                                colData = comparison[["samples"]],
                                design = as.formula(formula_parsed))
  mcols(dds)$basepairs<-gene_lengths
  dds <- estimateSizeFactors(dds)
  idx <- rowSums( counts(dds, normalized=TRUE) >= cpm ) >= min_samples
  dds <- dds[idx,]
  if (vst == TRUE) {
    dds<-varianceStabilizingTransformation(dds)
    print("running vst")
  }
  return(dds)
}

get_enriched_terms<-function(gene_list, mappings, return_sample_GOData=FALSE){
  # use the gene 2 GOterms mapping provided for D. incarnata
  geneID2GO<-mappings
  # the input genes form the input, use these to annotate all genes, 1 is present in input list, 0 is absent
  geneSel<-gene_list
  geneSel<-factor(as.integer(names(geneID2GO) %in% geneSel))
  names(geneSel)<-names(geneID2GO)
  
  # set up the topGO object
  sampleGOdata <- new("topGOdata",
                      ontology = "BP",
                      allGenes = geneSel, 
                      nodeSize = 10,
                      annot = annFUN.gene2GO,
                      gene2GO = geneID2GO)
  
  # run three tests, fisher, Kol-Smirn, and Kol-Smirn with elimination
  resultFisher <- runTest(sampleGOdata, algorithm = "weight01", statistic = "fisher")
  
  # generate summary tane and return it
  allRes <- GenTable(sampleGOdata, classicFisher = resultFisher,
                     orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 100,
                     numChar=1000 )
  #allRes<-GenTable(sampleGOdata, Fis = resultFisher, topNodes = 20)
  
  if (return_sample_GOData == TRUE){
    return(list(result=allRes, goData=sampleGOdata))
  } else {
    return(allRes)
  }

get_significant_genes_nofc<-function(results_object, pvalue=0.05, mappings_format=FALSE, directional=FALSE){
    if (directional == TRUE) {
      de_genes_up<-results_object$results %>% data.frame() %>% filter(log2FoldChange > !!fold_change & padj < !!pvalue) %>% rownames()
      de_genes_down<-results_object$results %>% data.frame() %>% filter(log2FoldChange < -!!fold_change & padj < !!pvalue) %>% rownames()
      if (mappings_format  == TRUE){
        de_genes_up %<>% str_remove(":cds")
        de_genes_down %<>% str_remove(":cds")
        de_genes_up %<>% str_remove("-RA")
        de_genes_up %<>% str_remove("-t1")
        de_genes_down %<>% str_remove("-RA")
        de_genes_down %<>% str_remove("-t1")
      }
      to_return_genes <- list(up=de_genes_up, down=de_genes_down)
      return(to_return_genes)
    } else if (directional == FALSE) {
      de_genes<-results_object$results %>% data.frame() %>% filter(abs(log2FoldChange) > !!fold_change & padj < !!pvalue) %>% rownames()
      to_return_genes <- de_genes
      return(de_genes)
    }
  }
}


get_significant_genes<-function(results_object, fold_change=2, pvalue=0.05, directional=FALSE){
  if (directional == TRUE) {
    de_genes_up<-results_object$results %>% data.frame() %>% filter(log2FoldChange > !!fold_change & padj < !!pvalue) %>% rownames()
    de_genes_down<-results_object$results %>% data.frame() %>% filter(log2FoldChange < -!!fold_change & padj < !!pvalue) %>% rownames()
    to_return_genes <- list(up=de_genes_up, down=de_genes_down)
    return(to_return_genes)
  } else if (directional == FALSE) {
    de_genes<-results_object$results %>% data.frame() %>% filter(abs(log2FoldChange) > !!fold_change & padj < !!pvalue) %>% rownames()
    to_return_genes <- de_genes
    return(de_genes)
  }
}


get_species_tree <- function(){
  species_tree<-ape::read.tree("/Users/katieemelianova/Desktop/Diospyros/TeernaDiospyrosTrees/short_names_378_hetfil.contree")
  tip_to_drop<-c(species_tree$tip.label[(startsWith(species_tree$tip.label, "cal"))][-1], 
                 species_tree$tip.label[(startsWith(species_tree$tip.label, "unk"))],
                 species_tree$tip.label[(startsWith(species_tree$tip.label, "ine"))][-1],
                 species_tree$tip.label[(startsWith(species_tree$tip.label, "che"))][-1],
                 species_tree$tip.label[(startsWith(species_tree$tip.label, "vei"))][-1],
                 species_tree$tip.label[(startsWith(species_tree$tip.label, "fas"))],
                 species_tree$tip.label[(startsWith(species_tree$tip.label, "sam"))],
                 species_tree$tip.label[(startsWith(species_tree$tip.label, "mac"))],
                 species_tree$tip.label[(startsWith(species_tree$tip.label, "oub"))],
                 species_tree$tip.label[(startsWith(species_tree$tip.label, "Mau"))],
                 species_tree$tip.label[(startsWith(species_tree$tip.label, "fer"))],
                 species_tree$tip.label[(startsWith(species_tree$tip.label, "gfer"))],
                 species_tree$tip.label[(startsWith(species_tree$tip.label, "cffer"))],
                 species_tree$tip.label[(startsWith(species_tree$tip.label, "hil"))],
                 species_tree$tip.label[(startsWith(species_tree$tip.label, "san"))][-1],
                 species_tree$tip.label[(startsWith(species_tree$tip.label, "fol"))],
                 species_tree$tip.label[(startsWith(species_tree$tip.label, "sub"))],
                 species_tree$tip.label[(startsWith(species_tree$tip.label, "vie"))][-1],
                 species_tree$tip.label[(startsWith(species_tree$tip.label, "fla"))][-1],
                 species_tree$tip.label[(startsWith(species_tree$tip.label, "umb"))][-1],
                 species_tree$tip.label[(startsWith(species_tree$tip.label, "lab"))][-1],
                 species_tree$tip.label[(startsWith(species_tree$tip.label, "tris"))][-1],
                 species_tree$tip.label[(startsWith(species_tree$tip.label, "par"))][-1],
                 species_tree$tip.label[(startsWith(species_tree$tip.label, "ptpar"))],
                 species_tree$tip.label[(startsWith(species_tree$tip.label, "gla"))][-1],
                 species_tree$tip.label[(startsWith(species_tree$tip.label, "ruf"))][-1],
                 species_tree$tip.label[(startsWith(species_tree$tip.label, "pan"))][-1],
                 species_tree$tip.label[(startsWith(species_tree$tip.label, "per"))][-1],
                 species_tree$tip.label[(startsWith(species_tree$tip.label, "yah"))][-1],
                 species_tree$tip.label[(startsWith(species_tree$tip.label, "eru"))][-1],
                 species_tree$tip.label[(startsWith(species_tree$tip.label, "heq"))][-1],
                 species_tree$tip.label[(startsWith(species_tree$tip.label, "rev"))][-1],
                 species_tree$tip.label[(startsWith(species_tree$tip.label, "cfpus"))],
                 species_tree$tip.label[(startsWith(species_tree$tip.label, "pus"))][-1],
                 species_tree$tip.label[(startsWith(species_tree$tip.label, "imp"))][-1],
                 species_tree$tip.label[(startsWith(species_tree$tip.label, "ctpar"))],
                 species_tree$tip.label[(startsWith(species_tree$tip.label, "min"))][-1],
                 species_tree$tip.label[(startsWith(species_tree$tip.label, "trid"))][-1],
                 species_tree$tip.label[(startsWith(species_tree$tip.label, "afmin"))],
                 species_tree$tip.label[(startsWith(species_tree$tip.label, "spn"))][-1],
                 species_tree$tip.label[(startsWith(species_tree$tip.label, "Mad"))])
  
  
  ######################################
  #          drop defined tips         #
  ######################################
  
  species_tree<-drop.tip(species_tree, tip_to_drop)
  
  
  #####################################################
  #         make tip labels full species names        #
  #####################################################
  
  species_tree$tip.label<-case_when(startsWith(species_tree$tip.label, "cal") ~"calciphila",
                                    startsWith(species_tree$tip.label, "unk") ~"unknown",
                                    startsWith(species_tree$tip.label, "ine") ~"inexplorata",
                                    startsWith(species_tree$tip.label, "che") ~"cherleria",
                                    startsWith(species_tree$tip.label, "vei") ~"veillonii",
                                    startsWith(species_tree$tip.label, "fas") ~"fasciculata",
                                    startsWith(species_tree$tip.label, "sam") ~"samoensis",
                                    startsWith(species_tree$tip.label, "mac") ~"macrocarpa",
                                    startsWith(species_tree$tip.label, "oub") ~"oubatchensis",
                                    startsWith(species_tree$tip.label, "Mau") ~"mauritius",
                                    startsWith(species_tree$tip.label, "fer") ~"ferrea",
                                    startsWith(species_tree$tip.label, "gfer") ~"ferrea",
                                    startsWith(species_tree$tip.label, "cffer") ~"ferrea",
                                    startsWith(species_tree$tip.label, "hil") ~"hillebrandii",
                                    startsWith(species_tree$tip.label, "san") ~"sandwicensis",
                                    startsWith(species_tree$tip.label, "fol") ~"foliosa",
                                    startsWith(species_tree$tip.label, "sub") ~"subsessilis",
                                    startsWith(species_tree$tip.label, "vie") ~"vieillardii",
                                    startsWith(species_tree$tip.label, "fla") ~"flavocarpa",
                                    startsWith(species_tree$tip.label, "umb") ~"umbrosa",
                                    startsWith(species_tree$tip.label, "lab") ~"labillardierei",
                                    startsWith(species_tree$tip.label, "tris") ~"trisulca",
                                    startsWith(species_tree$tip.label, "par") ~"parviflora",
                                    startsWith(species_tree$tip.label, "ptpar") ~"parviflora",
                                    startsWith(species_tree$tip.label, "gla") ~"glans",
                                    startsWith(species_tree$tip.label, "ruf") ~"rufutomentosa",
                                    startsWith(species_tree$tip.label, "pan") ~"pancheri",
                                    startsWith(species_tree$tip.label, "per") ~"perplexa",
                                    startsWith(species_tree$tip.label, "yah") ~"yahouensis",
                                    startsWith(species_tree$tip.label, "eru") ~"erudita",
                                    startsWith(species_tree$tip.label, "heq") ~"hequetiae",
                                    startsWith(species_tree$tip.label, "rev") ~"revolutissima",
                                    startsWith(species_tree$tip.label, "cfpus") ~"pustulata",
                                    startsWith(species_tree$tip.label, "pus") ~"pustulata",
                                    startsWith(species_tree$tip.label, "imp") ~"impolita",
                                    startsWith(species_tree$tip.label, "ctpar") ~"parviflora",
                                    startsWith(species_tree$tip.label, "min") ~"minimifolia",
                                    startsWith(species_tree$tip.label, "trid") ~"tridentata",
                                    startsWith(species_tree$tip.label, "afmin") ~"minimifolia",
                                    startsWith(species_tree$tip.label, "spn") ~"sp. Pic N'ga")
  
  species_tree <- root(species_tree, which(species_tree$tip.label == "sandwicensis"))
  
  # make a list where item name is species name and objects within are the tip labels belonging to that species
  groupInfo<-split(species_tree$tip.label, species_tree$tip.label)
  
  # use groupOTU to group the tips by species and plot
  species_tree<-groupOTU(species_tree, groupInfo, group_name = "species")
  
  return(species_tree)
}

