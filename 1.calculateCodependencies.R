

### load libraries for use
library(tidyverse)
library(dplyr)
library(tidyr)
library(org.Hs.eg.db)
library( annotables )
library(clusterProfiler)
library(ggplot2)
library(ggrepel)

### Get current directory
getwd()

### Set working directory
setwd('/Users/laura/Documents/PROJECTS/Broad_DepMap_20Q3 test/')


## load all data 
data <- read.delim("data/Achilles_gene_effect_CRISPR_AVANA_20Q3.csv", header=T, sep = ",")

## load DUB data
DUB_data <- read.delim("data/DUB_data_achilles_CRISPR_AVANA_20Q3.tsv", header=T, sep = "\t")

## Get DUB list
DUB_list <- as.data.frame(colnames(DUB_data)[2:95])

#DUB_list <- read.delim("data/DUB_list.txt", sep = "\t", header = F)


## convert numbers to numeric
DUB_data[, 2:96] <- sapply(DUB_data[, 2:96],as.numeric)

# set row names to be the cell line then transpose so each cell line is a column
row.names(DUB_data) <- DUB_data$CellLine
DUB_data_t <- as.data.frame(t(DUB_data[, 2:96]))


# count number of dependent cell lines (dependency score below -0.5)
DUB_data_t$depCount <- rowSums(DUB_data_t < -0.5, na.rm = T)
DUB_data_t$TestedCount <- rowSums(DUB_data_t < 10, na.rm = T)
DUB_data_t$fractionDep <- DUB_data_t$depCount / DUB_data_t$TestedCount

## make a new table with just three columns, then create a column for the DUB
DUB_data_fractionDep <- DUB_data_t[c("depCount", "TestedCount", "fractionDep")]
DUB_data_fractionDep$DUB <- row.names(DUB_data_fractionDep)

DUB_data_fractionDep <- subset(DUB_data_fractionDep, DUB_data_fractionDep$DUB %in% DUB_list$`colnames(DUB_data)[2:95]`)

## parse DUB name (get first word from string)
Names_temp <- strsplit(DUB_data_fractionDep$DUB, "[..]")
Names_edit = lapply(Names_temp, function(l) l[[1]])
DUB_data_fractionDep$DUB <- as.character(Names_edit)

#sort by the fraction of cell lines that are dependent
DUB_data_fractionDep <- DUB_data_fractionDep[order(DUB_data_fractionDep$fractionDep),]
DUB_data_fractionDep$rank <- c(1:length(DUB_data_fractionDep$DUB))

write.table(DUB_data_fractionDep, "CRISPR_depCellLineFraction/fraction_dependent_DepMap.txt", quote = F, row.names = F, sep = "\t")

#### plot the fraction of dependent cell lines for each DUB
DUB_data_fractionDep$geneLabels <- ""
DUB_data_fractionDep$geneLabels[54:94] <- DUB_data_fractionDep$DUB[54:94]

plot_title <- "DUB Essentiality in DepMap data"
  
ggplot(DUB_data_fractionDep) +
  geom_point(aes(x = rank, y = fractionDep)) +
  #scale_colour_manual(values = c( "#000000","#FF0000"))+
  #geom_point(alpha=0.7, size=1, aes(colour = threshold_drug)) +
  geom_text_repel(aes(x = rank, y = fractionDep, label = geneLabels), size =2.5,segment.size = 0.2, max.overlaps=Inf) +
  ggtitle(plot_title) +
  xlab("DUB") + 
  ylab("Fraction of Cell Lines Dependent") +
  #theme(legend.position = "none",
  #      plot.title = element_text(size = rel(1.5), hjust = 0.5),
  #      axis.title = element_text(size = rel(1.25))) +
  xlim(-5,110) +
  #ylim(0, 1.2) +
  scale_y_continuous(breaks=seq(0,1,0.1)) + # set axis ticks
  #geom_text_repel(aes(label = genelabels), size=2.5) +
  geom_hline(yintercept=0, size = 0.1) +
  geom_vline(xintercept=0, size = 0.1) +
  theme_bw()+
  theme(axis.text.x = element_text(size=10, angle = 0, hjust = 0.5 ),
        axis.text.y = element_text(size=10, angle = 0, hjust = 0.5 ),
        axis.title = element_text(size=10),
        panel.grid = element_blank(),
        plot.title = element_text(size=12,hjust = 0.5),
        axis.line = element_line(colour = "black"),
        aspect.ratio = 0.8,
        
        legend.position="none"
        #legend.text = element_text(size = 10),
        #legend.title = element_text(color = "#FFFFFF"),
        #legend.key.height = unit(0.5,"line")
  )

filesave = "CRISPR_depCellLineFraction/FractionDepCellLines_DUB_CRISPR_DepMap.png"
ggsave(filesave, height=4, width=4, dpi = 1200)




#### plot the fraction of dependent cell lines for each DUB as color coded bar plot

DUB_data_fractionDep$geneLabels[1:94] <- DUB_data_fractionDep$DUB[1:94]

DUB_data_fractionDep$class <- "Minimal Dependency"
DUB_data_fractionDep$class[DUB_data_fractionDep$fractionDep > 0.01] <- "Selectively Essential"
DUB_data_fractionDep$class[DUB_data_fractionDep$fractionDep > 0.9] <- "Commonly Essential"

# rearrange order
DUB_data_fractionDep$class <- factor(DUB_data_fractionDep$class, levels = c("Minimal Dependency", "Selectively Essential", "Commonly Essential"))

ggplot(DUB_data_fractionDep, aes(x=reorder(DUB, rank), y=fractionDep, fill = class)) + 
  geom_bar(stat="identity",position="dodge",) + 
  scale_fill_manual(values = c( "darkgrey","#f8c458" ,"#d17b76"))+
  #geom_point(alpha=0.7, size=1, aes(colour = threshold_drug)) +
  #geom_text_repel(aes(x = rank, y = fractionDep, label = geneLabels), size =2.5,segment.size = 0.2) +
  ggtitle(plot_title) +
  xlab("DUB") + 
  ylab("Fraction Cell Lines Dependent") +
  #theme(legend.position = "none",
  #      plot.title = element_text(size = rel(1.5), hjust = 0.5),
  #      axis.title = element_text(size = rel(1.25))) +
  #xlim(-5,100) +
  #ylim(0, 1.2) +
  scale_y_continuous(breaks=seq(0,1,0.1)) + # set axis ticks
  #geom_text_repel(aes(label = genelabels), size=2.5) +
  #geom_hline(yintercept=0, size = 0.1) +
  #geom_vline(xintercept=0, size = 0.1) +
  theme_bw()+
  theme(axis.text.x = element_text(size=7, angle = 90, vjust = 0.5 ),
        axis.text.y = element_text(size=7, angle = 0, hjust = 0.5 ),
        axis.title = element_text(size=8),
        panel.grid = element_blank(),
        plot.title = element_text(size=12,hjust = 0.5),
        axis.line = element_line(colour = "black"),
        #aspect.ratio = 0.3,
        
        legend.position="bottom",
        legend.text = element_text(size = 8),
        legend.title = element_text(color = "#FFFFFF"),
        #legend.key.height = unit(0.5,"line")
  )

filesave = "CRISPR_depCellLineFraction/FractionDepCellLines_DUB_CRISPR_DepMap_wide.png"
ggsave(filesave, height=3, width=8, dpi = 1200)



#### threshold cell lines dependent -- at least 3 cell lines dependent on the DUB

DUB_data_t <- DUB_data_t[order(DUB_data_t$depCount),]
  
#take at least 3 cell lines dependent
DUB_data_t_Dependent <- subset(DUB_data_t, DUB_data_t$depCount >= 3)

## convert numbers to numeric
data[, 2:18120] <- sapply(data[, 2:18120],as.numeric)

rownames(data) <- data[,1]

data <- data[,2:18120]


## calculate all codepedency correlations for all genes
cor <- as.data.frame(cor(data, use = "complete.obs"))

#cor_test <- cor.test(data[c("USP22..23326.")][,1], data[c("BRD9..65980.")][,1], use = 'pairwise.complete.obs')


## get matrix for DUB codependencies (pull out just the DUB codependencies)
cor_DUB <- cor[c(as.character(DUB_list$`colnames(DUB_data)`[1:96])),]


cor_DUB <- t(cor_DUB)

cor_DUB <- as.data.frame(cor_DUB)

#cor_DUB$
  
cor_DUB$gene <- rownames(cor_DUB)

#cor_DUB <- cor_DUB[1:94]

#turn wide table into long and skinny table
cor_DUB_l <- gather(cor_DUB, "gene2", "corr", ATXN3..4287.:ZUP1..221302.)

cor_DUB_l <- cor_DUB_l[c("gene","gene2","corr")]

## sort by p value
cor_DUB_l <- cor_DUB_l[order(cor_DUB_l$corr, decreasing = TRUE),]

## get rid of self correlations (e.g. USP1 : USP1 correlation = 1)
cor_DUB_l <- subset(cor_DUB_l, cor_DUB_l$corr < 1)


## pull out DUBs dependent on at least 3 cell lines
cor_DUB_l <- subset(cor_DUB_l, cor_DUB_l$gene2 %in% row.names(DUB_data_t_Dependent))


## parse gene2 name (get first word from string)
Names_temp <- strsplit(cor_DUB_l$gene2, "[..]")
Names_edit = lapply(Names_temp, function(l) l[[1]])

cor_DUB_l$gene2_HUGO <- as.character(Names_edit)



names(cor_DUB_l)[names(cor_DUB_l)=="gene"] <- "gene1"

## parse gene1 name (get first word from string)
Names_temp <- strsplit(cor_DUB_l$gene1, "[..]")
Names_edit = lapply(Names_temp, function(l) l[[1]])

cor_DUB_l$gene1_HUGO <- as.character(Names_edit)


### get codependencies by z score
#cor_DUB_sig <- subset(cor_DUB_l, z_score_sig == TRUE)


### count sig codependencies for each DUB
#codep_count <- cor_DUB_l %>% count(gene2_HUGO)



## save results
write.table(cor_DUB_l, "all_correlations/DUB_DepMap_correlations_dep3cellLines.tsv", sep = "\t", quote = F, row.names = F)  


# read in previously saved file
cor_DUB_l <- read.delim("all_correlations/DUB_DepMap_correlations_dep3cellLines.tsv", sep = "\t")


##run GO enrichment on top 7 codependencies for each DUB

MSigDB_GO_File <- "ClusterProfiler/geneSet_files_for_ClusterProfiler/c5.all.v7.0.entrez.gmt"
MSigDB_GO <- read.gmt(MSigDB_GO_File)


E2H <- annotables::grch38 %>% dplyr::select( HUGO = symbol, ENTREZ = entrez ) %>%
  filter( !duplicated( HUGO ) )


cumulative_sig_table_topResults <- data.frame()

compile_top_codependencies <- data.frame()

### one by one pull coDependency clusters from file to run GO enrichment
### for loop that will iterate through DUBs_in_data list of DUBs
for (i in unique(cor_DUB_l$gene2_HUGO)){
  
  ##for debugging or running one DUB in particular
  #i = "STAMBP"
  
  cluster <- subset(cor_DUB_l, gene2_HUGO == i)
  
  cluster$abs <- abs(cluster$corr)
  
  cluster <- cluster[order(cluster$abs, decreasing = T),]
  
  # take top 7
  cluster <- cluster[1:7,]
  
  compile_top_codependencies <- rbind(compile_top_codependencies, cluster)
  
  ## only take sig correlations
  #cluster <- subset(cluster, z_score > 3 | z_score < 3)
  
  ## only take negative correlations
  #cluster <- subset(cluster, corr < 0)
  
  ### set corr threshold
  #cluster <- cluster[c(abs((cluster$corr)) > 0.4 ),]
  
  cluster_genelist <- as.data.frame(cluster$gene1_HUGO)

  colnames(cluster_genelist)[colnames(cluster_genelist)=="cluster$gene1_HUGO"] <- "HUGO"
  cluster_genelist$HUGO <- as.character(cluster_genelist$HUGO)
  my_genes <- left_join( cluster_genelist, E2H )
  my_gene_names <- my_genes$ENTREZ
  my_gene_names <- my_genes$ENTREZ[!(is.na(my_genes$ENTREZ))]

  #sig_e_ENRICHER <- enricher(my_gene_names, TERM2GENE=GO_Kegg_React_Hall_combined, minGSSize = 0, pvalueCutoff=0.05, qvalueCutoff=0.05, pAdjustMethod = "bonferroni")

  if (length(my_gene_names) > 1){

    ######## take top network
    #my_genes$rank <- c(1:length(my_genes$HUGO))
    #my_genes <- my_genes[c(my_genes$rank <= 100),]

    sig_e_ENRICHER <- enricher(my_gene_names, TERM2GENE=MSigDB_GO, minGSSize = 0, pvalueCutoff=0.05, qvalueCutoff=0.05, pAdjustMethod = "bonferroni")


    #sig_e_Complex <- enrichGO(my_gene_names,  'org.Hs.eg.db', ont="CC", pvalueCutoff=0.05, minGSSize = 0, qvalueCutoff=0.05, readable = TRUE, pAdjustMethod = "bonferroni")
    #sig_e_GO <- enrichGO(my_gene_names,  'org.Hs.eg.db', ont="ALL", pvalueCutoff=0.05, minGSSize = 0,maxGSSize=1200, qvalueCutoff=0.05, readable = TRUE, pAdjustMethod = "bonferroni")

    #sig_e_All <- enrichGO(my_gene_names,  'org.Hs.eg.db', ont="ALL", pvalueCutoff=0.1, qvalueCutoff=0.25, readable = TRUE)

    ### make readable
    result_Readable <- setReadable(sig_e_ENRICHER, 'org.Hs.eg.db', 'ENTREZID')

    ### convert to dataframe
    sig_e <- data.frame(result_Readable)


    #e <- term_enrichment(my_genes, resources = resources, all_symbols = cached_coding_genes)

    results_filename = paste("ClusterProfiler/GO_Kegg_React_Hall_ClusterProfiler_R_top7/", i, "_Geneset_res_table.txt", sep="")

    write.table(sig_e, file=results_filename, sep="\t", quote=F, col.names=NA)

    #sig_e <- subset(e, e$q < 0.25)
    sig_e <- add_column(sig_e, DUB = i, .before = 1)
    #sig_e_cc <- add_column(sig_e_cc, DUB = i, .before = 1)
    #sig_e_go_all <- add_column(sig_e_go_all, DUB = i, .before = 1)

    ## have at least two genes
    sig_e <- filter(sig_e, sig_e$Count > 1)
    #sig_e_cc <- filter(sig_e_cc, sig_e_cc$Count > 1)
    #sig_e_go_all <- filter(sig_e_go_all, sig_e_go_all$Count > 1)

    cumulative_sig_table_topResults <- rbind(cumulative_sig_table_topResults, sig_e)

    #cumulative_sig_table_topResults_GO_CC <- rbind(cumulative_sig_table_topResults_GO_CC, sig_e_cc)
    #cumulative_sig_table_topResults_GO_ALL <- rbind(cumulative_sig_table_topResults_GO_ALL, sig_e_go_all)
    

  }
}

#### save cumulative table

write.table(cumulative_sig_table_topResults, file="ClusterProfiler/GO_Kegg_React_Hall_ClusterProfiler_R_top7/sigGO_DepMap_SigCodependencies.tsv", sep="\t", quote=F, col.names=NA)

write.table(compile_top_codependencies, "all_correlations/DUB_DepMap_top7_correlations.tsv", sep = "\t", quote = F, row.names = F)  



# read in previously saved file
cumulative_sig_table_topResults <- read.delim("ClusterProfiler/GO_Kegg_React_Hall_ClusterProfiler_R_top7/sigGO_DepMap_SigCodependencies.tsv", sep="\t")

compile_top_codependencies <- read.delim("all_correlations/DUB_DepMap_top7_correlations.tsv", sep = "\t")

## grab top term per DUB
cumulative_sig_table_topResults_sort <- cumulative_sig_table_topResults[order(cumulative_sig_table_topResults$qvalue),]

GO_compile_top <- cumulative_sig_table_topResults_sort[!duplicated(cumulative_sig_table_topResults_sort[c("DUB")]),]


#### cross reference codependencies with PPID

PPID_combined_file <- read.delim("CrossRef_Protein-Protein_Databases/PPID_combined_file/PPIDs_DUB_interactors_merged_long.tsv", sep = "\t")

PPID_combined_file$DUB_gene <- paste(PPID_combined_file$DUB, PPID_combined_file$interactors, sep = "_")

PPID_combined_file_w <- spread(PPID_combined_file, database, database)

cor_DUB_l$DUB_gene <- paste(cor_DUB_l$gene2_HUGO, cor_DUB_l$gene1_HUGO, sep = "_")
compile_top_codependencies$DUB_gene <- paste(compile_top_codependencies$gene2_HUGO, compile_top_codependencies$gene1_HUGO, sep = "_")

cor_DUB_l_PPID <- left_join(cor_DUB_l, PPID_combined_file_w, by = "DUB_gene")
cor_DUB_l_PPID$PPID_support[cor_DUB_l_PPID$Biogrid == "Biogrid" | cor_DUB_l_PPID$IntAct == "IntAct" | cor_DUB_l_PPID$PathwayCommons == "PathwayCommons" | cor_DUB_l_PPID$NURSA == "NURSA"] <- "yes"

cor_DUB_l_PPID$Biogrid <- gsub("Biogrid", "yes",cor_DUB_l_PPID$Biogrid)
cor_DUB_l_PPID$IntAct <- gsub("IntAct", "yes",cor_DUB_l_PPID$IntAct)
cor_DUB_l_PPID$PathwayCommons <- gsub("PathwayCommons", "yes",cor_DUB_l_PPID$PathwayCommons)
cor_DUB_l_PPID$NURSA <- gsub("NURSA", "yes",cor_DUB_l_PPID$NURSA)

cor_DUB_l_PPID$DUB <- cor_DUB_l_PPID$gene2_HUGO
cor_DUB_l_PPID$interactors <- cor_DUB_l_PPID$gene1_HUGO

top_DUB_codeps_PPID <- left_join(compile_top_codependencies, PPID_combined_file_w, by = "DUB_gene")
top_DUB_codeps_PPID$PPID_support[top_DUB_codeps_PPID$Biogrid == "Biogrid" | top_DUB_codeps_PPID$IntAct == "IntAct" | top_DUB_codeps_PPID$PathwayCommons == "PathwayCommons" | top_DUB_codeps_PPID$NURSA == "NURSA"] <- "yes"

top_DUB_codeps_PPID$Biogrid <- gsub("Biogrid", "yes",top_DUB_codeps_PPID$Biogrid)
top_DUB_codeps_PPID$IntAct <- gsub("IntAct", "yes",top_DUB_codeps_PPID$IntAct)
top_DUB_codeps_PPID$PathwayCommons <- gsub("PathwayCommons", "yes",top_DUB_codeps_PPID$PathwayCommons)
top_DUB_codeps_PPID$NURSA <- gsub("NURSA", "yes",top_DUB_codeps_PPID$NURSA)

top_DUB_codeps_PPID$DUB <- top_DUB_codeps_PPID$gene2_HUGO
top_DUB_codeps_PPID$interactors <- top_DUB_codeps_PPID$gene1_HUGO

write.table(cor_DUB_l_PPID, "CrossRef_Protein-Protein_Databases/DUB_DepMap_all_correlations_PPID.tsv", sep = "\t", quote = F, row.names = F)  
write.table(top_DUB_codeps_PPID, "CrossRef_Protein-Protein_Databases/DUB_DepMap_top7_correlations_PPID.tsv", sep = "\t", quote = F, row.names = F)  



top_DUB_codeps_PPID$DUB <- top_DUB_codeps_PPID$gene2_HUGO

top_DUB_codeps_PPID_GO <- left_join(top_DUB_codeps_PPID, GO_compile_top, by = "DUB")

DUB_data_fractionDep <- read.delim("CRISPR_depCellLineFraction/fraction_dependent_DepMap.txt", sep = "\t")

top_DUB_codeps_PPID_GO <- left_join(top_DUB_codeps_PPID_GO, DUB_data_fractionDep, by = "DUB")

write.table(top_DUB_codeps_PPID_GO, "tidy_combined_table/DUB_DepMap_top7_PPID_GO.tsv", sep = "\t", quote = F, row.names = F)  


top_DUB_codeps_PPID_GO_select <- subset(top_DUB_codeps_PPID_GO, top_DUB_codeps_PPID_GO$PPID_support == TRUE | top_DUB_codeps_PPID_GO$qvalue < 0.05 | top_DUB_codeps_PPID_GO$fractionDep > 0.2)

top_DUB_codeps_PPID_GO_drop <- subset(top_DUB_codeps_PPID_GO, !(top_DUB_codeps_PPID_GO$gene2_HUGO %in% top_DUB_codeps_PPID_GO_select$gene2_HUGO))


top_DUB_codeps_PPID_GO_select_table <- subset(top_DUB_codeps_PPID_GO, top_DUB_codeps_PPID_GO$gene2_HUGO %in% top_DUB_codeps_PPID_GO_select$gene2_HUGO)



write.table(top_DUB_codeps_PPID_GO_select_table, "tidy_combined_table/DUB_DepMap_top7_PPID_GO_forNetwork.tsv", sep = "\t", quote = F, row.names = F)  
write.table(top_DUB_codeps_PPID_GO_select_table, "network_figure/DUB_DepMap_top7_PPID_GO_forNetwork.tsv", sep = "\t", quote = F, row.names = F)  






#### Run GO enrichment on all DUB codependencies combined together ####
MSigDB_GO_File <- "ClusterProfiler/geneSet_files_for_ClusterProfiler/c5.all.v7.0.entrez.gmt"
MSigDB_GO <- read.gmt(MSigDB_GO_File)

E2H <- annotables::grch38 %>% dplyr::select( HUGO = symbol, ENTREZ = entrez ) %>%
  filter( !duplicated( HUGO ) )

## take top codependencies genes for all DUBs
cluster <- compile_top_codependencies

## only take positive correlations
#cluster <- subset(cluster, corr > 0)

### set corr threshold
#cluster <- cluster[c(abs((cluster$corr)) > 0.4 ),]

cluster_genelist <- as.data.frame(cluster$gene1_HUGO)

colnames(cluster_genelist)[colnames(cluster_genelist)=="cluster$gene1_HUGO"] <- "HUGO"
cluster_genelist$HUGO <- as.character(cluster_genelist$HUGO)
my_genes <- left_join( cluster_genelist, E2H )  
my_gene_names <- my_genes$ENTREZ

sig_e_ENRICHER <- enricher(my_gene_names, TERM2GENE=MSigDB_GO, minGSSize = 0, pvalueCutoff=0.05, qvalueCutoff=0.05, pAdjustMethod = "bonferroni")
sig_e_MF <- enrichGO(my_gene_names,  'org.Hs.eg.db', ont="MF", pvalueCutoff=0.05, minGSSize = 0, qvalueCutoff=0.05, readable = TRUE, pAdjustMethod = "bonferroni")

### make readable
result_Readable <- setReadable(sig_e_ENRICHER, 'org.Hs.eg.db', 'ENTREZID')

### convert to dataframe
sig_e <- data.frame(result_Readable)
sig_e_MF <- data.frame(sig_e_MF)

results_filename = paste("ClusterProfiler/GO_Kegg_React_Hall_ClusterProfiler_R_top7/depmap_top7_GO_allDUBsCombined_Geneset_res_table.txt", sep="")
write.table(sig_e, file=results_filename, sep="\t", quote=F, col.names=NA)

results_filename = paste("ClusterProfiler/GO_Kegg_React_Hall_ClusterProfiler_R_top7/depmap_top7_GO_allDUBsCombined_Geneset_res_table_GOMolecularFunction.txt", sep="")
write.table(sig_e_MF, file=results_filename, sep="\t", quote=F, col.names=NA)


sig_e_MF$order <- c(1:length(sig_e_MF$ID))
sig_e_MF$`-log10(q-value)` <- -log10(sig_e_MF$qvalue)

barplottitle = "GO Molecular Function Enrichment"

ggplot(sig_e_MF, aes(y=`-log10(q-value)`, x=reorder(Description, -order))) + 
  geom_bar(position="dodge", stat="identity", color = "black") +
  scale_fill_manual(values=c("red", "#ffbf00")) +
  #geom_text(aes(label=Description), hjust=1.1, color="black", size=2)+ 
  ggtitle(barplottitle) +
  xlab("")+
  theme_bw()+
  coord_flip() +
  theme(
    axis.text.x = element_text(size=8),
    axis.text.y =  element_text(size=8),
    axis.title = element_text(size=8),
    #panel.grid = element_blank(),
    plot.title = element_text(size=8,hjust = 0.5),
    axis.line = element_line(colour = "black"),
    #aspect.ratio = 0.5,
    legend.position="bottom"
    #legend.text = element_text(size = 10),
    #legend.title = element_text(color = "#FFFFFF"),
    #legend.key.height = unit(0.5,"line")
  )

filesave = "DUB_E3_network/DepMap_topMF.pdf"
ggsave(filesave, height=3, width=5, dpi = 1200)


