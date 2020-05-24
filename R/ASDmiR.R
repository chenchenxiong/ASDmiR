# Raw ASD datasets have matched ASD-related expression profiles of 172 miRNAs, 595 lncRNAs, 18114 mRNA, samples have been categorized as ASD (104 samples) and normal (82 samples).
raw_Autism <- load("ASD.RData") 

# To discover differentially expressed miRNAs, lncRNAs and mRNAs between ASD samples and normal samples, differential expression analysis can be used.
## we select top 100 miRNAs, top 300 lncRNAs and top 4000 mRNAs that are differentially expressed
library(miRLAB)
library(limma)
DiffExpAna_miR_lncR<-DiffExpAnalysis("t_miR_Exp_Autism.csv", "t_miR_Exp_Normal.csv", "t_lncR_Exp_Autism.csv", "t_lncR_Exp_Normal.csv", topkmiR = 100, topkmR = 300, p.miR = 1, p.mR = 1)
DiffExpAna_miR_mR<-DiffExpAnalysis("t_miR_Exp_Autism.csv", "t_miR_Exp_Normal.csv", "t_mR_Exp_Autism.csv", "t_mR_Exp_Normal.csv", topkmiR = 100, topkmR = 4000, p.miR = 1, p.mR = 1) ##186*4100


# Identification of miRNA-associated regulatory networks: miRNA-target regulatory network and miRNA sponge network.
## Predict miRNA targets using 12 computational methods to predict miRNA targets,
## comparison the performance of methods by using experimentally validated miRNA-target interactions in selecting the top 50, 100, 150, 200 miRNA-target interactions of each miRNA 
## visualization of results shows that best performing method among 12 methods,
## Identification of miRNA-target interactions by the best performing method.
library(miRLAB)
cause = 1:100 #column of 1:35 are miRNAs
effect_lncR = 101:400 #column of 101:400 are lncRNAs
effect_mR = 101:4100 #column of 101:4100 are mRNAs

##predict miRNA targets using Pearson correlation
DiffExp_ps_miR_lncR <- Pearson("DiffExp_miR_lncR.csv", cause, effect_lncR)
DiffExp_ps_miR_mR <- Pearson("DiffExp_miR_mR.csv", cause, effect_mR)

##predict miRNA targets using Spearman correlation
DiffExp_spearman_miR_lncR <- Spearman("DiffExp_miR_lncR.csv", cause, effect_lncR)
DiffExp_spearman_miR_mR <- Spearman("DiffExp_miR_mR.csv", cause, effect_mR)

##predict miRNA targets using Kendall correlation
DiffExp_kendall_miR_lncR<- Kendall("DiffExp_miR_lncR.csv", cause, effect_lncR)
DiffExp_kendall_miR_mR <- Kendall("DiffExp_miR_mR.csv", cause, effect_mR)

##predict miRNA targets using  Distance correlation
DiffExp_dcov_miR_lncR <- Dcov("DiffExp_miR_lncR.csv", cause, effect_lncR)
DiffExp_dcov_miR_mR <- Dcov("DiffExp_miR_mR.csv", cause, effect_mR)

##predict miRNA targets using Hoeffding's D measure
DiffExp_hoeffding_miR_lncR <- Hoeffding("DiffExp_miR_lncR.csv", cause, effect_lncR)
DiffExp_hoeffding_miR_mR <- Hoeffding("DiffExp_miR_mR.csv", cause, effect_mR)

##predict miRNA targets using Mutual Information
DiffExp_mi_miR_lncR <- MI("DiffExp_miR_lncR.csv", cause, effect_lncR)
DiffExp_mi_miR_mR <- MI("DiffExp_miR_mR.csv", cause, effect_mR)

##predict miRNA targets using IDA 
DiffExp_ida_miR_lncR <- IDA("DiffExp_miR_lncR.csv", cause, effect_lncR, pcmethod = "original", targetbinding = NA)
DiffExp_ida_miR_mR <- IDA("DiffExp_miR_mR.csv", cause, effect_mR, pcmethod = "original", targetbinding = NA)

##predict miRNA targets using RDC 
DiffExp_rdc_miR_lncR <- RDC("DiffExp_miR_lncR.csv", cause, effect_lncR)
DiffExp_rdc_miR_mR <- RDC("DiffExp_miR_mR.csv", cause, effect_mR)

##predict miRNA targets using Lasso
DiffExp_lasso_miR_lncR <- Lasso("DiffExp_miR_lncR.csv", cause, effect_lncR)
DiffExp_lasso_miR_mR <- Lasso("DiffExp_miR_mR.csv", cause, effect_mR)

##predict miRNA targets using Elastic
DiffExp_elastic_miR_lncR <- Elastic("DiffExp_miR_lncR.csv", cause, effect_lncR)
DiffExp_elastic_miR_mR <- Elastic("DiffExp_miR_mR.csv", cause, effect_mR)

##predict miRNA targets using Z-score
DiffExp_zscore_miR_lncR <- Zscore("DiffExp_miR_lncR.csv", cause, effect_lncR)
DiffExp_zscore_miR_mR <- Zscore("DiffExp_miR_mR.csv", cause, effect_mR)

##predict miRNA targets using ProMISe
DiffExp_promise_miR_lncR <- ProMISe("DiffExp_miR_lncR.csv", cause, effect_lncR)
DiffExp_promise_miR_mR <- ProMISe("DiffExp_miR_mR.csv", cause, effect_mR)

### Comparison study on the  number of confirmed miRNA-target interactions using 12 methods 
miR_lncR_12results<- list(DiffExp_ps_miR_lncR, DiffExp_spearman_miR_lncR, DiffExp_kendall_miR_lncR, DiffExp_dcov_miR_lncR, DiffExp_hoeffding_miR_lncR, DiffExp_mi_miR_lncR, 
                          DiffExp_ida_miR_lncR, DiffExp_rdc_miR_lncR, DiffExp_lasso_miR_lncR, DiffExp_elastic_miR_lncR, DiffExp_zscore_miR_lncR, DiffExp_promise_miR_lncR)
miR_mR_12results<- list(DiffExp_ps_miR_mR, DiffExp_spearman_miR_mR, DiffExp_kendall_miR_mR, DiffExp_dcov_miR_mR, DiffExp_hoeffding_miR_mR, DiffExp_mi_miR_mR, 
                        DiffExp_ida_miR_mR, DiffExp_rdc_miR_mR, DiffExp_lasso_miR_mR, DiffExp_elastic_miR_mR, DiffExp_zscore_miR_mR, DiffExp_promise_miR_mR)

source("zzz.R")
library(miRLAB)

## Validate the results of the top50 targets of each miRNA predicted by the 12 methods
miR_lncR_Experiment_top50 <- experiment_groundtruth(allmethods = miR_lncR_12results, topk = 50, Expgroundtruth = "miRNA_lncRNA_groundtruth_lncBase_v2.0+NPInter_v4.0.csv", LFC=1, downreg = TRUE)
miR_mR_Experiment_top50 <- experiment_groundtruth(allmethods = miR_mR_12results, topk = 50, Expgroundtruth = "miRNA_mRNA_groundtruth_miRTarBase_v8.0+TarBase_v8.0.csv", LFC=1, downreg = TRUE)

## Validate the results of the top100 targets of each miRNA predicted by the 12 methods
miR_lncR_Experiment_top100 <- experiment_groundtruth(allmethods = miR_lncR_12results, topk = 100, Expgroundtruth = "miRNA_lncRNA_groundtruth_lncBase_v2.0+NPInter_v4.0.csv", LFC=1, downreg = TRUE)
miR_mR_Experiment_top100 <- experiment_groundtruth(allmethods = miR_mR_12results, topk = 100, Expgroundtruth = "miRNA_mRNA_groundtruth_miRTarBase_v8.0+TarBase_v8.0.csv", LFC=1, downreg = TRUE)

## Validate the results of the top150 targets of each miRNA predicted by the 12 methods
miR_lncR_Experiment_top150 <- experiment_groundtruth(allmethods = miR_lncR_12results, topk = 150, Expgroundtruth = "miRNA_lncRNA_groundtruth_lncBase_v2.0+NPInter_v4.0.csv", LFC=1, downreg = TRUE)
miR_mR_Experiment_top150 <- experiment_groundtruth(allmethods = miR_mR_12results, topk = 150, Expgroundtruth = "miRNA_mRNA_groundtruth_miRTarBase_v8.0+TarBase_v8.0.csv", LFC=1, downreg = TRUE)

## Validate the results of the top200 targets of each miRNA predicted by the 12 methods
miR_lncR_Experiment_top200 <- experiment_groundtruth(allmethods = miR_lncR_12results, topk = 200, Expgroundtruth = "miRNA_lncRNA_groundtruth_lncBase_v2.0+NPInter_v4.0.csv", LFC=1, downreg = TRUE)
miR_mR_Experiment_top200 <- experiment_groundtruth(allmethods = miR_mR_12results, topk = 200, Expgroundtruth = "miRNA_mRNA_groundtruth_miRTarBase_v8.0+TarBase_v8.0.csv", LFC=1, downreg = TRUE)

# filter and Compare miRNA-lncRNA interactions
compare_miR_lncR_top50 <- filterAndCompare_groundtruth(miR_lncR_Experiment_top50, 1)
compare_miR_lncR_top100 <- filterAndCompare_groundtruth(miR_lncR_Experiment_top100, 1)
compare_miR_lncR_top150 <- filterAndCompare_groundtruth(miR_lncR_Experiment_top150, 1)
compare_miR_lncR_top200 <- filterAndCompare_groundtruth(miR_lncR_Experiment_top200, 1)
# filter and Compare miRNA-lncRNA interactions 
compare_miR_mR_top50 <- filterAndCompare_groundtruth(miR_mR_Experiment_top50, 1)
compare_miR_mR_top100 <- filterAndCompare_groundtruth(miR_mR_Experiment_top100, 1)
compare_miR_mR_top150 <- filterAndCompare_groundtruth(miR_mR_Experiment_top150, 1)
compare_miR_mR_top200 <- filterAndCompare_groundtruth(miR_mR_Experiment_top200, 1)


# Visualization of the results of comparison study on validated miRNA-target interactions
validate_methods<-c("Pearson","Spearman","Kendall", "Dcov", "Hoeffding","MI","IDA", "RDC", "Lasso","Elastic","Zscore", "ProMISe")
x_ <- c("Top 50", "Top 100", "Top 150", "Top 200")
## Integrate  comparing results of validated miRNA-lncRNA interactions in the case of top50, 100, 150, 200
miR_lncR_addtop50 <- data.frame(validate_num = compare_miR_lncR_top50[[2]], validate_methods, x_axis= x_[1], row.names = 1:12)
miR_lncR_addtop100 <- data.frame(validate_num =compare_miR_lncR_top100[[2]], validate_methods,x_axis= x_[2], row.names = 1:12)
miR_lncR_addtop150 <- data.frame(validate_num = compare_miR_lncR_top150[[2]], validate_methods, x_axis = x_[3], row.names = 1:12)
miR_lncR_addtop200 <- data.frame(validate_num = compare_miR_lncR_top200[[2]], validate_methods, x_axis = x_[4], row.names = 1:12)
miR_lncR <- rbind(miR_lncR_addtop50, miR_lncR_addtop100, miR_lncR_addtop150, miR_lncR_addtop200)

library(ggplot2)
col <- c("#FDDC9B","#E6E025","#E95658","#70C2BC","#5DB230","#9DC514","#99BC88","#F4A618","#BD77B0","#EA5997","#A8A6D3","#CA9855")
# Plot comparison results in terms of validated miRNA-lncRNA interactions
p1_lncR <- ggplot(miR_lncR, aes(x=x_axis, y=log2(validate_num)*20, fill = factor(validate_methods)))+geom_bar(stat = "identity", position ="stack")+geom_text(aes(label=validate_num),
                  stat = "identity",position = position_stack(vjust = 0.5),size = 10, colour = "black")+labs(x = NULL, y="Number of confirmed intearctions", 
                  title = "Comparison in terms of validated miRNA-lncRNA interactions")+theme(panel.background = element_blank(), axis.text.x = element_text(hjust = 0.5, vjust = 0.5,
                  face = "bold",size = 32), axis.text.y = element_text(hjust = 0.5, vjust = 0.5, face = "bold", size = 24), axis.title.y = element_text(hjust = 0.5, vjust = 0.5,
                  face = "bold", size = 32),axis.line = element_line(colour = "black", size = 0.6), plot.title = element_text(family = "myFont",hjust = 0.5,vjust = 12, size = 36, face = "bold"))
                  +coord_cartesian(expand = TRUE, clip = "off")+scale_y_continuous(expand = c(0,0))+guides(fill = FALSE)+scale_fill_manual(values = col[1:12])

## Integrate comparing results of validated miRNA-mRNA interactions in the case of top50, 100, 150, 200
miR_mR_addtop50 <- data.frame(validate_num = compare_miR_mR_top50[[2]], validate_methods, x_axis= x_[1], row.names = 1:12)
miR_mR_addtop100 <- data.frame(validate_num =compare_miR_mR_top100[[2]], validate_methods,x_axis= x_[2], row.names = 1:12)
miR_mR_addtop150 <- data.frame(validate_num = compare_miR_mR_top150[[2]], validate_methods, x_axis = x_[3], row.names = 1:12)
miR_mR_addtop200 <- data.frame(validate_num = compare_miR_mR_top200[[2]], validate_methods, x_axis = x_[4], row.names = 1:12)
miR_mR <- rbind(miR_mR_addtop50, miR_mR_addtop100, miR_mR_addtop150, miR_mR_addtop200)

library(ggplot2)
library(cowplot)
# Plot comparison results in terms of validated miRNA-mRNA interactions
windowsFonts(myFont = windowsFont("Times New Roman"))
p1_mR <- ggplot(miR_mR, aes(x=x_axis, y=log2(validate_num)*60, fill = validate_methods)) +geom_bar(stat = "identity", position ="stack")+labs(x = NULL, y="Number of confirmed interactions",
                title = "Comparison in terms of validated miRNA-mRNA interactions", face = "bold")+geom_text(aes(label = validate_num), vjust = 0.5, colour = "black",
                position = position_stack(vjust = 0.5), size = 10, face = "bold")+theme(panel.background = element_blank(),axis.text.x = element_text(hjust = 0.5, vjust = 0.5,face = "bold",
                size = 32), axis.text.y = element_text(hjust = 0.5, vjust = 0.5, face = "bold", size = 24),axis.title.y = element_text(hjust = 0.5, vjust = 0.5, face = "bold", size = 32),
                axis.line = element_line(colour = "black", size = 0.6), plot.title = element_text(family = "myFont",hjust = -0.5,vjust = 3, size = 36, face = "bold"),
                legend.title = element_text(size = 28),legend.text = element_text(colour = "black", size = 28))+coord_cartesian(expand = TRUE, clip = "off")+scale_y_continuous(expand = c(0,0))
                +scale_fill_manual(values = col[1:12])+guides(fill= guide_legend(title = "Methods",reverse = FALSE,size = 28, face = "bold"))

library(cowplot)
comparison_12_methods_figure <- plot_grid(p1_lncR, p1_mR, labels = c("A", "B"), label_size = 36,ncol = 2, hjust = -0.5, vjust = 1.25)

dev.new()
tiff(file="Comparison of number of confirmed miRNA target interactiions using 12 miRNA-target prediction methods.tif", width = 2400, height = 1200) # save as .tif format
library(ggplot2)
library(grid)
grid.newpage()
pushViewport(viewport(layout = grid.layout(1,1)))
vplayout <- function(x,y){viewport(layout.pos.row = x, layout.pos.col = y)}
print(comparison_12_methods_figure, vp=vplayout(1,1))
dev.off()
## visualization of results shows that ProMISe method performs the best among 12 methods.

### Identification of miRNA-target interactions by ProMISe
# Select top 200 targets(lncRNAs and miRNAs) of each miRNA
promise_100miR_200lncR = matrix(NA,nrow = 20000, ncol = 3)
for (m in 1:100) {
  temp <- as.matrix(bRank(DiffExp_promise_miR_lncR, m, 200, downreg = TRUE))
  promise_100miR_200lncR[(200*m-199):(200*m),1:3] =  temp[1:200,1:3]
}
colnames(promise_100miR_200lncR) = c("miRNA","gene","Correlation")

promise_100miR_200mR = matrix(NA,nrow = 20000, ncol = 3)
for (m in 1:100) {
  temp <- as.matrix(bRank(DiffExp_promise_miR_mR, m, 200, downreg = TRUE))
  promise_100miR_200mR[(200*m-199):(200*m),1:3] =  temp[1:200,1:3]
}
colnames(promise_100miR_200mR) = c("miRNA","gene","Correlation")
# Validation of miRNA-target interactions predicted by ProMISe method
library(miRLAB)
promise_miR_lncRtop200_Validated <- Validation(promise_100miR_200lncR, datacsv = "miRNA_lncRNA_groundtruth_lncBase_v2.0+NPInter_v4.0.csv")
promise_miR_mRtop200_Validated <- Validation(promise_100miR_200mR, datacsv = "miRNA_mRNA_groundtruth_miRTarBase_v8.0+TarBase_v8.0.csv")
## 241 validated miRNA-lncRNA interactions and 1438 validated miRNA-mRNA interactions are obtained by ProMISe method

promise_validated_miR_lncR_mR_1679<- rbind(promise_miR_lncRtop200_Validated[[1]],promise_miR_mRtop200_Validated[[1]])
##Integrate all validated miRNA--lncRNA interactions and miRNA-mRNA interactions as miRNA-target interactionsIn total, 1679 miRNA-target interactions have been identitied by ProMISe method


# Enrichment analysis of identified miRNA-target network
## Creat directed miRNA-target graph from miRNA-target interactions, and calculate node degree
library(igraph)
miR_target_graph <- make_graph(c(t(promise_validated_miR_lncR_mR_1679[,1:2])), directed = TRUE)
miR_target_graph_degree <- degree(graph_from_data_frame(as_data_frame(miR_target_graph),directed = TRUE))
miR_target_graph_outdegree <- degree(miR_target_graph, mode = "out")
# Due to the importance of hub miRNAs in miRNA target network, top 20% miRNAs with the largest degree are selected as hub miRNAs
hub_miRNA_target <- names(sort(miR_target_graph_outdegree[which(miR_target_graph_outdegree!=0)], decreasing = TRUE))[1:ceiling(0.2*length(which(miR_target_graph_outdegree!=0)))] 
## we have obtained 12 hub miRNAs(hsa-miR-195-5p, hsa-miR-15a-5p, hsa-miR-26b-5p, hsa-miR-23a-3p, hsa-miR-93-5p, hsa-miR-210-3p, hsa-miR-25-3p, 
## hsa-miR-30b-5p, hsa-miR-148b-3p, hsa-miR-149-5p, hsa-miR-200c-3p, hsa-miR-147a)  
 
## we integrate all targets of hub miRNAs, and conduct GO and KEGG enrichment analysis.
library(clusterProfiler)
# Integrate targets of hub miRNAs for enrichment analysis
miR_target_valid1679_DF <- as.data.frame(promise_validated_miR_lncR_mR_1679)
hub_target_inter <- data.frame()
for (n in 1:length(hub_target)) {
  temp_tar <- miR_target_valid1679_DF[which(hub_miRNA_target[n]==miR_target_valid1679_DF[1]),]
  hub_target_inter <- rbind(hub_target_inter, temp_tar)
}
entrezID_target <- bitr(hub_target_inter[,2], fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db") #Translate gene symbol into gene ID
hub_target_GO <- enrichGO(as.character(entrezID_target$ENTREZID), keyType = "ENTREZID", OrgDb = "org.Hs.eg.db", ont = "ALL", pvalueCutoff = 0.05, pAdjustMethod = "BH",qvalueCutoff = 0.1, readable = TRUE) 
hub_target_KEGG <- enrichKEGG(entrezID_target$ENTREZID, organism = "hsa",  pvalueCutoff = 0.05, pAdjustMethod = "BH",qvalueCutoff = 0.1,use_internal_data = TRUE, minGSSize = 10)
## As a result, 123 GO terms and 7 KEGG pathways are enriched for miRNA-target network


# Identification of miRNA sponge interaction network by Sensitivity Partial Pearson Correlation (SPPC)
## we select different SC cutoffs from 0.1 to 0.3 with a step of 0.05, to infer miRNA sponge interaction network. 
## Under different SC cutoffs, we use R square value to evaluate the goodness of power law degree distribution for the identified miRNA sponge interaction network.
## R Square of each miRNA sponge network can be obtained by using the Network Analyzer plugin in Cytoscape.

DiffExp_miR <- read.csv("DiffExp_miR.csv", header = FALSE, sep = ",")
DiffExp_lncR <- read.csv("DiffExp_lncR.csv", header = FALSE, sep = ",")
DiffExp_mR <- read.csv("DiffExp_mR.csv", header = FALSE, sep = ",")
DEA_miR_mR_lncR <- cbind(DiffExp_miR,DiffExp_lncR, DiffExp_mR)## 186*4400

library(miRspongeR)
# SC = 0.1
miR_sponge_sppc_0.1 <- spongeMethod(promise_validated_miR_lncR_mR_1679[,1:2], DEA_miR_mR_lncR, padjustvaluecutoff = 0.05, senscorcutoff = 0.1, method = "sppc") #R^2=0.729
# SC = 0.15
miR_sponge_sppc_0.15 <- spongeMethod(promise_validated_miR_lncR_mR_1679[,1:2], DEA_miR_mR_lncR, padjustvaluecutoff = 0.05, senscorcutoff = 0.15, method = "sppc") #R^2=0.797
# SC = 0.2
miR_sponge_sppc_0.2 <- spongeMethod(promise_validated_miR_lncR_mR_1679[,1:2], DEA_miR_mR_lncR, padjustvaluecutoff = 0.05, senscorcutoff = 0.2, method = "sppc") #R^2=0.812
# SC = 0.25
miR_sponge_sppc_0.25 <- spongeMethod(promise_validated_miR_lncR_mR_1679[,1:2], DEA_miR_mR_lncR, padjustvaluecutoff = 0.05, senscorcutoff = 0.25, method = "sppc") #R^2=0.815
# SC = 0.3
miR_sponge_sppc_0.3 <- spongeMethod(promise_validated_miR_lncR_mR_1679[,1:2], DEA_miR_mR_lncR, padjustvaluecutoff = 0.05, senscorcutoff = 0.3, method = "sppc") #R^2=0.770
## According to the principle of the largest R square value, we select the SC cutoff as 0.25 to infer miRNA sponge interaction network. And we know all miRNA sponge network are all scare-free.

## Due to the importance of hub genes in miRNA sponge network, top 20% target genes with the largest degree are selected as hub miRNA sponges
miR_sponge_graph <- make_graph(c(t(miR_sponge_sppc_0.25[,1:2])), directed = FALSE)
miR_sponge_graph_degree <- degree(graph_from_data_frame(as_data_frame(miR_sponge_graph), directed = FALSE))
hub_miR_sponge <- names(sort(miR_sponge_graph_degree[which(miR_sponge_graph_degree!=0)], decreasing = TRUE))[1:ceiling(0.2*length(which(miR_sponge_graph_degree!=0)))]
## We obtain 15 hub miRNA sponges (SLC38A2, SHOC2, DDX6, WSB1, PURB, DDX5, DLEU2, USP15, C6orf62, ADAM10, STK4, LBR, PNISR, ANKRD44, SERINC1)

## We integrate all hub miRNA sponges, and conduct GO and KEGG enrichment analysis.
library(clusterProfiler)
miR_sponge_sppc_0.25_DF <- as.data.frame(miR_sponge_sppc_0.25)
hub_sponge_inter <- data.frame()
for (i in 1:length(hub_gene)) {
  temp_sponge <- miR_sponge_sppc_0.25_DF[which(hub_miR_sponge[i] == miR_sponge_sppc_0.25_DF[1]),]
  hub_sponge_inter <- rbind(hub_sponge_inter,temp_sponge)
}
all_hubsponge <- union(hub_sponge_inter[,1],hub_sponge_inter[,2])
entrezID_sponge <- bitr(all_hubsponge, fromType = "SYMBOL", toType = "ENTREZID",OrgDb = "org.Hs.eg.db")
hub_sponge_GO <- enrichGO(entrezID_sponge$ENTREZID, keyType = "ENTREZID", OrgDb = "org.Hs.eg.db", ont = "ALL", pvalueCutoff = 0.05, pAdjustMethod = "BH",qvalueCutoff = 0.1, readable = TRUE)
hub_sponge_KEGG <- enrichKEGG(entrezID_sponge$ENTREZID, organism = "hsa", pvalueCutoff = 0.05, pAdjustMethod = "BH",qvalueCutoff = 0.1,minGSSize = 10)
## In total, 86 GO terms and 1 KEGG pathways are enriched in identified miRNA sponge network.


# Enrichment analysis of hub genes
# We use miEAA online tool to conduct enrichment analysis of hub miRNAs. 
# In the following, we conduct Enrichment analysis hub miRNA sponges.
# Disease enrichment analysis of miRNA target module
library(miRspongeR)
hub_miR_sponge_miRspongeR_DEA <- moduleDEA(list(hub_miR_sponge), ont = "DO", OrgDb = "org.Hs.eg.db", padjustvaluecutoff = 0.05, padjustedmethod = "BH")
# Functional enrichment analysis of miRNA target module
hub_miR_sponge_miRspongeR_FEA <- moduleFEA(list(hub_miR_sponge), ont = "ALL", padjustvaluecutoff = 0.05, padjustedmethod = "BH")
##As a result, 2 KEGG pathways are enriched for hub miRNA sponges.


# Identification of miRNA-assocaited modules
## Identification of miRNA target modules
library(biclique)
library(Rcpp)
# Obtain all bicliques
bi.format(filename = "promise_validated_miR_lncR_mR_1679.el", filetype = 0) # In total, 57 left vertices, 231 right vertices, 1679 edges 
# Maximal biclique enumeration in bicpartite graphs
miR_target_biclique <- bi.clique(filename = "promise_validated_miR_lncR_mR_1679.el", left_least = 3, right_least = 3, filetype = 0)
## In total, the number of bicliquesis 9265. We focus on conducting enrichment analysis of top 20 largest miRNA-target regulatory modules.

#Extract top 20 largest miRNA-target regulatory modules
biclique <- c()
n=1
for (i in 1:length(miR_target_biclique)){
  biclique$ID[n] <- i
  biclique$left[n] <- c(length(miR_target_biclique[[i]]$left))
  biclique$right[n] <- c(length(miR_target_biclique[[i]]$right))
  biclique$bicliques[n] <- biclique$left[n]*biclique$right[n]
  i = i + 1
  n = n + 1
}
biclique<- as.data.frame(biclique_0119)
biclique_top20 <- biclique[order(-biclique$bicliques),][1:20,]

p = 1
biclique_top20_left <- list()
biclique_top20_right <- list()
top20_bicliques <- list()
for (j in 1:length(biclique_top20$ID)) {
  top20_bicliques[[p]] <- miR_target_biclique[[as.numeric(biclique_top20$ID[j])]]
  p = p + 1
}
## Obtain targets of top20 miRNA sponge modules
top20_bicliques_targets <- list()
for (q in 1:20) {
  top20_biclique_temp  <- list(top20_bicliques[[q]]$right)
  q = q + 1
  top20_bicliques_targets <- c(top20_bicliques_targets, top20_biclique_temp)
}
## Disease and functional enrichment analysis of miRNA target module
library(miRspongeR)
miR_tarModule_DEA <- moduleDEA(top20_bicliques_targets, OrgDb = "org.Hs.eg.db", ont = "DO", padjustvaluecutoff = 0.05)
miR_tarModule_FEA <- moduleFEA(top20_bicliques_targets, ont = "ALL", padjustvaluecutoff = 0.05, padjustedmethod = "BH")

# Identification of miRNA sponge modules by Markov Cluster Algorithm (MCL) and enrichment analysis.
library(miRspongeR)
miR_sppc_MCL <- netModule(miR_sponge_sppc_0.25[,1:2], method = "MCL",modulesize = 3, save = TRUE)
## In total, 10 miRNA sponge modules have been identified.

## Disease and functional enrichment analysis of miRNA target module
library(miRspongeR)
module_DEA <- moduleDEA(miR_sppc_MCL, ont = "DO", OrgDb = "org.Hs.eg.db", padjustvaluecutoff = 0.05, padjustedmethod = "BH")
module_FEA <- moduleFEA(miR_sppc_MCL,ont = "ALL", padjustvaluecutoff = 0.05, padjustedmethod = "BH")






