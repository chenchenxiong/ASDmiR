# ASDmiR
ASDmiR: a step-wise method to uncover miRNA regulation related to autism spectrum disorder
# Introduction
MiRNAs (miRNAs) are involved in nervous system developmental, and have potential to cause ASD. However, the miRNA  regulation mechanism in ASD is largely unclear. In this work, we present a novel framework, ASDmiR, to identify miRNA-target networks and modules, miRNA sponge networks and modules for uncovering the pathogenesis of ASD, as well as conduct enrichment analysis.
# Description of each file
  DiffExp_lncR.csv: Differentially expressed lncRNAs.<br />
  DiffExp_miR.csv: Differentially expressed miRNAs.<br />
  DiffExp_mR.csv: Differentially expressed mRNAs.<br />
  DiffExp_miR_lncR.csv: Differentially expressed miRNAs and lncRNAs.<br />
  DiffExp_miR_mR.csv: Differentially expressed miRNAs and mRNAs.<br />
  miRNA_lncRNA_groundtruth_LncBase_v2.0+NPInter_v4.0.csv: Experimentally validated miRNA-lncRNA interactions from LncBase v2.0 and NPInter v4.0.<br />
  miRNA_mRNA_groundtruth_miRTarBase_v8.0+TarBase_v8.0.csv: Experimentally validated miRNA-mRNA interactions from miRTarBase v8.0 and TarBase v8.0.<br />
  promise_validated_miR_lncR_mR_1679.el: The format of el about miRNA-target interactions.<br />
  ASD.Rdata: ASD expression datasets.<br />
  t_lncR_Exp_Autism.csv: LncRNA expression profiles of ASD samples.<br />
  t_lncR_Exp_Normal.csv: LncRNA expression profiles of normal samples.<br />
  t_miR_Exp_Autism.csv: MiRNA expression profiles of ASD samples.<br />
  t_miR_Exp_Normal.csv: MiRNA expression profiles of normal samples.<br />
  t_mR_Exp_Autism.csv: MRNA expression profiles of ASD samples.<br />
  t_mR_Exp_Normal.csv: MRNA expression profiles of normal samples.<br />
# The usage of ASDmiR
Paste all files into a single folder (set the folder as the directory of R environment), the workflow of ASDmiR is implemented in ASDmiR.R.
#  Quick example to use ASDmiR
For uncovering miRNA regulation related to ASD, we prepare ASD-related miRNA, lncRNA and mRNA expression profiles. Paste the datasets, run script of quick example version (Quick_examples_ASDmiR.R) and source file (zzz.R) into a single folder (set the folder as the directory of R environment).
```R
# Load required R package
library(miRLAB)
library(limma)

# Load utility functions
source("zzz.R")

# Load prepared datasets
raw_Autism <- load("ASD.RData") 

# Differentially expressed analysis
DiffExpAna_miR_lncR<-DiffExpAnalysis("t_miR_Exp_Autism.csv", "t_miR_Exp_Normal.csv", "t_lncR_Exp_Autism.csv", "t_lncR_Exp_Normal.csv", topkmiR = 100, topkmR = 300, p.miR = 1, p.mR = 1)
DiffExpAna_miR_mR<-DiffExpAnalysis("t_miR_Exp_Autism.csv", "t_miR_Exp_Normal.csv", "t_mR_Exp_Autism.csv", "t_mR_Exp_Normal.csv", topkmiR = 100, topkmR = 4000, p.miR = 1, p.mR = 1)

# Identification of miRNA-associated regulatory networks by ProMISe
cause = 1:100 #column of 1:35 are miRNAs
effect_lncR = 101:400 #column of 101:400 are lncRNAs
effect_mR = 101:4100 #column of 101:4100 are mRNAs

##predict miRNA targets using ProMISe
DiffExp_promise_miR_lncR <- ProMISe("DiffExp_miR_lncR.csv", cause, effect_lncR)
DiffExp_promise_miR_mR <- ProMISe("DiffExp_miR_mR.csv", cause, effect_mR)

## Select top 200 targets(lncRNAs and miRNAs) of each miRNA
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

## Validation of miRNA-target interactions predicted by ProMISe method
library(miRLAB)
promise_miR_lncRtop200_Validated <- Validation(promise_100miR_200lncR, datacsv = "miRNA_lncRNA_groundtruth_LncBase_v2.0+NPInter_v4.0.csv")
promise_miR_mRtop200_Validated <- Validation(promise_100miR_200mR, datacsv = "miRNA_mRNA_groundtruth_miRTarBase_v8.0+TarBase_v8.0.csv")
##  miRNA-target interactions using ProMISe
promise_validated_miR_lncR_mR_1679 <- rbind(promise_miR_lncRtop200_Validated[[1]],promise_miR_mRtop200_Validated[[1]])

# Identification of miRNA target modules
library(biclique)
library(Rcpp)
bi.format(filename = "promise_validated_miR_lncR_mR_1679.el", filetype = 0) 
miR_target_biclique <- bi.clique(filename = "promise_validated_miR_lncR_mR_1679.el", left_least = 3, right_least = 3, filetype = 0)

# Identification of miRNA sponge interaction network by Sensitivity Partial Pearson Correlation (SPPC)
miR_sponge_sppc <- spongeMethod(promise_validated_miR_lncR_mR_1679[,1:2], DEA_miR_mR_lncR, padjustvaluecutoff = 0.05, senscorcutoff = 0.25, method = "sppc")

# Identification of miRNA sponge modules by Markov Cluster Algorithm (MCL) and enrichment analysis.
library(miRspongeR)
miR_sppc_MCL <- netModule(miR_sponge_sppc_0.25[,1:2], method = "MCL",modulesize = 3, save = TRUE)

# Disease and functional enrichment analysis of miRNA target module
library(miRspongeR)
module_DEA <- moduleDEA(miR_sppc_MCL, ont = "DO", OrgDb = "org.Hs.eg.db", padjustvaluecutoff = 0.05, padjustedmethod = "BH")
module_FEA <- moduleFEA(miR_sppc_MCL,ont = "ALL", padjustvaluecutoff = 0.05, padjustedmethod = "BH")
```
