ValidateAll_groundtruth <- function (CEmatrix, topk, groundtruth, LFC, downreg = TRUE) 
{
  if (!file.exists("logFC.imputed.rda")) 
    stop("Please download the transfection data from nugget.unisa.edu.au/Thuc/miRLAB/logFC.imputed.rda for validation. \n       If you would like to use experimentally confirmed data only, please use the Validation function")
  causes = 1:ncol(CEmatrix)
  NoExpConfirmed = c()
  names = c("miRNA", "mRNA", "Correlation")
  K1 = c()
  ResultK1 = matrix(ncol = 3)
  colnames(ResultK1) = names
  K2 = c()
  pE = c()
  S = nrow(CEmatrix)
  for (i in causes) {
    miRtopk = bRank(CEmatrix, i, topk, downreg)
    miRall = bRank(CEmatrix, i, S)
    temp1 = Validation(miRtopk, groundtruth)
    temp2 = Validation(miRall, groundtruth)
    pvalE = 1 - phyper((temp1[[2]] - 1), temp2[[2]], (S - temp2[[2]]), topk)
    pE = c(pE, pvalE)
    cat("EXPERIMENT:  ", colnames(CEmatrix)[i], ": ", temp1[[2]], "with ", temp2[[2]], " in the groundtruth. pvalue: ", pvalE, "\n")
    NoExpConfirmed = c(NoExpConfirmed, temp1[[2]])
    K1 = c(K1, temp2[[2]])
    if (temp1[[2]] > 0) 
      ResultK1 = rbind(ResultK1, temp1[[1]])
  }
  pvalueE = 1 - phyper((sum(NoExpConfirmed) - 1), sum(K1), (S * ncol(CEmatrix) - sum(K1)), topk * ncol(CEmatrix))
  cat("number of confirms by experiments: ", sum(NoExpConfirmed), ", pvalue: ", pvalueE, "\n")
  cat("groundtruths in experiments: ", sum(K1), "\n")
  result = list(ResultK1[2:nrow(ResultK1), ], cbind(K1,NoExpConfirmed, pE))
}

experiment_groundtruth <- function (allmethods, topk, Expgroundtruth, LFC, downreg) 
{
  psv = ValidateAll_groundtruth(allmethods[[1]], topk, Expgroundtruth, LFC, downreg)
  spv = ValidateAll_groundtruth(allmethods[[2]], topk, Expgroundtruth, LFC, downreg)
  kendallv = ValidateAll_groundtruth(allmethods[[3]], topk, Expgroundtruth, LFC, downreg)
  dcovv = ValidateAll_groundtruth(allmethods[[4]], topk, Expgroundtruth, LFC, downreg)
  hoeffdingv = ValidateAll_groundtruth(allmethods[[5]], topk, Expgroundtruth, LFC, downreg)
  miv = ValidateAll_groundtruth(allmethods[[6]], topk, Expgroundtruth, LFC, downreg)
  idav = ValidateAll_groundtruth(allmethods[[7]], topk, Expgroundtruth, LFC, downreg)
  rdcv = ValidateAll_groundtruth(allmethods[[8]], topk, Expgroundtruth, LFC, downreg)
  lassov = ValidateAll_groundtruth(allmethods[[9]], topk, Expgroundtruth, LFC, downreg)
  elasticv = ValidateAll_groundtruth(allmethods[[10]], topk, Expgroundtruth,LFC, downreg)
  zsv = ValidateAll_groundtruth(allmethods[[11]], topk, Expgroundtruth, LFC, downreg)
  promisev = ValidateAll_groundtruth(allmethods[[12]], topk, Expgroundtruth, LFC, downreg)
  miRs = colnames(allmethods[[1]])
  result1 = cbind(psv[[2]], spv[[2]][, 2:3], kendallv[[2]][, 2:3], dcovv[[2]][, 2:3], hoeffdingv[[2]][, 2:3], miv[[2]][, 2:3], idav[[2]][, 2:3], rdcv[[2]][, 2:3], lassov[[2]][, 2:3], elasticv[[2]][, 2:3], zsv[[2]][, 2:3], promisev[[2]][, 2:3])
  rownames(result1) = miRs
  return(result1)
}


filterAndCompare_groundtruth <- function (allresults, noVal) 
{
  ExpResult = allresults
  temp1 = apply(ExpResult, 1, function(x) all(x[c(2, 4, 6, 
                                                  8, 10, 12, 14, 16, 18, 20, 22, 24)] > (noVal - 1)))
  ExpResult = ExpResult[temp1, ]
  if (is.matrix(ExpResult)) {
    ExpResult = ExpResult[, c(2, 4, 6, 8, 10, 12, 14, 16, 
                              18, 20, 22, 24)]
    colnames(ExpResult) = c("Pearson", "Spearman", "Kendall", 
                            " Dcov", "Hoeffding", "MI", "IDA", "RDC", "Lasso", 
                            "Elastic", "Z-score", "ProMISe")
  }
  else {
    tt = as.matrix(ExpResult)
    ExpResult = t(tt)
    ExpResult = ExpResult[, c(2, 4, 6, 8, 10, 12, 14, 16, 
                              18, 20, 22, 24)]
    names(ExpResult) = c("Pearson", "Spearman", "Kendall", 
                         " Dcov", "Hoeffding", "MI", "IDA", "RDC", "Lasso", 
                         "Elastic", "Z-score", "ProMISe")
  }

  if (is.matrix(ExpResult)) {
    numberExpResult = apply(ExpResult, 2, sum)
    ranking = apply(ExpResult, 1, function(i) rank(i))
    ranking = t(ranking)
    rankExpResult = apply(ranking, 2, sum)
  }
  if (is.matrix(ExpResult)) {
    Exp = list(ExpResult, numberExpResult, rankExpResult)
  }
  else Exp = ExpResult
  result = Exp
}
