
# SERE SCORE

# data cleaning
obs.count <- read.table("Output_ceres.txt")
obs.count$contig <- NULL
obs.count$start <- NULL
obs.count$end <- NULL
obs.count$ids <- NULL
obs.count$perc_cov <- NULL 
obs.count$type <- NULL
sere.score <- function(obs.count, TH = 1) {
  num.samp <- ncol(obs.count);
  tot.count <- sum(obs.count);
  samp.count <- colSums(obs.count);
  row.count <- rowSums(obs.count);
  
# filter those genes for which total reads for all samples are > TH; where TH = 1
  idx <- which(row.count > 1);
  obs.count <- obs.count[idx, ];
  row.count <- row.count[idx];
  num.genes <- nrow(obs.count);
  
# expected counts
  expt.count <- data.frame()
  expt.count <- matrix(NA, nrow = num.genes, ncol = num.samp);
  for(gene.idx in seq(num.genes)) {
    expt.count[gene.idx, ] <- samp.count * row.count[gene.idx] / tot.count;
  }
  
# sere score
  disp.sum <- sum((obs.count - expt.count) ^ 2 / expt.count);
  sere <- sqrt(disp.sum / (num.genes * (num.samp - 1))); 
}

# SERE DENODROGRAM 
sere.dendro <- function(obs.count, TH = 1) {
  # number of samples
  num.samp <- ncol(obs.count);
  col.names <- names(obs.count);
  
# distance matrix
  dist.mat <- matrix(NA, nrow = num.samp, ncol = num.samp);
  rownames(dist.mat) <- col.names;
  colnames(dist.mat) <- col.names;	
  
# fill the distance matrix
  for(i in seq(num.samp)) {
    for(j in i : num.samp) {
      print(paste('Computing SERE score for : ', col.names[i], col.names[j]));
      dist.mat[i, j] <- sere.score(obs.count[, c(i, j)], TH);
      dist.mat[j, i] <- sere.score(obs.count[, c(i, j)], TH);
    }
  }
  
# make the dendrogram
  clust <- hclust(as.dist(dist.mat), method = "complete");
  plot(clust);
  
  return(clust);
}

# changing labels AFTER generating the dendrogram 
install.packages('dendextend')
library(dendextend)
clust2 <- clust
labels(clust2) <- c('S18_ISEQ_50ng_B2','S19_ISEQ_50ng_B2','S22_NS_50ng_B2','S23_NS_50ng_B2','S14_ISEQ_50ng_B2','S24_NS_50ng_B2','S25_NS_50ng_B2','S18_NS_50ng_B2','S21_NS_50ng_B2','S19_NS_50ng_B2','S20_NS_50ng_B2','S20_ISEQ_50ng_B2','S21_ISEQ_50ng_B2','S17_ISEQ_50ng_B2','S15_ISEQ_50ng_B2','S16_ISEQ_50ng_B2','S7_NS_10ng_B1','S8_NS_50ng_B1','S5_NS_10ng_B1','S6_NS_10ng_B1','S5_ISEQ_50ng_B1','S12_ISEQ_10ng_B2','S17_NS_10ng_B1','S13_ISEQ_10ng_B2','S2_ISEQ_10ng_B1','S9_NS_50ng_B1','S4_ISEQ_50ng_B1','S3_ISEQ_10ng_B1','S16_NS_10ng_B2','S1_ISEQ_10ng_B1')

# comparing the two dendrograms 
### par(mfrow = c(1,2)) [for small dendrograms]
plot(clust, main = "Cluster Dendrogram (Original)")
plot(clust2, main = "Cluster Dendrogram (Modified Labels)")
