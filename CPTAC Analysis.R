library(DEqMS)
library(mixtools)
library(mclust)
cptac <- read.csv("CPTAC_before_norm_iTRAQ.csv", row.names = "Gene")
#remove rows with NAs
cptac <- na.omit(cptac)

psm <- read.csv("iTRAQ_PSM_Counts.csv", row.names = "Gene")
#remove rows with NAs
psm <- na.omit(psm)

# Normalize as in paper
cptac_norm <- cptac
# row-wise zscore
saveRDS(t(scale(t(cptac_norm))), "ZscoreCPTAC.rds")

# R file from the docker image in PMID: 33212010
source("two_comp_normalize.R")
set.seed(123)
for (i in colnames(cptac_norm)) {
    cptac_norm[,i] <- two.comp.normalize(cptac_norm[,i])$norm.sample
}

cond <- ifelse(grepl("TN", colnames(cptac_norm)), "TN", "LuA")
cond <- as.factor(cond)
# The function model.matrix is used to generate the design matrix
design <- model.matrix(~0+cond) # 0 means no intercept for the linear model
colnames(design) <- gsub("cond","",colnames(design))

contrast <-  makeContrasts(TN - LuA,levels=design)
fit1 <- lmFit(cptac_norm, design)
fit2 <- contrasts.fit(fit1,contrasts = contrast)
fit3 <- eBayes(fit2)


psm.count.table <- data.frame(count = rowMins(as.matrix(psm)))
rownames(psm.count.table) <- rownames(psm)
fit3$count <- psm.count.table[rownames(fit3$coefficients),"count"]
fit4 <- spectraCounteBayes(fit3)


DEqMS.results <- outputResult(fit4,coef_col = 1)

prot_up <- DEqMS.results[DEqMS.results$sca.adj.pval < 0.05 & DEqMS.results$logFC >= 1,]
write.table(prot_up$gene, "CPTAC_UP.txt",col.names = F, row.names = F, quote = F)

