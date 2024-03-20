# 加载R???
library(dplyr)
library(ggplot2)
library(MutationalPatterns)
library(ggprism)
library(grid)
library(gridExtra)
library(ggsignif)
library(ggsci)
# 列出vcf的路???
vcf_files <- list.files(getwd())[grep(".vcf", list.files(getwd()))]
TYPE <- colnames(type_occurrences)[1:6]
GROUP <- c("Low grade", "High grade", "Cancer")
# 设置vcf文件对应的样本名???
sample_names <- c()
for (i in 1:length(vcf_files)){
  s <- vcf_files[i]
  sample_names <- c(sample_names, strsplit( s,"\\.")[[1]][1])
}


# 加载参考基因组
library(BSgenome.Hsapiens.UCSC.hg19)
ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
# 读取vcf文件
vcfs <- read_vcfs_as_granges(vcf_files, sample_names, ref_genome)
type_occurrences <- mut_type_occurrences(vcfs, ref_genome)

### Figure 2C
plot_SNV_distribution(type_occurrences)


### Figure 2D
mut_mat <- mut_matrix(vcf_list = vcfs, ref_genome = ref_genome)  
group <- c()
type <- c()
value <- c()
for (i in 1:length(sample_names)){
  s <- sample_names[i]
  group <- c(group, substr(s, nchar(s),nchar(s)))
}
group[group == "C"] <- "Cancer"
group[group == "H"] <- "High grade"
group[group == "L"] <- "Low grade"
mut_mat <- mut_mat[, group %in% c("High grade", "Cancer")]


plot_mutation_pattern(mut_mat, size = 4, colors = NA,  condensed = FALSE) 
nmf_res <- extract_signatures(mut_mat, rank = 5, nrun = 10)
colnames(nmf_res$signatures) <- paste("Signature", c(1:5))
rownames(nmf_res$contribution) <- paste("Signature", c(1:5))
plot_96_profile_M(nmf_res$signatures, size = c(0.08, 2.5, 0.01, 4, 2.5,1),condensed = TRUE)


plot_contribution_absolute_box(contribution = nmf_res$contribution, signatures = nmf_res$signature)
plot_contribution_relative_bar(contribution = nmf_res$contribution, signatures = nmf_res$signature)



