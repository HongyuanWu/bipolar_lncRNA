
# library(DESeq2)
# library(ggplot2)
# library(pheatmap)
# library(vsn)
# library(devtools)
# library(rafalib)
# library(BiocStyle)
# library(rmarkdown)
# library(geneplotter)
# library(plyr)
# library(LSD)
# library(gplots)
# library(RColorBrewer)
# library(stringr)
# library(topGO)
# library(genefilter)
# library(biomaRt)
# library(dplyr)
# library(EDASeq)
# library(fdrtool)
# library("sva")
library(limma)
library(edgeR)
library(qvalue)
library(tidyverse)
library(dplyr)
library(sva)
# library(gee)
# library(geepack)
# options(contrasts=c("contr.sum","contr.poly"))
# library(clusterProfiler)
# library(RColorBrewer)

# setwd("/data/NHLBI_BCB/Fayaz/34-lncRNAKB-BP-JHU/02-differential-expression-analysis/00-gene-level-redo-07222020/00-amygdala")
counts_data = read_tsv(file="/data/NHLBI_BCB/Fayaz/34-lncRNAKB-BP-JHU/01-lncRNAKB-gene-level-counts/amygdala.gene.featurecount.txt")

# setwd("/data/NHLBI_BCB/Fayaz/34-lncRNAKB-BP-JHU/02-differential-expression-analysis/00-gene-level-redo-07222020/00-sACC")
# counts_data = read_tsv(file="/data/NHLBI_BCB/Fayaz/34-lncRNAKB-BP-JHU/01-lncRNAKB-gene-level-counts/sacc.gene.featurecount.txt")

counts_data

phenodata = read_tsv(file="/data/NHLBI_BCB/Fayaz/34-lncRNAKB-BP-JHU/00-phenotype/amygdala_subjects_pheno.txt")
# phenodata = read_tsv(file="/data/NHLBI_BCB/Fayaz/34-lncRNAKB-BP-JHU/00-phenotype/sacc_subjects_pheno.txt")

phenodata

group1 = "Bipolar"
group2 = "Control"
feature = "gene"
# bregion = "sacc"
bregion = "amy"

phenodata = filter(phenodata, dx == group1 | dx == group2)
phenodata
counts_data = select(counts_data, c("gene",phenodata$SampleID))
counts_data

annot = read_csv(file="/data/NHLBI_BCB/Fayaz/34-lncRNAKB-BP-JHU/02-differential-expression-analysis/00-gene-level/annot_format_v7.csv")
annot

## filtering
#Keep lncRNAs with least 1 counts/reads in 90% of replicates/samples
counts_data_lncRNA = counts_data %>% semi_join(filter(annot, gene_type == "lncRNA"), by = c("gene" = "lnckbID"))
counts_data_lncRNA
counts_data_lncRNA = counts_data_lncRNA %>% filter(rowSums(counts_data_lncRNA>1) >= as.integer(0.90*(dim(counts_data_lncRNA)[2]-1)))
counts_data_lncRNA

#Keep protein coding genes with least 5 counts/reads in 90% of replicates/samples
counts_data_pc = counts_data %>% semi_join(filter(annot, gene_type == "protein_coding"), by = c("gene" = "lnckbID"))
counts_data_pc
counts_data_pc = counts_data_pc %>% filter(rowSums(counts_data_pc>5) >= as.integer(0.90*(dim(counts_data_pc)[2]-1)))
counts_data_pc

################################################################################################################
#
# Features were then filtered if they did not exceed the following mean expression thresholds: 
# 0.25 RPKM for genes
# 0.30 RPKM for exons
# 0.35 RP10M for junctions
# 0.40 TPM for transcripts 
# as determined with the expression_cutoff() function from the jaffelab R package (33). 
# Cutoffs were applied based on the samples from the amygdala and subgenual anterior cingulate cortex.
#
################################################################################################################
# pdf('test_expression_cutoff.pdf', width = 12)
# expression_cutoff(y2)
# dev.off()

################################################################################################################
#
#
# Graph of different threshold vs. number of lncRNA and protein coding <genes> compared to flagship count matrix 
#
#
################################################################################################################
# load("/data/NHLBI_BCB/Fayaz/31-leafcutter/lieber/andrew_rda/eQTL_expressed_rse_amygdala.rda")
load("/data/NHLBI_BCB/Fayaz/31-leafcutter/lieber/andrew_rda/eQTL_expressed_rse_sacc.rda")
genes_passed_threshold = as_tibble(rownames(rse_gene))
genes_passed_threshold
gencodev25annotation = read_tsv(file="/data/NHLBI_BCB/Fayaz/31-leafcutter/lieber/00-leafcutter/02-step4-run-fastQTL/01-annotate-sQTL-results-pz/gencode.v25.basic.annotation.sorted.bed", col_names=FALSE)
gencodev25annotation
genes_passed_threshold_annotated_gencodev25 = gencodev25annotation %>% semi_join(genes_passed_threshold, by = c("X5" = "value"))
count_of_types_of_genes = data.frame(table(genes_passed_threshold_annotated$X6))
colnames(count_of_types_of_genes) = c("gene_type","count")
count_of_types_of_genes

################################################################################################################
#
#
# PCA plot by brain region lncRNAKB GTEx v7
# All lncRNA 
#
#
################################################################################################################
subjecttobrainregions = read_tsv(file="/data/NHLBI_BCB/Fayaz/34-lncRNAKB-BP-JHU/GTEx_brain_samples_by_region.txt")
subjecttobrainregions

# pca = readRDS(file="/data/NHLBI_BCB/Fayaz/21-GTEx-data/12-gene-expression-boxplots-lncRNAKB-v6/log2_TPM_lncRNAsonly_lncRNAKBv6_pca.rds")
pca = readRDS(file="/data/NHLBI_BCB/Fayaz/21-GTEx-data/12-gene-expression-boxplots-lncRNAKB-v6/log2_TPM_pcgssonly_lncRNAKBv6_pca.rds")

pcs <- rownames_to_column(data.frame(pca$x), var = "subjectID") %>% as_tibble()
pcs

pcswbrainregion = pcs %>% inner_join(subjecttobrainregions, by = c("subjectID" = "SubjectID"))
pcswbrainregion

subjecttobrainregions = subjecttobrainregions %>% semi_join(pcswbrainregion, by = c("SubjectID" = "subjectID"))
subjecttobrainregions

library(randomcoloR)
n = length(levels(as.factor(subjecttobrainregions$Brain_Region)))
col = distinctColorPalette(k=n, runTsne=T)

BrainRegionColor=col[unclass(as.factor(subjecttobrainregions$Brain_Region))]
BrainRegionColor[1:10]

brainregions = unique(subjecttobrainregions$Brain_Region)
brainregions


samples_each_brainregion=c("SRR1078879",
"SRR1073755",
"SRR1075898",
"SRR1070986",
"SRR1071289",
"SRR1088389",
"SRR1072367",
"SRR1081741",
"SRR1069188",
"SRR1074164",
"SRR1072178",
"SRR1073143",
"SRR1077687"
)


# pcswtissue %>%
#	ggplot(aes(x=PC1,y=PC2)) +
#	geom_point(aes(colour=Tissuecolor), alpha=0.20) +
#	guides(colour=guide_legend(override.aes = list(alpha=1, size=1))) +
#	scale_color_discrete(name="Tissue", labels=tissues) +
#	theme(legend.text=element_text(size=7)) + 
#	labs(x=spc1, y=spc2, title="Principal Components_PC1.vs.PC2, lncRNAs, lncRNAKB, GTEx, v7") +
#	ggrepel::geom_text_repel(data = pcswtissue %>% filter(subjectID %in% samples_each_tissue), aes(label=Tissue), colour="black", size=2.5)
	

# dev.off()

eigs <- pca$sdev^2
e1=eigs[1] / sum(eigs) * 100
e2=eigs[2] / sum(eigs) * 100
e3=eigs[3] / sum(eigs) * 100
e4=eigs[4] / sum(eigs) * 100
e5=eigs[5] / sum(eigs) * 100
spc1=sprintf('PC1 (%% %.1f)',e1)
spc2=sprintf('PC2 (%% %.1f)',e2)
spc3=sprintf('PC3 (%% %.1f)',e3)
spc4=sprintf('PC4 (%% %.1f)',e4)
spc5=sprintf('PC5 (%% %.1f)',e5)

pcswbrainregion %>%
	ggplot(aes(x=PC1,y=PC2)) +
	geom_point(aes(colour=BrainRegionColor), alpha=0.9) +
	theme(legend.position="none") + 
	labs(x=spc1, y=spc2, title="Principal Components_PC1.vs.PC2, protein coding, lncRNAKB, GTEx, v7, brain") +
	ggrepel::geom_text_repel(data = pcswbrainregion %>% filter(subjectID %in% samples_each_brainregion), aes(label=Brain_Region), colour="black", size=3)
	
dev.off()

counts_data = bind_rows(counts_data_lncRNA, counts_data_pc)
counts_data

################################################################################################################
#
#
# How many lncRNA in GENCODE vs. Other annotations survive cut-offs
# 
#
#
################################################################################################################
counts_data_lncRNA_annot = counts_data %>% inner_join(filter(annot, gene_type == "lncRNA"), by = c("gene" = "lnckbID"))
counts_data_lncRNA_annot
data.frame(table(counts_data_lncRNA_annot$source))

counts_data = data.frame(counts_data[,2:dim(counts_data)[2]], row.names=counts_data$gene)
counts_data[1:10,1:10]
dge = DGEList(counts = counts_data)
dge = calcNormFactors(dge, method="TMM")

################################################################################################################
# covariates used by AJ
# modSep = model.matrix(~Dx + AgeDeath + Sex + snpPC1 + snpPC2 + snpPC3 + mitoRate + rRNA_rate + totalAssignedGene + RIN + ERCCsumLogErr, data=colData(rse_gene))
# age at death
# sex
# mitochondrial rate
# ribosomal RNA rate
# total gene assignment (EXCLUDED)
# RIN
# ERCC spike-in rates
# the top 3 principal components of the genotype data (calculated using ~100,000 LD-independent SNPs) (EXCLUDED)
# X quality surrogate variables (qSVs).  
library(SummarizedExperiment)
library(stringr)
# load("/data/NHLBI_BCB/Fayaz/31-leafcutter/lieber/andrew_rda/eQTL_expressed_rse_amygdala.rda")
load("/data/NHLBI_BCB/Fayaz/31-leafcutter/lieber/andrew_rda/eQTL_expressed_rse_sacc.rda")
all_covariates = as_tibble(colData(rse_gene))
all_covariates
covariates = select(all_covariates, c("PrimaryDx", "AgeDeath", "Sex", "mitoRate", "rRNA_rate", "RIN", "ERCCsumLogErr", "RNum", "BrNum"))
covariates
covariates = mutate(covariates, PrimaryDx = ifelse(PrimaryDx == "Other","Bipolar",as.character(PrimaryDx)))
covariates
covariates = covariates %>% semi_join(phenodata, by = c("RNum" = "SampleID"))
covariates
covariates = mutate(covariates, Dx = ifelse(PrimaryDx == "Bipolar",1,0))
covariates
load("/data/NHLBI_BCB/Fayaz/31-leafcutter/lieber/00-leafcutter/03-step5-run-differential-introl-excision-analysis/qSV_mat.rda")
qSV_mat = as_tibble(qSV_mat, rownames="SampleID")
qSV_mat
qSV_mat = qSV_mat %>% filter(str_detect(SampleID, "_sACC"))
# qSV_mat = qSV_mat %>% filter(str_detect(SampleID, "_Amygdala"))
qSV_mat
qSV_mat = qSV_mat %>% mutate(SampleID = str_replace_all(qSV_mat$SampleID, c("_Amygdala"="", "_sACC"="")))
qSV_mat
covariates = covariates %>% inner_join(qSV_mat, by = c("BrNum" = "SampleID"))
covariates
covariates = arrange(covariates, PrimaryDx)
covariates
#	load("/data/NHLBI_BCB/Fayaz/31-leafcutter/lieber/00-leafcutter/03-step5-run-differential-introl-excision-analysis/qSV_mat.rda")
#	rownames(qSV_mat) = str_replace_all(rownames(qSV_mat), c("_Amygdala"="", "_sACC"=""))
#	covs = merge(x=covariates, y=qSV_mat, by="row.names", all.x=T)
#	load("/data/NHLBI_BCB/Fayaz/31-leafcutter/lieber/andrew_rda/zandiHyde_bipolar_Genotypes_n511.rda")
#	snpPCs = mds[,1:3]	
#	covs2 = merge(x=covs, y=snpPCs, by.x="Row.names", by.y="row.names", all.x=T)
#	covs4 = data.frame(covs2[,2:dim(covs2)[2]], row.names=covs2$Row.names)
#	amygdala_subjects = read.table(file="/data/NHLBI_BCB/Fayaz/31-leafcutter/lieber/00-leafcutter/03-step5-run-differential-introl-excision-analysis/amy_transposed_sex-bp-5geno-pc_covs_07182019.txt", header=F, row.names=1)
#	covs5 = covs4[rownames(amygdala_subjects),]
#	covs5[,"Sex"] = str_replace_all(covs5[,"Sex"] , c("M"="1", "F"="0"))
#	covs5[,"Dx"] = str_replace_all(covs5[,"Dx"] , c("Bipolar"="1", "Control"="0", "Other"="0"))
#	write.table(covs5, file="/data/NHLBI_BCB/Fayaz/31-leafcutter/lieber/00-leafcutter/03-step5-run-differential-introl-excision-analysis/amy_transposed_covs_11212019.txt", quote=F, col.names=F)
# mod1 = model.matrix(~Dx, data=covariates)
# mod0 = cbind(mod1[,1])
# svseq = svaseq(as.matrix(counts_data),mod1,mod0,n.sv=5)$sv
# pdf(paste("bipseq_lncRNAKB_svaplot_plot_all_genes_pairs_plot_5svs_",feature,"_",bregion,"_",group1,"_",group2,".pdf", sep=""))
# pairs(svseq[,1:5], main=paste("scatter plot matrix of Surrogate Variable_5SVs_",feature,"_",bregion,"_",group1,"_",group2, sep=""), pch=21, bg=c('#e41a1c','#377eb8')[as.factor(phenodata$dx)], cex.main=0.75)
# dev.off()
#
################################################################################################################

# design = model.matrix(~Dx + AgeDeath + Sex + mitoRate + rRNA_rate + RIN + ERCCsumLogErr + svseq[,1:5], data=covariates)
design = model.matrix(~Dx + AgeDeath + Sex + mitoRate + rRNA_rate + RIN + ERCCsumLogErr + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18, data=covariates)
v = voom(dge, design) #Transform count data to log2-counts per million (logCPM+0.5), estimate the mean-variance relationship and use this to compute appropriate observation-level weights. The data are then ready for linear modelling.
fit = lmFit(v, design)
fit = eBayes(fit)
log2FC = fit$coefficients[, 2]
pvalue = fit$p.value[, 2]
coefs = fit$coefficients[, 2]
qvalue = qvalue(pvalue)$qvalues
gene_res = data.frame(log2FC, pvalue, qvalue)
gene_res = rownames_to_column(gene_res, var = "gene") %>% as_tibble()
gene_res
gene_res = gene_res %>% inner_join(annot, by = c("gene" = "lnckbID"))
gene_res
gene_res = arrange(gene_res, qvalue)
gene_res
write.csv(gene_res, paste("bipseq_lncRNAKB","_",feature,"_",bregion,"_",group1,"_",group2,".csv", sep=""), row.names=FALSE)

counts_data_normalized = as.data.frame(v$E) %>% as_tibble()
counts_data_normalized = rownames_to_column(as.data.frame(v$E), var = "gene") %>% as_tibble()
counts_data_normalized
# counts_data_lncRNA = counts_data %>% semi_join(filter(annot, gene_type == "lncRNA"), by = c("gene" = "lnckbID"))

# pca plot - using all genes (only filtered above for low counts)
pca = prcomp(t(counts_data_normalized %>% select(-one_of("gene"))), scale=F)
eigs <- pca$sdev^2
e1=eigs[1] / sum(eigs) * 100
e2=eigs[2] / sum(eigs) * 100
e3=eigs[3] / sum(eigs) * 100
e4=eigs[4] / sum(eigs) * 100
e5=eigs[5] / sum(eigs) * 100
spc1=sprintf('PC1 (%% %.1f)',e1)
spc2=sprintf('PC2 (%% %.1f)',e2)
spc3=sprintf('PC3 (%% %.1f)',e3)
spc4=sprintf('PC4 (%% %.1f)',e4)
spc5=sprintf('PC5 (%% %.1f)',e5)

pdf(paste("bipseq_lncRNAKB_pca_plot_all_genes_pairs_plot_first_five_PCs_TMM_CPM_",feature,"_",bregion,"_",group1,"_",group2,".pdf", sep=""))
pairs(pca$x[,1:5], main=paste("scatter plot matrix of Principal Components_first_five_PCs_",feature,"_",bregion,"_",group1,"_",group2, sep=""), pch=21, bg=c('#e41a1c','#377eb8')[as.factor(phenodata$dx)], cex.main=0.75)
dev.off()

##########################################################################################
#
# countfiltered >> log2CPM0.5 >> TMM >> lncRNA >> PCA
#
##########################################################################################
pca = prcomp(t(counts_data_normalized %>% semi_join(filter(gene_res, gene_type=="lncRNA"), by = c("gene" = "gene")) %>% select(-one_of("gene"))), scale=F)
eigs <- pca$sdev^2
e1=eigs[1] / sum(eigs) * 100
e2=eigs[2] / sum(eigs) * 100
e3=eigs[3] / sum(eigs) * 100
e4=eigs[4] / sum(eigs) * 100
e5=eigs[5] / sum(eigs) * 100
spc1=sprintf('PC1 (%% %.1f)',e1)
spc2=sprintf('PC2 (%% %.1f)',e2)
spc3=sprintf('PC3 (%% %.1f)',e3)
spc4=sprintf('PC4 (%% %.1f)',e4)
spc5=sprintf('PC5 (%% %.1f)',e5)

pdf(paste("bipseq_lncRNAKB_pca_plot_lncRNA_first_five_PCs_log2CPM0.5_TMM",feature,"_",bregion,"_",group1,"_",group2,".pdf", sep=""))
pairs(pca$x[,1:5], main=paste("scatter plot matrix of Principal Components_first_five_PCs_",feature,"_",bregion,"_",group1,"_",group2, sep=""), pch=21, bg=c('#e41a1c','#377eb8')[as.factor(phenodata$dx)], cex.main=0.75)
dev.off()

pdf(paste("bipseq_lncRNAKB_pca_plot_lncRNA_PC1_PC2_log2CPM0.5_TMM_",feature,"_",bregion,"_",group1,"_",group2,".pdf", sep=""))
# Add extra space to right of plot area; change clipping to figure
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(pca$x[,1], pca$x[,2], main=paste("scatter plot matrix of Principal Components_PC1.vs.PC2_",feature,"_",bregion,"_",group1,"_",group2, sep=""), pch=16, col=c('#e41a1c','#377eb8')[as.factor(phenodata$dx)], cex=1.20, cex.main=0.80, xlab=spc1, ylab=spc2)
# text(pca$x[,1], pca$x[,2], labels=paste(rownames(phenodata),phenodata$dx,sep="_"), pos=3, cex=0.30)
legend("topright", inset=c(-0.23,0), legend=c(group1,group2), pch=c(16,16), col=rev(levels(as.factor(c('#e41a1c','#377eb8')[as.factor(phenodata$dx)]))), title="Group", cex=0.60)
dev.off()

##########################################################################################
#
# countfiltered >> log2CPM0.5 >> TMM >> lncRNA >> qvalue<0.05 >> PCA
#
##########################################################################################
pca = prcomp(t(counts_data_normalized %>% semi_join(filter(gene_res, qvalue < 0.05, gene_type=="lncRNA"), by = c("gene" = "gene")) %>% select(-one_of("gene"))), scale=F)
eigs <- pca$sdev^2
e1=eigs[1] / sum(eigs) * 100
e2=eigs[2] / sum(eigs) * 100
e3=eigs[3] / sum(eigs) * 100
e4=eigs[4] / sum(eigs) * 100
e5=eigs[5] / sum(eigs) * 100
spc1=sprintf('PC1 (%% %.1f)',e1)
spc2=sprintf('PC2 (%% %.1f)',e2)
spc3=sprintf('PC3 (%% %.1f)',e3)
spc4=sprintf('PC4 (%% %.1f)',e4)
spc5=sprintf('PC5 (%% %.1f)',e5)

pdf(paste("bipseq_lncRNAKB_pca_plot_lncRNA_qvalue0.05_first_five_PCs_log2CPM0.5_TMM",feature,"_",bregion,"_",group1,"_",group2,".pdf", sep=""))
pairs(pca$x[,1:5], main=paste("scatter plot matrix of Principal Components_first_five_PCs_",feature,"_",bregion,"_",group1,"_",group2, sep=""), pch=21, bg=c('#e41a1c','#377eb8')[as.factor(phenodata$dx)], cex.main=0.75)
dev.off()

pdf(paste("bipseq_lncRNAKB_pca_plot_lncRNA_qvalue0.05_PC1_PC2_log2CPM0.5_TMM_",feature,"_",bregion,"_",group1,"_",group2,".pdf", sep=""))
# Add extra space to right of plot area; change clipping to figure
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(pca$x[,1], pca$x[,2], main=paste("scatter plot matrix of Principal Components_PC1.vs.PC2_",feature,"_",bregion,"_",group1,"_",group2, sep=""), pch=16, col=c('#e41a1c','#377eb8')[as.factor(phenodata$dx)], cex=1.20, cex.main=0.80, xlab=spc1, ylab=spc2)
# text(pca$x[,1], pca$x[,2], labels=paste(rownames(phenodata),phenodata$dx,sep="_"), pos=3, cex=0.30)
legend("topright", inset=c(-0.23,0), legend=c(group1,group2), pch=c(16,16), col=rev(levels(as.factor(c('#e41a1c','#377eb8')[as.factor(phenodata$dx)]))), title="Group", cex=0.60)
dev.off()

# pca plot - using results with q-value < 0.05 or any other subset based on significant DE for e.g. fc, log2fc
pca = prcomp(t(counts_data_normalized %>% semi_join(filter(gene_res, log2FC < -0.50 | log2FC > 0.50, qvalue < 0.05), by = c("gene" = "gene")) %>% select(-one_of("gene"))), scale=F)
eigs <- pca$sdev^2
e1=eigs[1] / sum(eigs) * 100
e2=eigs[2] / sum(eigs) * 100
e3=eigs[3] / sum(eigs) * 100
e4=eigs[4] / sum(eigs) * 100
e5=eigs[5] / sum(eigs) * 100
spc1=sprintf('PC1 (%% %.1f)',e1)
spc2=sprintf('PC2 (%% %.1f)',e2)
spc3=sprintf('PC3 (%% %.1f)',e3)
spc4=sprintf('PC4 (%% %.1f)',e4)
spc5=sprintf('PC5 (%% %.1f)',e5)

pdf(paste("bipseq_lncRNAKB_pca_plot_genes_q_less_than_0.05_foldchange_1_pairs_plot_first_five_PCs_TMM_CPM_",feature,"_",bregion,"_",group1,"_",group2,".pdf", sep=""))
pairs(pca$x[,1:5], main=paste("scatter plot matrix of Principal Components_first_five_PCs_",feature,"_",bregion,"_",group1,"_",group2, sep=""), pch=21, bg=c('#e41a1c','#377eb8')[as.factor(phenodata$dx)], cex.main=0.75)
dev.off()

# pca plot - using results with q-value < 0.05 or any other subset based on significant DE for e.g. fc, log2fc
pca = prcomp(t(counts_data_normalized %>% semi_join(filter(gene_res, log2FC < -0.50 | log2FC > 0.50, qvalue < 0.05, gene_type=="protein_coding"), by = c("gene" = "gene")) %>% select(-one_of("gene"))), scale=F)
eigs <- pca$sdev^2
e1=eigs[1] / sum(eigs) * 100
e2=eigs[2] / sum(eigs) * 100
e3=eigs[3] / sum(eigs) * 100
e4=eigs[4] / sum(eigs) * 100
e5=eigs[5] / sum(eigs) * 100
spc1=sprintf('PC1 (%% %.1f)',e1)
spc2=sprintf('PC2 (%% %.1f)',e2)
spc3=sprintf('PC3 (%% %.1f)',e3)
spc4=sprintf('PC4 (%% %.1f)',e4)
spc5=sprintf('PC5 (%% %.1f)',e5)

pdf(paste("bipseq_lncRNAKB_pca_plot_proteincoding_genes_q_less_than_0.05_foldchange_1_pairs_plot_first_five_PCs_TMM_CPM_",feature,"_",bregion,"_",group1,"_",group2,".pdf", sep=""))
pairs(pca$x[,1:5], main=paste("scatter plot matrix of Principal Components_first_five_PCs_",feature,"_",bregion,"_",group1,"_",group2, sep=""), pch=21, bg=c('#e41a1c','#377eb8')[as.factor(phenodata$dx)], cex.main=0.75)
dev.off()

# pca plot - using results with q-value < 0.05 or any other subset based on significant DE for e.g. fc, log2fc
pca = prcomp(t(counts_data_normalized %>% semi_join(filter(gene_res, log2FC < -0.50 | log2FC > 0.50, qvalue < 0.05, gene_type=="lncRNA"), by = c("gene" = "gene")) %>% select(-one_of("gene"))), scale=F)
eigs <- pca$sdev^2
e1=eigs[1] / sum(eigs) * 100
e2=eigs[2] / sum(eigs) * 100
e3=eigs[3] / sum(eigs) * 100
e4=eigs[4] / sum(eigs) * 100
e5=eigs[5] / sum(eigs) * 100
spc1=sprintf('PC1 (%% %.1f)',e1)
spc2=sprintf('PC2 (%% %.1f)',e2)
spc3=sprintf('PC3 (%% %.1f)',e3)
spc4=sprintf('PC4 (%% %.1f)',e4)
spc5=sprintf('PC5 (%% %.1f)',e5)

pdf(paste("bipseq_lncRNAKB_pca_plot_lncRNA_genes_q_less_than_0.05_foldchange_1_pairs_plot_first_five_PCs_TMM_CPM_",feature,"_",bregion,"_",group1,"_",group2,".pdf", sep=""))
pairs(pca$x[,1:5], main=paste("scatter plot matrix of Principal Components_first_five_PCs_",feature,"_",bregion,"_",group1,"_",group2, sep=""), pch=21, bg=c('#e41a1c','#377eb8')[as.factor(phenodata$dx)], cex.main=0.75)
dev.off()


# pdf(paste("pca_plot_genes_q_less_than_0.05_pc1_vs_pc2_cpm_TMM_CPM_",feature,"_",group1,"_",group2,".pdf", sep=""))
# Add extra space to right of plot area; change clipping to figure
# par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
# plot(pca$x[,1], pca$x[,2], main=paste("scatter plot matrix of Principal Components_PC1.vs.PC2_",feature,"_",group1,"_",group2, sep=""), pch=16, col=c('blue','orange')[as.factor(phenodata$dx)], cex=1.20, cex.main=0.80, xlab=spc1, ylab=spc2)
# text(pca$x[,1], pca$x[,2], labels=paste(rownames(phenodata),phenodata$dx,sep="_"), pos=3, cex=0.30)
# legend("topright", inset=c(-0.23,0), legend=c(group1,group2), pch=c(16,16), col=levels(as.factor(c('#e41a1c','#377eb8')[as.factor(phenodata$dx)])), title="Group", cex=0.60)
# dev.off()

## Plotting abundance distribution, after applying filter, normalization
pdf(file=paste("bipseq_lncRNAKB_log2_TMM_CPM_",feature,"_",bregion,"_",group1,"_",group2,".pdf", sep=""))
boxplot(counts_data_normalized %>% select(-one_of("gene")) %>% rename_all(funs(str_replace(.,covariates$RNum,as.character(covariates$Dx)))), col=c('#e41a1c','#377eb8')[as.factor(phenodata$dx)], las=2, ylab='log2(cpm(counts)+0.5)', main=paste("cpm:",feature,"_",bregion,"_",group1,"_",group2), xlab="sample", cex.main=0.80, cex.axis=0.50)
dev.off()

## Plotting volcano plot, Differential Expression Analysis-qval
q.trans = -1*log(gene_res$qvalue, base = 10)
logfold = gene_res$log2FC
pdf(file=paste("bipseq_lncRNAKB_volcano_qval_",feature,"_",bregion,"_",group1,"_",group2,".pdf", sep=""))
plot(logfold,q.trans,type="n",ylab="-1*log10(q-value)",xlab="log2(fold change)",main=paste("cpm:",feature,"_",bregion,"_",group1,"_",group2), cex.main=0.80, xlim=range(logfold), ylim=range(q.trans))
points(logfold,q.trans,col="black",cex=0.65)
points(logfold[(q.trans>1.3&logfold>0.50)],q.trans[(q.trans>1.3&logfold>0.50)],col="red",pch=16,cex=0.65)
points(logfold[(q.trans>1.3&logfold<(-0.50))],q.trans[(q.trans>1.3&logfold<(-0.50))],col="green",pch=16,cex=0.65)
# text(logfold[(q.trans>1.3&logfold>1.0)],q.trans[(q.trans>1.3&logfold>1.0)],labels=as.character(data$genesymbol[(q.trans>1.3&logfold>1.0)]), pos=3, cex=0.40)
# text(logfold[(q.trans>1.3&logfold<(-1.0))],q.trans[(q.trans>1.3&logfold<(-1.0))],labels=as.character(data$genesymbol[(q.trans>1.3&logfold<(-1.0))]), pos=3, cex=0.40)
abline(h=1.3)
abline(v=-0.50)
abline(v=0.50)
dev.off()

gene_res = gene_res %>%
  mutate(
    gene_type_significant = case_when(
      qvalue < 0.05 & gene_type == "lncRNA" ~ "lncRNA",
      qvalue < 0.05 & gene_type == "protein_coding" ~ "protein-coding",
      TRUE                      ~ "not_significant"
    )
  )

ggplot(gene_res, aes(log2FC, -log10(qvalue))) + geom_point(aes(col=gene_type_significant)) + scale_color_manual(values=c("green","lightblue","red")) +
coord_cartesian(xlim = c(-1.0, 1.0), expand = TRUE) + 
geom_hline(aes(yintercept = 1.3), color="black", linetype="dashed", size=0.5)
ggsave(file=paste("bipseq_lncRNAKB_volcano_plot_qvalue",feature,"_",bregion,"_",group1,"_",group2,".pdf", sep=""))
dev.off()

## Plotting volcano plot, Differential Expression Analysis-pval
p.trans = -1*log(gene_res$pvalue, base = 10)
logfold = gene_res$log2FC
pdf(file=paste("bipseq_lncRNAKB_volcano_pval_",feature,"_",bregion,"_",group1,"_",group2,".pdf", sep=""))
plot(logfold,p.trans,type="n",ylab="-1*log10(p-value)",xlab="log2(fold change)",main=paste("LimmaVoom:",feature,"_",bregion,"_",group1,"_",group2), cex.main=0.80, xlim=c(-1.0,1.0), ylim=range(p.trans))
points(logfold,p.trans,col="black",cex=0.65)
points(logfold[(p.trans>1.3)], p.trans[(p.trans>1.3)],col="lightblue",pch=16,cex=0.65)
# points(logfold[(p.trans>1.3&logfold>1.0)],p.trans[(p.trans>1.3&logfold>1.0)],col="red",pch=16,cex=0.65)
# points(logfold[(p.trans>1.3&logfold<(-1.0))],p.trans[(p.trans>1.3&logfold<(-1.0))],col="green",pch=16,cex=0.65)
# text(logfold[(p.trans>1.3&logfold>1.0)],p.trans[(p.trans>1.3&logfold>1.0)],labels=as.character(data$genesymbol[(p.trans>1.3&logfold>1.0)]), pos=3, cex=0.40)
# text(logfold[(p.trans>1.3&logfold<(-1.0))],p.trans[(p.trans>1.3&logfold<(-1.0))],labels=as.character(data$genesymbol[(p.trans>1.3&logfold<(-1.0))]), pos=3, cex=0.40)
abline(h=1.3)
# abline(v=-1.0)
# abline(v=1.0)
dev.off()

pdf(file=paste("bipseq_lncRNAKB_qval_distribution_",feature,"_",bregion,"_",group1,"_",group2,".pdf", sep=""))
hist(gene_res$qvalue, col="grey", border="white", xlab="qval", ylab="", main=paste("cpm:",feature,"_",bregion,"_",group1,"_",group2))
dev.off()
pdf(file=paste("bipseq_lncRNAKB_pval_distribution_",feature,"_",bregion,"_",group1,"_",group2,".pdf", sep=""))
hist(gene_res$pvalue, col="grey", border="white", xlab="pval", ylab="", main=paste("cpm:",feature,"_",bregion,"_",group1,"_",group2))
dev.off()

gene_res = read_csv(file="/data/NHLBI_BCB/Fayaz/34-lncRNAKB-BP-JHU/02-differential-expression-analysis/00-gene-level/00-sacc/bipseq_lncRNAKB_gene_sacc_Bipolar_Control.csv")
# /data/NHLBI_BCB/Fayaz/34-lncRNAKB-BP-JHU/02-differential-expression-analysis/00-gene-level/00-amygdala/bipseq_lncRNAKB_gene_amy_Bipolar_Control.csv
# /data/NHLBI_BCB/Fayaz/34-lncRNAKB-BP-JHU/02-differential-expression-analysis/00-gene-level/00-sacc/bipseq_lncRNAKB_gene_sacc_Bipolar_Control.csv

filter(gene_res, log2FC < -0.50 | log2FC > 0.50, qvalue < 0.05)


#######################################################################################################################
#
#
#
#######################################################################################################################
gene_res_amy = read_csv(file="/data/NHLBI_BCB/Fayaz/34-lncRNAKB-BP-JHU/02-differential-expression-analysis/00-gene-level-redo-07222020/00-amygdala/bipseq_lncRNAKB_gene_amy_Bipolar_Control.csv")
gene_res_amy

gene_res_amy_sig = filter(gene_res_amy, qvalue < 0.05, gene_type=="lncRNA")
gene_res_amy_sig = filter(gene_res_amy, qvalue < 0.05)

# gene_res_amy_sig = filter(gene_res_amy, log2FC < -0.50 | log2FC > 0.50, qvalue < 0.05)
gene_res_amy_sig

gene_res_sacc = read_csv(file="/data/NHLBI_BCB/Fayaz/34-lncRNAKB-BP-JHU/02-differential-expression-analysis/00-gene-level-redo-07222020/00-sACC/bipseq_lncRNAKB_gene_sacc_Bipolar_Control.csv")
gene_res_sacc

gene_res_sacc_sig = filter(gene_res_sacc, qvalue < 0.05, gene_type=="lncRNA")
gene_res_sacc_sig = filter(gene_res_sacc, qvalue < 0.05)
# gene_res_sacc_sig = filter(gene_res_sacc, log2FC < -0.50 | log2FC > 0.50, qvalue < 0.05)
gene_res_sacc_sig

gene_res_sig_common_amy_sacc = gene_res_amy_sig %>% inner_join(gene_res_sacc_sig, by = c("gene" = "gene"))
gene_res_sig_common_amy_sacc

gene_res_amy_uniq = anti_join(gene_res_amy_sig, gene_res_sacc_sig, by = c("gene" = "gene"))
gene_res_amy_uniq

gene_res_sacc_uniq = anti_join(gene_res_sacc_sig, gene_res_amy_sig, by = c("gene" = "gene"))
gene_res_sacc_uniq

#######################################################################################################################


###############################################################################################################################
#
# Plot Expression by bin vs. significance -log10(pvalue)
#
###############################################################################################################################
dge = DGEList(counts = counts_data)
dge = calcNormFactors(dge, method="TMM")

counts_data_cpm = cpm(dge, normalized.lib.sizes=TRUE)
counts_data_cpm = as_tibble(counts_data_cpm, rownames="gene")
counts_data_cpm
df = counts_data_cpm %>% select(!("gene")) %>% mutate(mean_exp_tissues = rowMeans(.))
summary(df$mean_exp_tissues)
df
df = df %>% add_column(counts_data_cpm[,'gene'], .before = 1) %>% select(gene,mean_exp_tissues)
df 

gene_res_p0.05 = filter(gene_res, pvalue < 0.05) %>% select(gene,pvalue,gene_type)
gene_res_p0.05 = as_tibble(data.frame(gene=gene_res_p0.05$gene, pvalue=gene_res_p0.05$pvalue, log10pvalue=-1*log10(gene_res_p0.05$pvalue), gene_type=gene_res_p0.05$gene_type))
gene_res_p0.05

df1 = gene_res_p0.05 %>% inner_join(df, by = "gene")
df1

df = df1
df
# df = df %>% filter(gene_type == "lncRNA")
df = df %>% filter(gene_type == "protein_coding")
df

# df = as_tibble(mean_expression_across_all_tissue_tau_pem)

# ggplot(data = df, mapping = aes(x=mean_exp_tissues)) + 
#  geom_histogram(aes(y=..density..),fill="bisque",color="white",alpha=0.7,binwidth=0.05) + 
#  geom_density() +
#  geom_rug() +
#  labs(x='mean expression per tissue') +
#  theme_minimal()

# mutate(dat, dist = ifelse(speed == 4, dist * 100, dist)

# set up cut-off values 
breaks = c(seq(from=0, to=0.1, by = 0.01), seq(from=0.2, to=0.4, by=0.1), seq(from=0.5, to=15, by=0.50), seq(from=20, to=1000, by=100), seq(from=1000, to=12000, by=1000), max(df$mean_exp_tissues))
# specify interval/bin labels
tags <- c("[0.0-0.01)","[0.01-0.02)","[0.02-0.03)","[0.03-0.04)","[0.04-0.05)","[0.05-0.06)","[0.06-0.07)","[0.07-0.08)","[0.08-0.09)","[0.09-0.1)","[0.1-0.2)","[0.2-0.3)","[0.3-0.4)","[0.4-0.5)","[0.5-1.0)","[1.0-1.5)","[1.5-2.0)","[2.0-2.5)","[2.5-3.0)","[3.0-3.5)","[3.5-4.0)","[4.0-4.5)","[4.5-5.0)","[5.0-5.5)","[5.5-6.0)","[6.0-6.5)", "[6.5-7.0)","[7.0-7.5)","[7.5-8.0)","[8.0-8.5)","[8.5-9.0)","[9.0-9.5)","[9.5-10.0)","[10.0-10.5)","[10.5-11.0)","[11.0-11.5)","[11.5-12.0)","[12.0-12.5)","[12.5- 13.0)","[13.0-13.5)","[13.5-14.0)","[14.0-14.5)","[14.5-15.0)","[15.0-20.0)","[20.0-120.0)", "[120.0-220.0)", "[220.0-320.0)", "[320.0-420.0)", "[420.0-520.0)", "[520.0-620.0)", "[620.0-720.0)", "[720.0-820.0)", "[820.0-920.0)","[920.0-1000.0)","[1000-2000)","[2000-3000)","[3000-4000)","[4000-5000)","[5000-6000)","[6000-7000)","[7000-8000)","[8000-9000)","[9000-10000)","[10000-11000)","[11000-12000)","[>12000]")
# bucketing values into bins
group_tags <- cut(df$mean_exp_tissues, 
                  breaks=breaks, 
                  include.lowest=TRUE, 
                  right=FALSE, 
                  labels=tags)
# inspect bins
summary(group_tags)

vgroup <- df %>% 
  mutate(tag = case_when(
mean_exp_tissues < 0.01 ~ tags[1],
mean_exp_tissues >= 0.01 & mean_exp_tissues < 0.02 ~ tags[2],
mean_exp_tissues >= 0.02 & mean_exp_tissues < 0.03 ~ tags[3],
mean_exp_tissues >= 0.03 & mean_exp_tissues < 0.04 ~ tags[4],
mean_exp_tissues >= 0.04 & mean_exp_tissues < 0.05 ~ tags[5],
mean_exp_tissues >= 0.05 & mean_exp_tissues < 0.06 ~ tags[6],
mean_exp_tissues >= 0.06 & mean_exp_tissues < 0.07 ~ tags[7],
mean_exp_tissues >= 0.07 & mean_exp_tissues < 0.08 ~ tags[8],
mean_exp_tissues >= 0.08 & mean_exp_tissues < 0.09 ~ tags[9],
mean_exp_tissues >= 0.09 & mean_exp_tissues < 0.1 ~ tags[10],
mean_exp_tissues >= 0.1 & mean_exp_tissues < 0.2 ~ tags[11],
mean_exp_tissues >= 0.2 & mean_exp_tissues < 0.3 ~ tags[12],
mean_exp_tissues >= 0.3 & mean_exp_tissues < 0.4 ~ tags[13],
mean_exp_tissues >= 0.4 & mean_exp_tissues < 0.5 ~ tags[14],
mean_exp_tissues >= 0.5 & mean_exp_tissues < 1.0 ~ tags[15],
mean_exp_tissues >= 1.0 & mean_exp_tissues < 1.5 ~ tags[16],
mean_exp_tissues >= 1.5 & mean_exp_tissues < 2.0 ~ tags[17],
mean_exp_tissues >= 2.0 & mean_exp_tissues < 2.5 ~ tags[18],
mean_exp_tissues >= 2.5 & mean_exp_tissues < 3.0 ~ tags[19],
mean_exp_tissues >= 3.0 & mean_exp_tissues < 3.5 ~ tags[20],
mean_exp_tissues >= 3.5 & mean_exp_tissues < 4.0 ~ tags[21],
mean_exp_tissues >= 4.0 & mean_exp_tissues < 4.5 ~ tags[22],
mean_exp_tissues >= 4.5 & mean_exp_tissues < 5.0 ~ tags[23],
mean_exp_tissues >= 5.0 & mean_exp_tissues < 5.5 ~ tags[24],
mean_exp_tissues >= 5.5 & mean_exp_tissues < 6.0 ~ tags[25],
mean_exp_tissues >= 6.0 & mean_exp_tissues < 6.5 ~ tags[26],
mean_exp_tissues >= 6.5 & mean_exp_tissues < 7.0 ~ tags[27],
mean_exp_tissues >= 7.0 & mean_exp_tissues < 7.5 ~ tags[28],
mean_exp_tissues >= 7.5 & mean_exp_tissues < 8.0 ~ tags[29],
mean_exp_tissues >= 8.0 & mean_exp_tissues < 8.5 ~ tags[30],
mean_exp_tissues >= 8.5 & mean_exp_tissues < 9.0 ~ tags[31],
mean_exp_tissues >= 9.0 & mean_exp_tissues < 9.5 ~ tags[32],
mean_exp_tissues >= 9.5 & mean_exp_tissues < 10.0 ~ tags[33],
mean_exp_tissues >= 10.0 & mean_exp_tissues < 10.5 ~ tags[34],
mean_exp_tissues >= 10.5 & mean_exp_tissues < 11.0 ~ tags[35],
mean_exp_tissues >= 11.0 & mean_exp_tissues < 11.5 ~ tags[36],
mean_exp_tissues >= 11.5 & mean_exp_tissues < 12.0 ~ tags[37],
mean_exp_tissues >= 12.0 & mean_exp_tissues < 12.5 ~ tags[38],
mean_exp_tissues >= 12.5 & mean_exp_tissues < 13.0 ~ tags[39],
mean_exp_tissues >= 13.0 & mean_exp_tissues < 13.5 ~ tags[40],
mean_exp_tissues >= 13.5 & mean_exp_tissues < 14.0 ~ tags[41],
mean_exp_tissues >= 14.0 & mean_exp_tissues < 14.5 ~ tags[42],
mean_exp_tissues >= 14.5 & mean_exp_tissues < 15.0 ~ tags[43],
mean_exp_tissues >= 15.0 & mean_exp_tissues < 20.0 ~ tags[44],
mean_exp_tissues >= 20.0 & mean_exp_tissues < 120.0 ~ tags[45],
mean_exp_tissues >= 120.0 & mean_exp_tissues < 220.0 ~ tags[46],
mean_exp_tissues >= 220.0 & mean_exp_tissues < 320.0 ~ tags[47],
mean_exp_tissues >= 320.0 & mean_exp_tissues < 420.0 ~ tags[48],
mean_exp_tissues >= 420.0 & mean_exp_tissues < 520.0 ~ tags[49],
mean_exp_tissues >= 520.0 & mean_exp_tissues < 620.0 ~ tags[50],
mean_exp_tissues >= 620.0 & mean_exp_tissues < 720.0 ~ tags[51],
mean_exp_tissues >= 720.0 & mean_exp_tissues < 820.0 ~ tags[52],
mean_exp_tissues >= 820.0 & mean_exp_tissues < 920.0 ~ tags[53],
mean_exp_tissues >= 920.0 & mean_exp_tissues < 1000.0 ~ tags[54],
mean_exp_tissues >= 1000 & mean_exp_tissues < 2000 ~ tags[55],
mean_exp_tissues >= 2000 & mean_exp_tissues < 3000 ~ tags[56],
mean_exp_tissues >= 3000 & mean_exp_tissues < 4000 ~ tags[57],
mean_exp_tissues >= 4000 & mean_exp_tissues < 5000 ~ tags[58],
mean_exp_tissues >= 5000 & mean_exp_tissues < 6000 ~ tags[59],
mean_exp_tissues >= 6000 & mean_exp_tissues < 7000 ~ tags[60],
mean_exp_tissues >= 7000 & mean_exp_tissues < 8000 ~ tags[61],
mean_exp_tissues >= 8000 & mean_exp_tissues < 9000 ~ tags[62],
mean_exp_tissues >= 9000 & mean_exp_tissues < 10000 ~ tags[63],
mean_exp_tissues >= 10000 & mean_exp_tissues < 11000 ~ tags[64],
mean_exp_tissues >= 11000 & mean_exp_tissues < 12000 ~ tags[65],
    ))
summary(vgroup)

vgroup$tag <- factor(vgroup$tag,
                       levels = tags,
                       ordered = FALSE)
summary(vgroup$tag)

ggplot(data = vgroup, mapping = aes(x=tag,y=log10pvalue)) + 
  geom_jitter(aes(color='blue'),alpha=0.2) +
  geom_boxplot(fill="bisque",color="black",alpha=0.3) + 
  labs(title="Plot of -log10pvalue by mean CPM expression across all samples-lncRNA", x='mean CPM expression across all samples') +
  guides(color=FALSE) +
  theme_minimal() +
  theme(axis.text.x = element_text(face="bold",size=5, angle=90)) +
  scale_y_continuous(breaks=seq(from=0,to=10, by=1))

dev.off()

#######################################################################################################################################
#
#
#
# WGCNA - CEMiTools
#
#
#######################################################################################################################################

library("CEMiTool")
library("AnnotationDbi")
library("org.Hs.eg.db")

# setwd("/data/NHLBI_BCB/Fayaz/34-lncRNAKB-BP-JHU/03-WGCNA-analysis-CEMiTools/00-amygdala")
# counts_data = read_tsv(file="/data/NHLBI_BCB/Fayaz/34-lncRNAKB-BP-JHU/01-lncRNAKB-gene-level-counts/amygdala.gene.featurecount.txt")

setwd("/data/NHLBI_BCB/Fayaz/34-lncRNAKB-BP-JHU/03-WGCNA-analysis-CEMiTools/00-sACC")
counts_data = read_tsv(file="/data/NHLBI_BCB/Fayaz/34-lncRNAKB-BP-JHU/01-lncRNAKB-gene-level-counts/sacc.gene.featurecount.txt")

counts_data

# phenodata = read_tsv(file="/data/NHLBI_BCB/Fayaz/34-lncRNAKB-BP-JHU/00-phenotype/amygdala_subjects_pheno.txt")
phenodata = read_tsv(file="/data/NHLBI_BCB/Fayaz/34-lncRNAKB-BP-JHU/00-phenotype/sacc_subjects_pheno.txt")

phenodata

group1 = "Bipolar"
group2 = "Control"
feature = "gene"
bregion = "sacc"
# bregion = "amy"

phenodata = filter(phenodata, dx == group1 | dx == group2)
phenodata
counts_data = dplyr::select(counts_data, c("gene",phenodata$SampleID))
counts_data

annot = read_csv(file="/data/NHLBI_BCB/Fayaz/34-lncRNAKB-BP-JHU/02-differential-expression-analysis/00-gene-level/annot_format_v7.csv")
annot

counts_data_unfiltered = counts_data

###################################################################################################################### 
#
# Calculate TPM 
# total lncRNAs vs diff. expressed lncRNAs Within vs. outside
# Test enrichment of differentially lncRNAs within modules
# Axis1 - Differentially expressed, not expressed
# Axis2 - in module, out module
#
# Wilcoxon rank sum test to test the enrichment of eigengenes with case-control status
#
######################################################################################################################
counts_data = counts_data %>% inner_join(annot, by = c("gene" = "lnckbID"))
counts_data
gene_lenghts_kb = (counts_data$end - counts_data$start)/1000
counts_data = data.frame(counts_data[,2:269], row.names=counts_data$gene) #change column number to (# of samples)+1
counts_data_scaled_by_gene_length = sapply(1:nrow(counts_data), function(i) counts_data[i,]/gene_lenghts_kb[i])
counts_data_scaled_by_gene_length_df = data.frame(matrix(unlist(counts_data_scaled_by_gene_length), nrow=dim(counts_data)[2], byrow=F), stringsAsFactors=FALSE)
counts_data_scaled_by_gene_length_df = t(counts_data_scaled_by_gene_length_df)
dim(counts_data_scaled_by_gene_length_df)
lib_sizes = colSums(counts_data_scaled_by_gene_length_df)
scaling_factor = lib_sizes/1000000
counts_data_scaled_by_gene_length_and_scaling_factor_df = sapply(1:ncol(counts_data_scaled_by_gene_length_df), function(i) counts_data_scaled_by_gene_length_df[,i]/scaling_factor[i])
counts_data_scaled_by_gene_length_and_scaling_factor_df = data.frame(counts_data_scaled_by_gene_length_and_scaling_factor_df, row.names=rownames(counts_data))
colnames(counts_data_scaled_by_gene_length_and_scaling_factor_df) = colnames(counts_data)
dim(counts_data_scaled_by_gene_length_and_scaling_factor_df)
tpm = counts_data_scaled_by_gene_length_and_scaling_factor_df
tpm = as_tibble(tpm, rownames="gene")
write.table(tpm, file="/data/NHLBI_BCB/Fayaz/34-lncRNAKB-BP-JHU/01-lncRNAKB-gene-level-counts/sacc.gene.featurecount.TPM.txt")
######################################################################################################################

## filtering-level1 by TPM
#Keep genes with TPM>0.50 in at least 20% of samples

# tpm = read.table(file="/data/NHLBI_BCB/Fayaz/34-lncRNAKB-BP-JHU/01-lncRNAKB-gene-level-counts/amygdala.gene.featurecount.TPM.txt", header=T)
tpm = read.table(file="/data/NHLBI_BCB/Fayaz/34-lncRNAKB-BP-JHU/01-lncRNAKB-gene-level-counts/sacc.gene.featurecount.TPM.txt", header=T)
tpm[1:10,1:10]
tpm = as_tibble(tpm)
tpm

tpm_filtered = tpm %>% filter(rowSums(tpm>0.50) >= as.integer(0.20*(dim(tpm)[2]-1)))
tpm_filtered

## filtering-level2 by raw counts
#Keep lncRNAs with least 1 counts/reads in 90% of replicates/samples
counts_data = counts_data %>% semi_join(tpm_filtered, by="gene")
counts_data

counts_data_lncRNA = counts_data %>% semi_join(filter(annot, gene_type == "lncRNA"), by = c("gene" = "lnckbID"))
counts_data_lncRNA
counts_data_lncRNA = counts_data_lncRNA %>% filter(rowSums(counts_data_lncRNA>1) >= as.integer(0.90*(dim(counts_data_lncRNA)[2]-1)))
counts_data_lncRNA

#Keep protein coding genes with least 5 counts/reads in 90% of replicates/samples
counts_data_pc = counts_data %>% semi_join(filter(annot, gene_type == "protein_coding"), by = c("gene" = "lnckbID"))
counts_data_pc
counts_data_pc = counts_data_pc %>% filter(rowSums(counts_data_pc>5) >= as.integer(0.90*(dim(counts_data_pc)[2]-1)))
counts_data_pc

counts_data = bind_rows(counts_data_lncRNA, counts_data_pc)
counts_data

df = data.frame(counts_data)
df1 = data.frame(df[,2:dim(df)[2]], row.names=df$gene)

counts_data_cpm = as_tibble(cpm(df1), rownames="gene")
counts_data_cpm

## filtering-level3 by CPM
#Keep lncRNAs with least 1 counts/reads in 90% of replicates/samples
counts_data_cpm_lncRNA = counts_data_cpm %>% semi_join(filter(annot, gene_type == "lncRNA"), by = c("gene" = "lnckbID"))
counts_data_cpm_lncRNA
counts_data_cpm_lncRNA = counts_data_cpm_lncRNA %>% filter(rowSums(counts_data_cpm_lncRNA>1) >= as.integer(0.90*(dim(counts_data_cpm_lncRNA)[2]-1)))
counts_data_cpm_lncRNA

#Keep protein coding genes with least 5 counts/reads in 90% of replicates/samples
counts_data_cpm_pc = counts_data_cpm %>% semi_join(filter(annot, gene_type == "protein_coding"), by = c("gene" = "lnckbID"))
counts_data_cpm_pc
counts_data_cpm_pc = counts_data_cpm_pc %>% filter(rowSums(counts_data_cpm_pc>2) >= as.integer(0.90*(dim(counts_data_cpm_pc)[2]-1)))
counts_data_cpm_pc

counts_data_cpm_filt = bind_rows(counts_data_cpm_lncRNA, counts_data_cpm_pc)
counts_data_cpm_filt

# annotv2 = read_csv(file="/data/NHLBI_BCB/Fayaz/34-lncRNAKB-BP-JHU/03-WGCNA-analysis-CEMiTools/annot_format_v3.csv")

counts_data_cpm_filt = counts_data_cpm_filt %>% inner_join(annot, by = c("gene" = "lnckbID"))
counts_data_cpm_filt

counts_data_cpm_filt_expr = dplyr::select(counts_data_cpm_filt, -c("gene", "lnckb", "chr","start","end","source","gene_type","info","chr_st_en","Gene_stable_ID","Gene_name","Gene_type","NCBI_gene_formerly_Entrezgene_ID","HGNC_ID"))
counts_data_cpm_filt_expr

df = data.frame(counts_data_cpm_filt_expr)
df1 = data.frame(df[,1:dim(df)[2]-1], row.names=make.unique(df$idents))

counts_data_cpm_filt_df = df1
dim(counts_data_cpm_filt_df)
counts_data_cpm_filt_df[1:10,1:10]

######################################################################################################################
#
#
#
######################################################################################################################

lncRNA = counts_data_cpm_filt %>% filter(gene_type == "lncRNA")
lncRNA = distinct(lncRNA %>% dplyr::select(idents))
lncRNA

pc = counts_data_cpm_filt %>% filter(gene_type == "protein_coding")
pc = distinct(pc %>% dplyr::select(idents))
pc

gene_res = read_csv(file="/data/NHLBI_BCB/Fayaz/34-lncRNAKB-BP-JHU/02-differential-expression-analysis/00-gene-level-redo-07222020/00-sACC/bipseq_lncRNAKB_gene_sacc_Bipolar_Control.csv")
gene_res

# gene_res = read_csv(file="/data/NHLBI_BCB/Fayaz/34-lncRNAKB-BP-JHU/02-differential-expression-analysis/00-gene-level-redo-07222020/00-amygdala/bipseq_lncRNAKB_gene_amy_Bipolar_Control.csv")
# gene_res

gene_res = gene_res %>% semi_join(counts_data_cpm_filt, by = "gene")
gene_res

gene_res = gene_res %>% inner_join(annot, by = c("gene" = "lnckbID"))
gene_res


####################################################################
#
# pvalue < 0.05
#
####################################################################

lncRNA_p0.05 = filter(gene_res, pvalue < 0.05) %>% filter(gene_type.y == "lncRNA") %>% dplyr::select(idents)
lncRNA_p0.05

lncRNA_pgreater0.05 = filter(gene_res, pvalue > 0.05) %>% filter(gene_type.y == "lncRNA") %>% dplyr::select(idents)
lncRNA_pgreater0.05

lncRNA_p0.05_up = filter(gene_res, pvalue < 0.05) %>% filter(gene_type.y == "lncRNA") %>% filter(log2FC > 0) %>% dplyr::select(idents)
lncRNA_p0.05_up

lncRNA_p0.05_down = filter(gene_res, pvalue < 0.05) %>% filter(gene_type.y == "lncRNA") %>% filter(log2FC < 0) %>% dplyr::select(idents)
lncRNA_p0.05_down

pc_p0.05 = filter(gene_res, pvalue < 0.05) %>% filter(gene_type.y == "protein_coding") %>% dplyr::select(idents)
pc_p0.05

pc_pgreater0.05 = filter(gene_res, pvalue > 0.05) %>% filter(gene_type.y == "protein_coding") %>% dplyr::select(idents)
pc_pgreater0.05

pc_p0.05_up = filter(gene_res, pvalue < 0.05) %>% filter(gene_type.y == "protein_coding") %>% filter(log2FC > 0) %>% dplyr::select(idents)
pc_p0.05_up

pc_p0.05_down = filter(gene_res, pvalue < 0.05) %>% filter(gene_type.y == "protein_coding") %>% filter(log2FC < 0) %>% dplyr::select(idents)
pc_p0.05_down


####################################################################
#
# qvalue < 0.05
#
####################################################################

lncRNA_q0.05 = filter(gene_res, qvalue < 0.05) %>% filter(gene_type.y == "lncRNA") %>% dplyr::select(idents)
lncRNA_q0.05

lncRNA_qgreater0.05 = filter(gene_res, qvalue > 0.05) %>% filter(gene_type.y == "lncRNA") %>% dplyr::select(idents)
lncRNA_qgreater0.05

lncRNA_q0.05_up = filter(gene_res, qvalue < 0.05) %>% filter(gene_type.y == "lncRNA") %>% filter(log2FC > 0) %>% dplyr::select(idents)
lncRNA_q0.05_up

lncRNA_q0.05_down = filter(gene_res, qvalue < 0.05) %>% filter(gene_type.y == "lncRNA") %>% filter(log2FC < 0) %>% dplyr::select(idents)
lncRNA_q0.05_down

pc_q0.05 = filter(gene_res, qvalue < 0.05) %>% filter(gene_type.y == "protein_coding") %>% dplyr::select(idents)
pc_q0.05

pc_q0.05_up = filter(gene_res, qvalue < 0.05) %>% filter(gene_type.y == "protein_coding") %>% filter(log2FC > 0) %>% dplyr::select(idents)
pc_q0.05_up

pc_q0.05_down = filter(gene_res, qvalue < 0.05) %>% filter(gene_type.y == "protein_coding") %>% filter(log2FC < 0) %>% dplyr::select(idents)
pc_q0.05_down

#1. lncRNA
#2. lncRNA_p0.05
#3. lncRNA_pgreater0.05
#4. lncRNA_upreg
#5. lncRNA_downreg
#6. protein_coding
#7. protein_coding_p0.05
#8. protein_coding_upreg
#9. protein_coding_downreg
#10. lncRNA_q0.05
#11. lncRNA_qgreater0.05
#12. lncRNA_upreg
#13. lncRNA_downreg
#14. protein_coding_q0.05
#15. protein_coding_upreg
#16. protein_coding_downreg

######################################################################################################################

df = data.frame(phenodata)
colnames(df) = c("SampleName","Class")
phenodata = df
phenodata[1:10,]

gmt.cemitool = read_gmt("/data/NHLBI_BCB/Larochelle_Andre/05_Notch_Hypoxia_RNAseq/08-WGCNA-Cemitools/c5.all.v7.0.symbols.gmt")

cemres = cemitool(expr=counts_data_cpm_filt_df, 
	annot=phenodata, 
	gmt=gmt.cemitool, 
	verbose=T, 
	cor_method="pearson", 
	apply_vst=F, 
	filter=F, 
	ora_pval=0.05, 
	network_type="unsigned", 
	diss_thresh=0.90)

cem <- mod_gsea(cemres)
cem <- plot_gsea(cemres)
show_plot(cem, "gsea")
dev.off()

eigengene = mod_summary(cemres, method="eigengene")
eigengenedf = data.frame(eigengene[,2:dim(eigengene)[2]], row.names=eigengene$modules)

eigengenedf_tr = t(eigengenedf)
eig = as_tibble(eigengenedf_tr, rownames = "subjectID")
# eig = eig %>% dplyr::select(-c("Not.Correlated"))
eig

###
# What to do if there are binary traits? Examples of binary traits:Sex (male/female) Disease status (healthy/diseased)

# Pearsonand Spearmancorrelation (function cor): will work 
# Biweightmidcorrelation(function bicor) needs modification- restrict number of values that will be treated as outliers - turn off robust treatment for the binary variable
# bicor(moduleEigengenes, datTraits, use=“p",robustY=FALSE, maxPOutliers=0.1)

phenodata = as_tibble(phenodata)

phenodata = phenodata %>%
  mutate(
    dx = case_when(
      Class == "Bipolar" ~ 1,
      Class == "Control" ~ 0,
      TRUE ~ NA_real_
    )
  )
phenodata

x <- correlate(x = eig, y = phenodata %>% dplyr::select(c("dx"))) %>%    # Create correlation data frame (cor_df)
     rearrange() %>%  # rearrange by correlations
     shave() # Shave off the upper triangle for a clean result
#> 

fashion(x)
rplot(x)

correlations = data.frame(cor(x = data.frame(eig), y = data.frame(phenodata %>% dplyr::select(c("dx")))))
correlations = as_tibble(data.frame(modules = rownames(correlations), corr = correlations$dx))
correlations

correlations = correlations %>%
  mutate(
    corrd = case_when(
      corr > 0 ~ "positive",
      corr < 0 ~ "negative",
      TRUE ~ "other"
    )
  )
correlations

# correlations <- correlations %>% group_by(corrd)
# correlations <- correlations %>% arrange(desc(corr))
# correlations

p<-ggplot(data=correlations, aes(x=reorder(modules,corr), y=corr, fill=corrd)) +
  geom_bar(stat="identity") +
  xlab("Modules") + ylab("Pearsons Correlation with Phenotype (Bipolar/Control)") +
  theme(axis.text.x = element_text(size = 7.0, angle = 90))
p
dev.off()   


cor(x = eigengenedf_tr, y = y)

symbols(x = 1.5, 
	y = 1.5, 
	circles=0.30, 
	inches = FALSE, 
	add = FALSE,
    fg = par("col"), 
    bg = NA,
    xlab = NULL, 
    ylab = NULL, 
    main = NULL,
    xlim = c(1,5), 
    ylim = c(1,5))

text(x=1.5, 
	y = 1.5, 
	labels = c("1178/1305"), 
	adj = NULL,
    pos = NULL, 
    offset = 0.5, 
    vfont = NULL,
    cex = 1, 
    col = NULL, 
    font = NULL)

symbols(x=1.5, 
	y = 2.5, 
	circles=0.30, 
	inches = FALSE, 
	add = TRUE,
    fg = par("col"), 
    bg = NA,
    xlab = NULL, 
    ylab = NULL, 
    main = NULL,
    xlim = c(1,5), 
    ylim = c(1,5))

text(x=1.5, 
	y = 2.5, 
	labels = c("872/664"), 
	adj = NULL,
    pos = NULL, 
    offset = 0.5, 
    vfont = NULL,
    cex = 1, 
    col = NULL, 
    font = NULL)


# Plot soft-threshold beta and r2
cemres = plot_beta_r2(cemres)
show_plot(cemres, "beta_r2")

# Plot scale-free model fit as a function of the soft-thresholding beta parameter choice
cemres <- plot_mean_k(cemres)
# Check resulting plot
show_plot(cemres, "mean_k")
dev.off()

module_names = mod_names(cemres)
hubs_top20 = get_hubs(cem=cemres, n=20) # Returns n genes in each module with high connectivity, Default: "adjacency"
hubs = get_hubs(cem=cemres, n="all", method = "kME")
moduleinfo = matrix(nrow=length(module_names), ncol=39, byrow=T)
for(i in 1:length(module_names)) {
	mgenes = module_genes(cemres, module=module_names[i])
	mgenesvec = as_tibble(mgenes$genes)
	mgenesvec

	eigen = eig %>% dplyr::select(c("subjectID",module_names[i]))
	eigen = eigen %>% inner_join(phenodata, by = c("subjectID" = "SampleID"))
	eigen
	res = wilcox.test(eval(parse(text=module_names[i])) ~ dx, data = eigen, exact = FALSE)
	u = res$p.value

	kME = hubs %>% pluck(module_names[i]) %>% as_tibble(rownames="Genes")

	a = lncRNA %>% semi_join(mgenesvec, by = c("idents" = "value"))
	a

	e = pc %>% semi_join(mgenesvec, by = c("idents" = "value"))
	e

	####################################################################
	#
	# pvalue < 0.05
	#
	####################################################################

	b = lncRNA_p0.05 %>% semi_join(mgenesvec, by = c("idents" = "value"))
	b
	b_out = dim(lncRNA_p0.05)[1] - dim(b)[1]

	k = lncRNA_pgreater0.05 %>% semi_join(mgenesvec, by = c("idents" = "value"))
	k
	k_out = dim(lncRNA_pgreater0.05)[1] - dim(k)[1]

	c = lncRNA_p0.05_up %>% semi_join(mgenesvec, by = c("idents" = "value"))
	c

	d = lncRNA_p0.05_down %>% semi_join(mgenesvec, by = c("idents" = "value"))
	d

	f = pc_p0.05 %>% semi_join(mgenesvec, by = c("idents" = "value"))
	f
	f_out = dim(pc_p0.05)[1] - dim(f)[1]

	w = pc_pgreater0.05 %>% semi_join(mgenesvec, by = c("idents" = "value"))
	w
	w_out = dim(pc_pgreater0.05)[1] - dim(w)[1]

	g = pc_p0.05_up %>% semi_join(mgenesvec, by = c("idents" = "value"))
	g

	h = pc_p0.05_down %>% semi_join(mgenesvec, by = c("idents" = "value"))
	h

	####################################################################
	#
	# qvalue < 0.05
	#
	####################################################################

	l = lncRNA_q0.05 %>% semi_join(mgenesvec, by = c("idents" = "value"))
	l
	l_out = dim(lncRNA_q0.05)[1] - dim(l)[1]

	m = lncRNA_qgreater0.05 %>% semi_join(mgenesvec, by = c("idents" = "value"))
	m
	m_out = dim(lncRNA_pgreater0.05)[1] - dim(m)[1]

	n = lncRNA_q0.05_up %>% semi_join(mgenesvec, by = c("idents" = "value"))
	n

	o = lncRNA_q0.05_down %>% semi_join(mgenesvec, by = c("idents" = "value"))
	o

	q = pc_q0.05 %>% semi_join(mgenesvec, by = c("idents" = "value"))
	q

	r = pc_q0.05_up %>% semi_join(mgenesvec, by = c("idents" = "value"))
	r

	s = pc_q0.05_down %>% semi_join(mgenesvec, by = c("idents" = "value"))
	s

	t = fisher.test(matrix(c(dim(b)[1], b_out, dim(k)[1], k_out),2,2), alternative='greater')$p.value

	v = fisher.test(matrix(c(dim(f)[1], f_out, dim(w)[1], w_out),2,2), alternative='greater')$p.value

	x = fisher.test(matrix(c(dim(f)[1]+dim(b)[1], f_out+b_out, dim(w)[1]+dim(k)[1], w_out+k_out),2,2), alternative='greater')$p.value

	####################################################################
	#
	# Mean module memberships (KME) of gene sets/list of genes
	#
	####################################################################

	kME_a = a %>% inner_join(kME, by = c("idents" = "Genes"))
	mu_kME_a = mean(kME_a$value)
	mu_kME_a

	kME_b = b %>% inner_join(kME, by = c("idents" = "Genes"))
	mu_kME_b = mean(kME_b$value)
	mu_kME_b

	kME_c = c %>% inner_join(kME, by = c("idents" = "Genes"))
	mu_kME_c = mean(kME_c$value)
	mu_kME_c

	kME_d = d %>% inner_join(kME, by = c("idents" = "Genes"))
	mu_kME_d = mean(kME_d$value)
	mu_kME_d

	kME_e = e %>% inner_join(kME, by = c("idents" = "Genes"))
	mu_kME_e = mean(kME_e$value)
	mu_kME_e

	kME_f = f %>% inner_join(kME, by = c("idents" = "Genes"))
	mu_kME_f = mean(kME_f$value)
	mu_kME_f

	kME_g = g %>% inner_join(kME, by = c("idents" = "Genes"))
	mu_kME_g = mean(kME_g$value)
	mu_kME_g

	kME_h = h %>% inner_join(kME, by = c("idents" = "Genes"))
	mu_kME_h = mean(kME_h$value)
	mu_kME_h

	moduleinfo[i,1] = bregion # brain region name
	moduleinfo[i,2] = module_names[i] # module names
	moduleinfo[i,3] = dim(mgenesvec)[1] # number of genes in the module
	moduleinfo[i,4] = dim(a)[1] # lncRNA in the module
	moduleinfo[i,5] = dim(b)[1] # lncRNA_p0.05 in the module
	moduleinfo[i,6] = dim(c)[1] # lncRNA_p0.05_up in the module
	moduleinfo[i,7] = dim(d)[1] # lncRNA_p0.05_down in the module
	moduleinfo[i,8] = b_out # lncRNA_p0.05 outside the module
	moduleinfo[i,9] = dim(k)[1] # lncRNA_pgreater0.05 in the module
	moduleinfo[i,10] = k_out # lncRNA_pgreater0.05 outside the module
	moduleinfo[i,11] = t # fisher test for overrepresentation of lncRNA_p0.05 in the module
	moduleinfo[i,12] = v # fisher test for overrepresentation of pc_p0.05 in the module
	moduleinfo[i,13] = x # fisher test for overrepresentation of pc_p0.05+lncRNA_p0.05 in the module
	moduleinfo[i,14] = u # wilcoxon rank sum test for testing the difference in median of eigen gene expression between cases and controls for each module 
	moduleinfo[i,15] = dim(e)[1] # pc in the module
	moduleinfo[i,16] = dim(f)[1] # pc_p0.05 in the module
	moduleinfo[i,17] = dim(g)[1] # pc_p0.05_up in the module
	moduleinfo[i,18] = dim(h)[1] # pc_p0.05_down in the module
	moduleinfo[i,19] = f_out # pc_p0.05 outside the module
	moduleinfo[i,20] = dim(w)[1] # pc_pgreater0.05 in the module
	moduleinfo[i,21] = w_out # pc_pgreater0.05 outside the module
	moduleinfo[i,22] = dim(l)[1] # lncRNA_q0.05 in the module
	moduleinfo[i,23] = dim(n)[1] # lncRNA_q0.05_up in the module
	moduleinfo[i,24] = dim(o)[1] # lncRNA_q0.05_down in the module
	moduleinfo[i,25] = l_out # lncRNA_q0.05 outside the module
	moduleinfo[i,26] = dim(m)[1] # lncRNA_qgreater0.05 in the module
	moduleinfo[i,27] = m_out # lncRNA_qgreater0.05 outside the module
	moduleinfo[i,28] = dim(q)[1] # pc_q0.05 in the module
	moduleinfo[i,29] = dim(r)[1] # pc_q0.05_up in the module
	moduleinfo[i,30] = dim(s)[1] # pc_q0.05_down in the module
	moduleinfo[i,31] = mu_kME_a # mean module membership kME lncRNA in the module
	moduleinfo[i,32] = mu_kME_b # mean module membership kME lncRNA_p0.05 in the module
	moduleinfo[i,33] = mu_kME_c # mean module membership kME lncRNA_p0.05_up in the module
	moduleinfo[i,34] = mu_kME_d # mean module membership kME lncRNA_p0.05_down in the module
	moduleinfo[i,35] = mu_kME_e # mean module membership kME pc in the module
	moduleinfo[i,36] = mu_kME_f # mean module membership kME pc_p0.05 in the module
	moduleinfo[i,37] = mu_kME_g # mean module membership kME pc_p0.05_up in the module
	moduleinfo[i,38] = mu_kME_h # mean module membership kME pc_p0.05_down in the module
	moduleinfo[i,39] = paste(names(hubs_top20[[module_names[i]]]), collapse="/") #top 20 most connected genes in each module by adjacency method
}

moduleinfo = as.data.frame(moduleinfo)
colnames(moduleinfo) = c("Brain_Region",
		"ModuleName",
		"number_genes_in_module",
		"lncRNA",
		"lncRNA_p0.05",
		"lncRNA_p0.05_up",
		"lncRNA_p0.05_down",
		"lncRNA_p0.05_out",
		"lncRNA_pgreater0.05",
		"lncRNA_pgreater0.05_out",
		"lncRNA_p0.05_fisher_test_pvalue",
		"pc_p0.05_fisher_test_pvalue",
		"lncRNA_pc_p0.05_fisher_test_pvalue",
		"wilcoxon_test_pvalue",
		"pc",
		"pc_p0.05",
		"pc_p0.05_up",
		"pc_p0.05_down",
		"pc_p0.05_out",
		"pc_pgreater0.05",
		"pc_pgreater0.05_out",
		"lncRNA_q0.05",
		"lncRNA_q0.05_up",
		"lncRNA_q0.05_down",
		"lncRNA_q0.05_out",
		"lncRNA_qgreater0.05",
		"lncRNA_qgreater0.05_out",
		"pc_q0.05",
		"pc_q0.05_up",
		"pc_q0.05_down",
		"mu_kME_lncRNA",
		"mu_kME_lncRNA_p0.05",
		"mu_kME_lncRNA_p0.05_up",
		"mu_kME_lncRNA_p0.05_down",
		"mu_kME_pc",
		"mu_kME_pc_p0.05",
		"mu_kME_pc_p0.05_up",
		"mu_kME_pc_p0.05_down",
		"top20_hub_genes")
write.csv(moduleinfo, file="Modules_Identified.csv", row.names=F)

orares = ora_data(cem=cemres) # Retrieve over representation analysis (ORA) results
moduleenrichmentinfo = data.frame(head(orares[orares$Module %in% module_names[1],], n=3))
for(i in 2:length(module_names)) {
	module_name = module_names[i]
	moduleenrichmentinfo = rbind(moduleenrichmentinfo, head(orares[orares$Module %in% module_name,], n=3))
}
colnames(moduleenrichmentinfo) = c("Module","ID","Description","GeneRatio","BgRatio","pvalue","p.adjust","qvalue","geneID","Count")
moduleenrichmentinfo = as.data.frame(moduleenrichmentinfo, colClasses=c("character","character","character","numeric","numeric","numeric","numeric","numeric","character","numeric"), stringsAsFactors=FALSE, make.names=T)
moduleenrichmentinfo_moduleinfo_merged = merge(x=moduleenrichmentinfo, y=moduleinfo, by.x="Module", by.y="ModuleName", all.x=T)
write.csv(moduleenrichmentinfo, file="Modules_Identified_Enrichment_top3_pathways.csv", row.names=F)
write.csv(moduleenrichmentinfo_moduleinfo_merged, file="Modules_Identified_Enrichment_top3_pathways_with_module_info.csv", row.names=F)
write.csv(orares, file="Modules_Identified_Enrichment_all_pathways.csv")

orares_moduleinfo_merged = merge(x=orares, y=moduleinfo, by.x="Module", by.y="ModuleName", all.x=T)
write.csv(orares_moduleinfo_merged, file="Modules_Identified_Enrichment_all_pathways_with_module_info.csv")
############################
moduleinfo = as_tibble(moduleinfo)
moduleinfo
moduleinfo_stats = dplyr::select(moduleinfo, c("ModuleName","lncRNA_p0.05_fisher_test_pvalue","pc_p0.05_fisher_test_pvalue","lncRNA_pc_p0.05_fisher_test_pvalue","wilcoxon_test_pvalue"))
moduleinfo_stats
moduleinfo_stats = pivot_longer(moduleinfo_stats, cols=!ModuleName, names_to="enrichment_test", values_to="pvalue") %>% mutate("-log10pvalue" = -log10(as.numeric(pvalue)))
moduleinfo_stats
moduleinfo_stats = filter(moduleinfo_stats, ModuleName != "Not.Correlated")
moduleinfo_stats
moduleinfo_stats = moduleinfo_stats %>% mutate("ModuleNameMod" = as.numeric(str_sub(ModuleName, 2)))
moduleinfo_stats

ggplot(moduleinfo_stats, aes(x = fct_reorder(ModuleName,ModuleNameMod, .desc=F), `-log10pvalue`)) + geom_point(aes(shape=enrichment_test, col=enrichment_test)) +
geom_hline(aes(yintercept = 1.3), color="black", linetype="dashed", size=0.5) +
theme(axis.text.x = element_text(size = 5, angle = 90), plot.title = element_text(size = 7)) +
xlab("ModuleName") +
ggsave(file=paste("Module_enrichment",feature,"_",bregion,"_",group1,"_",group2,".pdf", sep=""))
dev.off()

moduleinfo = as_tibble(moduleinfo)
moduleinfo
moduleinfo_stats = dplyr::select(moduleinfo, c("ModuleName","lncRNA","pc"))
moduleinfo_stats
moduleinfo_stats = pivot_longer(moduleinfo_stats, cols=!ModuleName, names_to="gene_type", values_to="count")
moduleinfo_stats
moduleinfo_stats = filter(moduleinfo_stats, ModuleName != "Not.Correlated")
moduleinfo_stats

ggplot(moduleinfo_stats, aes(x = reorder(ModuleName, -as.numeric(count)), y = as.numeric(count))) + 
 	geom_bar(aes(fill = gene_type), stat="identity", width = 0.7) +
 	xlab("ModuleName") +
 	ylab("Number of Genes in Module") +
 	theme(axis.text.x = element_text(size = 5, angle = 90), plot.title = element_text(size = 7)) +
    ggsave(file=paste("Module_lncRNA_vs_pc",feature,"_",bregion,"_",group1,"_",group2,".pdf", sep=""))
    dev.off()

moduleinfo = as_tibble(moduleinfo)
moduleinfo
moduleinfo_stats = dplyr::select(moduleinfo, c("ModuleName","number_genes_in_module","lncRNA","pc","lncRNA_p0.05","pc_p0.05"))
moduleinfo_stats
moduleinfo_stats = moduleinfo_stats %>% mutate("%lncRNA" = (as.numeric(lncRNA)-as.numeric(lncRNA_p0.05))/as.numeric(number_genes_in_module), 
	"%lncRNA_p0.05" = (as.numeric(lncRNA_p0.05)/as.numeric(number_genes_in_module)),
	"%pc" = (as.numeric(pc)-as.numeric(pc_p0.05))/as.numeric(number_genes_in_module),
	"%pc_p0.05" = (as.numeric(pc_p0.05)/as.numeric(number_genes_in_module)),
	)
moduleinfo_stats = moduleinfo_stats %>% dplyr::select(c("ModuleName","%lncRNA","%lncRNA_p0.05","%pc","%pc_p0.05"))
moduleinfo_stats
moduleinfo_stats = pivot_longer(moduleinfo_stats, cols=!ModuleName, names_to="gene_type", values_to="proportions")
moduleinfo_stats
moduleinfo_stats = filter(moduleinfo_stats, ModuleName != "Not.Correlated")
moduleinfo_stats
moduleinfo_stats = moduleinfo_stats %>% mutate("ModuleNameMod" = as.numeric(str_sub(ModuleName, 2)))
moduleinfo_stats

ggplot(moduleinfo_stats, aes(x = fct_reorder(ModuleName,ModuleNameMod, .desc=F), y = proportions)) + 
 	geom_bar(aes(fill = gene_type), stat="identity", width = 0.7) +
 	xlab("ModuleName") +
 	ylab("Proportions of Genes in Module") +
 	theme(axis.text.x = element_text(size = 5, angle = 90), plot.title = element_text(size = 7)) +
    ggsave(file=paste("Module_lncRNA_vs_pc_p0.05",feature,"_",bregion,"_",group1,"_",group2,".pdf", sep=""))
    dev.off()

moduleinfo = as_tibble(moduleinfo) 	
moduleinfo
moduleinfo_stats = dplyr::select(moduleinfo, c("ModuleName","lncRNA_p0.05_up","lncRNA_p0.05_down"))
moduleinfo_stats
moduleinfo_stats = pivot_longer(moduleinfo_stats, cols=!ModuleName, names_to="gene_type", values_to="count")
moduleinfo_stats
moduleinfo_stats = filter(moduleinfo_stats, ModuleName != "Not.Correlated")
moduleinfo_stats
moduleinfo_stats = moduleinfo_stats %>% mutate("ModuleNameMod" = as.numeric(str_sub(ModuleName, 2)))
moduleinfo_stats

ggplot(moduleinfo_stats, aes(x = fct_reorder(ModuleName,ModuleNameMod, .desc=F)), y = as.numeric(count)) + 
 	geom_bar(aes(y = as.numeric(count), fill = gene_type), stat="identity", width = 0.7) +
 	xlab("ModuleName") +
 	ylab("Number of lncRNA in Module") +
 	theme(axis.text.x = element_text(size = 5, angle = 90), plot.title = element_text(size = 7)) +
    ggsave(file=paste("Module_lncRNA_p0.05_up_down",feature,"_",bregion,"_",group1,"_",group2,".pdf", sep=""))
    dev.off()

moduleinfo = as_tibble(moduleinfo) 	
moduleinfo
moduleinfo_stats = dplyr::select(moduleinfo, c("ModuleName","pc_p0.05_up","pc_p0.05_down"))
moduleinfo_stats
moduleinfo_stats = pivot_longer(moduleinfo_stats, cols=!ModuleName, names_to="gene_type", values_to="count")
moduleinfo_stats
moduleinfo_stats = filter(moduleinfo_stats, ModuleName != "Not.Correlated")
moduleinfo_stats
moduleinfo_stats = moduleinfo_stats %>% mutate("ModuleNameMod" = as.numeric(str_sub(ModuleName, 2)))
moduleinfo_stats

ggplot(moduleinfo_stats, aes(x = fct_reorder(ModuleName,ModuleNameMod, .desc=F)), y = as.numeric(count)) + 
 	geom_bar(aes(y = as.numeric(count), fill = gene_type), stat="identity", width = 0.7) +
 	xlab("ModuleName") +
 	ylab("Number of protein_coding in Module") +
 	theme(axis.text.x = element_text(size = 5, angle = 90), plot.title = element_text(size = 7)) +
    ggsave(file=paste("Module_pc_p0.05_up_down",feature,"_",bregion,"_",group1,"_",group2,".pdf", sep=""))
    dev.off()

moduleinfo = as_tibble(moduleinfo)
moduleinfo
moduleinfo_stats = dplyr::select(moduleinfo, c("ModuleName","mu_kME_lncRNA","mu_kME_lncRNA_p0.05","mu_kME_lncRNA_p0.05_up","mu_kME_lncRNA_p0.05_down"))
moduleinfo_stats
moduleinfo_stats = pivot_longer(moduleinfo_stats, cols=!ModuleName, names_to="mu_kME_lncRNA", values_to="mu_kME")
moduleinfo_stats
moduleinfo_stats = filter(moduleinfo_stats, ModuleName != "Not.Correlated")
moduleinfo_stats
moduleinfo_stats = moduleinfo_stats %>% mutate("ModuleNameMod" = as.numeric(str_sub(ModuleName, 2)))
moduleinfo_stats

ggplot(moduleinfo_stats, aes(x = fct_reorder(ModuleName,ModuleNameMod, .desc=F), as.numeric(mu_kME))) + geom_point(aes(shape=mu_kME_lncRNA, col=mu_kME_lncRNA)) +
theme(axis.text.x = element_text(size = 5, angle = 90), plot.title = element_text(size = 7)) +
xlab("ModuleName") +
ylab("Mean Module Membership") +
scale_y_continuous(limits = c(0,1), breaks = breaks_extended(5)) +
ggsave(file=paste("Module_membership_lcnRNA",feature,"_",bregion,"_",group1,"_",group2,".pdf", sep=""))
dev.off()

moduleinfo = as_tibble(moduleinfo)
moduleinfo
moduleinfo_stats = dplyr::select(moduleinfo, c("ModuleName","mu_kME_pc","mu_kME_pc_p0.05","mu_kME_pc_p0.05_up","mu_kME_pc_p0.05_down"))
moduleinfo_stats
moduleinfo_stats = pivot_longer(moduleinfo_stats, cols=!ModuleName, names_to="mu_kME_pc", values_to="mu_kME")
moduleinfo_stats
moduleinfo_stats = filter(moduleinfo_stats, ModuleName != "Not.Correlated")
moduleinfo_stats
moduleinfo_stats = moduleinfo_stats %>% mutate("ModuleNameMod" = as.numeric(str_sub(ModuleName, 2)))
moduleinfo_stats

ggplot(moduleinfo_stats, aes(x = fct_reorder(ModuleName,ModuleNameMod, .desc=F), as.numeric(mu_kME))) + geom_point(aes(shape=mu_kME_pc, col=mu_kME_pc)) +
theme(axis.text.x = element_text(size = 5, angle = 90), plot.title = element_text(size = 7)) +
xlab("ModuleName") +
ylab("Mean Module Membership") +
scale_y_continuous(limits = c(0,1), breaks = breaks_extended(5)) +
ggsave(file=paste("Module_membership_pc",feature,"_",bregion,"_",group1,"_",group2,".pdf", sep=""))
dev.off()

#######################################################################################################################################
#######################################################################################################################################


###GAGE GO pathway analysis
# from http://bioconductor.org/packages/release/bioc/vignettes/gage/inst/doc/RNA-seqWorkflow.pdf

library(gage)
library(gageData)
library(pheatmap)
library("AnnotationDbi")
library("org.Mm.eg.db")

data(go.sets.mm)
data(go.subs.mm)

library(dplyr)
library(DESeq2)
library(topGO)
library(genefilter)
library(fdrtool)
library(xlsx)

res = read.csv(file="/data/NHLBI_BCB/Fayaz/25-limma-liu-Jie-RNAseq/02-0_vs_1_females_only/stat_results_gene_counts_0_1.csv", header=T, row.names=1)
res$symbol = mapIds(org.Mm.eg.db, keys=as.character(res$genesymbol), column="SYMBOL",   keytype="ALIAS", multiVals="first")
res$ensmbl = mapIds(org.Mm.eg.db, keys=as.character(res$genesymbol), column="ENSEMBL",  keytype="ALIAS", multiVals="first")
res$entrez = mapIds(org.Mm.eg.db, keys=as.character(res$genesymbol), column="ENTREZID", keytype="ALIAS", multiVals="first")
res$name =   mapIds(org.Mm.eg.db, keys=as.character(res$genesymbol), column="GENENAME", keytype="ALIAS", multiVals="first") 

foldchanges = res$log2FC
names(foldchanges) = res$entrez
head(foldchanges)

gobpsets = go.sets.mm[go.subs.mm$BP]
gobpres = gage(foldchanges, gsets=gobpsets, same.dir=TRUE)
lapply(gobpres, head)
head(gobpres$greater[, 2:5])
head(gobpres$less[, 2:5])
go_header=c("GO term","stat.mean","p.val","q.val","set.size")
write.table(t(go_header), file="GO-BP.up.sig.tsv", sep="\t", col.names = FALSE, row.names=F)
write.table(t(go_header), file="GO-BP.down.sig.tsv", sep="\t", col.names = FALSE, row.names=F)
write.table(na.omit(gobpres$greater[, 2:5]), file = "GO-BP.up.sig.tsv",  sep = "\t", col.names = FALSE,append = TRUE )
write.table(na.omit(gobpres$less[, 2:5]),    file = "GO-BP.down.sig.tsv",  sep = "\t", col.names = FALSE,append = TRUE )

nametag="BP"
pdf(file= paste0(nametag,"_up_and_down_pathways_pval.pdf"))
par(mar=c(5,20,2,2))
x1=gobpres$greater[20:1,"p.val"]
x2=gobpres$less[1:20,"p.val"]
xx=c(log10(x2),-log10(x1))
head(xx)
barplot(xx,horiz = T,names.arg = sub("GO:[0-9]*","",names(xx)),col=c(rep("green",20),rep("red",20)),las=2, xlab ="log10(pvalue)",cex.names = 0.50,xlim = c(-8,30))
 
dev.off()

pdf(file = paste0(nametag,"_up_and_down_pathways_qval.pdf"))
par(mar=c(5,20,2,2))
x1=gobpres$greater[20:1,"q.val"]
x2=gobpres$less[1:20,"q.val"]
xx=c(log10(x2),-log10(x1))
head(xx)
barplot(xx,horiz = T,names.arg = sub("GO:[0-9]*","",names(xx)),col=c(rep("green",20),rep("red",20)),las=2, xlab ="log10(qvalue)",cex.names = 0.50,xlim = c(-8,30))
 
dev.off()

gobpsets = go.sets.mm[go.subs.mm$MF]
gobpres = gage(foldchanges, gsets=gobpsets, same.dir=TRUE)
lapply(gobpres, head)
head(gobpres$greater[, 2:5])
head(gobpres$less[, 2:5])
go_header=c("GO term","stat.mean","p.val","q.val","set.size")
write.table(t(go_header), file="GO-MF.up.sig.tsv", sep="\t", col.names = FALSE, row.names=F)
write.table(t(go_header), file="GO-MF.down.sig.tsv", sep="\t", col.names = FALSE, row.names=F)
write.table(na.omit(gobpres$greater[, 2:5]), file = "GO-MF.up.sig.tsv",  sep = "\t", col.names = FALSE,append = TRUE )
write.table(na.omit(gobpres$less[, 2:5]),    file = "GO-MF.down.sig.tsv",  sep = "\t", col.names = FALSE,append = TRUE )

nametag="MF"
pdf(file= paste0(nametag,"_up_and_down_pathways_pval.pdf"))
par(mar=c(5,20,2,2))
x1=gobpres$greater[20:1,"p.val"]
x2=gobpres$less[1:20,"p.val"]
xx=c(log10(x2),-log10(x1))
head(xx)
barplot(xx,horiz = T,names.arg = sub("GO:[0-9]*","",names(xx)),col=c(rep("green",20),rep("red",20)),las=2, xlab ="log10(pvalue)",cex.names = 0.50,xlim = c(-8,30))
 
dev.off()

pdf(file = paste0(nametag,"_up_and_down_pathways_qval.pdf"))
par(mar=c(5,20,2,2))
x1=gobpres$greater[20:1,"q.val"]
x2=gobpres$less[1:20,"q.val"]
xx=c(log10(x2),-log10(x1))
head(xx)
barplot(xx,horiz = T,names.arg = sub("GO:[0-9]*","",names(xx)),col=c(rep("green",20),rep("red",20)),las=2, xlab ="log10(qvalue)",cex.names = 0.50,xlim = c(-8,30))
 
dev.off()

gobpsets = go.sets.mm[go.subs.mm$CC]
gobpres = gage(foldchanges, gsets=gobpsets, same.dir=TRUE)
lapply(gobpres, head)
head(gobpres$greater[, 2:5])
head(gobpres$less[, 2:5])
go_header=c("GO term","stat.mean","p.val","q.val","set.size")
write.table(t(go_header), file="GO-CC.up.sig.tsv", sep="\t", col.names = FALSE, row.names=F)
write.table(t(go_header), file="GO-CC.down.sig.tsv", sep="\t", col.names = FALSE, row.names=F)
write.table(na.omit(gobpres$greater[, 2:5]), file = "GO-CC.up.sig.tsv",  sep = "\t", col.names = FALSE,append = TRUE )
write.table(na.omit(gobpres$less[, 2:5]),    file = "GO-CC.down.sig.tsv",  sep = "\t", col.names = FALSE,append = TRUE )

nametag="CC"
pdf(file= paste0(nametag,"_up_and_down_pathways_pval.pdf"))
par(mar=c(5,20,2,2))
x1=gobpres$greater[20:1,"p.val"]
x2=gobpres$less[1:20,"p.val"]
xx=c(log10(x2),-log10(x1))
head(xx)
barplot(xx,horiz = T,names.arg = sub("GO:[0-9]*","",names(xx)),col=c(rep("green",20),rep("red",20)),las=2,xlab ="log10(pvalue)",cex.names = 0.50,xlim = c(-8,30))
 
dev.off()

pdf(file = paste0(nametag,"_up_and_down_pathways_qval.pdf"))
par(mar=c(5,20,2,2))
x1=gobpres$greater[20:1,"q.val"]
x2=gobpres$less[1:20,"q.val"]
xx=c(log10(x2),-log10(x1))
head(xx)
barplot(xx,horiz = T,names.arg = sub("GO:[0-9]*","",names(xx)),col=c(rep("green",20),rep("red",20)),las=2,xlab ="log10(qvalue)",cex.names = 0.50,xlim = c(-8,30))
 
dev.off()

# get the gene names for top 10
#go=go.gsets(species = "mouse", id.type = "EG")
#go.gs=go$go.sets
#head(go.gs)
#countData <- read.table("10.5_14.5_p_mm10_coding_cleaned_entrez.tsv", header=T, sep="\t", row.names=1)
#for (gs in rownames(gobpres$greater)[1:10]) {
#  outname = gsub(" |:|/", "_", substr(gs, 12, 100))
#  geneData(genes = go.gs, exprs = countData, ref = NULL, 
#           samp = NULL, outname = outname, txt = T, heatmap = T,
#           limit = 3, scatterplot = T)
#  }

# http://genomespot.blogspot.com.au/2014/09/data-analysis-step-8-pathway-analysis.html

xlsx::write.xlsx(read.table("GO-BP.up.sig.tsv", header=T, sep="\t"),file='GO.xlsx',row.names=F,sheetName = 'BP.UP')
xlsx::write.xlsx(read.table("GO-BP.down.sig.tsv", header=T, sep="\t"),file='GO.xlsx',row.names=F,sheetName = 'BP.DOWN', append=T)
xlsx::write.xlsx(read.table("GO-MF.up.sig.tsv", header=T, sep="\t"),file='GO.xlsx',row.names=F,sheetName = 'MF.UP', append=T)
xlsx::write.xlsx(read.table("GO-MF.down.sig.tsv", header=T, sep="\t"),file='GO.xlsx',row.names=F,sheetName = 'MF.DOWN', append=T)
xlsx::write.xlsx(read.table("GO-CC.up.sig.tsv", header=T, sep="\t"),file='GO.xlsx',row.names=F,sheetName = 'CC.UP', append=T)
xlsx::write.xlsx(read.table("GO-CC.down.sig.tsv", header=T, sep="\t"),file='GO.xlsx',row.names=F,sheetName = 'CC.DOWN', append=T)

pathw = c("GO:0002449 lymphocyte mediated immunity","GO:0002443 leukocyte mediated immunity","GO:0042110 T cell activation","GO:0035710 CD4-positive, alpha-beta T cell activation","GO:0050856 regulation of T cell receptor signaling pathway","GO:0030217 T cell differentiation","GO:0045622 regulation of T-helper cell differentiation","GO:0043367 CD4-positive, alpha-beta T cell differentiation","GO:0070741 response to interleukin-6","GO:0071354 cellular response to interleukin-6","GO:0051056 regulation of small GTPase mediated signal transduction","GO:0007264 small GTPase mediated signal transduction")

for (i in 1:length(pathw))
{
  	print(pathw[i])
  	ee = unlist(go.sets.hs[pathw[i]])
  	eenum = as.numeric(ee)		
  	#length(ee)
#	res = res.unique #only for transcript level analyses
  	res_subset_hmap = res[as.numeric(res$entrez) %in% eenum,]
  	#dim(pathway_genes)
  	#res_subset=arrange(res_subset,fc,desc(fc))
  	res_subset_hmap = res_subset_hmap[order(abs(as.numeric(res_subset_hmap$log2fc))),]
	pathway = sub(':', "_", pathw[i])
	write.csv(res_subset_hmap,file=paste("pathway_genes_counts","_",pathway,".csv",sep=""), row.names = F)
  
	#use top 20% of genes (upregulated and downregulated) if the number of genes in pathway is greater than or equal to 50
	if((dim(res_subset_hmap)[1]) >= 50){
		res_subset_hmap_updown_top25 = res_subset_hmap[(dim(res_subset_hmap)[1]-(0.20 * dim(res_subset_hmap)[1])):dim(res_subset_hmap)[1],]
	} else{
		res_subset_hmap_updown_top25 = res_subset_hmap	
	}
#  	res_subset_hmap_updown_top25 = res_subset_hmap[(dim(res_subset_hmap)[1]-(0.20 * dim(res_subset_hmap)[1])):dim(res_subset_hmap)[1],]
  	res_subset_hmap_updown_top25 = res_subset_hmap_updown_top25[order(as.numeric(res_subset_hmap_updown_top25$fc)),]
	#xmat=res_subset[,8:13]
#  	xmat = subset(res_subset_hmap_updown_top25, select = grep("FPKM", names(res_subset_hmap)))
#  	rownames(xmat) = res_subset_hmap_updown_top25$gene_name
  		
	##Heat Map for top results based on criteria sent by Kim
	res_subset_hmap = subset(res, (log2FC < -1.0 | log2FC > 1.0) & pvalue < 0.05)
#	res_subset_hmap = subset(res.unique, (as.numeric(res.unique$log2fc) < -0.50 | as.numeric(res.unique$log2fc) > 0.50) & as.numeric(res.unique$pval) < 0.05) #only for transcript level analyses
  	res_subset_hmap = res_subset_hmap[order(abs(as.numeric(res_subset_hmap$log2FC))),]

	#use top 20% of genes (upregulated and downregulated) if the number of genes in pathway is greater than or equal to 50
	if ((dim(res_subset_hmap)[1]) >= 50) {
		res_subset_hmap_updown_top25 = res_subset_hmap[(dim(res_subset_hmap)[1]-(0.20 * dim(res_subset_hmap)[1])):dim(res_subset_hmap)[1],]
	} else {
		res_subset_hmap_updown_top25 = res_subset_hmap	
	}

	scaleforheatmap = function(x){
		indexes = which(x==0)
		x[indexes] = mean(x)
		newx = log2(x/mean(x))
		meannewx = mean(newx)
		newxx = newx - meannewx
	}
#	xmat = subset(res_subset_hmap_updown_top25, select = grep("S", names(res_subset_hmap)))
#	colnames(xmat) = c("S37_F","S39_F","S41_F","S43_F","S45_F","S38_R","S40_R","S42_R","S44_R","S46_R")
	res_subset_hmap_updown_top25 = res_subset_hmap_updown_top25[order(as.numeric(res_subset_hmap_updown_top25$log2FC)),]
	xmat = subset(res_subset_hmap_updown_top25, select = grep("S*.*[.]1$", names(res_subset_hmap)))
  	rownames(xmat) = make.names(res_subset_hmap_updown_top25$genesymbol, unique=T)
#	xmat = log2(cpm(xmat)+1)
  	x = apply(xmat, 1, function(x) scaleforheatmap(x))	
	newxmat = cbind(t(x))
#	colnames(newxmat) = c("S37_F","S39_F","S41_F","S43_F","S45_F","S38_R","S40_R","S42_R","S44_R","S46_R")
#	colnames(newxmat) = colnames(xmat)$
#	colnames(newxmat) = pheno_data_comparison$SampleName_group
	pheatmap(newxmat,filename = paste("pheatmap","_","top_genes_count_based_analysis_pvalue_0.05_log2fc_0.50updown_0_1",".pdf",sep=""),fontsize_row = 5,fontsize_col = 5,fontsize = 3,show_rownames = T, scale="row", cluster_rows=F, cluster_cols=F, col=colorRampPalette( c("green", "black", "red"), space="rgb")(64), cellwidth=5, cellheight=5, main='Top_Genes_count_based_analysis_0_1', border_color=NA)
#	pheatmap(newxmat,filename = paste("pheatmap","_",pathway,".pdf",sep=""),fontsize_row = 6,fontsize_col = 6,scale="row", cluster_rows=F, cluster_cols=F, col=colorRampPalette( c("green", "black", "red"), space="rgb")(64), cellwidth=10, cellheight=10, fontsize = 4.3,show_rownames = T, main=pathway)
	dev.off()
}
############
# Heat map for snoRNAs only
counts_data = read.table(file="/data/NHLBI_BCB/Larochelle_Andre/05_Notch_Hypoxia_RNAseq/07-limma-DGE/snoRNAs_gencodev25_notch_hypoxia_rnaseq.txt", header=T, row.names=1)
snoRNAs_IDs = as.character(counts_data$snoRNAID)
counts_data = data.frame(counts_data[,1:dim(counts_data)[2]-1])
dim(counts_data)
counts_data[1:10,]
length(snoRNAs_IDs)
snoRNAs_IDs[1:10]

## filtering
#Keep genes with least 5 counts/reads in at least 6 replicates/samples
isexpr = rowSums(counts_data>5) >= 4
table(isexpr)
counts_data = counts_data[isexpr,]
dim(counts_data)
counts_data[1:10,]
counts_data_cpm = cpm(y=counts_data, log=T, normalized.lib.sizes=T)
dim(counts_data_cpm)
counts_data_cpm[1:10,]
snoRNAs_isexpr = snoRNAs_IDs[isexpr]
length(snoRNAs_isexpr)
snoRNAs_isexpr[1:10]

phenodata = read.csv(file="/data/NHLBI_BCB/Larochelle_Andre/05_Notch_Hypoxia_RNAseq/07-limma-DGE/phenodata.csv", header=T, row.names=1, check.names=F)
dim(phenodata)
phenodata
design = model.matrix(~phenodata$dx)

counts_data_cpm_adj = removeBatchEffect(x=counts_data_cpm, batch=phenodata$pair, design=design)
counts_data_cpm = scale(x=counts_data_cpm_adj)
dim(counts_data_cpm)
counts_data_cpm[1:10,]

paletteLength = 150
col = colorRampPalette( c("green", "black", "red"), space="rgb")(paletteLength)
metadata <- data.frame(c(rep("hypoxia_only",6), rep("notch_only",6), rep("notch_hypoxia",6)),row.names=phenodata$SampleNames)
colnames(metadata) <- c("condition")
pheatmap(counts_data_cpm,filename = paste("pheatmap","_","notch_hypoxia_notch_only_hypoxia_only_snoRNAs_v4",".pdf",sep=""), fontsize_col=5, fontsize=5, show_rownames=T, scale="row", cluster_rows=T, cluster_cols=T, col=col, cellwidth=10, cellheight=1, main='notch_hypoxia_notch_only_hypoxia_only', border_color=NA, cex=1, annotation_col=metadata, fontsize_row=1, labels_row=snoRNAs_isexpr)


############
# Heat map for select pathways

“HSC-R_FDR_0.05_Probe_List”
“Notch_signaling_pathway”
“GRAHAM_NORMAL_QUIESCENT_VS_NORMAL_DIVIDING_UP”




# out <- pheatmap(data, 
#      show_rownames=F, cluster_cols=T, cluster_rows=T, scale="row",
#      cex=1, clustering_distance_rows="euclidean", cex=1,
#      clustering_distance_cols="euclidean", clustering_method="complete", border_color=FALSE,
#      annotation_col=metadata,
#      annotation_row=metadata_gene)

# annotation_row = data.frame(gene_direction_ibrutinib = as.vector(statres_gsigs$Direction.of.change))
# ann_colors = list(gene_direction_ibrutinib=c(UP="yellow", DN="blue"), timepoint=c(pre="orange", sixm="purple"))
# rownames(annotation_row) = make.names(statres_gsigs$genesymbol, unique=T)

# metadata_gene <- data.frame(
#      c(rep("Tcell", nrow(data)/2), rep("Bcell", nrow(data)/2)),
#      row.names=rownames(data))
# colnames(metadata_gene) <- c("Cell")

# dev.off()

#=====================================================================================
#
#  Code chunk 1 - heat maps for interesting pathways
#
#=====================================================================================

setwd("/data/NHLBI_BCB/Larochelle_Andre/05_Notch_Hypoxia_RNAseq/07-limma-DGE/03-run-WebGesstaltR")
# setwd("/gpfs/gsfs5/users/NHLBI_BCB/Larochelle_Andre/05_Notch_Hypoxia_RNAseq/07-limma-DGE/02-GSEA/02-intersect_3_conditions")

# needed for reading GMT files
library(clusterProfiler)
library(pheatmap)
library(edgeR)
library(DESeq2)

hypoxia_Vs_notch_hypoxia_gsea_results = read.table(file="/data/NHLBI_BCB/Larochelle_Andre/05_Notch_Hypoxia_RNAseq/07-limma-DGE/03-run-WebGesstaltR/Project_hypoxia_Vs_notch_hypoxia_GSEA_custom_pathways/enrichment_results_hypoxia_Vs_notch_hypoxia_GSEA_custom_pathways.txt", header=T, row.names=1)
notch_Vs_notch_hypoxia_gsea_results = read.table(file="/data/NHLBI_BCB/Larochelle_Andre/05_Notch_Hypoxia_RNAseq/07-limma-DGE/03-run-WebGesstaltR/Project_notch_Vs_notch_hypoxia_GSEA_custom_pathways/enrichment_results_notch_Vs_notch_hypoxia_GSEA_custom_pathways.txt", header=T, row.names=1)
gsea_results = merge(x=hypoxia_Vs_notch_hypoxia_gsea_results, y=notch_Vs_notch_hypoxia_gsea_results, by="row.names")

# pathways_of_interest = read.gmt("/data/NHLBI_BCB/Larochelle_Andre/05_Notch_Hypoxia_RNAseq/08-WGCNA-Cemitools/c5.bp.v7.0.symbols.pathways.of.interest.gmt")
# pathways_gmt = pathways_of_interest
# pathways = as.character(unique(pathways_of_interest$ont))

counts_data = read.table(file="/data/NHLBI_BCB/Larochelle_Andre/05_Notch_Hypoxia_RNAseq/06-featureCounts/gene.featurecount.txt", header=T, row.names=1, check.names=F)
phenodata = read.csv(file="/data/NHLBI_BCB/Larochelle_Andre/05_Notch_Hypoxia_RNAseq/07-limma-DGE/phenodata.csv", header=T, row.names=1, check.names=F)
## filtering
#Keep genes with least 5 counts/reads in at least 4 replicates/samples
isexpr = rowSums(counts_data>5) >= 4
table(isexpr)
counts_data = counts_data[isexpr,]
# dge = DGEList(counts = counts_data)
# dge = calcNormFactors(dge, method="TMM")
# design = model.matrix(~ phenodata$dx + phenodata$pair)
# v = voom(dge, design) #Transform count data to log2-counts per million (logCPM+0.5), estimate the mean-variance relationship and use this to compute appropriate observation-level weights. The data are then ready for linear modelling.
annot = read.table(file="/data/NHLBI_BCB/Fayaz/02-RNA-seq-transcript-level/data_source/GENCODE_primary_assembly_Mouse_Human/gencode.v25.primary_assembly.annotation.txt", header=1, row.names=1)
design = model.matrix(~ phenodata$pair + phenodata$dx) # batch first, condition second for DESeq2
dds <- DESeqDataSetFromMatrix(countData = counts_data, colData = phenodata, design = design)
rlog_counts_data = rlog(dds,  blind=FALSE)
rlog_counts_data = assay(rlog_counts_data)
expressiondf = merge(x=rlog_counts_data, y=annot, by="row.names")

notch_only_Vs_notch_hypoxia = read.csv(file="/data/NHLBI_BCB/Larochelle_Andre/05_Notch_Hypoxia_RNAseq/07-limma-DGE/01-notch_only_vs_notch_hypoxia_paired_analysis/stat_results_gene_counts_notch_only_notch_hypoxia.csv", header=T, row.names=1)
hypoxia_only_Vs_notch_hypoxia = read.csv(file="/data/NHLBI_BCB/Larochelle_Andre/05_Notch_Hypoxia_RNAseq/07-limma-DGE/00-hypoxia_only_vs_notch_hypoxia_paired_analysis/stat_results_gene_counts_hypoxia_only_notch_hypoxia.csv", header=T, row.names=1)

expressiondf_notch = merge(x=expressiondf, y=notch_only_Vs_notch_hypoxia, by.x="Row.names", by.y="row.names")
expressiondf_notch = data.frame(expressiondf_notch[,2:19], genesymbol=expressiondf_notch$genesymbol.x, log2FC_n_Vs_nh=expressiondf_notch$log2FC, row.names=expressiondf_notch$Row.names)

expressiondf_hypoxia = merge(x=expressiondf, y=hypoxia_only_Vs_notch_hypoxia, by.x="Row.names", by.y="row.names")
expressiondf_hypoxia = data.frame(expressiondf_hypoxia[,2:19], genesymbol=expressiondf_hypoxia$genesymbol.x, log2FC_h_Vs_nh=expressiondf_hypoxia$log2FC, row.names=expressiondf_hypoxia$Row.names)

expressiondf1 = merge(x=expressiondf_notch, y=expressiondf_hypoxia, by="row.names")
expressiondf2 = data.frame(expressiondf1[,2:19], expressiondf1$genesymbol.x, expressiondf1$log2FC_n_Vs_nh, expressiondf1$log2FC_h_Vs_nh, row.names=expressiondf1$Row.names)
colnames(expressiondf2) = c("D1","D2","D3","D4","D5","D6","D7","D8","D9","D10","D11","D12","D13","D14","D15","D16","D17","D18","genesymbol","log2FC_n_Vs_nh","log2FC_h_Vs_nh")

design = model.matrix(~ phenodata$dx) # used for removing subject specific effects (removeBatchEffect)
# pairedttest_stats = matrix(nrow=length(pathways), ncol=5, byrow=T);
# pairedttest_stats = matrix(nrow=dim(gsea_results)[1], ncol=5, byrow=T);

# for (i in 1:length(pathways)) {
#	print(pathways[i])
# for (i in 1:dim(gsea_results)[1]) {
#	print(gsea_results[i,1])	
#	leadingedge_genes_hypoxia = gsea_results[i,"userId.x"]
#	leadingedge_genes_hypoxia = strsplit(as.character(leadingedge_genes_hypoxia), split=";")[[1]]
#	leadingedge_genes_notch = gsea_results[i,"userId.y"]
#	leadingedge_genes_notch = strsplit(as.character(leadingedge_genes_notch), split=";")[[1]]
#	genesinpathway = unique(c(leadingedge_genes_hypoxia,leadingedge_genes_notch))
#	genesinpathway = intersect(x=leadingedge_genes_hypoxia,y=leadingedge_genes_notch)
#	genesinpathway = pathways_gmt[pathways_gmt$ont==pathways[i],]$gene
#	print(length(genesinpathway))

#	pathway_expressiondf = expressiondf2[which(expressiondf2$genesymbol %in% genesinpathway),]
#	pathway_expressiondf = pathway_expressiondf[order(-abs(pathway_expressiondf$log2FC_h_Vs_nh)),]
#	head(pathway_expressiondf)

#	heatmap_expressiondf = data.frame(pathway_expressiondf[,1:18], row.names=make.names(pathway_expressiondf$genesymbol,unique=TRUE), check.names=FALSE, check.rows=FALSE)
#	heatmap_expressiondf_subject_effect_removed = removeBatchEffect(x=heatmap_expressiondf, batch=phenodata$pair, design=design)
#	heatmap_expressiondf_subject_effect_removed_scaled = scale(x=heatmap_expressiondf_subject_effect_removed)

#	paletteLength = 150
#	col = colorRampPalette( c("green", "black", "red"), space="rgb")(paletteLength)
#	metadata <- data.frame(c(rep("hypoxia_only",6), rep("notch_only",6), rep("notch_hypoxia",6)),row.names=phenodata$SampleNames)
#	colnames(metadata) <- c("condition")
#	pheatmap(heatmap_expressiondf_subject_effect_removed,filename = paste("pheatmap","_",pathways[i],".pdf",sep=""), fontsize_col=5, fontsize=5, show_rownames=F, scale="row", cluster_rows=T, cluster_cols=F, main=pathways[i], col=col, border_color=NA, cex=1, annotation_col=metadata, fontsize_row=2, width=4, height=4)
#	pheatmap(heatmap_expressiondf_subject_effect_removed,filename = paste("pheatmap","_",gsea_results[i,1],".v2.pdf",sep=""), fontsize_col=5, fontsize=5, show_rownames=F, scale="row", cluster_rows=F, cluster_cols=F, main=gsea_results[i,1], border_color=NA, cex=1, col=col, annotation_col=metadata, fontsize_row=2, width=4, height=4)

#	hypoxia_only = heatmap_expressiondf[,1:6]
#	notch_only = heatmap_expressiondf[,7:12]
#	notch_hypoxia = heatmap_expressiondf[,13:18]

#	hypoxia_only_pathway_mean = apply(hypoxia_only,1,mean)
#	notch_only_pathway_mean = apply(notch_only,1,mean)
#	notch_hypoxia_pathway_mean = apply(notch_hypoxia,1,mean)

#	hypoxia_only_Vs_notch_hypoxia = t.test(x=notch_hypoxia_pathway_mean, y=hypoxia_only_pathway_mean, paired=T, alternative="two.sided")
#	hypoxia_only_Vs_notch_hypoxia_tstat = as.vector(hypoxia_only_Vs_notch_hypoxia$statistic)
#	hypoxia_only_Vs_notch_hypoxia_pvalue = hypoxia_only_Vs_notch_hypoxia$p.value

#	notch_only_Vs_notch_hypoxia = t.test(x=notch_hypoxia_pathway_mean, y=notch_only_pathway_mean, paired=T, alternative="two.sided")
#	notch_only_Vs_notch_hypoxia_tstat = as.vector(notch_only_Vs_notch_hypoxia$statistic)
#	notch_only_Vs_notch_hypoxia_pvalue = notch_only_Vs_notch_hypoxia$p.value

#	pairedttest_stats[i,1] =  pathways[i]
#	pairedttest_stats[i,1] =  gsea_results[i,1]
#	pairedttest_stats[i,2] =  hypoxia_only_Vs_notch_hypoxia_tstat
#	pairedttest_stats[i,3] =  hypoxia_only_Vs_notch_hypoxia_pvalue
#	pairedttest_stats[i,4] =  notch_only_Vs_notch_hypoxia_tstat
#	pairedttest_stats[i,5] =  notch_only_Vs_notch_hypoxia_pvalue
# }

# colnames(pairedttest_stats) = c("pathway","h_Vs_nh_tstat","h_Vs_nh_pvalue","n_Vs_nh_tstat","n_Vs_nh_pvalue")
# head(pairedttest_stats)
# write.csv(pairedttest_stats, file="pairedttest_stats_GOBP.csv", quote=F, row.names=F)

hscs = which(gsea_results[,1] %in% c("DOULATOV_HSC","FARES_HSC","EPPERT_HSC"))
for (i in 1:length(hscs)){
	if(gsea_results[i,1]=="DOULATOV_HSC"){
		leadingedge_genes_hypoxia = gsea_results[i,"userId.x"]
		leadingedge_genes_hypoxia = strsplit(as.character(leadingedge_genes_hypoxia), split=";")[[1]]
		leadingedge_genes_notch = gsea_results[i,"userId.y"]
		leadingedge_genes_notch = strsplit(as.character(leadingedge_genes_notch), split=";")[[1]]
		a = unique(c(leadingedge_genes_hypoxia,leadingedge_genes_notch))
#		print(a)
	}

	if(gsea_results[i,1]=="EPPERT_HSC"){
		leadingedge_genes_hypoxia = gsea_results[i,"userId.x"]
		leadingedge_genes_hypoxia = strsplit(as.character(leadingedge_genes_hypoxia), split=";")[[1]]
		leadingedge_genes_notch = gsea_results[i,"userId.y"]
		leadingedge_genes_notch = strsplit(as.character(leadingedge_genes_notch), split=";")[[1]]
		b = unique(c(leadingedge_genes_hypoxia,leadingedge_genes_notch))
#		print(b)	
	}
		
	if(gsea_results[i,1]=="FARES_HSC"){
		leadingedge_genes_hypoxia = gsea_results[i,"userId.x"]
		leadingedge_genes_hypoxia = strsplit(as.character(leadingedge_genes_hypoxia), split=";")[[1]]
		leadingedge_genes_notch = gsea_results[i,"userId.y"]
		leadingedge_genes_notch = strsplit(as.character(leadingedge_genes_notch), split=";")[[1]]
		c = unique(c(leadingedge_genes_hypoxia,leadingedge_genes_notch))
#		print(c)	
	}
}

print(a)
print(b)
print(c)
Reduce(intersect, list(a,b,c))

#=====================================================================================
#
#  Code chunk 2 - run WebGestaltR for the HSC,notch signaling custom pathways
#
#=====================================================================================

library("WebGestaltR")

# rankfile = read.table(file="/data/NHLBI_BCB/Larochelle_Andre/05_Notch_Hypoxia_RNAseq/07-limma-DGE/03-run-WebGesstaltR/hypoxia_Vs_notch_hypoxia.rnk", header=T)
rankfile = read.table(file="/data/NHLBI_BCB/Larochelle_Andre/05_Notch_Hypoxia_RNAseq/07-limma-DGE/03-run-WebGesstaltR/notch_Vs_notch_hypoxia.rnk", header=T)

WebGestaltR(enrichMethod="GSEA",
	organism="hsapiens",
	enrichDatabaseFile="/data/NHLBI_BCB/Larochelle_Andre/05_Notch_Hypoxia_RNAseq/07-limma-DGE/03-run-WebGesstaltR/HSC.gmt",
	enrichDatabaseType="genesymbol",
	interestGene=rankfile,
	interestGeneType="genesymbol",
	collapseMethod="mean",
	minNum=10,
	maxNum=5000,
	sigMethod="fdr",
	fdrMethod="BH",
	fdrThr=1,
	topThr=10,
	reportNum=20,
	perNum=1000,
	isOutput=TRUE,
	outputDirectory="/data/NHLBI_BCB/Larochelle_Andre/05_Notch_Hypoxia_RNAseq/07-limma-DGE/03-run-WebGesstaltR",
	projectName="notch_Vs_notch_hypoxia_GSEA_custom_pathways",
	saveRawGseaResult=TRUE,
	gseaPlotFormat="png",
	nThreads=4
)

pdf("lncRNAKB_pca_pairs_plot_15PCs_lncRNAs_tissueorigin.pdf")
pairs(pcs[,9:15], main="scatter_plot_matrix_of_Principal_Components_PCs9toPCs15", pch=21, bg=c('#e41a1c','#377eb8','#4daf4a')[as.factor(pcs$Origin)], cex.main=0.75)
dev.off()

pdf("lncRNAKB_pca_pairs_plot_5PCs_lncRNAs_tissueorigin.pdf")
pairs(pcs[,1:5], main="scatter_plot_matrix_of_Principal_Components_PCs1toPCs5", pch=21, bg=c('#e41a1c','#377eb8','#4daf4a')[as.factor(pcs$Origin)], cex.main=0.75)
dev.off()


pdf("pca_plot_pc9_vs_pc12.pdf")
# Add extra space to right of plot area; change clipping to figure
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(pcs[,9], pcs[,12], main="scatter plot matrix of Principal Components_PC9.vs.PC12", pch=16, col=c('#e41a1c','#377eb8','#4daf4a')[as.factor(pcs$Origin)], cex=1.20, cex.main=0.80, xlab="PC9", ylab="PC12")
text(pcs[,9], pcs[,12], labels=pcs$Origin, pos=3, cex=0.30)
legend("topright", inset=c(-0.23,0), legend=c('Endoderm','Mesoderm','Ectoderm'), pch=c(16,16), col=levels(as.factor(c('#e41a1c','#377eb8','#4daf4a')[as.factor(pcs$Origin)])), title="Origin", cex=0.60)
dev.off()


