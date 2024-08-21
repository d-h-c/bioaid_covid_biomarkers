######################################################################
# Normalisation and identification of COVID-19 signatures from BioAid A and B
# Written in R 4.0.3
######################################################################


source(file = 'functions.R')
gene.names <- read.csv("gene_names_ensembl89.csv")
colnames(gene.names) <- c("ensembl_gene_id","external_gene_name","biotype")

load.packages()
library(PCAtools)
library(stringr)
library(ggplot2)
library(DESeq2)
library(sva)
library(biomaRt)
library(ashr)
library(nnet)
library(caTools)
library(viridis)
library(tidyverse)



################################################################################################################
# read in batch A and remove the plate effects
################################################################################################################

####### Read in Batch A 
# ArrayExpress E-MTAB-10527
bioaid_a_exp <- read.csv(gzfile("E-MTAB-10527_gene_count_matrix.csv.gz"), row.names = 1, check.names =F)
bioaid_a_meta <- read.csv("batch_a_metadata.csv", row.names = 1 )

outliers <-  c("I-P-0255", "I-P-0277", "I-P-0328", "I-P-1092", "I-P-1195", "I-P-1222")
keepsams <- rownames(bioaid_a_meta)[!rownames(bioaid_a_meta) %in% outliers]

bioaid_a_meta <- bioaid_a_meta[match(colnames(bioaid_a_exp), bioaid_a_meta$Patient.ID.E),]
bioaid_a_meta <- bioaid_a_meta[keepsams,]
bioaid_a_exp <- bioaid_a_exp[,keepsams]

# quantify how many genes are completely missing
a_zero <- rowSums(bioaid_a_exp)
a_zero <- a_zero[a_zero==0]
# remove those completely missing
bioaid_a_exp <- bioaid_a_exp[!(rownames(bioaid_a_exp) %in% names(a_zero)),]

# run DEA to identify genes SDE between plates but only those SDE from bacterial and viral combined
plate_cor_meta <- bioaid_a_meta[bioaid_a_meta$diagnosis.E %in% c("Gram Negative BC","Gram Positive BC","Viral positive"),]
plate_cor_meta[plate_cor_meta$diagnosis.E %in% c("Gram Negative BC","Gram Positive BC"),"diagnosis.E"] <- "Bacteremia"

plate_cor_exp <- bioaid_a_exp[,match(plate_cor_meta$Patient.ID.E, colnames(bioaid_a_exp))]

DES_plate <- DESeqDataSetFromMatrix(countData = plate_cor_exp,
                                    colData = plate_cor_meta, 
                                    design = ~ 1)

keep <- rowSums(counts(DES_plate) >= 10 ) >= 3
DES_plate <- DES_plate[keep,]
DES_plate <- estimateSizeFactors(DES_plate)

mm <- model.matrix(~ 0 + factor(Plate.E) + factor(diagnosis.E) + factor(sex), data = data.frame(colData(DES_plate)))
colnames(mm) <- c("plate 13", 'plate 14', 'viral', 'male')

design(DES_plate) <- mm

DES_plate <- DESeq(DES_plate, betaPrior = FALSE)
DES_plate <- DES_plate[which(mcols(DES_plate)$betaConv),]

sde_plate <- results(DES_plate, 
                     independentFiltering=TRUE,
                     alpha=0.05, 
                     pAdjustMethod="BH", 
                     parallel=FALSE, 
                     contrast=list(c("plate.13"), c("plate.14")))

sde_plate <- deseq.results(sde_plate, gene.names, '../Results/SDE_Plate')
sde_plate <- sde_plate[sde_plate$padj < 0.001,]
# 5734 genes SDE 

# remove the genes identified as SDE between plate 13 and 14 
bioaid_a_exp <- bioaid_a_exp[!(rownames(bioaid_a_exp) %in% rownames(sde_plate)),]

# plot PCA of the raw counts
bioaid_a_meta <- bioaid_a_meta[match(colnames(bioaid_a_exp), rownames(bioaid_a_meta)),]

pc_a <- pca(log2(bioaid_a_exp+1), metadata = bioaid_a_meta, scale = F)

biplot(pc_a, 
       lab = NULL,
       colby = 'Plate.E', 
       legendPosition = 'right', 
       colLegendTitle = "plate")

biplot(pc_a, 
       lab = NULL,
       colby = 'disease', 
       legendPosition = 'right', 
       colLegendTitle = "group")

pdf(file = 'figures/eigen_cor_a.pdf', height=7, width =10)
eigencorplot(pc_a, metavars = c('disease', 'gend', 'Plate.E', 'age'))
dev.off()

################################################################################################################
# read in batch B
################################################################################################################
# E-MTAB-13307
bioaid_b_exp <- read.delim(file = 'batch_P200157_gene_count_matrix.txt', check.names =F, row.names = 1)
bioaid_b_meta <- read.csv(file = 'Metadata_Updated_Merged.csv', header = T)


bioaid_b_meta <- bioaid_b_meta[match(colnames(bioaid_b_exp), as.character(bioaid_b_meta$X)),]
rownames(bioaid_b_meta) <- bioaid_b_meta$X

# plot PCA of the raw counts
pc_b <- pca(log2(bioaid_b_exp+1), metadata = bioaid_b_meta, scale = F)

biplot(pc_b, 
       lab = NULL,
       colby = 'Group', 
       legendPosition = 'right', 
       colLegendTitle = "Group")

pdf(file = 'figures/eigen_cor_b.pdf', height=7, width =10)
eigencorplot(pc_b, metavars = c('Group', 'Age', "Sex", 'COVID_CATEGORY'))
dev.off()

################################################################################################################
# merge datasets, remove ribosomal genes, combat to remove batch effects
################################################################################################################

#### before merging the datasets, quantify how many genes have total of zero in each 
a_zero <- rowSums(bioaid_a_exp)
a_zero <- a_zero[a_zero==0]

b_zero <- rowSums(bioaid_b_exp)
b_zero <- b_zero[b_zero==0]

length(a_zero)
length(b_zero)



# remove 12529 genes completely missing in bioaid B
bioaid_b_exp <- bioaid_b_exp[!(rownames(bioaid_b_exp) %in% names(b_zero)),]

# merge datasets
genes_match <- intersect(rownames(bioaid_a_exp), rownames(bioaid_b_exp))
# 35167 genes are present in both 
bioaid_b_exp <- bioaid_b_exp[match(genes_match, rownames(bioaid_b_exp)),]
bioaid_a_exp <- bioaid_a_exp[match(genes_match, rownames(bioaid_a_exp)),]

merged_df <- cbind(bioaid_a_exp, bioaid_b_exp)

bioaid_a_meta <- bioaid_a_meta[match(colnames(bioaid_a_exp), bioaid_a_meta$Patient.ID.E),]
bioaid_b_meta <- bioaid_b_meta[match(colnames(bioaid_b_exp), bioaid_b_meta$X),]

# load in premade merged metadata file 
merged_meta <- read.csv(file = 'merged_meta.csv', check.names = F)
merged_meta <- merged_meta[match(colnames(merged_df), merged_meta$ID),]
rownames(merged_meta) <- merged_meta$ID

# check whether XIST matches with sex on the database 
xist <- data.frame(expression = t(merged_df['ENSG00000229807',]), 
                   sex = merged_meta$Sex)

ggplot(xist, aes(x = sex, y = ENSG00000229807))+geom_boxplot()+theme_bw()+geom_point()


# remove these 
remove_id <- read.csv("genes_remove.csv", row.names = 1)[,1]
merged_df <- merged_df[!(rownames(merged_df) %in% remove_id),]

# also identify genes with RPL-, and RPS-
remove_g <- unique(c(unlist(gene.names$ensembl_gene_id[str_detect(gene.names$external_gene_name, 'RPL')]), 
                     unlist(gene.names$ensembl_gene_id[str_detect(gene.names$external_gene_name, 'RPS')])))


# remove the ribosomal genes
merged_df <- merged_df[!(rownames(merged_df) %in% remove_g),]
dim(merged_df)

### visualise the merged dataframe before batch effect correction 
merged_meta$Group2 <- merged_meta$Group
merged_meta$Group2[merged_meta$COVID_CATEGORY %in% c("READMITTED PREVIOUS COVID", "EXCLUDE", "1ST TEST NEGATIVE", "MIXED INFECTION")] <- "COVID_Other"

pc_merge <- pca(log2(merged_df+1), metadata = merged_meta, scale = F)

pdf(file = 'figures/PCAs_Normalisation/Merged_Raw.pdf', height = 8, width = 10)
biplot(pc_merge, 
       lab = NULL,
       colby = 'Group2', 
       legendPosition = 'right', 
       colLegendTitle = "Group")
biplot(pc_merge, 
       lab = NULL,
       colby = 'Sex', 
       legendPosition = 'right', 
       colLegendTitle = "Sex")
biplot(pc_merge, 
       lab = NULL,
       colby = 'Cohort', 
       legendPosition = 'right', 
       colLegendTitle = "Cohort")
dev.off()

# Need to correct the major batch effects
merged_meta$Plate.Batch <- merged_meta$Plate
merged_meta$Plate.Batch[merged_meta$Cohort == 'New'] <- 3
merged_meta$Plate.Batch[merged_meta$Plate.Batch == "Plate 13"] <- 1
merged_meta$Plate.Batch[merged_meta$Plate.Batch == "Plate 14"] <- 2
merged_meta$Plate.Batch <- as.factor(merged_meta$Plate.Batch)

# remove those with all zero again 
zero <- names(rowSums(merged_df)[rowSums(merged_df)==0])
merged_df <- merged_df[!(rownames(merged_df) %in% zero),]

# run ComBat
combat_merge <- ComBat_seq(as.matrix(merged_df),
                           batch = merged_meta$Plate.Batch)

merged_meta <- merged_meta[match(colnames(combat_merge), merged_meta$ID),]
rownames(merged_meta) <- merged_meta$ID

pc.combat.merge <- pca(log2(combat_merge+1), meta = merged_meta, scale = F)

pdf(file = 'figures/PCAs_Normalisation/Merged_Combat.pdf', height = 5, width =6)
biplot(pc.combat.merge, 
       lab = NULL, 
       colby = 'Plate.Batch', 
       legendPosition = 'right', 
       colLegendTitle = "Batch", 
       pointSize = 2)

biplot(pc.combat.merge, 
       lab = NULL, 
       colby = 'Group2', 
       legendPosition = 'right', 
       colLegendTitle = "Batch", 
       pointSize = 2)

biplot(pc.combat.merge, 
       lab = NULL, 
       colby = 'Cohort', 
       legendPosition = 'right', 
       colLegendTitle = "Batch", 
       pointSize = 2)

dev.off()

eigencorplot(pc.combat.merge, metavars = c('Plate.Batch', "Cohort", "Group", "Sex", "Age", "Plate"))

merged_meta$Group[merged_meta$Group %in% c("Gram Negative BC", "Gram Positive BC")] <- 'DB'
merged_meta$Group[merged_meta$Group == "Viral positive"] <- 'DV'


################################################################################################################
# DESeq2 normalisation
################################################################################################################

DES_obj <- DESeqDataSetFromMatrix(countData = combat_merge+1,
                                  colData = merged_meta, 
                                  design = ~1)

keep <- rowSums(counts(DES_obj) >= 10 ) >= 3
DES_obj <- DES_obj[keep,]
DES_obj <- estimateSizeFactors(DES_obj)

pdf(file = "figures/Histogram_SizeFactors.pdf", height = 6, width = 6)
hist(sizeFactors(DES_obj))
dev.off()

normalized <- counts(DES_obj, normalized=TRUE)
dim(normalized)

pc_deseq <- pca(log2(normalized+1), metadata = merged_meta, scale = F)
biplot(pc_deseq, 
       lab = NULL, 
       colby = 'Cohort', 
       legendPosition = 'right', 
       colLegendTitle = "Cohort")
biplot(pc_deseq, 
            lab = NULL, 
            colby = 'Group', 
            legendPosition = 'right', 
            colLegendTitle = "Group")

# look at totals 
sums_norm <- data.frame(totals = colSums(normalized), 
                        cohort = merged_meta$Cohort, 
                        group = merged_meta$Group)

sums_norm$group <- factor(sums_norm$group, levels = c('HC', 'Control', 'DB', 'DV', "COVID19"))

pdf(file= 'figures/Sums_Norm.pdf', height = 5, width =6)
ggplot(sums_norm, aes(x = group, y = totals, color = cohort))+
  geom_boxplot()+
  theme_bw()+
  geom_point(position = position_jitterdodge(jitter.width = 0.1))
dev.off()


###### COVID data phenotypes 
ex <- colData(DES_obj)$ID[!(colData(DES_obj)$COVID_CATEGORY == "DEFINITE") &
                            colData(DES_obj)$Group == "COVID19"]
ex <- ex[!(is.na(ex))]
samples <- colData(DES_obj)$ID[!(colData(DES_obj)$ID %in% ex)]

# remove the COVID exclude samples
DES_obj <- DES_obj[,samples]

DES_obj$admission_severity[!(DES_obj$Group == "COVID19")] <- 0



# create a column for severity which also includes sex but for the 3 categories of covid
# mild = 1, 2, moderate = 3,4, severe = 5,6,7
DES_obj$severity_grouped <- "Not_COVID"
DES_obj$severity_grouped[DES_obj$admission_severity %in% c(1,2)] <- 'Mild'
DES_obj$severity_grouped[DES_obj$admission_severity %in% c(3,4)] <- 'Moderate'
DES_obj$severity_grouped[DES_obj$admission_severity %in% c(5,6)] <- 'Severe'
DES_obj$severity_grouped <- as.factor(DES_obj$severity_grouped)

colData(DES_obj)$maximum_severity <- as.factor(colData(DES_obj)$maximum_severity)
colData(DES_obj)$admission_severity <- as.factor(colData(DES_obj)$admission_severity)

save(file = 'DESeq_object.RData', DES_obj, merged_meta)
### this data is used for subsequent analyses and also used for the sex differences analysis 

################################################################################################################
# DESeq analysis
################################################################################################################

load(file = 'DESeq_object.RData', verbose = T)

# remove the probable viral samples 
DES.obj.group <- DES_obj[,colData(DES_obj)$ID[!(colData(DES_obj)$Group == "PV")]]
colData(DES.obj.group)$BV <- "N"
colData(DES.obj.group)$BV[colData(DES.obj.group)$Group %in% c("DB", "DV")] <- "Y"
DES.obj.group$Group <- as.factor(DES.obj.group$Group)

# model matrix 
mm <- model.matrix(~0 + Group + Group + factor(Sex) + Age + Plate.Batch, data = colData(DES.obj.group))
colnames(mm) <- c("Control", "COVID19", "DB", "DV", "HC", "M","Age", "Batch_2", "Batch_3")

# run model 
design(DES.obj.group) <- mm
DES.obj.group <- DESeq(DES.obj.group, betaPrior = FALSE, fitType='local')
DES.obj.group <- DES.obj.group[which(mcols(DES.obj.group)$betaConv),]

# PCA 
counts_norm <- counts(DES.obj.group, normalized = TRUE)




DES.obj.group$group_pca <- as.character(DES.obj.group$Group)
DES.obj.group$group_pca[DES.obj.group$group_pca=="Control"] <- "Non-infected controls"
DES.obj.group$group_pca[DES.obj.group$group_pca=="COVID19"] <- "COVID-19"
DES.obj.group$group_pca[DES.obj.group$group_pca=="DB"] <- "Bacterial"
DES.obj.group$group_pca[DES.obj.group$group_pca=="DV"] <- "Viral"
DES.obj.group$group_pca[DES.obj.group$group_pca=="HC"] <- "Healthy controls"
DES.obj.group$sex_pca <- as.character(DES.obj.group$Sex)
DES.obj.group$sex_pca[DES.obj.group$sex_pca=="M"] <- "Male"
DES.obj.group$sex_pca[DES.obj.group$sex_pca=="F"] <- "Female"
DES.obj.group$batch_pca <- DES.obj.group$Cohort
DES.obj.group$batch_pca[DES.obj.group$batch_pca=="New"] <- 'B'
DES.obj.group$group_pca <- fct_relevel(as.factor(DES.obj.group$group_pca), 
                                       'COVID-19', 'Bacterial', 'Viral', "Non-infected controls", 'Healthy controls')

pc_norm <- pca(log2(counts_norm+1), metadata = data.frame(colData(DES.obj.group)), scale = F)

# save PCA plots 
pal <- viridis(5)
pc_group <- biplot(pc_norm, 
       lab = NULL,
       colkey=c('COVID-19' = pal[1], 
                'Bacterial' = pal[2], 
                'Viral' = pal[3],
                "Non-infected controls" = pal[4], 
                "Healthy controls" = pal[5]),
       colby = 'group_pca',
       legendPosition = 'bottom', 
       colLegendTitle = "Disease group")

pc_group_3_4 <- biplot(pc_norm, 
                   lab = NULL,
                   colkey=c('COVID-19' = pal[1], 
                            'Bacterial' = pal[2], 
                            'Viral' = pal[3],
                            "Non-infected controls" = pal[4], 
                            "Healthy controls" = pal[5]),
                   x = "PC3",
                   y = "PC4",
                   colby = 'group_pca', 
                   legendPosition = 'bottom', 
                   colLegendTitle = "Disease group")

pc_batch <- biplot(pc_norm, 
                   lab = NULL,
                   colkey = c("A" = '#c51b7d', 
                              'B' = '#7fbc41'),
                   colby = 'batch_pca', 
                   legendPosition = 'bottom', 
                   colLegendTitle = "Sequencing batch")

pc_batch_3_4 <- biplot(pc_norm, 
                   lab = NULL,
                   x = "PC3", 
                   y = "PC4",
                   colkey = c("A" = '#c51b7d', 
                              'B' = '#7fbc41'),
                   colby = 'batch_pca', 
                   legendPosition = 'bottom', 
                   colLegendTitle = "Sequencing batch")

pc_age <- biplot(pc_norm, 
                   lab = NULL,
                 x = "PC1",
                 y = "PC2",
                   colby = 'Age', 
                   legendPosition = 'bottom', 
                   colLegendTitle = "Age")

pc_age_3_4 <- biplot(pc_norm, 
                 lab = NULL,
                 x = "PC3",
                 y = "PC4",
                 colby = 'Age', 
                 legendPosition = 'bottom', 
                 colLegendTitle = "Age")

pc_sex <- biplot(pc_norm, 
                 lab = NULL,
                 colkey = c('Male' = '#a63603', 
                            'Female' = '#fd8d3c'),
                 colby = 'sex_pca', 
                 legendPosition = 'bottom', 
                 colLegendTitle = "Sex")

pc_sex_3_4 <- biplot(pc_norm, 
                 lab = NULL,
                 x="PC3",
                 y='PC4',
                 colkey = c('Male' = '#a63603', 
                            'Female' = '#fd8d3c'),
                 colby = 'sex_pca', 
                 legendPosition = 'bottom', 
                 colLegendTitle = "Sex")

pdf(file = 'figures/PCA_plots.pdf', height = 22, width = 8)
grid.arrange(pc_group, pc_group_3_4,
             pc_batch, pc_batch_3_4, 
             pc_age, pc_age_3_4, 
             pc_sex, pc_sex_3_4, ncol = 2)
dev.off()
pdf(file = 'figures/justgroup_pca.pdf', height = 10, width = 11)
pc_group
dev.off()


biplot(pc_norm, 
       lab = NULL,
       x = "PC1",
       y = "PC2",
       colby = 'severity_grouped', 
       legendPosition = 'bottom', 
       colLegendTitle = "Severity")


# take out the results for all of the comparisons 

#  take out the results from DESeq, comparing COVID19 vs DB
COVID.DB.res <- results(DES.obj.group, 
                        independentFiltering=TRUE,
                        alpha=0.05, 
                        pAdjustMethod="BH", 
                        parallel=FALSE, 
                        contrast=list(c("COVID19"), c("DB")))

COVID.DB.res <- lfcShrink(DES.obj.group, 
                          res = COVID.DB.res, 
                          type = 'ashr')

# COVID19 vs DV
COVID.DV.res <- results(DES.obj.group, 
                        independentFiltering=TRUE,
                        alpha=0.05, 
                        pAdjustMethod="BH", 
                        parallel=FALSE, 
                        contrast=list(c("COVID19"), c("DV")))

COVID.DV.res <- lfcShrink(DES.obj.group, 
                          res = COVID.DV.res, 
                          type = 'ashr')

# COVID19 vs control
COVID.Control.res <- results(DES.obj.group, 
                             independentFiltering=TRUE,
                             alpha=0.05, 
                             pAdjustMethod="BH", 
                             parallel=FALSE, 
                             contrast=list(c("COVID19"), c("Control")))

COVID.Control.res <- lfcShrink(DES.obj.group, 
                               res = COVID.Control.res, 
                               type= 'ashr')

# COVID vs DB/DV combined 
COVID.BV.res <- results(DES.obj.group, 
                        independentFiltering=TRUE,
                        alpha=0.05, 
                        pAdjustMethod="BH", 
                        parallel=FALSE, 
                        contrast=list(c("COVID19"), c("DB", "DV")),
                        listValues=c(1, -1/2))

COVID.BV.res <- lfcShrink(DES.obj.group, 
                          res = COVID.BV.res, 
                          type= 'ashr')


# COVID vs DB/DV/control combined 
COVID.C.BV.res <- results(DES.obj.group, 
                          independentFiltering=TRUE,
                          alpha=0.05, 
                          pAdjustMethod="BH", 
                          parallel=FALSE, 
                          contrast=list(c("COVID19"), c("DB", "DV", "Control")),
                          listValues=c(1, -1/3))

COVID.C.BV.res <- lfcShrink(DES.obj.group, 
                            res = COVID.C.BV.res, 
                            type= 'ashr')

# process the results and save to file 
COVID.DB.res <- deseq.results(COVID.DB.res, gene.names, 'Results/DEGs_COVID_DB')
COVID.DV.res <- deseq.results(COVID.DV.res, gene.names, 'Results/DEGs_COVID_DV')
COVID.Con.res <- deseq.results(COVID.Control.res, gene.names, 'Results/DEGs_COVID_Control')
COVID.BV.res <- deseq.results(COVID.BV.res, gene.names, 'Results/DEGs_COVID_BV')
COVID.C.BV.res <- deseq.results(COVID.C.BV.res, gene.names, 'Results/DEGs_COVID_Control_BV')

# infected vs not infected
infected.not.infected.res <- results(DES.obj.group, 
                          independentFiltering=TRUE,
                          alpha=0.05, 
                          pAdjustMethod="BH", 
                          parallel=FALSE, 
                          contrast=list(c("Control"), c("DB", "DV", "COVID19")),
                          listValues=c(1, -1/3))

infected.not.infected.res <- lfcShrink(DES.obj.group, 
                            res = infected.not.infected.res, 
                            type= 'ashr')
infected.not.infected.res <- deseq.results(infected.not.infected.res, gene.names, 'Results/DEGs_Inf_Not_Inf')

table(infected.not.infected.res$Significant)
dim(infected.not.infected.res[infected.not.infected.res$Significant==1 & infected.not.infected.res$log2FoldChange>0,])
dim(infected.not.infected.res[infected.not.infected.res$Significant==1 & infected.not.infected.res$log2FoldChange<0,])

top_hits <- infected.not.infected.res[infected.not.infected.res$padj < 0.0001 & 
                                abs(infected.not.infected.res$log2FoldChange) > 1.5,]

write.csv(file = 'figures/Not_infected_vs_Infected_genes.csv', top_hits)

dim(top_hits[top_hits$log2FoldChange > 0,])
dim(top_hits[top_hits$log2FoldChange < 0,])

top_hits[top_hits$gene %in% c("SLC24A2", "CD1E"),]

# Create and save volcano plot 
# COVID vs DB
pdf(file = 'figures/gVolcano_DB_COVID.pdf', height = 10, width = 10)
ggplot(COVID.DB.res, aes(x = log2FoldChange, y = -log10(padj), color = Color))+
  geom_point(size = 2)+
  theme_bw()+
  scale_x_continuous(limits = c(-5.9, 5.9))+
  scale_color_manual(values = c("black", 'gold', 'red'), 
                     labels = c('NS',
                                "P-value", 
                                expression('P-value and log'[2]*'FC')))+
  labs(y= "-log10 adjusted p-value", 
       color = "")+
  ggtitle("COVID-19 vs Bacterial")+
  geom_text_repel(
    data = subset(COVID.DB.res, -log10(padj) > 25 |
                    abs(log2FoldChange) > 5),
    aes(label = gene), 
    color = 'black', 
    size = 3,
    box.padding = unit(0.25, "lines"),
    point.padding = unit(0.25, "lines"))+
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed')+
  geom_vline(xintercept = 1, linetype = 'dashed')+
  geom_vline(xintercept = -1, linetype = 'dashed')+
  theme(legend.text.align = 0, 
        legend.position = 'bottom', 
        plot.title = element_text(hjust = 0.5))
dev.off()

table(COVID.DB.res$Significant)
dim(COVID.DB.res[COVID.DB.res$log2FoldChange > 0 & COVID.DB.res$Significant == 1,])
dim(COVID.DB.res[COVID.DB.res$log2FoldChange < 0 & COVID.DB.res$Significant == 1,])

# COVID vs DV
pdf(file = 'figures/gVolcano_DV_COVID.pdf', height = 10, width = 10)
ggplot(COVID.DV.res, aes(x = log2FoldChange, y = -log10(padj), color = Color))+
  geom_point(size = 2)+
  theme_bw()+
  scale_x_continuous(limits = c(-3.5, 3.5))+
  scale_color_manual(values = c("black", 'gold', 'red'), 
                     labels = c('NS',
                                "P-value", 
                                expression('P-value and log'[2]*'FC')))+
  labs(y= "-log10 adjusted p-value", 
       color = "")+
  ggtitle("COVID-19 vs viral")+
  geom_text_repel(
    data = subset(COVID.DV.res, -log10(padj) > 4 & abs(log2FoldChange) > 1.5 |
                    abs(log2FoldChange) > 2 |
                    log2FoldChange > 1.5 & -log10(padj) > 5),
    aes(label = gene), 
    color = 'black', 
    size = 3,
    box.padding = unit(0.25, "lines"),
    point.padding = unit(0.25, "lines"))+
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed')+
  geom_vline(xintercept = 1, linetype = 'dashed')+
  geom_vline(xintercept = -1, linetype = 'dashed')+
  theme(legend.text.align = 0, 
        legend.position = 'bottom', 
        plot.title = element_text(hjust = 0.5))
dev.off()

table(COVID.DV.res$Significant)
dim(COVID.DV.res[COVID.DV.res$log2FoldChange > 0 & COVID.DV.res$Significant == 1,])
dim(COVID.DV.res[COVID.DV.res$log2FoldChange < 0 & COVID.DV.res$Significant == 1,])

# COVID vs Control
pdf(file = 'figures/gVolcano_Control_COVID.pdf', height = 10, width = 10)
ggplot(COVID.Con.res, aes(x = log2FoldChange, y = -log10(padj), color = Color))+
  geom_point(size = 2)+
  theme_bw()+
  scale_x_continuous(limits = c(-5.6, 5.6))+
  scale_color_manual(values = c("black", 'gold', 'red'), 
                     labels = c('NS',
                                "P-value", 
                                expression('P-value and log'[2]*'FC')))+
  labs(y= "-log10 adjusted p-value", 
       color = "")+
  ggtitle("COVID-19 vs Control")+
  geom_text_repel(
    data = subset(COVID.Con.res, -log10(padj) > 20 & abs(log2FoldChange) > 2 |
                    abs(log2FoldChange) > 3.5),
    aes(label = gene), 
    color = 'black', 
    size = 3,
    box.padding = unit(0.25, "lines"),
    point.padding = unit(0.25, "lines"))+
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed')+
  geom_vline(xintercept = 1, linetype = 'dashed')+
  geom_vline(xintercept = -1, linetype = 'dashed')+
  theme(legend.text.align = 0, 
        legend.position = 'bottom', 
        plot.title = element_text(hjust = 0.5))
dev.off()

table(COVID.Con.res$Significant)
dim(COVID.Con.res[COVID.Con.res$log2FoldChange > 0 & COVID.Con.res$Significant == 1,])
dim(COVID.Con.res[COVID.Con.res$log2FoldChange < 0 & COVID.Con.res$Significant == 1,])

# COVID vs DB+DV
pdf(file = 'figures/gVolcano_BV_COVID.pdf', height = 10, width = 10)
ggplot(COVID.BV.res, aes(x = log2FoldChange, y = -log10(padj), color = Color))+
  geom_point(size = 2)+
  theme_bw()+
  scale_x_continuous(limits = c(-3,3))+
  scale_color_manual(values = c("black", 'gold', 'red'), 
                     labels = c('NS',
                                "P-value", 
                                expression('P-value and log'[2]*'FC')))+
  labs(y= "-log10 adjusted p-value", 
       color = "")+
  ggtitle("COVID-19 vs bacterial + viral")+
  geom_text_repel(
    data = subset(COVID.BV.res, -log10(padj) > 10 & abs(log2FoldChange)>1.5|
                    -log10(padj) > 15 ),
    aes(label = gene), 
    color = 'black', 
    size = 3,
    box.padding = unit(0.25, "lines"),
    point.padding = unit(0.25, "lines"))+
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed')+
  geom_vline(xintercept = 1, linetype = 'dashed')+
  geom_vline(xintercept = -1, linetype = 'dashed')+
  theme(legend.text.align = 0, 
        legend.position = 'bottom', 
        plot.title = element_text(hjust = 0.5))
dev.off()

table(COVID.BV.res$Significant)
dim(COVID.BV.res[COVID.BV.res$log2FoldChange > 0 & COVID.BV.res$Significant == 1,])
dim(COVID.BV.res[COVID.BV.res$log2FoldChange < 0 & COVID.BV.res$Significant == 1,])


# COVID vs DB+DV+Control
pdf(file = 'figures/gVolcano_BV_Control_COVID.pdf', height = 10, width = 10)
ggplot(COVID.C.BV.res, aes(x = log2FoldChange, y = -log10(padj), color = Color))+
  geom_point(size = 2)+
  theme_bw()+
  scale_x_continuous(limits = c(-4, 4))+
  scale_color_manual(values = c("black", 'gold', 'red'), 
                     labels = c('NS',
                                "P-value", 
                                expression('P-value and log'[2]*'FC')))+
  labs(y= "-log10 adjusted p-value", 
       color = "")+
  ggtitle("COVID-19 vs bacterial + viral + control")+
  geom_text_repel(
    data = subset(COVID.C.BV.res, -log10(padj) > 10 & abs(log2FoldChange) > 2|
                    abs(log2FoldChange) > 3 ),
    aes(label = gene), 
    color = 'black', 
    size = 3,
    box.padding = unit(0.25, "lines"),
    point.padding = unit(0.25, "lines"))+
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed')+
  geom_vline(xintercept = 1, linetype = 'dashed')+
  geom_vline(xintercept = -1, linetype = 'dashed')+
  theme(legend.text.align = 0, 
        legend.position = 'bottom', 
        plot.title = element_text(hjust = 0.5))
dev.off()

table(COVID.C.BV.res$Significant)
dim(COVID.C.BV.res[COVID.C.BV.res$log2FoldChange > 0 & COVID.C.BV.res$Significant == 1,])
dim(COVID.C.BV.res[COVID.C.BV.res$log2FoldChange < 0 & COVID.C.BV.res$Significant == 1,])

# 



######## make boxplots of the 10 genes in the signature 
genes_signature <- unique(c("OASL", "UPB1", 'ZNF684', 'IL1RN', 'ENTPD7', 'NFKBIE', 'CDKN1C',
                     'CD44', 'OTOF', 'MSR1', 'UPB1'))

genes_signature <- data.frame(gene_names = genes_signature,
                              id = gene.names$ensembl_gene_id[match(genes_signature,
                                                                    gene.names$external_gene_name)])

counts_sig <- counts(DES.obj.group, normalized = T)
counts_sig <- counts_sig[genes_signature$id,]
counts_sig <- data.frame(t(counts_sig))
counts_sig$group <- DES.obj.group$Group[match(rownames(counts_sig), 
                                              DES.obj.group$ID)]

counts_sig$group <- as.character(counts_sig$group)
counts_sig$group[counts_sig$group=="DB"] <- "Bacterial"
counts_sig$group[counts_sig$group=="DV"] <- "Viral"
counts_sig$group[counts_sig$group=="HC"] <- "Healthy controls"
counts_sig$group[counts_sig$group=="Control"] <- "Not infected controls"

comparisons <- list(c('COVID19', 'Bacterial'), 
                    c('COVID19', 'Viral'), 
                    c('COVID19', 'Not infected controls'))

for(i in 1:10){
  df <- data.frame(group = counts_sig$group, 
                   gene = counts_sig[,i])

  name <- genes_signature$gene_names[genes_signature$id == colnames(counts_sig)[i]]
  path <- paste('figures//', name, '.pdf', sep = '')
  p <- ggplot(df, aes(x = fct_relevel(group, "COVID19", "Bacterial", "Viral", "Not infected controls"), 
                      y = log2(gene+1), 
                      fill = fct_relevel(group, "COVID19", "Bacterial", "Viral", "Not infected controls")))+
    theme_bw()+
    geom_boxplot(outlier.shape = NA)+
    geom_jitter(position = position_jitterdodge())+
    ggtitle(name)+
    theme(legend.position = 'none')+
    scale_fill_manual(values = c("#440154FF", "#3B528BFF", "#21908CFF", "#5DC863FF", 'darkgrey'))+
    labs(x = "", y = paste("Log2 transformed, normalised counts - ", name, sep = ''))+
    stat_compare_means(comparisons = comparisons, label = 'p.signif')+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  pdf(file = path, height = 6, width = 5)
  plot(p)
  dev.off()
}

#########################################################################################################
# FSPLS 
#########################################################################################################

# use list of genes from DESeq 
# only need one set of normalized data 
exp.values <- counts(DES.obj.group, normalized = TRUE)

# filter so samples included are DB, DV, COVID, Control
exp.values <- exp.values[,colnames(exp.values) %in% DES.obj.group$ID[DES.obj.group$Group %in% c("COVID19", 'DB', 'DV', 'Control')]]
group <- droplevels(DES.obj.group$Group[match(colnames(exp.values), DES.obj.group$ID)])


norm.bv <- log2(exp.values+1)
merged.subset <- merged_meta[match(colnames(norm.bv), merged_meta$ID),]

# check PCA 
pca.new <- pca(norm.bv, metadata = merged.subset)

biplot(pca.new, 
       lab = NULL, 
       colby = 'Plate.Batch', 
       legendPosition = 'right', 
       colLegendTitle = "Plate.Batch")

pdf(file = "figures/gFSPLS_PCA.pdf", height = 7, width = 8)
biplot(pca.new, 
       lab = NULL, 
       colby = 'Group', 
       legendPosition = 'right', 
       colLegendTitle = "Group")
dev.off()

path.fspls <- paste("fspls_lars_multi_meta.R")

##########################################################################################
# FSPLS for COVID vs DB
##########################################################################################

# filter so only DB DV and COVID and also the genes of interest 
normalized.db <- norm.bv[rownames(norm.bv) %in% rownames(COVID.DB.res)[COVID.DB.res$padj < 0.01 & 
                                                                         abs(COVID.DB.res$log2FoldChange) > 1 & 
                                                                         COVID.DB.res$baseMean > 150],]
normalized.db <- normalized.db[,colnames(normalized.db) %in% 
                                 merged.subset$ID[merged.subset$Group %in% 
                                                    c("COVID19", "DB")]]

response.db.covid <- merged.subset$Group[match(colnames(normalized.db), 
                                               merged.subset$ID)]

response.db.covid <- class.ind(response.db.covid)[,1]

# run FS-PLS 
fspls_db_once <- fspls.binomial(path.fspls = path.fspls, 
                                 expression = t(normalized.db), 
                                 response = response.db.covid, 
                                 p.thresh = 0.001, 
                                 max = 2, 
                                 beam = 10, 
                                 split = F)


top_db_binomial <- fspls_db_once$`3`[[1]]$beta_l[which(fspls_db_once$`3`[[1]]$traineval[,'auc'] == max(fspls_db_once$`3`[[1]]$traineval[,'auc']))]

fspls_db_once$`3`[[1]]$traineval[which(fspls_db_once$`3`[[1]]$traineval[,'auc'] == max(fspls_db_once$`3`[[1]]$traineval[,'auc'])),]

gene.names$external_gene_name[match(rownames(top_db_binomial[[1]]$train), gene.names$ensembl_gene_id)]

##########################################################################################
# FSPLS for COVID vs DV
##########################################################################################

# filter so only DB DV and COVID and also the genes of interest
normalized.dv <- norm.bv[rownames(norm.bv) %in% rownames(COVID.DV.res)[COVID.DV.res$padj < 0.01 &
                                                                         abs(COVID.DV.res$log2FoldChange) > 1 & 
                                                                         COVID.DV.res$baseMean > 150],]
normalized.dv <- normalized.dv[,colnames(normalized.dv) %in% 
                                 merged.subset$ID[merged.subset$Group %in% 
                                                    c("COVID19", "DV")]]

response.dv.covid <- merged.subset$Group[match(colnames(normalized.dv), 
                                               merged.subset$ID)]

response.dv.covid <- class.ind(response.dv.covid)[,1]


# run FS-PLS 
fspls_dv_once <- fspls.binomial(path.fspls = path.fspls, 
                                expression = t(normalized.dv), 
                                response = response.dv.covid, 
                                p.thresh = 0.001, 
                                max = 2, 
                                beam = 10, 
                                split = F)

top_dv_binomial <- fspls_dv_once$`3`[[1]]$beta_l[which(fspls_dv_once$`3`[[1]]$traineval[,'auc'] == max(fspls_dv_once$`3`[[1]]$traineval[,'auc']))]

fspls_dv_once$`3`[[1]]$traineval[which(fspls_dv_once$`3`[[1]]$traineval[,'auc'] == max(fspls_dv_once$`3`[[1]]$traineval[,'auc'])),]

gene.names$external_gene_name[match(rownames(top_dv_binomial[[1]]$train), gene.names$ensembl_gene_id)]

########################################################
# COVID vs DB/DV combined
########################################################

# filter so only DB DV and COVID and also the genes of interest
normalized.bv <- norm.bv[rownames(norm.bv) %in% rownames(COVID.BV.res)[COVID.BV.res$padj < 0.01 & 
                                                                         abs(COVID.BV.res$log2FoldChange) > 1 & 
                                                                         COVID.BV.res$baseMean > 150],]
normalized.bv <- normalized.bv[,colnames(normalized.bv) %in% 
                                 merged.subset$ID[merged.subset$Group %in% 
                                                    c("COVID19", "DB", "DV")]]

response.bv.covid <- merged.subset$Group[match(colnames(normalized.bv), 
                                               merged.subset$ID)]

response.bv.covid[response.bv.covid %in% c("DB", "DV")] <- "BV"

response.bv.covid <- class.ind(response.bv.covid)[,2]

# run FSPLS
fspls_bv_once <- fspls.binomial(path.fspls = path.fspls, 
                               expression = t(normalized.bv), 
                               response = response.bv.covid, 
                               p.thresh = 0.001, 
                               max = 4, 
                               beam = 20, 
                               split = F)

top_bv_binomial <- fspls_bv_once$`5`[[1]]$beta_l[which(fspls_bv_once$`5`[[1]]$traineval[,'auc'] == max(fspls_bv_once$`5`[[1]]$traineval[,'auc']))]

gene.names$external_gene_name[match(rownames(top_bv_binomial[[1]]$train), gene.names$ensembl_gene_id)]

fspls_bv_once$`5`[[1]]$traineval[which(fspls_bv_once$`5`[[1]]$traineval[,'auc'] == max(fspls_bv_once$`5`[[1]]$traineval[,'auc'])),]

COVID.BV.res$log2FoldChange[match(rownames(top_bv_binomial[[1]]$train), rownames(COVID.BV.res))]


##########################################################################################
# FSPLS for COVID vs Control + DB + DV
##########################################################################################


# filter so only DB DV CONTROL and COVID and also the genes of interest
normalized.c.bv <- norm.bv[rownames(norm.bv) %in% rownames(COVID.C.BV.res)[COVID.C.BV.res$padj < 0.01 & 
                                                                             abs(COVID.C.BV.res$log2FoldChange) > 1 & 
                                                                             COVID.C.BV.res$baseMean > 150],]
normalized.c.bv <- normalized.c.bv[,colnames(normalized.c.bv) %in% 
                                     merged.subset$ID[merged.subset$Group %in% 
                                                        c("COVID19", "Control", "DB", "DV")]]

response.c.bv.covid <- merged.subset$Group[match(colnames(normalized.c.bv), 
                                                 merged.subset$ID)]

response.c.bv.covid <- class.ind(response.c.bv.covid)[,2]


# run FS-PLS once to check if it is the same signature
fspls_bvc_once <- fspls.binomial(path.fspls = path.fspls, 
                                expression = t(normalized.c.bv), 
                                response = response.c.bv.covid, 
                                p.thresh = 0.001, 
                                max = 4, 
                                beam = 10, 
                                split = F)

top_cbv_binomial <- fspls_bvc_once$`5`[[1]]$beta_l[which(fspls_bvc_once$`5`[[1]]$traineval[,'auc'] == max(fspls_bvc_once$`5`[[1]]$traineval[,'auc']))]
gene.names$external_gene_name[match(rownames(top_cbv_binomial[[1]]$train), gene.names$ensembl_gene_id)]
fspls_bvc_once$`5`[[1]]$traineval[which(fspls_bvc_once$`5`[[1]]$traineval[,'auc'] == max(fspls_bvc_once$`5`[[1]]$traineval[,'auc'])),]

COVID.C.BV.res$log2FoldChange[match(rownames(top_cbv_binomial[[1]]$train), rownames(COVID.C.BV.res))]


##########################################################################################
# Visualisation and evaluation of signatures
##########################################################################################


# create a signature combining the top signatures from the single runs of FS-PLS 
# COVID vs DB

# for LFC threshold of 0.5 there are 2 signatures equally good (one has higher upper AUC, select that)

db_sig <- fspls_db_once$`3`[[1]]$beta_l[which(fspls_db_once$`3`[[1]]$traineval[,'auc'] == max(fspls_db_once$`3`[[1]]$traineval[,'auc']))][[1]]$train
# COVID vs DV
dv_sig <- fspls_dv_once$`3`[[1]]$beta_l[which(fspls_dv_once$`3`[[1]]$traineval[,'auc'] == max(fspls_dv_once$`3`[[1]]$traineval[,'auc']))][[1]]$train
# COVID vs DB+DV 
bv_sig <- fspls_bv_once$`5`[[1]]$beta_l[which(fspls_bv_once$`5`[[1]]$traineval[,'auc'] == max(fspls_bv_once$`5`[[1]]$traineval[,'auc']))][[1]]$train
# COVID vs DB+DV+C
bvc_sig <- fspls_bvc_once$`5`[[1]]$beta_l[which(fspls_bvc_once$`5`[[1]]$traineval[,'auc'] == max(fspls_bvc_once$`5`[[1]]$traineval[,'auc']))][[1]]$train

# plot each signature separetely 
db_exp_sig <- data.frame(Group = merged.subset$Group[match(colnames(norm.bv), merged.subset$ID)],
                                          t(norm.bv[rownames(db_sig),]))
db_exp_sig <- db_exp_sig[db_exp_sig$Group %in% c("COVID19", 'DB'),]

dv_exp_sig <- data.frame(Group = merged.subset$Group[match(colnames(norm.bv), merged.subset$ID)],
                         t(norm.bv[rownames(dv_sig),]))
dv_exp_sig <- dv_exp_sig[dv_exp_sig$Group %in% c("COVID19", 'DV'),]

bv_exp_sig <- data.frame(Group = merged.subset$Group[match(colnames(norm.bv), merged.subset$ID)],
                         t(norm.bv[rownames(bv_sig),]))
bv_exp_sig <- bv_exp_sig[bv_exp_sig$Group %in% c("COVID19", 'DB', 'DV'),]

bvc_exp_sig <- data.frame(Group = merged.subset$Group[match(colnames(norm.bv), merged.subset$ID)],
                         t(norm.bv[rownames(bvc_sig),]))
bvc_exp_sig <- bvc_exp_sig[bvc_exp_sig$Group %in% c("COVID19", 'DB', 'DV', 'Control'),]


db_exp_sig$resp <-  as.numeric(db_exp_sig$Group  == "COVID19")
db_2 <- glm(resp ~ ., data = db_exp_sig[,2:4], family = "binomial")

dv_exp_sig$resp <-  as.numeric(dv_exp_sig$Group  == "COVID19")
dv_2 <- glm(resp ~ ., data = dv_exp_sig[,2:4], family = "binomial")

bv_exp_sig$resp <-  as.numeric(bv_exp_sig$Group  == "COVID19")
bv_4 <- glm(resp ~ ., data = bv_exp_sig[,2:6], family = "binomial")

bvc_exp_sig$resp <-  as.numeric(bvc_exp_sig$Group  == "COVID19")
bvc_4 <- glm(resp ~ ., data = bvc_exp_sig[,2:6], family = "binomial")

# calculate DRS
db_exp_sig$DRS <- predict(db_2, newdata  = db_exp_sig , type = "response")
dv_exp_sig$DRS <- predict(dv_2, newdata  = dv_exp_sig , type = "response")
bv_exp_sig$DRS <- predict(bv_4, newdata  = bv_exp_sig , type = "response")
bvc_exp_sig$DRS <- predict(bvc_4, newdata  = bvc_exp_sig , type = "response")

db_roc <- roc(db_exp_sig$resp, db_exp_sig$DRS, ci = T, print.auc = T)
dv_roc <- roc(dv_exp_sig$resp, dv_exp_sig$DRS, ci = T, print.auc = T)
bv_roc <- roc(bv_exp_sig$resp, bv_exp_sig$DRS, ci = T, print.auc = T)
bvc_roc <- roc(bvc_exp_sig$resp, bvc_exp_sig$DRS, ci = T, print.auc = T)


db_roc_plot <- roc_plot(db_roc, "COVID-19 vs. bacterial", '#EBB424')
dv_roc_plot <- roc_plot(dv_roc, "COVID-19 vs. viral", '#49C39E')
bv_roc_plot <- roc_plot(bv_roc, "COVID-19 vs. bacterial + viral", '#E51670')
bvc_roc_plot <- roc_plot(bvc_roc, "COVID-19 vs. bacterial, viral + non-infected", '#3F0F3F')

pdf(file = 'figures/discovery_ROC_plots.pdf', height = 8, width = 8)
grid.arrange(bvc_roc_plot, 
             bv_roc_plot,
             db_roc_plot, 
             dv_roc_plot, 
             ncol = 2)
dev.off()





# extract expression values for these genes 
expression_pooled_signature <- data.frame(Group = merged.subset$Group[match(colnames(norm.bv), merged.subset$ID)],
                                          t(norm.bv[unique(c(rownames(db_sig),
                                                rownames(dv_sig),
                                                rownames(bv_sig),
                                                rownames(bvc_sig))),]))

gene.names$external_gene_name[match(colnames(expression_pooled_signature)[2:ncol(expression_pooled_signature)], gene.names$ensembl_gene_id)]


expression_pooled_signature$Response_COVID <- 0
expression_pooled_signature$Response_COVID[expression_pooled_signature$Group=="COVID19"] <- 1
genes.model <-  c("ENSG00000117010", "ENSG00000115155", "ENSG00000136689",
                  "ENSG00000146232", "ENSG00000038945", "ENSG00000198018",
                  "ENSG00000129757", "ENSG00000026508", "ENSG00000135114",
                  "ENSG00000100024")

model_10gene <- glm(Response_COVID ~., data = expression_pooled_signature[,c("Response_COVID",genes.model)], family ="binomial")
beta_rest <- summary(model_10gene)$coefficients
beta_rest <- data.frame(beta_rest[order(beta_rest[,4]),][,'Estimate'])

check_combinations <- function(bfmod, num){
  com <- combn(genes.model,num)
  for(i in 1:ncol(com)){
    upgenes <- com[,i]
    downgenes <- genes.model[!genes.model %in% upgenes]
    if(length(upgenes)==1){
      vals <- expression_pooled_signature[,upgenes]
    }else if(length(upgenes)>1){
      vals <- rowSums(expression_pooled_signature[,upgenes])
    }else{
      vals <- rep(0, nrow(expression_pooled_signature))
    }
    if(length(downgenes)==1){
      vals <- vals - expression_pooled_signature[,downgenes]
    }else if(length(downgenes)>1){
      vals <- vals - rowSums(expression_pooled_signature[,downgenes])
    }else{
      vals <- vals - rep(0, nrow(expression_pooled_signature))
    }
    aucv <- auc(roc(expression_pooled_signature$Response_COVID, vals, quiet =T))
    if(aucv > bfmod$a){
      bfmod$genes <- upgenes
      bfmod$a <- aucv
    }
  }
  bfmod
}

bfmod <- list(genes = c(), a = 0.5)
for(i in 1:10){
  print(i)
  bfmod <- check_combinations(bfmod,i)
}

beta_rest[bfmod$genes,"direction"] <- "Up"
beta_rest[genes.model[!genes.model %in% bfmod$genes],"direction"] <- "Down"

# save the betas 
beta_rest$ID <- gene.names$external_gene_name[match(rownames(beta_rest), gene.names$ensembl_gene_id)]
write.csv(file = 'betas_discovery.csv', beta_rest)

# evaluate performance of 10 genes combined 
expression_pooled_signature$DRS <- predict(model_10gene , expression_pooled_signature[,genes.model], type = "response")

ten_roc_all <- roc(expression_pooled_signature$Response_COVID, 
    expression_pooled_signature$DRS, plot = T, ci = T)

pdf(file = 'figures/10gene_discovery.pdf', height = 4, width = 4)
db_roc_plot <- roc_plot(ten_roc_all, "10-gene signature performance", 'darkgrey')
dev.off()


# COVID vs rest AUC 
covid_rest <- roc(expression_pooled_signature$Response_COVID, expression_pooled_signature$DRS, plot = TRUE, print.auc = TRUE, ci = TRUE, 
    main = "10-gene signature tested in COVID-19 vs. combined")

covid_db <- roc(expression_pooled_signature$Response_COVID[expression_pooled_signature$Group%in% c('COVID19', 'DB')],
    expression_pooled_signature$DRS[expression_pooled_signature$Group %in% c('COVID19', 'DB')],
    plot = TRUE, print.auc = TRUE, ci = TRUE,
    main = "10-gene signature tested in COVID-19 vs DB")

covid_dv <- roc(expression_pooled_signature$Response_COVID[expression_pooled_signature$Group %in% c('COVID19', 'DV')],
                expression_pooled_signature$DRS[expression_pooled_signature$Group %in% c('COVID19', 'DV')],
                plot = TRUE, print.auc = TRUE, ci = TRUE,
                main = "10-gene signature tested in COVID-19 vs DV")

covid_bv <- roc(expression_pooled_signature$Response_COVID[expression_pooled_signature$Group%in% c('COVID19', 'DV', "DB")],
                expression_pooled_signature$DRS[expression_pooled_signature$Group %in% c('COVID19', 'DV', "DB")],
                plot = TRUE, print.auc = TRUE, ci = TRUE,
                main = "10-gene signature tested in COVID-19 vs DB+DV")

covid_c <- roc(expression_pooled_signature$Response_COVID[expression_pooled_signature$Group %in% c('COVID19', 'Control')],
                expression_pooled_signature$DRS[expression_pooled_signature$Group %in% c('COVID19', 'Control')],
                plot = TRUE, print.auc = TRUE, ci = TRUE,
                main = "10-gene signature tested in COVID-19 vs Control")

# overlaid plot 
pdf(file = 'figures/AUC_separate_groups.pdf', height = 6, width = 6)
plot(covid_rest)
plot(ci.sp(covid_rest, sensitivities=seq(0, 1, .01)), type="shape", col = alpha('grey', 0.1))


text(x = 0.4, y = 0.6, 
     paste('COVID-19 vs. combined\n AUC: ',
           round(auc(covid_rest), 3), ' (95% CI: ', 
           round(ci(covid_rest)[[1]], 3), "-",
           round(ci(covid_rest)[[3]], 3), ")", sep=''))
lines(covid_db, col = '#FF495C')
text(x = 0.4, y = 0.45, col = '#FF495C',
     paste('COVID-19 vs. DB\n AUC: ',
           round(auc(covid_db), 3), ' (95% CI: ', 
           round(ci(covid_db)[[1]], 3), "-",
           round(ci(covid_db)[[3]], 3), ")", sep=''))
lines(covid_dv, col = '#256EFF')
text(x = 0.4, y = 0.3, col = '#256EFF',
     paste('COVID-19 vs. DV\n AUC: ',
           round(auc(covid_dv), 3), ' (95% CI: ', 
           round(ci(covid_dv)[[1]], 3), "-",
           round(ci(covid_dv)[[3]], 3), ")", sep=''))
lines(covid_c, col = '#3DDC97')
text(x = 0.4, y = 0.15, col = '#3DDC97',
     paste('COVID-19 vs. Control\n AUC: ',
           round(auc(covid_c), 3), ' (95% CI: ', 
           round(ci(covid_c)[[1]], 3), "-",
           round(ci(covid_c)[[3]], 3), ")", sep=''))
dev.off()


 


