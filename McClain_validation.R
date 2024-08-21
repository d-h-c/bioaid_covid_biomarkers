########################################################################
# 25th Nov 2021
# using the McClain et al 2021 dataset to validate the COVID19 diagnostic gene signature
########################################################################


library(DESeq2)
library(PCAtools)
library(biomaRt)
load.packages()


# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE161731
count_mat <- read.csv(file = gzfile('CountData_Clean.csv.gz'), header = T, row.names = 1)
meta <- read.csv(file = 'Metadata_gendered_clean.csv', header = T, row.names =1)

gene_names <- read.csv("gene_names_ensembl89.csv")
ind <-  str_detect(gene_names$symbol,"RPL") | str_detect(gene_names$symbol,"RPS")


# remove the ribosomal genes
count_mat <- count_mat[!(rownames(count_mat) %in% gene_names[ind,1]),]
dim(count_mat)
# 59k genes 

count_mat <- count_mat[!rownames(count_mat) %in% names(rowSums(count_mat)[rowSums(count_mat) < 10]),]
# 41k genes
meta <- meta[match(colnames(count_mat), rownames(meta)),]
meta$Batch <- as.factor(meta$Batch)

pc <- PCAtools::pca(log2(count_mat+1), metadata = meta, scale = F)

biplot(pc, 
       lab = NULL,
       colby = 'Batch', 
       legendPosition = 'right', 
       colLegendTitle = "Batch")

biplot(pc, 
       lab = NULL,
       colby = 'cohort', 
       showLoadings = T,
       legendPosition = 'right', 
       colLegendTitle = "Disease group")

biplot(pc, 
       lab = NULL,
       shape = 'gender', 
       colby  = 'cohort',
       legendPosition = 'right')

biplot(pc, 
       lab = NULL,
       colby = 'time_since_onset', 
       legendPosition = 'right', 
       colLegendTitle = "time_since_onset")

eigencorplot(pc, metavars = c('cohort', 'Batch', 
                              'gender', 'AGE', 'hospitalized', 
                              'time_since_onset', "RACE"))


pc_b <- pcaMethods::pca(t(log2(count_mat+1)), nPcs = 10, scale = 'none')
el <- pcaMethods:::simpleEllipse(pc_b@scores[,1], pc_b@scores[,2], alfa = 0.99, nrow(pc_b@scores))
table(point.in.polygon(pc_b@scores[,1], pc_b@scores[,2], el[,1], el[,2]))

ggplot(data.frame(pc_b@scores, 
                  group = meta$cohort[match(rownames(pc_b@scores), rownames(meta))]), 
       aes(x = PC1, y = PC2, color = group))+geom_point()+
  annotate("path",
           x=el[,1],
           y=el[,2])

out_b <- meta[rownames(meta) %in% rownames(pc_b@scores)[point.in.polygon(pc_b@scores[,1], pc_b@scores[,2], el[,1], el[,2]) == 0],]

loadings <- pc_b@loadings[,
loadings <- loadings[order(abs(loadings), decreasing = T)]
loadings <- loadings[1:20]
names(loadings) <- gene_list$external_gene_name[match(names(loadings), gene_list$ensembl_gene_id)]

# normalise using DESeq
dds_val <- DESeqDataSetFromMatrix(countData = count_mat,
                                  colData = meta, 
                                  design = ~1)

keep <- rowSums(counts(dds_val) >= 10 ) >= 3
dds_val <- dds_val[keep,]
dds_val <- estimateSizeFactors(dds_val)

normalised_counts <- counts(dds_val, normalized=TRUE)
dim(normalised_counts)
# 31k

pc <- PCAtools::pca(log2(normalised_counts+1), metadata = meta, scale = F)

biplot(pc, 
       lab = NULL,
       showLoadings = T,
       colby = 'time_since_onset', 
       legendPosition = 'right', 
       colLegendTitle = "time_since_onset")

eigencorplot(pc, metavars = c('cohort', 'Batch', 
                              'gender', 'AGE', 'hospitalized', 
                              'time_since_onset', "RACE"))

############################################################################

table(meta$cohort)
genes <- gene_names[gene_names$symbol %in% c("CD44", "OTOF", "MSR1", "UPB1", 
                                                       "ENTPD7", "NFKBIE", "CDKN1C", "OASL", 
                                                       "IL1RN", "ZNF684"),]
expression_genes <- data.frame(group = meta$cohort, 
                               t(log(normalised_counts[genes$gene,]+1)))


colnames(expression_genes)[2:ncol(expression_genes)] <- genes$symbol[match(colnames(expression_genes)[2:ncol(expression_genes)],genes$gene)]
expression_genes$time <- meta$time_since_onset[match(rownames(expression_genes), rownames(meta))]
expression_genes <- expression_genes[!expression_genes$group == 'healthy',]
expression_genes <- expression_genes[!expression_genes$time %in% c('late'),]

expression_genes$covid <- 0
expression_genes$covid[expression_genes$group=="COVID-19"] <- 1

# read in original model weights 
original_weights <- read.csv(file = 'betas_discovery.csv', row.names = 1)
original_weights <- original_weights[rownames(original_weights)!="(Intercept)",]

table(meta$time_since_onset)
table(expression_genes$time)

####################################################################################
# test the performance of the full 10 gene signature 
# add on the simple disease risk score

intercept <- 14.08092
response <- intercept + as.vector(as.matrix(expression_genes[,match(original_weights$ID,colnames(expression_genes))]) %*% as.matrix(original_weights$beta_rest.order.beta_rest...4.........Estimate..))
expression_genes$original_drs  <- exp(response) / (1 + exp(response))


expression_genes$simple_drs <- rowSums(expression_genes[colnames(expression_genes) %in% original_weights$ID[original_weights$direction == 'Up']]) - rowSums(expression_genes[colnames(expression_genes) %in% original_weights$ID[original_weights$direction == 'Down']])



######### look at the early samples first 
expression_genes$time[!expression_genes$time %in% c("early", "middle")] <- 'not_covid'
expression_early <- subset(expression_genes, time != 'middle') 

# performance with simple DRS
early_covid_rest_simple <- roc(expression_early$covid, expression_early$simple_drs, ci  = T)

early_covid_bacterial_simple <- roc(expression_early$covid[expression_early$group %in% c("COVID-19", "Bacterial")], 
                                    expression_early$simple_drs[expression_early$group %in% c("COVID-19", "Bacterial")], ci  = T)

early_covid_viral_simple <- roc(expression_early$covid[expression_early$group %in% c("COVID-19", "CoV other", "Influenza")], 
                                    expression_early$simple_drs[expression_early$group %in% c("COVID-19", "CoV other", "Influenza")], ci  = T)


# retrain weights 
model.mcclain <- glm(covid~., data = expression_early[,c(2:11, 13)], family ="binomial")

# calculate retrained model 
expression_early$retrained_drs <- predict(model.mcclain,expression_early[,2:11], type ="response")

early_covid_rest_retrained <- roc(expression_early$covid, expression_early$retrained_drs, ci  = T)

early_covid_bacterial_retrained <- roc(expression_early$covid[expression_early$group %in% c("COVID-19", "Bacterial")], 
                                    expression_early$retrained_drs[expression_early$group %in% c("COVID-19", "Bacterial")], ci  = T)

early_covid_viral_retrained <- roc(expression_early$covid[expression_early$group %in% c("COVID-19", "CoV other", "Influenza")], 
                                expression_early$retrained_drs[expression_early$group %in% c("COVID-19", "CoV other", "Influenza")], ci  = T)

# roc plots
# covid vs all 
early_c_rest_simple_plot <- roc_plot(early_covid_rest_simple, "Simple disease risk score", 'darkgrey')
early_c_rest_retrained_plot <- roc_plot(early_covid_rest_retrained, "Retrained weights", '#3F0F3F')

covid_rest_plots <- ggarrange(early_c_rest_simple_plot, early_c_rest_retrained_plot, ncol=2, nrow=1, common.legend = TRUE,legend="bottom")
covid_rest_plots <- annotate_figure(covid_rest_plots, top = text_grob("Early COVID-19 vs. bacterial + viral"))

# covid vs bacterial
early_c_bac_simple_plot <- roc_plot(early_covid_bacterial_simple, "Simple disease risk score", 'darkgrey')
early_c_bac_retrained_plot <- roc_plot(early_covid_bacterial_retrained, "Retrained weights", '#EBB424')

covid_bac_plots <- ggarrange(early_c_bac_simple_plot, early_c_bac_retrained_plot, ncol=2, nrow=1, common.legend = TRUE,legend="bottom")
covid_bac_plots <- annotate_figure(covid_bac_plots, top = text_grob("Early COVID-19 vs. bacterial"))

# covid vs viral
early_c_vir_simple_plot <- roc_plot(early_covid_viral_simple, "Simple disease risk score", 'darkgrey')
early_c_vir_retrained_plot <- roc_plot(early_covid_viral_retrained, "Retrained weights", '#49C39E')

covid_vir_plots <- ggarrange(early_c_vir_simple_plot, early_c_vir_retrained_plot, ncol=2, nrow=1, common.legend = TRUE,legend="bottom")
covid_vir_plots <- annotate_figure(covid_vir_plots, top = text_grob("Early COVID-19 vs. viral"))

pdf(file = 'McClain_10_genes_roc.pdf', height = 12, width = 7.5)
grid.arrange(covid_rest_plots, 
             covid_bac_plots, 
             covid_vir_plots, 
             ncol = 1)
dev.off()


############ 10 gene performance in the middle cohort as well. 
expression_genes$time[!expression_genes$time %in% c("early", "middle")] <- 'not_covid'

# performance with simple DRS
all_covid_rest_simple <- roc(expression_genes$covid, expression_genes$simple_drs, ci  = T)

all_covid_bacterial_simple <- roc(expression_genes$covid[expression_genes$group %in% c("COVID-19", "Bacterial")], 
                                    expression_genes$simple_drs[expression_genes$group %in% c("COVID-19", "Bacterial")], ci  = T)

all_covid_viral_simple <- roc(expression_genes$covid[expression_genes$group %in% c("COVID-19", "CoV other", "Influenza")], 
                                expression_genes$simple_drs[expression_genes$group %in% c("COVID-19", "CoV other", "Influenza")], ci  = T)


# retrain weights 
mod_all <- glm(covid~., data = expression_genes[,c(2:11, 13)], family = "binomial")

# calculate retrained model 
expression_genes$retrained_drs <- predict(mod_all,expression_genes[,2:11], type = "response")

all_covid_rest_retrained <- roc(expression_genes$covid, expression_genes$retrained_drs, ci  = T)

all_covid_bacterial_retrained <- roc(expression_genes$covid[expression_genes$group %in% c("COVID-19", "Bacterial")], 
                                     expression_genes$retrained_drs[expression_genes$group %in% c("COVID-19", "Bacterial")], ci  = T)

all_covid_viral_retrained <- roc(expression_genes$covid[expression_genes$group %in% c("COVID-19", "CoV other", "Influenza")], 
                                 expression_genes$retrained_drs[expression_genes$group %in% c("COVID-19", "CoV other", "Influenza")], ci  = T)

# roc plots
# covid vs all 
all_c_rest_simple_plot <- roc_plot(all_covid_rest_simple, "Simple disease risk score", 'darkgrey')
all_c_rest_retrained_plot <- roc_plot(all_covid_rest_retrained, "Retrained weights", '#3F0F3F')

all_covid_rest_plots <- ggarrange(all_c_rest_simple_plot, all_c_rest_retrained_plot, ncol=2, nrow=1, common.legend = TRUE,legend="bottom")
all_covid_rest_plots <- annotate_figure(all_covid_rest_plots, top = text_grob("Early + middle COVID-19 vs. bacterial + viral"))

# covid vs bacterial
all_c_bac_simple_plot <- roc_plot(all_covid_bacterial_simple, "Simple disease risk score", 'darkgrey')
all_c_bac_retrained_plot <- roc_plot(all_covid_bacterial_retrained, "Retrained weights", '#EBB424')

all_covid_bac_plots <- ggarrange(all_c_bac_simple_plot, all_c_bac_retrained_plot, ncol=2, nrow=1, common.legend = TRUE,legend="bottom")
all_covid_bac_plots <- annotate_figure(all_covid_bac_plots, top = text_grob("Early + middle COVID-19 vs. bacterial"))

# covid vs viral
all_c_vir_simple_plot <- roc_plot(all_covid_viral_simple, "Simple disease risk score", 'darkgrey')
all_c_vir_retrained_plot <- roc_plot(all_covid_viral_retrained, "Retrained weights", '#49C39E')

all_covid_vir_plots <- ggarrange(all_c_vir_simple_plot, all_c_vir_retrained_plot, ncol=2, nrow=1, common.legend = TRUE,legend="bottom")
all_covid_vir_plots <- annotate_figure(all_covid_vir_plots, top = text_grob("Early + middle COVID-19 vs. viral"))

pdf(file = 'McClain_10_genes_roc_middle.pdf', height = 12, width = 7.5)
grid.arrange(all_covid_rest_plots, 
             all_covid_bac_plots, 
             all_covid_vir_plots, 
             ncol = 1)
dev.off()

