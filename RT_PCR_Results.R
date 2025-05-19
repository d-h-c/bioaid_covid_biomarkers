###########################################################################################
# 17/08/2021
# Analysis of RT-PCR LAMP results for validation of genes in COVID signature
###########################################################################################
library(tidyverse)
library(viridis)
library(caret)
library(pROC)
library(PCAtools)
library(data.table)
library(ggpubr)
source("functions.R")
library(grid)
library(gridExtra)

########################################################################################

# first load in metadata
meta <- read.csv(file = 'data/validation/RT-PCR_validation_phenotypes.csv', header = T)


###############

# read in CT values 
ct_values <- read.csv(file = 'data/validation/CT_data.csv', header = T)

ct_values$group <- meta$group_short[match(ct_values$Sample, meta$ID)]
ct_values$CT[ct_values$CT==999] <- 40

# what is the range of GAPDH
ct_values$CT[ct_values$Gene=="GAPDH_ref"]
samples.remove  <- read.csv("data/validation/samples_no_gapdh.csv", header = F)[,1]

# some have no GADPH so remove those 
ct_values <- ct_values[!ct_values$Sample %in% samples.remove,]

even_indexes <- seq(2,nrow(ct_values),2)
odd_indexes <- seq(1,nrow(ct_values)-1,2)
ct_values$replicate <- NA
ct_values$replicate[odd_indexes] <- 1
ct_values$replicate[even_indexes] <- 2

ct_values$replicate_gene <- paste(ct_values$Gene, ct_values$replicate, sep = "_")

# make into wide matrix format
ct_wide <- pivot_wider(ct_values, names_from = replicate_gene, values_from = CT, id_cols = Sample)

pdf(file = 'figures/correlation_plots.pdf', height = 5, width = 5)
for(i in 1:length(levels(as.factor(ct_values$Gene)))){
  gene <- levels(as.factor(ct_values$Gene))[i]
 plot(ct_wide[,unique(ct_values$replicate_gene[ct_values$Gene==gene])],
     # col = as.factor(meta$group_short[match(ct_wide$Sample, meta$ID)]),
      ylim = c(min(ct_wide[,unique(ct_values$replicate_gene[ct_values$Gene==gene])]),
               max(ct_wide[,unique(ct_values$replicate_gene[ct_values$Gene==gene])])),
      xlim = c(min(ct_wide[,unique(ct_values$replicate_gene[ct_values$Gene==gene])]),
               max(ct_wide[,unique(ct_values$replicate_gene[ct_values$Gene==gene])])),
      main = paste(levels(as.factor(ct_values$Gene))[i],
                   round(cor(unlist(ct_wide[,unique(ct_values$replicate_gene[ct_values$Gene==gene])][,1]),
                      unlist(ct_wide[,unique(ct_values$replicate_gene[ct_values$Gene==gene])][,2])), 5),
                   sep = ": "))
 cor(ct_wide[,unique(ct_values$replicate_gene[ct_values$Gene==gene])])
# calculate how many samples have 40 as their value
 apply(ct_wide[,unique(ct_values$replicate_gene[ct_values$Gene==gene])], 2, function(y){
   # print(length(y[y==40]))
   #print(sd(y))
 })
 missing_both <- data.frame(missing = apply(ct_wide[,unique(ct_values$replicate_gene[ct_values$Gene==gene])], 1, function(p){
   if(sum(p) == 80){
     return('missing-both')
   }
   else{
     return('not-both')
   }
 }))
 missing_both$group <- ct_values$group[match(ct_wide$Sample, ct_values$Sample)]
 print(table(missing_both$missing))
 print(table(missing_both$missing, missing_both$group))
 }
dev.off()

for(i in 1:length(levels(as.factor(ct_values$Gene)))){
  gene <- levels(as.factor(ct_values$Gene))[i]
  plot(ct_wide[,unique(ct_values$replicate_gene[ct_values$Gene==gene])], 
       ylim = c(min(na.rm = TRUE, ct_wide[,unique(ct_values$replicate_gene[ct_values$Gene==gene])]), 
                max(na.rm = TRUE, ct_wide[,unique(ct_values$replicate_gene[ct_values$Gene==gene])])), 
       xlim = c(min(na.rm = TRUE, ct_wide[,unique(ct_values$replicate_gene[ct_values$Gene==gene])]), 
                max(na.rm = TRUE, ct_wide[,unique(ct_values$replicate_gene[ct_values$Gene==gene])])), 
       main = cor(x = unlist(ct_wide[,unique(ct_values$replicate_gene[ct_values$Gene==gene])][,1]), 
                  y = unlist(ct_wide[,unique(ct_values$replicate_gene[ct_values$Gene==gene])][,2]), 
                  use = "complete.obs"))
  cor(ct_wide[,unique(ct_values$replicate_gene[ct_values$Gene==gene])])
  #calculate how many samples have 40 as their value
  apply(ct_wide[,unique(ct_values$replicate_gene[ct_values$Gene==gene])], 2, function(y){
    print(length(y[y==40]))
  })
}


######################################################################
# logistic regression 
# COVID vs rest 

# with averaged values 

# replace measurements
# if a sample has 2 40s -> NA
# if a sample has 1 40 -> take the other replicate
# if a sample has 0 40s -> take the average
ct_wide <- data.frame(ct_wide)

new_expression <- data.frame(t(apply(ct_wide[,2:ncol(ct_wide)], 1, function(x){
  values_sample <- vector(length=11)
  for(i in 1:length(levels(as.factor(ct_values$Gene)))){
    gene_cols <- unique(ct_values$replicate_gene[ct_values$Gene == levels(as.factor(ct_values$Gene))[i]])
    x_gene <- c(x[gene_cols[1]], x[gene_cols[2]])
    if(x_gene[1]  == 40 & x_gene[2] == 40){
      value <- 40
    }
    else if(x_gene[1] == 40 | x_gene[2] == 40){
     value <- x_gene[!x_gene==40]
    }
    else if(x_gene[1] != 40 & x_gene[2] != 40){
      value <- mean(x_gene)
    }
    values_sample[i] <- value
  }
  return(values_sample)
})))

colnames(new_expression) <- levels(as.factor(ct_values$Gene))
rownames(new_expression) <- ct_wide$Sample

new_expression <- new_expression[,c(1:3, 5:ncol(new_expression))]
ct_values$group <- fct_relevel(as.factor(ct_values$group), "COVID19", "DB", "PB", "DV", "PV", "Other_Infection", "Not_Infected", "Inconclusive")

# boxplot of each 
pal <- viridis(5)

ct_values$comparator_groups <- as.character( ct_values$group )
ct_values$comparator_groups[ct_values$Sample == '2310'] <- "Definite Infections"

ct_values$comparator_groups[ct_values$group %in% c("DB", 'DV', 'Other_Infection', "Mixed")] <- "Definite Infections"
ct_values$comparator_groups[ct_values$comparator_groups %in% c("Inconclusive", "PB", "PV")] <- "Probable Infections"
ct_values$comparator_groups[ct_values$comparator_groups == "Not_Infected"] <- 'Not Infected'

comparisons <- list(c('COVID19', 'Definite Infections'), 
                    c('COVID19', 'Probable Infections'), 
                    c('COVID19', 'Not Infected'))

rownames(new_expression) <- str_remove_all(rownames(new_expression), 'IP\\-')

for(i in 1:10){
  df <- data.frame(group = ct_values$comparator_groups[match(rownames(new_expression), ct_values$Sample)], 
                   gene = new_expression[,i])
  name <- unlist(str_split(colnames(new_expression)[i], "\\_"))[1]
  path <- paste('figures/', name, '.pdf', sep = '')
  p <- ggplot(df, aes(x = group, y = gene, fill = group))+
    theme_bw()+
    geom_boxplot(outlier.shape = NA)+
    geom_jitter(position = position_jitterdodge())+
    ggtitle(name)+
    theme(legend.position = 'none')+
    scale_fill_viridis_d()+
    labs(x = "", y = name)+
    stat_compare_means(comparisons = comparisons, label = 'p.signif')+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  pdf(file = path, height = 6, width = 5)
  plot(p)
  dev.off()
}

comparisons_detail <- list(c("COVID19\nn=37", "Bacterial Infections\nn=27"), 
                           c("COVID19\nn=37", "Viral Infections\nn=4"),
                           c("COVID19\nn=37", "Not Infected\nn=20"), 
                           c("COVID19\nn=37", "Probable+Indeterminate\nn=113"))

ct_values$group_detail <- ct_values$comparator_groups
ct_values$group_detail[ct_values$group == 'DB'] <- "Bacterial Infections"
ct_values$group_detail[ct_values$group == 'DV'] <- "Viral Infections"
ct_values$group_detail[ct_values$group == 'Other_Infection'] <- "Other Infections"
ct_values$group_detail[ct_values$group %in% c('Inconclusive','PV','PB')] <- "Probable+Indeterminate"

table(ct_values$group_detail[match(rownames(new_expression), ct_values$Sample)])
ct_values$group_detail[ct_values$group_detail == 'Bacterial Infections'] <- "Bacterial Infections\nn=27"
ct_values$group_detail[ct_values$group_detail == 'COVID19'] <- "COVID19\nn=37"
ct_values$group_detail[ct_values$group_detail == 'Not Infected'] <- "Not Infected\nn=20"
ct_values$group_detail[ct_values$group_detail == 'Viral Infections'] <- "Viral Infections\nn=4"
ct_values$group_detail[ct_values$group_detail == 'Probable+Indeterminate'] <- "Probable+Indeterminate\nn=113"
ct_values$group_detail[ct_values$group_detail == 'Other Infections'] <- "Other Infections\nn=3"


plotl <- list()
for(i in 1:10){
  df <- data.frame(group = ct_values$group_detail[match(rownames(new_expression), ct_values$Sample)], 
                   gene = new_expression[,i])
  name <- unlist(str_split(colnames(new_expression)[i], "\\_"))[1]
  path <- paste('figures/', name, '.pdf', sep = '')
  p <- ggplot(df, aes(x = fct_relevel(group, "COVID19\nn=37", 
                                      "Bacterial Infections\nn=27",
                                      "Viral Infections\nn=4"), 
                      y = gene, 
                      fill = fct_relevel(group, "COVID19\nn=37", 
                                         "Bacterial Infections\nn=27",
                                         "Viral Infections\nn=4")))+
    theme_bw()+
    geom_boxplot(outlier.shape = NA)+
    geom_jitter(position = position_jitterdodge())+
    ggtitle(name)+
    theme(legend.position = 'none')+
    scale_fill_viridis_d()+
    labs(x = "", y = name)+
    stat_compare_means(comparisons = comparisons_detail, label = 'p.signif')+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  pdf(file = path, height = 5, width = 6)
  plot(p)
  dev.off()
  plotl[[i]] <- p
}

pdf(file = "figures/s6.pdf", height = 25, width = 12)
ggarrange(plotl[[1]], plotl[[2]],plotl[[3]],
          plotl[[4]],plotl[[5]],plotl[[6]],
          plotl[[7]],plotl[[8]],plotl[[9]],plotl[[10]],
          ncol = 2, nrow = 5)
dev.off()



######### first test the performance using the original model weights 
original_weights <- read.csv(file = 'data/discovery/betas_discovery.csv', row.names =1)
original_weights <- original_weights[rownames(original_weights)!="(Intercept)",]
ids <- c('NFKBIE_4590', 'ZNF684_1585', 'OASL_8261', 'UPB1_9684', 
                         'OTOF_0334', 'CDKN1C_6036y', 'CD44_8274', 'MSR1_6209', 
                         'IL1RN_1686', 'ENTPD7_2397')
for(i in ids){original_weights[original_weights$ID == gsub("_.*","",i),"ID"] <- i}

original_weights$beta_reversed <- -1 *original_weights$beta

resp <- as.vector(as.matrix(new_expression[,match(original_weights$ID, colnames(new_expression))])  %*% as.matrix(original_weights$beta_reversed))
new_expression$original_weights_drs <- resp # exp(resp) / (1 + exp(resp)) 


new_expression$group <- ct_values$group[match(rownames(new_expression), ct_values$Sample)]
new_expression$group <- as.character(new_expression$group)
new_expression$group[rownames(new_expression) == '2310'] <- 'Mixed'
new_expression$group <- as.factor(new_expression$group)

new_expression$covid_group <- 0
new_expression$covid_group[new_expression$group=="COVID19"] <- 1
new_expression$covid_group <- as.numeric(new_expression$covid_group)

rownames(new_expression) <- paste("IP-", rownames(new_expression), sep = '')

# first look at the performance in COVID vs definite infections 
covid_definites <- new_expression[new_expression$group %in% c('COVID19', 'Mixed', 'DB', 'DV', 'Other_Infection'),]
covid_definites_roc_original <- roc(covid_definites$covid_group, covid_definites$original_weights_drs, ci = T)

# now look at the performance in COVID vs probable infections
covid_probable <- new_expression[new_expression$group %in% c('COVID19', 'Inconclusive', 'PB', 'PV'),]
covid_probable_roc_original <- roc(covid_probable$covid_group, covid_probable$original_weights_drs, ci = T)

# now look at the performance in COVID vs not infected
covid_not_infected <- new_expression[new_expression$group %in% c('COVID19', 'Not_Infected'),]
covid_not_infected_roc_original <- roc(covid_not_infected$covid_group, covid_not_infected$original_weights_drs, ci = T)

# finally COVID vs all groups
covid_all <- new_expression
covid_all_roc_original <- roc(covid_all$covid_group, covid_all$original_weights_drs, ci = T)

# COVID vs probables and definites 
covid_def_prob <- new_expression[new_expression$group %in% c('COVID19', 'Mixed', 'DB', 'DV', 'Inconclusive', 'PB', 'PV'),]
covid_def_prob_roc_original <- roc(covid_def_prob$covid_group, covid_def_prob$original_weights_drs, ci = T)

##### now calculate simple DRS for all comparisons

# test out the performance if we do a simple DRS including adding up those that increase (decrease here) 
#original_weights$direction <- 'Up'
#original_weights$direction[original_weights$beta_reversed < 0] <- "Down"

# COVID vs definites 
covid_definites$simple_drs <- (rowSums(covid_definites[,colnames(covid_definites) %in% original_weights$ID[original_weights$direction=="Up"]], na.rm = T)-
                                rowSums(covid_definites[,colnames(covid_definites) %in% original_weights$ID[original_weights$direction=="Down"]], na.rm = T))

covid_definites_roc_simple <- roc(covid_definites$covid_group, covid_definites$simple_drs, ci = T)

# COVID vs probables
covid_probable$simple_drs <- (rowSums(covid_probable[,colnames(covid_probable) %in% original_weights$ID[original_weights$direction=="Up"]], na.rm = T)-
                                 rowSums(covid_probable[,colnames(covid_probable) %in% original_weights$ID[original_weights$direction=="Down"]], na.rm = T))

covid_probables_roc_simple <- roc(covid_probable$covid_group, covid_probable$simple_drs, ci = T)

# COVID vs not infected
covid_not_infected$simple_drs <- (rowSums(covid_not_infected[,colnames(covid_not_infected) %in% original_weights$ID[original_weights$direction=="Up"]], na.rm = T)-
                                rowSums(covid_not_infected[,colnames(covid_not_infected) %in% original_weights$ID[original_weights$direction=="Down"]], na.rm = T))

covid_not_infected_roc_simple <- roc(covid_not_infected$covid_group, covid_not_infected$simple_drs, ci = T)

# COVID vs all 
covid_all$simple_drs <- (rowSums(covid_all[,colnames(covid_all) %in% original_weights$ID[original_weights$direction=="Up"]], na.rm = T)-
                                    rowSums(covid_all[,colnames(covid_all) %in% original_weights$ID[original_weights$direction=="Down"]], na.rm = T))

covid_all_roc_simple <- roc(covid_all$covid_group, covid_all$simple_drs, ci = T)

# COVID vs probables and definites 
covid_def_prob$simple_drs <- (rowSums(covid_def_prob[,colnames(covid_def_prob) %in% original_weights$ID[original_weights$direction=="Up"]], na.rm = T)-
                           rowSums(covid_def_prob[,colnames(covid_def_prob) %in% original_weights$ID[original_weights$direction=="Down"]], na.rm = T))

covid_def_prob_roc_simple <- roc(covid_def_prob$covid_group, covid_def_prob$simple_drs, ci = T)

########## retrain the model weights, contrasting COVID to all groups 
model <- glm(covid_group ~., data = new_expression[new_expression$group %in% c("COVID19", "DB", "DV", "Not_Infected", "Mixed"),c(1:10,13)], family = "binomial")
new_weights <- summary(model)$coefficients
new_weights <- new_weights[2:nrow(new_weights),1]
new_intercept <- summary(model)$coefficients["(Intercept)","Estimate"]



# COVID vs definites 
covid_definites[is.na(covid_definites)] <- 0
response <- new_intercept + as.vector(as.matrix(covid_definites[,match(names(new_weights), 
                                                                            colnames(covid_definites))]) %*% as.matrix(new_weights))
covid_definites$retrained_drs <- exp(response) / (1 + exp(response))
covid_definites_roc_retrained <- roc(covid_definites$covid_group, covid_definites$retrained_drs, ci = T)

# COVID vs probables 
covid_probable[is.na(covid_probable)] <- 0
response <- new_intercept + as.vector(as.matrix(covid_probable[,match(names(new_weights), 
                                                                            colnames(covid_probable))]) %*% as.matrix(new_weights))
covid_probable$retrained_drs <- exp(response) / (1 + exp(response))
covid_probable_roc_retrained <- roc(covid_probable$covid_group, covid_probable$retrained_drs, ci = T)

# COVID vs not infected 
covid_not_infected[is.na(covid_not_infected)] <- 0
response <- new_intercept + as.vector(as.matrix(covid_not_infected[,match(names(new_weights), 
                                                                          colnames(covid_not_infected))]) %*% as.matrix(new_weights))
covid_not_infected$retrained_drs <- exp(response) / (1 + exp(response))
covid_not_infected_roc_retrained <- roc(covid_not_infected$covid_group, covid_not_infected$retrained_drs, ci = T)

# COVID vs all
covid_all[is.na(covid_all)] <- 0
response <- new_intercept + as.vector(as.matrix(covid_all[,match(names(new_weights), 
                                                                                  colnames(covid_all))]) %*% as.matrix(new_weights))
covid_all$retrained_drs <- exp(response) / (1 + exp(response))
covid_all_roc_retrained <- roc(covid_all$covid_group, covid_all$retrained_drs, ci = T)

# covid vs definites and probables 
# COVID vs all
covid_def_prob[is.na(covid_def_prob)] <- 0
response <- as.vector(as.matrix(covid_def_prob[,match(names(new_weights), 
                                                                colnames(covid_def_prob))]) %*% as.matrix(new_weights))
covid_def_prob$retrained_drs <- exp(response) / (1 + exp(response))
covid_def_prob_roc_retrained <- roc(covid_def_prob$covid_group, covid_def_prob$retrained_drs, ci = T)


# generate roc curves for each comparison with simple DRS and retrained 
# COVID vs definites
covid_def_roc_simple <- roc_plot(covid_definites_roc_simple, "Simple disease risk score", 'darkgrey')
covid_def_roc_retrained <- roc_plot(covid_definites_roc_retrained, "Retrained weights", '#118ab2')

covid_def_plots <- ggarrange(covid_def_roc_simple,covid_def_roc_retrained, ncol=2, nrow=1, common.legend = TRUE,legend="bottom")
covid_def_plots <- annotate_figure(covid_def_plots, top = text_grob("COVID-19 vs. definite infections"))

# COVID vs probables
covid_prob_roc_simple <- roc_plot(covid_probables_roc_simple, "Simple disease risk score", 'darkgrey')
covid_prob_roc_retrained <- roc_plot(covid_probable_roc_retrained, "Retrained weights", '#ff5400')

covid_prob_plots <- ggarrange(covid_prob_roc_simple,covid_prob_roc_retrained, ncol=2, nrow=1, common.legend = TRUE,legend="bottom")
covid_prob_plots <- annotate_figure(covid_prob_plots, top = text_grob("COVID-19 vs. probable infections"))

# COVID vs not infected 
covid_ni_roc_simple <- roc_plot(covid_not_infected_roc_simple, "Simple disease risk score", 'darkgrey')
covid_ni_roc_retrained <- roc_plot(covid_not_infected_roc_retrained, "Retrained weights", '#adff02')

covid_ni_plots <- ggarrange(covid_ni_roc_simple,covid_ni_roc_retrained, ncol=2, nrow=1, common.legend = TRUE,legend="bottom")
covid_ni_plots <- annotate_figure(covid_ni_plots, top = text_grob("COVID-19 vs. not infected"))

# COVID vs all
covid_all_plot_simple <- roc_plot(covid_all_roc_simple, "Simple disease risk score", 'darkgrey')
covid_all_plot_retrained <- roc_plot(covid_all_roc_retrained, "Retrained weights", '#3F0F3F')

covid_all_plots <- ggarrange(covid_all_plot_simple,covid_all_plot_retrained, ncol=2, nrow=1, common.legend = TRUE,legend="bottom")
covid_all_plots <- annotate_figure(covid_all_plots, top = text_grob("COVID-19 vs. all comparator groups"))


pdf(file = '../remade_plots/RT_PCR_Validation_Rocs.pdf', 
    height = 15.5, width = 7.5)
grid.arrange(covid_all_plots, 
             covid_def_plots, 
             covid_prob_plots, 
             covid_ni_plots, nrow = 4)
dev.off()


############ compare to WCC and CRP

crp_wcc <- read.csv(file = 'data/validation/crpwcc_validation.csv', header = T)
crp_wcc$UIN <- paste("IP-", crp_wcc$UIN, sep = "")

#crp_wcc$GROUP <- new_expression$group[match(crp_wcc$UIN, rownames(new_expression))]

crp_wcc <- crp_wcc[crp_wcc$UIN %in% rownames(new_expression),]

crp_wcc$covid <- 0
crp_wcc$covid[crp_wcc$GROUP == "COVID19"] <- 1

# all
roc_all_wcc <- roc(crp_wcc$covid, crp_wcc$WCC, ci = T)
roc_all_crp <- roc(crp_wcc$covid, as.numeric(crp_wcc$CRP), ci = T)

# definites 
roc_def_wcc <- roc(crp_wcc$covid[crp_wcc$GROUP %in% c("COVID19", 'DB', 'DV', 'Mixed', "Other_Infection")],
    crp_wcc$WCC[crp_wcc$GROUP %in% c("COVID19", 'DB', 'DV', 'Mixed', "Other_Infection")], ci = T)

roc_def_crp <- roc(crp_wcc$covid[crp_wcc$GROUP %in% c("COVID19", 'DB', 'DV', 'Mixed', "Other_Infection")],
    as.numeric(crp_wcc$CRP[crp_wcc$GROUP %in% c("COVID19", 'DB', 'DV', 'Mixed', "Other_Infection")]), ci = T)

# probable
roc_prob_wcc <- roc(crp_wcc$covid[crp_wcc$GROUP %in% c("COVID19", 'Inconclusive', 'PB', 'PV')],
    crp_wcc$WCC[crp_wcc$GROUP %in% c("COVID19", 'Inconclusive', 'PB', 'PV')], ci = T)

roc_prob_crp <- roc(crp_wcc$covid[crp_wcc$GROUP %in% c("COVID19", 'Inconclusive', 'PB', 'PV')],
    as.numeric(crp_wcc$CRP[crp_wcc$GROUP %in% c("COVID19", 'Inconclusive', 'PB', 'PV')]), ci = T)


# not infected
roc_ni_wcc <- roc(crp_wcc$covid[crp_wcc$GROUP %in% c("COVID19", 'Not_Infected')],
    crp_wcc$WCC[crp_wcc$GROUP %in% c("COVID19", 'Not_Infected')], ci = T)

roc_ni_crp <- roc(crp_wcc$covid[crp_wcc$GROUP %in% c("COVID19", 'Not_Infected')],
    as.numeric(crp_wcc$CRP[crp_wcc$GROUP %in% c("COVID19", 'Not_Infected')]), ci = T)


### plot the WCC and CRP ROC curves
roc_all_wcc_p <- roc_plot(roc_all_wcc, "WCC", '#3F0F3F')
roc_all_crp_p <- roc_plot(roc_all_crp, "CRP", '#3F0F3F')

roc_def_wcc_p <- roc_plot(roc_def_wcc, "WCC", '#118ab2')
roc_def_crp_p <- roc_plot(roc_def_crp, "CRP", '#118ab2')

roc_prob_wcc_p <- roc_plot(roc_prob_wcc, "WCC", '#ff5400')
roc_prob_crp_p <- roc_plot(roc_prob_crp, "CRP", '#ff5400')

roc_ni_wcc_p <- roc_plot(roc_ni_wcc, "WCC", '#adff02')
roc_ni_crp_p <- roc_plot(roc_ni_crp, "CRP", '#adff02')



covid_all_plots_wcc_crp <- ggarrange(roc_all_wcc_p,roc_all_crp_p, ncol=2, nrow=1, common.legend = TRUE,legend="bottom")
covid_all_plots_wcc_crp <- annotate_figure(covid_all_plots_wcc_crp, top = text_grob("COVID-19 vs. all comparator groups"))

covid_def_plots_wcc_crp <- ggarrange(roc_def_wcc_p,roc_def_crp_p, ncol=2, nrow=1, common.legend = TRUE,legend="bottom")
covid_def_plots_wcc_crp <- annotate_figure(covid_def_plots_wcc_crp, top = text_grob("COVID-19 vs. definite infections"))

covid_prob_plots_wcc_crp <- ggarrange(roc_prob_wcc_p,roc_prob_crp_p, ncol=2, nrow=1, common.legend = TRUE,legend="bottom")
covid_prob_plots_wcc_crp <- annotate_figure(covid_prob_plots_wcc_crp, top = text_grob("COVID-19 vs. probable infections"))

covid_ni_plots_wcc_crp <- ggarrange(roc_ni_wcc_p,roc_ni_crp_p, ncol=2, nrow=1, common.legend = TRUE,legend="bottom")
covid_ni_plots_wcc_crp <- annotate_figure(covid_ni_plots_wcc_crp, top = text_grob("COVID-19 vs. not infected"))

png(file = 'figures/CRP_WCC_rocs.png', 
    res = 700, units = 'in', 
    height = 15.5, width = 7.5)
grid.arrange(covid_all_plots_wcc_crp, 
             covid_def_plots_wcc_crp, 
             covid_prob_plots_wcc_crp, 
             covid_ni_plots_wcc_crp, nrow = 4)
dev.off()


# COVID vs definites and probables, for the supplementary 
covid_def_prob_plot_simple <- roc_plot(covid_def_prob_roc_simple, "Simple disease risk score", 'darkgrey')
covid_def_prob_plot_retrained <- roc_plot(covid_def_prob_roc_retrained, "Retrained weights", 'darkorchid1')

covid_def_prob_plots <- ggarrange(covid_def_prob_plot_simple,covid_def_prob_plot_retrained, ncol=2, nrow=1, common.legend = TRUE,legend="bottom")
covid_def_prob_plots <- annotate_figure(covid_def_prob_plots, top = text_grob("COVID-19 vs. definite and probable infections"))

pdf(file = 'figures/COVID_probable_definites_rt_pcr.pdf', height = 6, width = 11)
covid_def_prob_plots
dev.off()


