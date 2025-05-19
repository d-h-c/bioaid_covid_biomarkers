# functions for package - omicsHJ
load.packages <- function(){
  library('ggplot2')
  library(viridis)
  library(ggpubr)
  library(glmnet)
  library(pROC)
  library(limma)
  library(ROCR)
  library(dplyr)
  library(grid)
  library(gridExtra)
  library(stringr)
  library(clipr)
}

# boxplot function
# data must have samples in rows and features in columns
# var.name is the quantitative variable, e.g. protein, rna, to compare
# groups is the categorical variable between which the quantatitive variable is compared
# comps should contain all the comparisons you want to make in the following format:
#   list(c("group1", "group2"), c("group1", "group3"), c("group3, "group2"))
boxplots.disease <- function(data, var.name, groups, compare_means, comps){
  # subset data to contain column with feature expression values and groups of interest
  df <- data.frame(data[,groups], data[,var.name])
  print(paste(dim(df)[1], "samples;", sep = " ", length(levels(df[,1])), "groups:"))
  print(paste(levels(df[,1])))
  colnames(df) <- c("Groups", "Feature")
  # boxplot
  plot <- ggplot(df, aes(x = Groups, y = Feature, fill = Groups)) +
    geom_boxplot() +
    scale_fill_viridis_d()
  if(isTRUE(compare_means)){plot <- plot + stat_compare_means(comparisons = comps)}
  else(plot <- plot)
  plot <- plot + ggtitle(paste("Boxplots comparing expression of", var.name, "between groups", sep = " "))
  print(plot)
}


# ensure that the groups are in a column named 'groups'
# features must be in format like this: f.list <- c("protein_2", "protein_1", "protein_6")
# comps is in case there are more than 2 groups
# roc won't really work for more than 2 groups
# comps should contain all the comparisons you want to make in the following format:
#   list(c("group1", "group2"), c("group1", "group3"), c("group3, "group2"))
# the output from this function is drs
disease.risk.score <- function(data, features, groups, comps){
  # subset so the selected features are taken out of dataframe
  exp <- data[,features]
  group.df <- data[,groups]
  df <- data.frame(group.df, exp)
  colnames(df)[1] <- 'groups'
  formula.f <-  as.formula(paste(groups, paste(features, collapse=" + "), sep=" ~ "))
  # make linear model
  glm.res <- glm(formula.f, data = df, family=binomial(logit))
  # extract beta coefficients from glm results
  coefs <- glm.res$coefficients[2:length(glm.res$coefficients)]
  drs <- data.frame(as.matrix(df[,2:ncol(df)]) %*% as.matrix(coefs))
  drs <- data.frame(drs, df$groups)
  colnames(drs) <- c('drs', 'groups')
  # roc
  # note the roc will only work with 2 classes
  roc(groups ~ drs, data = drs, ci = TRUE,
      plot = TRUE, print.auc = TRUE, main = "ROC curve from DRS")
  # plot drs boxplots
  drs.plot <- ggplot(drs, aes(x = groups, y = drs, fill = groups)) + geom_boxplot()+
    ylab("Disease risk score")+xlab("Disease group")+
    scale_fill_viridis_d()+
    labs(fill = "Disease group")+
    ggtitle("Disease risk score (DRS) per disease group")
  if(length(levels(drs$groups))==2){
    drs.plot <- drs.plot+stat_compare_means(label.x = 1.4)}
  else if (length(levels(drs$groups))> 2){
    drs.plot <- drs.plot + stat_compare_means(comparisons = comps)}
  print(drs.plot)
  return(drs)
}


# Limma function
limma.fun <- function(data, p.thresh, comparisons, start.col, end.col, n.features, model.mat){
  # construct linear models
  l.m <- lmFit(t(data[,start.col:end.col]), model.mat)
  # for pairwise comparisons between groups, make a contrast matrix
  constrasts.design <- paste(unlist(comparisons), collapse = "-")
  constrasts.mat <- makeContrasts(contrasts = constrasts.design, levels = colnames(model.mat))
  # linear model for age, gender, disease + plate
  cb.fit <- eBayes(contrasts.fit(l.m, constrasts.mat))
  # list of DE features
  top.features <- topTable(cb.fit, coef = 1, adjust = "BH", number = n.features, p.value = p.thresh)
  print(top.features)
  return(top.features)
  }


# density plot function
# input to this function is the DRS and the groups of the samples
# ensure that column 1 of drs is the group and column 2 is drs
overlapping.density <- function(drs.df){
  colnames(drs.df) <- c("group", "drs")
  den.group.1 <- density(drs.df$drs[drs.df$group==levels(drs.df$group)[1]])
  den.group.2 <- density(drs.df$drs[drs.df$group==levels(drs.df$group)[2]])
  min.x <- min(min(den.group.1$x), min(den.group.2$x))-0.1
  max.x <- max(max(den.group.1$x), max(den.group.2$x))+0.1
  max.y <- max(max(den.group.1$y), max(den.group.2$y))+0.2
  plot(den.group.1, xlim = c(min.x, max.x), ylim = c(0,max.y), col = viridis(2)[1], lwd = 3,
       main = "Kernel Density Estimation plot from the DRS")
  lines(den.group.2, col = viridis(2)[2], lwd = 3)
  legend("topright", legend = levels(drs.df$group), fill = viridis(2), cex = 2)
  # calculate overlapping portition
  den.1 <- density(drs.df$drs[drs.df$group==levels(drs.df$group)[1]], na.rm = T)
  den.2 <- density(drs.df$drs[drs.df$group==levels(drs.df$group)[2]], na.rm = T)
  X <- den.1$x
  Y.gp1 <- den.1$y
  Y.gp2 <- den.2$y
  overlap <- pmin(Y.gp1, Y.gp2)
  # this colours in the overlap
  polygon(c(X, X[ 1 ]), c (overlap, overlap[1]),
          lwd = 2, col = 'darkgreen', density = 20)
  # this calualtes the number
  Total <- trapz(X, Y.gp1) + trapz(X, Y.gp2)
  (Surface <- trapz ( X, overlap ) / Total)
  SText <- paste ( sprintf ( "%.3f", 100*Surface ), "%" )
  # this adds the number to the plot
  text( X[which.max (overlap)], 1.2 * max(overlap), SText)
}


# plotly volcano plot 
# limma output must include ALL features in the dataset, irrespective of p value
volcano.plotly <- function(limma.output, p.thresh, lfc.thresh, username, key, name){
  limma.output$significance <- NA
  limma.output$significance[limma.output$adj.P.Val < p.thresh] <- "Significant"
  limma.output$significance[limma.output$adj.P.Val > p.thresh] <- "Not significant"
  limma.output$significance[limma.output$adj.P.Val > p.thresh && abs(limma.output$logFC) > lfc.thresh] <- "Not significant, high LFC"
  limma.output$significance[limma.output$adj.P.Val < p.thresh & abs(limma.output$logFC) > lfc.thresh] <- "Significant, high LFC"
  limma.output$name <- rownames(limma.output)
  # log into plotly 
  Sys.setenv("plotly_username"=username)
  Sys.setenv("plotly_api_key"=key)
  # create plot 
  plot <- plot_ly()
  plot <- add_trace(plot,
                    data = limma.output,
                    x=~logFC,
                    y=-log(limma.output$P.Value),
                    type = 'scatter', 
                    text = ~name,
                    color = ~significance,
                    mode = 'markers') %>% 
    layout(title = paste("Volcano plot for: ", name, sep = ""), 
           xaxis = list(
             range = c(-3,3)
             ))
  print(plot)
  api_create(plot, filename = 'Limma: KD vs DB-DV (untargeted)')
}
# example call: 
#untargeted.kd <- volcano.plotly(limma.output = kd.bv.limma, p.thresh = 0.05, lfc.thresh = 1, username = "hj4817", key = "NTPQPlYxPaOwjspCopKM", name = "KD vs DB+DV (untargeted)")

# fspls packages
fspls.packages <- function(){
  library(nnet)
  library(dplyr)
  library(pROC)
  library(ROCR)
  library(glmnet)
}



# FSPLS function
# this function assumes that the version of fspls running is: fspls_lars_multi_meta.R, for me this is in my githib account 
# the code must be modified, if not using the git version: line 2892 
# FROM:		all_results[[global_index]] = do.call(c, list(all_results[[global_index]], l2))
# TO: 		all_results[[global_index]] = do.call(list, list(all_results[[global_index]], l2))
# expression must be data frame with n x p (samples x variables)
# response must have the same length as the number of rows of expression and be in numeric format
fspls.binomial <- function(path.fspls, expression, response, p.thresh, beam, split, split.ratio, max){
  source(path.fspls)

  fspls.packages()
  ##FOLLOWING JUST SETS DEFAULT VALUES, can change with options if you want to
  # params =readJSONParams(getOption("fspls.default_params"))
  options("fspls.family"= "binomial") ## binomial method 
  options("method" = "fspls")
  options("fspls.lambda" = 0)       ## if 0, no shrinkage. IF specified, then uses that shrinkage, if NULL, then uses CV
  options("fspls.lambda1" = 1)      ## for pvalue adjustment
  options("fspls.debug" = "off")     ## a lot of debugging information printed
  options("fspls.log" = NULL)
  options("fspls.beam" = beam)
  options("fspls.pv_thresh" = p.thresh)
  options("fspls.max" = max)
  elapse_start = proc.time()
  # create input for fspls
  # check if test training split is wanted 
  if(isTRUE(split)){
    print(paste("Splitting by: ", split.ratio, sep = ""))
    ind <- sample.split(response, SplitRatio = split.ratio)
    train.exp <- expression[ind == TRUE,]
    test.exp <- expression[ind == FALSE,]
    train.res <- response[ind == TRUE]
    test.res <- response[ind == FALSE]
    print(train.res)
    print(test.res)
    # create test and training data items 
    train.obj <- list(data = train.exp, y = data.frame(train.res), 
                      weights = data.frame(rep(1, length(train.res))))
    test.obj <- list(data = test.exp, y = data.frame(test.res), 
                     weights = data.frame(rep(1, length(test.res))))
    train.obj <- list("train" = train.obj)
    test.obj <- list("test" = test.obj)
  }
  else if(isFALSE(split)){
    weights = rep(1, length(response))
    obj <- list(data = expression, y = data.frame(response), weights = data.frame(response = weights))
    train.obj <- list("train" = obj)
    test.obj <- list("test" = obj)
  }
  #FOLLOWING RESTRICTS TO VARIABLES WITH NA
  options("fspls.pheno" = dimnames(train.obj[[1]]$y)[[2]][1])
  print("LOADING DATA ELPASED TIME:")
  print(proc.time()-elapse_start)
  ## p-value threshold will determine the number of selected features
  elapse_start = proc.time()
  # run fs-pls
  
  model = trainModel(trainOriginal_l1 = train.obj,testOriginal_l1 = test.obj)#, pv_thresh = p.thresh, max = max)
}


#### fspls with multiple iterations returning top signaure, auc, acc for each iteration 
fspls.iterate <- function(n.iterations, path.fspls, expression, response, max, p.thresh, beam, split, split.ratio, seed){
  accuracy <- list()
  signature <- list()
  aucs <- list()
  aucs_train <- list()
  models <- list()
  df.names <- data.frame()
  for(i in 1:n.iterations){
    seed.i <- seed+i
    set.seed(seed.i)
    print(paste("Iteration",i, sep = " "))
    model <- fspls.binomial(path.fspls = path.fspls,
                            expression = expression, 
                            response = response, p.thresh = p.thresh, 
                            beam = beam, max = max, 
                            split = split, split.ratio = split.ratio)
    # extract highest AUC 
    idx <- length(model)
    ind <- which(model[[idx]][[1]][['testeval']][,'auc'] == max (model[[idx]][[1]][['testeval']][,'auc']))
    print(ind)
    idx_genes <- rownames(model[[idx]][[1]][['testeval']])[ind]
    print(paste("ID:", idx_genes, sep = " "))
    if (length(idx_genes) > 1) {
      idx_sig <- which(model[[idx]][[1]][["testeval"]][idx_genes,'acc']==max(model[[idx]][[1]][["testeval"]][idx_genes,'acc']))[1]
      idx_sig <- idx_genes[idx_sig]
    } else idx_sig <- idx_genes
    print(idx_sig)
    acc <- model[[idx]][[1]]$testeval[rownames(model[[idx]][[1]]$testeval) == idx_sig,"acc"]
    auc <- model[[idx]][[1]]$testeval[rownames(model[[idx]][[1]]$testeval) == idx_sig, 'auc']
    auc_train <- model[[idx]][[1]]$traineval[rownames(model[[idx]][[1]]$testeval) == idx_sig, 'auc']
    names.frame <- model[[idx]][[1]]$variablesl[rownames(model[[idx]][[1]]$testeval) == idx_sig]
    df <- c(table(names.frame), auc = auc)
    pros <- str_split(idx_sig, "\\_\\_")
    pros <- c(unlist(pros), "AUC")
    pros <- str_remove_all(pros, "\\.test")
    df <- t(data.frame(df))
    colnames(df) <- c(pros)
    accuracy[[i]] <- acc
    signature[[i]] <- idx_sig    
    aucs[[i]] <- auc
    models[[i]] <- model
    aucs_train[[i]] <- auc_train
    #print(df.names)
    df <- data.frame(df)
    df.names <- bind_rows(df.names, df)
  }
  return(list(acc = accuracy, sig = signature, auc = aucs, model = models, table = df.names, training_aucs = aucs_train))
}




# uniprot IDs from somascan IDs 
uniprot.df <- function(top.pro){
  uni.idx <- lapply(rownames(top.pro), function(x){
    uni <- proteinData$UniProt[grep(x, proteinData$SeqId)]
    return(uni)
  })
  uni.idx <- unlist(as.matrix(uni.idx))
  return(uni.idx)
}


# names from somascan IDs 
names.df <- function(top.pro){
  uni.idx <- lapply(rownames(top.pro), function(x){
    uni <- proteinData$Target[grep(x, proteinData$SeqId)]
    return(uni)
  })
  uni.idx <- unlist(as.matrix(uni.idx))
  return(uni.idx)
}

# entrez gene id from somascan IDs 
entrez.df <- function(top.pro){
  uni.idx <- lapply(rownames(top.pro), function(x){
    ent <- proteinData$EntrezGeneID[grep(x, proteinData$SeqId)]
    return(ent)
  })
  uni.idx <- unlist(as.matrix(uni.idx))
  return(uni.idx)
}

# PCA plot function 
# data must have samples in ROWS and features in COLUMNS 
pca.fun <- function(data, colour_by, dim1, dim2){
  pc <- pca(data, nPcs = 10, center = TRUE, scale = "uv", method = 'ppca')
  pc <- data.frame(pc@scores[,dim1:dim2], col_by = colour_by)
  plot <- ggplot(pc, aes(x = pc[,1], y = pc[,2], color = col_by))+geom_point()#+
    #geom_text(label=rownames(pc))
    #stat_ellipse(geom = "polygon", alpha = 1/5, aes(fill = col_by))
  print(plot)
  return(plot)
}


# upset plot, output from fs-pls 
# ensure matrix is in correct format: last column corresponding to AUC (capitals)
upset.plot <- function(frequency.table, path){
  viridis.pal <- viridis(4)
  upset.plot <- upset(data.frame(frequency.table), nset = (ncol(frequency.table)-1), nintersects = nrow(frequency.table), order.by = "freq", 
                      text.scale = 0.8, boxplot.summary = c("AUC"),mb.ratio = c(0.3, 0.7), line.size = 0.5,
                      matrix.color = viridis.pal[1], main.bar.color = viridis.pal[2], sets.bar.color = viridis.pal[3], att.color = viridis.pal[4])
  pdf(file = path, height = 16, width = 14)
  print(upset.plot)
  dev.off()
  print(upset.plot)
}


# disco plot function
# if create = Y, make sure system settings have been loaded so file can be saved in the right place 
disco.fun <- function(data.set, x, y, colour, name.list, comparison, create){
  plot <- plot_ly()
  plot <- add_trace(plot, data = data.set, x = ~data.set[,x], y = ~data.set[,y], type = 'scatter', mode = 'markers', text = name.list, color = ~data.set[,colour]) %>%
    layout(title = comparison, xaxis = list(
      title = x, 
      range = c(-4,4)),
      yaxis = list(
        title = y, 
        range= c(-4,4)
      ))
  if(create == "Y"){
    api_create(plot, filename = comparison)}
  else(return(plot))
}


# function to extract DE proteins according to group of interst 
de.fun <- function(limma.list, up.down, soma, ms){
  if(up.down == "up"){
    list <- limma.list[limma.list$logFC > 0,]
  } else if(up.down == "down"){
    list <- limma.list[limma.list$logFC < 0,]
  }
  if(soma == TRUE){
    list$entrez <- entrez.df(list)
  }
  if(ms == TRUE){
    list$uniprot <- str_remove_all(rownames(list), "sp\\|")
    list$uniprot <- str_remove_all(list$uniprot, "\\_HUMAN.{1,50}")
    list$uniprot <- str_remove_all(list$uniprot, "\\_HUMAN")
    list$uniprot <- str_remove_all(list$uniprot, "tr\\|")
    list$uniprot <- str_remove_all(list$uniprot, "\\|.{1,50}")
  }
  list <- list[!(is.na(list$logFC)),]
  return(list)
}


#######
# iterative k means clustering 
kmeans.iterate <- function(n.it, dataset, ks){
  seed = 123
  total.ss <- vector()
  k.means.models <- list()
  for(i in 1:n.it){
    seed = seed+i
    set.seed(seed)
    k.means.res <- kmeans(dataset, centers = ks)
    print(k.means.res$tot.withinss)
    total.ss[[i]] <- k.means.res$tot.withinss
    k.means.models[[i]] <- k.means.res
  }
  top.mod <- k.means.models[which(total.ss == max(total.ss))]
  return(top.mod)
}

# coefficient of variation function 
CV <- function(mean, sd){
  (sd/mean)*100
}


pairwise.ratios <- function(x, prefix="probeset", char=":"){
  n <- ncol(x)
  cn <- colnames(x)
  if(length(cn) == 0){
    cn <- gsub(" ", "0", formatC(seq.int(n), width=nchar(n)))
    cn <- paste(prefix, cn, sep="")
  }
  cmb <- combn(n, 2)
  r1 <- apply(cmb, 2, function(j) x[, j[1]]/x[, j[2]])
  r2 <- apply(cmb, 2, function(j) x[, j[2]]/x[, j[1]])
  colnames(r1) <- apply(cmb, 2, function(j) paste(cn[j], collapse=char))
  colnames(r2) <- apply(cmb, 2, function(j) paste(cn[rev(j)], collapse=char))
  cbind(r1, r2)[, order(c(colnames(r1), colnames(r2)))]
}


roc_plot <- function(roc, title, colour){
  ci <- ci.se(roc, specificities=seq(0, 1, l=25))
  ci <- data.frame(x = as.numeric(rownames(ci)),
                   lower = ci[, 1],
                   upper = ci[, 3])
  p <- ggroc(roc) + theme_bw() + geom_abline(slope=1, intercept = 1, linetype = "dashed", alpha=0.7, color = "grey")+
    labs(x = 'Specificity', y = "Sensitivity")+
    geom_ribbon(
      data = ci,
      aes(x = x, ymin = lower, ymax = upper),
      alpha = 0.3,
      fill = colour,
      inherit.aes = F) +
    annotate("text", x=0.4, y=0.1, label= paste("AUC: ", 
                                                round(auc(roc),    3)*100, "% (95% CI: ",
                                                round( ci(roc)[1], 3)*100,'%-', 
                                                round( ci(roc)[3], 3)*100, "%)",
                                                sep = ''))+
    theme(plot.title = element_text(size=12))+
    ggtitle(title)
  print(p)
}



recalculate_auc <- function(expression, directions, n_features, group){
  rocs <- list()
  for(i in 1:n_features){
    if(i == 1){
      drs <- expression[,i]
      roc = roc(response = group, predictor = drs)
      rocs[[i]] <- roc
    }
    else{
      drs <- calculate_drs(as.data.frame(expression[,1:i]), directions[1:i])
      roc <- roc(response = group, predictor = drs$DRS)
      rocs[[i]] <- roc
    }
  }
  return(rocs)
}

calculate_drs <- function(expression, weights){
  drs <- data.frame(DRS = as.matrix(expression) %*% as.matrix(weights))
  return(drs)
}


deseq.results <- function(res, gene.names, outfile){
  res <- res[order(res$padj),]
  res.df <- data.frame(res)
  rownames(res.df) <- rownames(res)
  
  gene.names <- gene.names[match(rownames(res.df), gene.names$ensembl_gene_id),]
  res.df$gene <- gene.names$external_gene_name
  res.df$Significant <- 0
  res.df$Significant[res.df$padj < 0.05] <- 1
  
  
  # create volcano plot column
  res.df$Color <- 'black'
  res.df$Color[abs(res.df$log2FoldChange) > 1] <- 'limegreen'
  res.df$Color[res.df$padj < 0.05] <- 'gold'
  res.df$Color[res.df$padj < 0.05 & 
                 abs(res.df$log2FoldChange) > 1] <- 'red'
  
  res.df$Color <- factor(res.df$Color, levels = c("black", 'limegreen', "gold", 'red'))
  
  res.df <- res.df[!(is.na(res.df$log2FoldChange)),]
  res.df <- res.df[!(is.na(res.df$padj)),]
  print(range(res.df$log2FoldChange))
  write.csv(file = paste(outfile, '.csv', sep = ""), res.df)
  return(res.df)
}
