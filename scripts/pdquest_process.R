library(xlsx)
library(reshape)
library(dplyr)
library(limma)
library(ggplot2)
library(ggfortify)
library(tidyr)
library("languageR")
library(nlme)
library(gplots)
library(devtools)
library(ComplexHeatmap)
library(dendextend)
library(circlize)
library(remotes)


#install_version("rjson", "0.2.15")
#install_version("circlize", "0.4.14")

#install_github("jokergoo/ComplexHeatmap")

#BiocManager::install("limma")


roots_raw_dat <- read.xlsx('PQquest_data/rice_roots_exported.xlsx', sheetIndex = 1)
redone_roots_all <- roots_raw_dat %>%  dplyr::select("SSP","Gel.Name", "Quantity" ) %>% 
  reshape(idvar = "SSP", timevar = "Gel.Name", direction = "wide")
new_colnames=c()

inds_col=c('1','2','3','1','2','3','4','5','1','2','3','4','5')
i=1
for (name in colnames(redone_roots_all)[2:14]){
  split1<-strsplit(name,split='.', fixed=TRUE)[[1]][2]
  split2<-strsplit(split1,split='_', fixed=TRUE)[[1]][1]
  new_colnames <- c(new_colnames,paste(split2,inds_col[i], sep= '_'))
  i<-i+1
}

#select samples
redone_roots <- redone_roots_all[,c(1,2,3,4,5,6,7,10,11,12)]
colnames(redone_roots)[2:10] <- c('A1','A2','A3','C1','C2','C3','R1','R2','R3')
new_roots = read.csv('PQquest_data/new_spots_roots.csv', sep=";", na.strings=c("","NA"), header = TRUE)
redone_roots <- rbind(redone_roots, new_roots)


#take logarithm   
redone_roots_log <- cbind(redone_roots[1], log2(redone_roots[2:10]))
redone_roots_log_reshaped <- gather(redone_roots_log[2:10], condition, measurement, A1:R3, factor_key=TRUE)
redone_roots_log_reshaped$factor_col <- 'condition'
redone_roots_log_reshaped[grepl("A", redone_roots_log_reshaped$condition), 3] <- 'Anoxia'
redone_roots_log_reshaped[grepl("C", redone_roots_log_reshaped$condition), 3] <- 'Control'
redone_roots_log_reshaped[grepl("R", redone_roots_log_reshaped$condition), 3] <- 'Re-aeration'

redone_roots_log_reshaped$factor_col <- factor(redone_roots_log_reshaped$factor_col, 
                                              levels=c('Control', 'Anoxia','Re-aeration'))

redone_roots_log_reshaped$condition <- factor(redone_roots_log_reshaped$condition ,
                                              levels= c('C1','C2','C3','A1','A2','A3','R1','R2','R3'))
#non-?normalized box-plots
ggplot(redone_roots_log_reshaped, aes(y=measurement, x=condition, fill=factor_col))+
  geom_boxplot(alpha=0.9)+theme_bw()+
  scale_fill_manual(values = c("#678667","#375DB6","#A5432B"))+
  xlab('Sample') + ylab('Optical density')+
  theme(axis.text.y = element_text(color='black', 
                                   size=16),
        axis.title.y=element_text(color="black", 
                                  size=18),
        panel.background = element_blank(), 
        axis.line = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(color='black', 
                                   size=16),
        axis.title.x = element_text( color="black", 
                                     size=18),
        legend.title=element_text(size=18), 
        legend.text=element_text(size=16)
  )+ 
  guides(fill = guide_legend(override.aes = list(size = 1),title="Condition")) 


#mean OD values
anoxia = redone_roots_log_reshaped[redone_roots_log_reshaped$factor_col=='Anoxia', 2] #19.8538 - anoxia
control = redone_roots_log_reshaped[redone_roots_log_reshaped$factor_col=='Control', 2] #20.76237 - control
reaeration <- redone_roots_log_reshaped[redone_roots_log_reshaped$factor_col=='Re-aeration', 2] #21.6658 - re-aeration
mean_OD_df <- data.frame(Density = c(anoxia, control, reaeration), 
                         Condition=c(rep('anoxia',279),rep('control',279),rep('reaeration',279)))

t_test_df=data.frame(comparision=c("Anoxia_vs_Control","Anoxia_vs_Re-aeration", 'Control_vs_Re-aeration'), p.value=1)
t_test_df[1, 2] <- t.test(anoxia, control)$p.value
t_test_df[2, 2] <- t.test(anoxia, reaeration)$p.value
t_test_df[3, 2] <- t.test(control, reaeration)$p.value
t_test_df$p.adj <- p.adjust(t_test_df$p.value, method = "fdr")
write.xlsx(t_test_df,'t_test_OD_roots.xlsx')


#Mean OD violins
ggplot(redone_roots_log_reshaped, aes(y=measurement, x=factor_col, fill=factor_col))+
  geom_violin(alpha=0.9, color='black')+theme_bw()+
  scale_fill_manual(values = c("#678667","#375DB6","#A5432B"))+
  geom_hline(yintercept=19.8538, size=2, linetype='dashed', color='#375DB6')+
  geom_hline(yintercept=20.76237, size=2, linetype='dashed', color='#678667')+
  geom_hline(yintercept=21.6658, size=2, linetype='dashed', color='#A5432B')+
  xlab('Condition') + ylab('Optical density')+
  theme(axis.text.y = element_text(color='black', 
                                   size=16),
        axis.title.y=element_text(color="black", 
                                  size=18),
        panel.background = element_blank(), 
        axis.line = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(color='black', 
                                   size=16),
        axis.title.x = element_text( color="black", 
                                     size=18),
        legend.title=element_text(size=18), 
        legend.text=element_text(size=16)
  )+ 
  guides(fill = guide_legend(override.aes = list(size = 1),title="Condition"))


#normalize data
redone_roots_normalized <- cbind(redone_roots_log[1],  normalizeQuantiles(redone_roots_log[2:10]))
redone_roots_normalized_reshaped <- gather(redone_roots_normalized[2:10], condition, measurement, A1:R3, factor_key=TRUE)
redone_roots_normalized_reshaped$factor_col <- 'condition'
redone_roots_normalized_reshaped[grepl("A", redone_roots_normalized_reshaped$condition), 3] <- 'Anoxia'
redone_roots_normalized_reshaped[grepl("C", redone_roots_normalized_reshaped$condition), 3] <- 'Control'
redone_roots_normalized_reshaped[grepl("R", redone_roots_normalized_reshaped$condition), 3] <- 'Re-aeration'


redone_roots_normalized_reshaped$factor_col <- factor(redone_roots_normalized_reshaped$factor_col, 
                                               levels=c('Control', 'Anoxia','Re-aeration'))

redone_roots_normalized_reshaped$condition <- factor(redone_roots_normalized_reshaped$condition ,
                                              levels= c('C1','C2','C3','A1','A2','A3','R1','R2','R3'))

#normalized box-plots
ggplot(redone_roots_normalized_reshaped, aes(y=measurement, x=condition, fill=factor_col))+
  geom_boxplot(alpha=0.9)+theme_bw()+
  scale_fill_manual(values = c("#678667","#375DB6","#A5432B"))+
  xlab('Sample') + ylab('Optical density')+
  theme(axis.text.y = element_text(color='black', 
                                   size=16),
        axis.title.y=element_text(color="black", 
                                  size=18),
        panel.background = element_blank(), 
        axis.line = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(color='black', 
                                   size=16),
        axis.title.x = element_text( color="black", 
                                     size=18),
        legend.title=element_text(size=18), 
        legend.text=element_text(size=16)
  )+ 
  guides(fill = guide_legend(override.aes = list(size = 1),title="Condition")) 


#annotate points
annotation_roots_df <-  read.csv('PQquest_data/PD_quest_root_annot.csv', sep=";", na.strings=c("","NA"), header = TRUE,
                                 stringsAsFactors = F)
annotation_roots_df_spots <-annotation_roots_df[annotation_roots_df$Accession!='-',1]
annotation_roots <- redone_roots_normalized[redone_roots_normalized$SSP %in% annotation_roots_df_spots,-1]



redone_roots_normalized$spot_num <- '0'
redone_roots_normalized$spot_acc <- '0'
redone_roots_normalized$spot_annot <- '0'

for(ind in 1:nrow(annotation_roots_df)){
  spot=annotation_roots_df[ind, 1]
  redone_roots_normalized[redone_roots_normalized$SSP==spot, 11] <- annotation_roots_df[ind,3]
  redone_roots_normalized[redone_roots_normalized$SSP==spot, 12] <- annotation_roots_df[ind,2]
  redone_roots_normalized[redone_roots_normalized$SSP==spot, 13] <- annotation_roots_df[ind,4]
}
rownames(redone_roots_normalized) <- redone_roots_normalized$spot_num




annotation_roots_annot_with_accessions <- redone_roots_normalized[redone_roots_normalized$SSP %in% annotation_roots_df_spots,-c(1,11,13)]
annotation_roots_annot_with_accessions$fixed_name <- 'name'

for (ind in 1:nrow(annotation_roots_annot_with_accessions)){
  acc_str=annotation_roots_annot_with_accessions[ind,10]
  acc_split=strsplit(acc_str,split='[', fixed=TRUE)[[1]][1]
  annotation_roots_annot_with_accessions[ind,11] <- acc_split
}

annotation_roots_annot_with_accessions_means <- data.frame(control=1:41, anoxia=1:41, reaeration=1:41)
annotation_roots_annot_with_accessions_means$anoxia <- rowMeans(annotation_roots_annot_with_accessions[,c(1,2,3)])
annotation_roots_annot_with_accessions_means$control <- rowMeans(annotation_roots_annot_with_accessions[,c(4,5,6)])
annotation_roots_annot_with_accessions_means$reaeration <- rowMeans(annotation_roots_annot_with_accessions[,c(7,8,9)])
annotation_roots_annot_with_accessions_means$accessions <- annotation_roots_annot_with_accessions$fixed_name

annotation_roots_annot_with_accessions_means <- annotation_roots_annot_with_accessions_means %>% 
  group_by(accessions) %>% 
  summarize(OD_control=mean(control), OD_anoxia=mean(anoxia), OD_reaeration =mean(reaeration)) %>% as.data.frame()
annotation_roots_annot_with_accessions_means <- annotation_roots_annot_with_accessions_means[-c(3,7,9,10,11,13,15),]

write.xlsx(redone_roots_normalized,'roots_with_new_log.xlsx')


#All conditions k-means
#find optimal clusters
sc_roots<-data.frame(t(redone_roots_normalized[2:10]))

#elbow method visualization
wss <- (nrow(sc_roots)-1)*sum(apply(sc_roots,2,var))
for (i in 1:4) wss[i] <- sum(kmeans(sc_roots,
                                    centers=i)$withinss,nstart=25,iter.max=1000)

wss_gf=data.frame(clust=c(1:4), wss=wss)

ggplot(wss_gf, aes(clust,wss))+
  geom_line()+theme_bw()+ xlab('Number of clusters') + ylab('WSS')+
  geom_vline(xintercept =2, linetype="dashed")+
  theme(axis.text.y = element_text(color='black', 
                                   size=12),
        axis.title.y=element_text(color="black", 
                                  size=16),
        panel.background = element_blank(), 
        axis.line = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(color='black', 
                                   size=12),
        axis.title.x = element_text( color="black", 
                                     size=16)
  ) 

#k-means clustering
roots_kmeans <- kmeans(sc_roots,centers=2)
roots_kmeans$cluster
sc_roots$name=rownames(sc_roots)

#k-means visualization
autoplot(roots_kmeans, 
         data=sc_roots, label = TRUE,frame = TRUE,
         frame.type = 'norm',label.size = 6,label.col='black',alpha=0)+
  scale_colour_manual(values = c('black','black','black'))+
  scale_fill_manual(values = c("#A9A9A9","#A60B0B",'red'))+
  theme_bw()+
  theme(axis.text.y = element_text(color='black', 
                                   size=24),
        axis.title.y=element_text(color="black", 
                                  size=28),
        panel.background = element_blank(), 
        axis.line = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(color='black', 
                                   size=24),
        axis.title.x = element_text(color="black", 
                                    size=28),
        legend.position = 'none'
  ) + 
  scale_x_continuous(expand = c(.1, .1)) 
scale_y_continuous(expand = c(.1, .1))



#significance test for normalized data

#kruskell test function
k_p_val <- function(x) { vals_df <- x %>% gather(condition, measurement, A1:R3, factor_key=TRUE)
   vals_df$factor_col <- 'factor'
   vals_df[grepl("A", vals_df$condition), 3] <- 'Anoxia'
   vals_df[grepl("C", vals_df$condition), 3] <- 'Control'
   vals_df[grepl("R", vals_df$condition), 3] <- 'Re-aeration'
  tryCatch(kruskal.test(measurement ~ factor_col, data = vals_df)$p.value,
           error = function(e) NA)
}

#mixed linear model function
lme_p_val <- function(x) { vals_df <- x %>% gather(condition, measurement, A1:R3, factor_key=TRUE)
    vals_df$factor_col <- 'factor'
    vals_df[grepl("A", vals_df$condition), 3] <- 'Anoxia'
    vals_df[grepl("C", vals_df$condition), 3] <- 'Control'
    vals_df[grepl("R", vals_df$condition), 3] <- 'Re-aeration'
    fit.lme <- lme(measurement ~ factor_col, random=~1|condition, data=vals_df)
    ret_vec <- c(as.data.frame(anova(fit.lme))$`p-value`[2])
    fit_sum <- as.data.frame(summary(fit.lme)$tTable)
    ret_vec <- c(ret_vec, fit_sum$`p-value`[c(2,3)])
    tryCatch(ret_vec,
         error = function(e) NA)
}

pvals_root_kraskell <- c()
pvals_root_lme <- data.frame(all_pval=c(1), cont_pval = c(1), re_pval = c(1))


#find p-vals for tests
for (i in 1:nrow(redone_roots_normalized)){
  pvals_root_kraskell <- c(pvals_root_kraskell, k_p_val(redone_roots_normalized[i,2:10]))
  pvals_root_lme <- rbind(pvals_root_lme, lme_p_val(redone_roots_normalized[i,2:10]))
}
pvals_root_lme <- pvals_root_lme[-1,]

#adjust p-values
pval_df=data_frame(spot=redone_roots_normalized[,1], kraskell_raw=pvals_root_kraskell, 
                   kraskell_adj = p.adjust(pvals_root_kraskell, method = "fdr"),
                   lme_all_raw = pvals_root_lme[,1],
                   lme_all_adj = p.adjust(pvals_root_lme[,1], method = "fdr"),
                   lme_cont_raw = pvals_root_lme[,2],
                   lme_cont_adj = p.adjust(pvals_root_lme[,2], method = "fdr"),
                   lme_re_raw = pvals_root_lme[,3],
                   lme_re_adj = p.adjust(pvals_root_lme[,3], method = "fdr"))



#adjustment results
#FDR
sum(pval_df$kraskell_raw <= 0.05, na.rm = TRUE) #3 proteins with significant 
sum(pval_df$lme_all_raw <= 0.05, na.rm = TRUE) #20 proteins with significant
sum(pval_df$lme_cont_raw <= 0.05, na.rm = TRUE) #21 proteins with significant 
sum(pval_df$lme_re_raw <= 0.05, na.rm = TRUE) #7 proteins with significant
sum(pval_df$kraskell_adj <= 0.05, na.rm = TRUE) #0 proteins with significant 
sum(pval_df$lme_all_adj <= 0.05, na.rm = TRUE) #0 proteins with significant
sum(pval_df$lme_cont_adj <= 0.05, na.rm = TRUE) #0 proteins with significant 
sum(pval_df$lme_re_adj <= 0.05, na.rm = TRUE) #0 proteins with significant

#holm
sum(pval_df$kraskell_adj <= 0.05, na.rm = TRUE) #0 proteins with significant 
sum(pval_df$lme_all_adj <= 0.05, na.rm = TRUE) #0 proteins with significant
sum(pval_df$lme_cont_adj <= 0.05, na.rm = TRUE) #0 proteins with significant 
sum(pval_df$lme_re_adj <= 0.05, na.rm = TRUE) #0 proteins with significant

#BH
sum(pval_df$kraskell_adj <= 0.05, na.rm = TRUE) #0 proteins with significant 
sum(pval_df$lme_all_adj <= 0.05, na.rm = TRUE) #0 proteins with significant
sum(pval_df$lme_cont_adj <= 0.05, na.rm = TRUE) #0 proteins with significant 
sum(pval_df$lme_re_adj <= 0.05, na.rm = TRUE) #0 proteins with significant

#Modded_linear_model
root.fac <- data.frame(Condition=c('anoxia','anoxia','anoxia','control','control','control', 're-aeration', 're-aeration', 're-aeration'))
rownames(root.fac) <- colnames( redone_roots_normalized[,2:10])
rownames(redone_roots_normalized) <- gsub(" ", "_", redone_roots_normalized[,1])

X <- model.matrix(~ Condition, data = root.fac)
fit <- lmFit(redone_roots_normalized[,2:10], design = X, method = "robust", maxit = 10000)
efit <- eBayes(fit)

#modelled_p-values adjustment
num_spots <- nrow(redone_roots_normalized)
pval_df$eBayes_cont_raw <- topTable(efit, coef = 2, number = num_spots,
                      sort.by = "none")$P.Value
pval_df$eBayes_cont_adj <- topTable(efit, coef = 2, number = num_spots,
                                    sort.by = "none")$adj.P.Val
pval_df$eBayes_re_raw <- topTable(efit, coef = 3, number = num_spots,
                                    sort.by = "none")$P.Value
pval_df$eBayes_re_adj <- topTable(efit, coef = 3, number = num_spots,
                                    sort.by = "none")$adj.P.Val

#the number of signifficant proteins
sum(pval_df$eBayes_cont_adj<=0.05) #26
sum(pval_df$eBayes_cont_raw<=0.05) #37
sum(pval_df$eBayes_re_adj<=0.05) #5
sum(pval_df$eBayes_re_raw<=0.05) #19

write.xlsx(pval_df,'dots_pvals_roots.xlsx')

#extract significant spots
signif_pval_roots <- redone_roots_normalized[pval_df$eBayes_cont_adj<=0.05 | pval_df$eBayes_re_adj<=0.05,]
rownames(signif_pval_roots) <-redone_roots_normalized[pval_df$eBayes_cont_adj<=0.05 | pval_df$eBayes_re_adj<=0.05,11]

write.xlsx(signif_pval_roots,'rice_roots_signif_spots_annot.xlsx')


#Signif conditions k-means
#Clusters for significant points
sc_roots<-data.frame(t(signif_pval_roots[, 2:10]))


#elbow method visualization
wss <- (nrow(sc_roots)-1)*sum(apply(sc_roots,2,var))
for (i in 1:4) wss[i] <- sum(kmeans(sc_roots,
                                    centers=i)$withinss,nstart=25,iter.max=1000)

wss_gf=data.frame(clust=c(1:4), wss=wss)

ggplot(wss_gf, aes(clust,wss))+
  geom_line()+theme_bw()+ xlab('Number of clusters') + ylab('WSS')+
  geom_vline(xintercept =2, linetype="dashed")+
  theme(axis.text.y = element_text(color='black', 
                                   size=12),
        axis.title.y=element_text(color="black", 
                                  size=16),
        panel.background = element_blank(), 
        axis.line = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(color='black', 
                                   size=12),
        axis.title.x = element_text( color="black", 
                                     size=16)
  ) 

#k-means clustering
roots_kmeans <- kmeans(sc_roots,centers=2)
roots_kmeans$cluster
sc_roots$name=rownames(sc_roots)

#k-means visualization
autoplot(roots_kmeans, 
         data=sc_roots, label = TRUE,frame = TRUE,
         frame.type = 'norm',label.size = 6,label.col='black',alpha=0)+
  scale_colour_manual(values = c('black','black','black'))+
  scale_fill_manual(values = c("#A9A9A9","#A60B0B",'red'))+
  theme_bw()+
  theme(axis.text.y = element_text(color='black', 
                                   size=24),
        axis.title.y=element_text(color="black", 
                                  size=28),
        panel.background = element_blank(), 
        axis.line = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(color='black', 
                                   size=24),
        axis.title.x = element_text(color="black", 
                                    size=28),
        legend.position = 'none'
  ) + 
  scale_x_continuous(expand = c(.1, .1)) 


#annotated k-means

#Clusters for significant points
sc_roots<-data.frame(t(annotation_roots))


#elbow method visualization
wss <- (nrow(sc_roots)-1)*sum(apply(sc_roots,2,var))
for (i in 1:4) wss[i] <- sum(kmeans(sc_roots,
                                    centers=i)$withinss,nstart=25,iter.max=1000)

wss_gf=data.frame(clust=c(1:4), wss=wss)

ggplot(wss_gf, aes(clust,wss))+
  geom_line()+theme_bw()+ xlab('Number of clusters') + ylab('WSS')+
  geom_vline(xintercept =3, linetype="dashed")+
  theme(axis.text.y = element_text(color='black', 
                                   size=12),
        axis.title.y=element_text(color="black", 
                                  size=16),
        panel.background = element_blank(), 
        axis.line = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(color='black', 
                                   size=12),
        axis.title.x = element_text( color="black", 
                                     size=16)
  ) 

#k-means clustering
roots_kmeans <- kmeans(sc_roots,centers=3)
roots_kmeans$cluster
sc_roots$name=rownames(sc_roots)

#k-means visualization
autoplot(roots_kmeans, 
         data=sc_roots, label = TRUE,frame = TRUE,
         frame.type = 'norm',label.size = 6,label.col='black',alpha=0)+
  scale_colour_manual(values = c('black','black','black'))+
  scale_fill_manual(values = c("#A9A9A9","#A60B0B",'red'))+
  theme_bw()+
  theme(axis.text.y = element_text(color='black', 
                                   size=24),
        axis.title.y=element_text(color="black", 
                                  size=28),
        panel.background = element_blank(), 
        axis.line = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(color='black', 
                                   size=24),
        axis.title.x = element_text(color="black", 
                                    size=28),
        legend.position = 'none'
  ) + 
  scale_x_continuous(expand = c(.1, .1)) 


#annotated + significant k-means

annotated_sigif_roots <- signif_pval_roots[signif_pval_roots$SSP %in% annotation_roots_df_spots,]

#Clusters for significant points
sc_roots<-data.frame(t(annotated_sigif_roots[2:10]))


#elbow method visualization
wss <- (nrow(sc_roots)-1)*sum(apply(sc_roots,2,var))
for (i in 1:4) wss[i] <- sum(kmeans(sc_roots,
                                    centers=i)$withinss,nstart=25,iter.max=1000)

wss_gf=data.frame(clust=c(1:4), wss=wss)

ggplot(wss_gf, aes(clust,wss))+
  geom_line()+theme_bw()+ xlab('Number of clusters') + ylab('WSS')+
  geom_vline(xintercept =2, linetype="dashed")+
  theme(axis.text.y = element_text(color='black', 
                                   size=12),
        axis.title.y=element_text(color="black", 
                                  size=16),
        panel.background = element_blank(), 
        axis.line = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(color='black', 
                                   size=12),
        axis.title.x = element_text( color="black", 
                                     size=16)
  ) 

#k-means clustering
roots_kmeans <- kmeans(sc_roots,centers=2)
roots_kmeans$cluster
sc_roots$name=rownames(sc_roots)

#k-means visualization
autoplot(roots_kmeans, 
         data=sc_roots, label = TRUE,frame = TRUE,
         frame.type = 'norm',label.size = 6,label.col='black',alpha=0)+
  scale_colour_manual(values = c('black','black','black'))+
  scale_fill_manual(values = c("#A9A9A9","#A60B0B",'red'))+
  theme_bw()+
  theme(axis.text.y = element_text(color='black', 
                                   size=24),
        axis.title.y=element_text(color="black", 
                                  size=28),
        panel.background = element_blank(), 
        axis.line = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(color='black', 
                                   size=24),
        axis.title.x = element_text(color="black", 
                                    size=28),
        legend.position = 'none'
  ) + 
  scale_x_continuous(expand = c(.1, .1)) 



## New heatmap All
rownames(redone_roots_normalized) <- redone_roots_normalized$spot_num
row_dend = as.dendrogram(hclust(dist((redone_roots_normalized[,2:10]))))
cluster_dend <- hclust(dist((redone_roots_normalized[,2:10])))
plot(cluster_dend)

clusterCut <- cutree(cluster_dend, k = 3)
names_frame <- data.frame(keyName=names(clusterCut), value=clusterCut, row.names=NULL) %>% arrange(value)
colnames(names_frame) <- c('Spot', 'Cluster')
names_frame$Control <- 0
names_frame$Anoxia <- 0
names_frame$Reaeration <- 0
names_frame$Annotation <- 0
names_frame$Accession <- 0

for (ind in 1:nrow(names_frame)){
  normalized_row <- c()
  normalized_row <- redone_roots_normalized[redone_roots_normalized$spot_num==names_frame[ind, 1], -1]
  names_frame[ind, 3] <- mean(as.numeric(normalized_row[4:6]))
  names_frame[ind, 4] <- mean(as.numeric(normalized_row[1:3]))
  names_frame[ind, 5] <- mean(as.numeric(normalized_row[7:9]))
  names_frame[ind, 6] <- normalized_row[11]
  names_frame[ind, 7] <- normalized_row[12]
}

write.xlsx(names_frame,'h_clusters_root_all.xlsx')

col_fun = colorRamp2(c(min(redone_roots_normalized[,2:10]), 
                       max(redone_roots_normalized[,2:10])-min(redone_roots_normalized[,2:10])+4,
                       max(redone_roots_normalized[,2:10])), c("#5B3D56", "#F8F5EF", "#D8B03E"))
col_fun(seq(min(redone_roots_normalized[,2:10]),  max(redone_roots_normalized[,2:10])))

Heatmap(redone_roots_normalized[,2:10], name = "mat",
                 heatmap_width  = unit(19, "cm"), heatmap_height = unit(27, "cm"),
                 cluster_rows = row_dend, 
        row_names_gp = gpar(fontsize = 9),
        column_names_gp = gpar(fontsize = 18),
        col = col_fun)

##Mean_acc heatmap

row.names(annotation_roots_annot_with_accessions_means) <- annotation_roots_annot_with_accessions_means$accessions

row_dend = as.dendrogram(hclust(dist((annotation_roots_annot_with_accessions_means[2:4]))))



col_fun = colorRamp2(c(min(annotation_roots_annot_with_accessions_means[2:4])+3, 
                       max(annotation_roots_annot_with_accessions_means[2:4])-min(annotation_roots_annot_with_accessions_means[2:4])+10,
                       max(annotation_roots_annot_with_accessions_means[2:4])), c("#F8F5EF","#5B3D56",  "#D8B03E"))# FCC119
col_fun(seq(min(annotation_roots_annot_with_accessions_means[2:4]),  max(annotation_roots_annot_with_accessions_means[2:4])))


Heatmap(annotation_roots_annot_with_accessions_means[2:4], name = "mat",
        heatmap_width  = unit(15, "cm"), heatmap_height = unit(15, "cm"),
        cluster_rows = row_dend, 
        row_names_gp = gpar(fontsize = 12),
        column_names_gp = gpar(fontsize = 18),
        col = col_fun,
        cluster_columns = FALSE)



##Signif heatmap
row_dend = as.dendrogram(hclust(dist((signif_pval_roots[2:10]))))
cluster_dend <- hclust(dist((signif_pval_roots[2:10])))
plot(cluster_dend)
#cutree(cluster_dend, k = 3)

clusterCut <- cutree(cluster_dend, h=14)
names_frame <- data.frame(keyName=names(clusterCut), value=clusterCut, row.names=NULL) %>% arrange(value)
colnames(names_frame) <- c('Spot', 'Cluster')
names_frame$Control <- 0
names_frame$Anoxia <- 0
names_frame$Reaeration <- 0
names_frame$Annotation <- 0
names_frame$Accession <- 0

for (ind in 1:nrow(names_frame)){
  normalized_row <- c()
  normalized_row <- signif_pval_roots[signif_pval_roots$spot_num==names_frame[ind, 1], -1]
  names_frame[ind, 3] <- mean(as.numeric(normalized_row[4:6]))
  names_frame[ind, 4] <- mean(as.numeric(normalized_row[1:3]))
  names_frame[ind, 5] <- mean(as.numeric(normalized_row[7:9]))
  names_frame[ind, 6] <- normalized_row[11]
  names_frame[ind, 7] <- normalized_row[12]
}

write.xlsx(names_frame,'h_clusters_roots_signif.xlsx')



col_fun = colorRamp2(c(min(signif_pval_roots[2:10]), 
                       max(signif_pval_roots[2:10])-min(signif_pval_roots[2:10])+4,
                       max(signif_pval_roots[2:10])), c("#5B3D56", "#F8F5EF", "#D8B03E"))
col_fun(seq(min(signif_pval_roots[2:10]),  max(signif_pval_roots[2:10])))

Heatmap(signif_pval_roots[2:10], name = "mat",
        heatmap_width  = unit(15, "cm"), heatmap_height = unit(15, "cm"),
        cluster_rows = row_dend, 
        row_names_gp = gpar(fontsize = 12),
        column_names_gp = gpar(fontsize = 18),
        col = col_fun,
        cluster_columns = FALSE)

##Annot heatmap

annotation_roots_annot <- redone_roots_normalized[redone_roots_normalized$SSP %in% annotation_roots_df_spots,]

row_dend = as.dendrogram(hclust(dist((annotation_roots_annot[2:10]))))
cluster_dend <- hclust(dist((annotation_roots_annot[2:10])))
plot(cluster_dend)
#cutree(cluster_dend, k = 3)

clusterCut <- cutree(cluster_dend, h=14)
names_frame <- data.frame(keyName=names(clusterCut), value=clusterCut, row.names=NULL) %>% arrange(value)
colnames(names_frame) <- c('Spot', 'Cluster')
names_frame$Control <- 0
names_frame$Anoxia <- 0
names_frame$Reaeration <- 0
names_frame$Annotation <- 0
names_frame$Accession <- 0

for (ind in 1:nrow(names_frame)){
  normalized_row <- c()
  normalized_row <- annotation_roots_annot[annotation_roots_annot$spot_num==names_frame[ind, 1], -1]
  names_frame[ind, 3] <- mean(as.numeric(normalized_row[4:6]))
  names_frame[ind, 4] <- mean(as.numeric(normalized_row[1:3]))
  names_frame[ind, 5] <- mean(as.numeric(normalized_row[7:9]))
  names_frame[ind, 6] <- normalized_row[11]
  names_frame[ind, 7] <- normalized_row[12]
}

write.xlsx(names_frame,'h_clusters_roots_annot.xlsx')


col_fun = colorRamp2(c(min(annotation_roots_annot[2:10]), 
                       max(annotation_roots_annot[2:10])-min(annotation_roots_annot[2:10])+4,
                       max(annotation_roots_annot[2:10])), c("#5B3D56", "#F8F5EF", "#D8B03E"))
col_fun(seq(min(annotation_roots_annot[2:10]),  max(annotation_roots_annot[2:10])))

Heatmap(annotation_roots_annot[2:10], name = "mat",
        heatmap_width  = unit(15, "cm"), heatmap_height = unit(15, "cm"),
        cluster_rows = row_dend, 
        row_names_gp = gpar(fontsize = 12),
        column_names_gp = gpar(fontsize = 18),
        col = col_fun)

##Annot signif heatmap
row_dend = as.dendrogram(hclust(dist((annotated_sigif_roots[2:10]))))
cluster_dend <- hclust(dist((annotated_sigif_roots[2:10])))
plot(cluster_dend)
#cutree(cluster_dend, k = 3)

clusterCut <- cutree(cluster_dend, h=13)
names_frame <- data.frame(keyName=names(clusterCut), value=clusterCut, row.names=NULL) %>% arrange(value)
colnames(names_frame) <- c('Spot', 'Cluster')
names_frame$Control <- 0
names_frame$Anoxia <- 0
names_frame$Reaeration <- 0
names_frame$Annotation <- 0
names_frame$Accession <- 0

for (ind in 1:nrow(names_frame)){
  normalized_row <- c()
  normalized_row <- annotated_sigif_roots[annotated_sigif_roots$spot_num==names_frame[ind, 1], -1]
  names_frame[ind, 3] <- mean(as.numeric(normalized_row[4:6]))
  names_frame[ind, 4] <- mean(as.numeric(normalized_row[1:3]))
  names_frame[ind, 5] <- mean(as.numeric(normalized_row[7:9]))
  names_frame[ind, 6] <- normalized_row[11]
  names_frame[ind, 7] <- normalized_row[12]
}

write.xlsx(names_frame,'h_clusters_roots_signif_annot.xlsx')



col_fun = colorRamp2(c(min(annotated_sigif_roots[2:10]), 
                       max(annotated_sigif_roots[2:10])-min(annotated_sigif_roots[2:10])+4,
                       max(annotated_sigif_roots[2:10])), c("#5B3D56", "#F8F5EF", "#D8B03E"))
col_fun(seq(min(annotated_sigif_roots[2:10]),  max(annotated_sigif_roots[2:10])))

Heatmap(annotated_sigif_roots[2:10], name = "mat",
        heatmap_width  = unit(15, "cm"), heatmap_height = unit(15, "cm"),
        cluster_rows = row_dend, 
        row_names_gp = gpar(fontsize = 12),
        column_names_gp = gpar(fontsize = 18),
        col = col_fun)


#K-means All spots 
sc_roots<-data.frame(redone_roots_normalized[2:10])


#elbow method visualization
wss <- (nrow(sc_roots)-1)*sum(apply(sc_roots,2,var))
for (i in 1:4) wss[i] <- sum(kmeans(sc_roots,
                                    centers=i)$withinss,nstart=25,iter.max=1000)

wss_gf=data.frame(clust=c(1:4), wss=wss)

ggplot(wss_gf, aes(clust,wss))+
  geom_line()+theme_bw()+ xlab('Number of clusters') + ylab('WSS')+
  geom_vline(xintercept =2, linetype="dashed")+
  theme(axis.text.y = element_text(color='black', 
                                   size=12),
        axis.title.y=element_text(color="black", 
                                  size=16),
        panel.background = element_blank(), 
        axis.line = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(color='black', 
                                   size=12),
        axis.title.x = element_text( color="black", 
                                     size=16)
  ) 

#k-means clustering
roots_kmeans <- kmeans(sc_roots,centers=2)
roots_kmeans$cluster
sc_roots$name=rownames(sc_roots)

#k-means visualization
autoplot(roots_kmeans, 
         data=sc_roots, label = TRUE,frame = TRUE,
         frame.type = 'norm',label.size = 6,label.col='black',alpha=0)+
  scale_colour_manual(values = c('black','black','black'))+
  scale_fill_manual(values = c("#A9A9A9","#A60B0B",'red'))+
  theme_bw()+
  theme(axis.text.y = element_text(color='black', 
                                   size=24),
        axis.title.y=element_text(color="black", 
                                  size=28),
        panel.background = element_blank(), 
        axis.line = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(color='black', 
                                   size=24),
        axis.title.x = element_text(color="black", 
                                    size=28),
        legend.position = 'none'
  ) + 
  scale_x_continuous(expand = c(.1, .1)) 



clusterCut <- roots_kmeans$cluster
names_frame <- data.frame(keyName=names(clusterCut), value=clusterCut, row.names=NULL) %>% arrange(value)
colnames(names_frame) <- c('Spot', 'Cluster')
names_frame$Control <- 0
names_frame$Anoxia <- 0
names_frame$Reaeration <- 0
names_frame$Annotation <- 0
names_frame$Accession <- 0

for (ind in 1:nrow(names_frame)){
  normalized_row <- c()
  normalized_row <- redone_roots_normalized[redone_roots_normalized$spot_num==names_frame[ind, 1], -1]
  names_frame[ind, 3] <- mean(as.numeric(normalized_row[4:6]))
  names_frame[ind, 4] <- mean(as.numeric(normalized_row[1:3]))
  names_frame[ind, 5] <- mean(as.numeric(normalized_row[7:9]))
  names_frame[ind, 6] <- normalized_row[11]
  names_frame[ind, 7] <- normalized_row[12]
}

write.xlsx(names_frame,'kmeans_clusters_root_all.xlsx')


#K-means signif spots
sc_roots<-data.frame(signif_pval_roots[2:10])


#elbow method visualization
wss <- (nrow(sc_roots)-1)*sum(apply(sc_roots,2,var))
for (i in 1:4) wss[i] <- sum(kmeans(sc_roots,
                                    centers=i)$withinss,nstart=25,iter.max=1000)

wss_gf=data.frame(clust=c(1:4), wss=wss)

ggplot(wss_gf, aes(clust,wss))+
  geom_line()+theme_bw()+ xlab('Number of clusters') + ylab('WSS')+
  geom_vline(xintercept =2, linetype="dashed")+
  theme(axis.text.y = element_text(color='black', 
                                   size=12),
        axis.title.y=element_text(color="black", 
                                  size=16),
        panel.background = element_blank(), 
        axis.line = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(color='black', 
                                   size=12),
        axis.title.x = element_text( color="black", 
                                     size=16)
  ) 

#k-means clustering
roots_kmeans <- kmeans(sc_roots,centers=2)
roots_kmeans$cluster
sc_roots$name=rownames(sc_roots)

#k-means visualization
autoplot(roots_kmeans, 
         data=sc_roots, label = TRUE,frame = TRUE,
         frame.type = 'norm',label.size = 6,label.col='black',alpha=0)+
  scale_colour_manual(values = c('black','black','black'))+
  scale_fill_manual(values = c("#A9A9A9","#A60B0B",'red'))+
  theme_bw()+
  theme(axis.text.y = element_text(color='black', 
                                   size=24),
        axis.title.y=element_text(color="black", 
                                  size=28),
        panel.background = element_blank(), 
        axis.line = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(color='black', 
                                   size=24),
        axis.title.x = element_text(color="black", 
                                    size=28),
        legend.position = 'none'
  ) + 
  scale_x_continuous(expand = c(.1, .1)) 



clusterCut <- roots_kmeans$cluster
names_frame <- data.frame(keyName=names(clusterCut), value=clusterCut, row.names=NULL) %>% arrange(value)
colnames(names_frame) <- c('Spot', 'Cluster')
names_frame$Control <- 0
names_frame$Anoxia <- 0
names_frame$Reaeration <- 0
names_frame$Annotation <- 0
names_frame$Accession <- 0

for (ind in 1:nrow(names_frame)){
  normalized_row <- c()
  normalized_row <- signif_pval_roots[signif_pval_roots$spot_num==names_frame[ind, 1], -1]
  names_frame[ind, 3] <- mean(as.numeric(normalized_row[4:6]))
  names_frame[ind, 4] <- mean(as.numeric(normalized_row[1:3]))
  names_frame[ind, 5] <- mean(as.numeric(normalized_row[7:9]))
  names_frame[ind, 6] <- normalized_row[11]
  names_frame[ind, 7] <- normalized_row[12]
}

write.xlsx(names_frame,'kmeans_clusters_root_signif.xlsx')


#K-means annot spots
sc_roots<-data.frame(annotation_roots_annot[2:10])


#elbow method visualization
wss <- (nrow(sc_roots)-1)*sum(apply(sc_roots,2,var))
for (i in 1:4) wss[i] <- sum(kmeans(sc_roots,
                                    centers=i)$withinss,nstart=25,iter.max=1000)

wss_gf=data.frame(clust=c(1:4), wss=wss)

ggplot(wss_gf, aes(clust,wss))+
  geom_line()+theme_bw()+ xlab('Number of clusters') + ylab('WSS')+
  geom_vline(xintercept =2, linetype="dashed")+
  theme(axis.text.y = element_text(color='black', 
                                   size=12),
        axis.title.y=element_text(color="black", 
                                  size=16),
        panel.background = element_blank(), 
        axis.line = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(color='black', 
                                   size=12),
        axis.title.x = element_text( color="black", 
                                     size=16)
  ) 

#k-means clustering
roots_kmeans <- kmeans(sc_roots,centers=2)
roots_kmeans$cluster
sc_roots$name=rownames(sc_roots)

#k-means visualization
autoplot(roots_kmeans, 
         data=sc_roots, label = TRUE,frame = TRUE,
         frame.type = 'norm',label.size = 6,label.col='black',alpha=0)+
  scale_colour_manual(values = c('black','black','black'))+
  scale_fill_manual(values = c("#A9A9A9","#A60B0B",'red'))+
  theme_bw()+
  theme(axis.text.y = element_text(color='black', 
                                   size=24),
        axis.title.y=element_text(color="black", 
                                  size=28),
        panel.background = element_blank(), 
        axis.line = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(color='black', 
                                   size=24),
        axis.title.x = element_text(color="black", 
                                    size=28),
        legend.position = 'none'
  ) + 
  scale_x_continuous(expand = c(.1, .1)) 



clusterCut <- roots_kmeans$cluster
names_frame <- data.frame(keyName=names(clusterCut), value=clusterCut, row.names=NULL) %>% arrange(value)
colnames(names_frame) <- c('Spot', 'Cluster')
names_frame$Control <- 0
names_frame$Anoxia <- 0
names_frame$Reaeration <- 0
names_frame$Annotation <- 0
names_frame$Accession <- 0

for (ind in 1:nrow(names_frame)){
  normalized_row <- c()
  normalized_row <- annotation_roots_annot[annotation_roots_annot$spot_num==names_frame[ind, 1], -1]
  names_frame[ind, 3] <- mean(as.numeric(normalized_row[4:6]))
  names_frame[ind, 4] <- mean(as.numeric(normalized_row[1:3]))
  names_frame[ind, 5] <- mean(as.numeric(normalized_row[7:9]))
  names_frame[ind, 6] <- normalized_row[11]
  names_frame[ind, 7] <- normalized_row[12]
}

write.xlsx(names_frame,'kmeans_clusters_root_annot.xlsx')


#K-means signif annot spots
sc_roots<-data.frame(annotated_sigif_roots[2:10])


#elbow method visualization
wss <- (nrow(sc_roots)-1)*sum(apply(sc_roots,2,var))
for (i in 1:4) wss[i] <- sum(kmeans(sc_roots,
                                    centers=i)$withinss,nstart=25,iter.max=1000)

wss_gf=data.frame(clust=c(1:4), wss=wss)

ggplot(wss_gf, aes(clust,wss))+
  geom_line()+theme_bw()+ xlab('Number of clusters') + ylab('WSS')+
  geom_vline(xintercept =2, linetype="dashed")+
  theme(axis.text.y = element_text(color='black', 
                                   size=12),
        axis.title.y=element_text(color="black", 
                                  size=16),
        panel.background = element_blank(), 
        axis.line = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(color='black', 
                                   size=12),
        axis.title.x = element_text( color="black", 
                                     size=16)
  ) 

#k-means clustering
roots_kmeans <- kmeans(sc_roots,centers=2)
roots_kmeans$cluster
sc_roots$name=rownames(sc_roots)

#k-means visualization
autoplot(roots_kmeans, 
         data=sc_roots, label = TRUE,frame = TRUE,
         frame.type = 'norm',label.size = 6,label.col='black',alpha=0)+
  scale_colour_manual(values = c('black','black','black'))+
  scale_fill_manual(values = c("#A9A9A9","#A60B0B",'red'))+
  theme_bw()+
  theme(axis.text.y = element_text(color='black', 
                                   size=24),
        axis.title.y=element_text(color="black", 
                                  size=28),
        panel.background = element_blank(), 
        axis.line = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(color='black', 
                                   size=24),
        axis.title.x = element_text(color="black", 
                                    size=28),
        legend.position = 'none'
  ) + 
  scale_x_continuous(expand = c(.1, .1)) 



clusterCut <- roots_kmeans$cluster
names_frame <- data.frame(keyName=names(clusterCut), value=clusterCut, row.names=NULL) %>% arrange(value)
colnames(names_frame) <- c('Spot', 'Cluster')
names_frame$Control <- 0
names_frame$Anoxia <- 0
names_frame$Reaeration <- 0
names_frame$Annotation <- 0
names_frame$Accession <- 0

for (ind in 1:nrow(names_frame)){
  normalized_row <- c()
  normalized_row <- annotated_sigif_roots[annotated_sigif_roots$spot_num==names_frame[ind, 1], -1]
  names_frame[ind, 3] <- mean(as.numeric(normalized_row[4:6]))
  names_frame[ind, 4] <- mean(as.numeric(normalized_row[1:3]))
  names_frame[ind, 5] <- mean(as.numeric(normalized_row[7:9]))
  names_frame[ind, 6] <- normalized_row[11]
  names_frame[ind, 7] <- normalized_row[12]
}

write.xlsx(names_frame,'kmeans_clusters_root_signif_annot.xlsx')



#____________________________ shoots _____________________

#read shoots data
shoots_raw_dat <- read.xlsx('PQquest_data/rice_shoots_exported.xlsx', sheetIndex = 1)
redone_shoots <- shoots_raw_dat %>%  dplyr::select("SSP","Gel.Name", "Quantity" ) %>% reshape(idvar = "SSP", timevar = "Gel.Name", direction = "wide")
old_shoots_names <- colnames(redone_shoots)
new_colnames=c()

inds_col=c('1','2','3','4','5','6','7','1','2','3','4','5','6','7','8','1','2','3','4','5','6','7')
i=1
for (name in colnames(redone_shoots)[2:23]){
  split1<-strsplit(name,split='.', fixed=TRUE)[[1]][2]
  split2<-strsplit(split1,split='_', fixed=TRUE)[[1]][1]
  new_colnames <- c(new_colnames,paste(split2,inds_col[i], sep= '_'))
  i<-i+1
}
#extract values
colnames(redone_shoots)[2:23] <- new_colnames
colnames(redone_shoots)[17] <- 'reaeration_1'
redone_shoots <- redone_shoots[,c(1,6,7,8,12,13,17,21)]
colnames(redone_shoots)[2:8] <- c('A1','A2','A3','C1','C2','R1','R2')
#selected gels: A 5-8, C 9,12,13 R 17,18,21

#add manual spots
new_shoots = read.csv('PQquest_data/new_spots_shoots.csv', sep=";", na.strings=c("","NA"), header = TRUE)
redone_shoots <- rbind(redone_shoots, new_shoots)

#take logarithm
redone_shoots_log <- cbind(redone_shoots[1], log2(redone_shoots[2:8]))
redone_shoots_log_reshaped <- gather(redone_shoots_log[2:8], condition, measurement, A1:R2, factor_key=TRUE)
redone_shoots_log_reshaped$factor_col <- 'condition'
redone_shoots_log_reshaped[grepl("A", redone_shoots_log_reshaped$condition), 3] <- 'Anoxia'
redone_shoots_log_reshaped[grepl("C", redone_shoots_log_reshaped$condition), 3] <- 'Control'
redone_shoots_log_reshaped[grepl("R", redone_shoots_log_reshaped$condition), 3] <- 'Re-aeration'


#NON Rubusco
for (ind in 1:nrow(redone_shoots)){
  print(ind)
  print( redone_shoots[ind, 1])
}

#37, 38, 73, 74
redone_shoots_log_no_rubisco <- cbind(redone_shoots[-c(37, 38, 73, 74),1], log2(redone_shoots[-c(37, 38, 73, 74),2:8]))
redone_shoots_log_no_rubisco_reshaped <- gather(redone_shoots_log_no_rubisco[2:8], condition, measurement, A1:R2, factor_key=TRUE)
redone_shoots_log_no_rubisco_reshaped$factor_col <- 'condition'
redone_shoots_log_no_rubisco_reshaped[grepl("A", redone_shoots_log_no_rubisco_reshaped$condition), 3] <- 'Anoxia'
redone_shoots_log_no_rubisco_reshaped[grepl("C", redone_shoots_log_no_rubisco_reshaped$condition), 3] <- 'Control'
redone_shoots_log_no_rubisco_reshaped[grepl("R", redone_shoots_log_no_rubisco_reshaped$condition), 3] <- 'Re-aeration'


redone_shoots_log_no_rubisco_reshaped$factor_col <- factor(redone_shoots_log_no_rubisco_reshaped$factor_col, 
                                                      levels=c('Control', 'Anoxia','Re-aeration'))

redone_shoots_log_no_rubisco_reshaped$condition <- factor(redone_shoots_log_no_rubisco_reshaped$condition ,
                                                     levels= c('C1','C2','C3','A1','A2','A3','R1','R2'))

#non-normalized box-plots
ggplot(redone_shoots_log_no_rubisco_reshaped, aes(y=measurement, x=condition, fill=factor_col))+
  geom_boxplot(alpha=0.9)+theme_bw()+
  scale_fill_manual(values = c("#678667","#375DB6","#A5432B"))+
  xlab('Sample') + ylab('Optical density')+
  theme(axis.text.y = element_text(color='black', 
                                   size=16),
        axis.title.y=element_text(color="black", 
                                  size=18),
        panel.background = element_blank(), 
        axis.line = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(color='black', 
                                   size=16),
        axis.title.x = element_text( color="black", 
                                     size=18),
        legend.title=element_text(size=18), 
        legend.text=element_text(size=16)
  )+ 
  guides(fill = guide_legend(override.aes = list(size = 1),title="Condition")) 


#mean OD values
anoxia = redone_shoots_log_reshaped[redone_shoots_log_reshaped$factor_col=='Anoxia', 2]  #19.11415 - anoxia
control = redone_shoots_log_reshaped[redone_shoots_log_reshaped$factor_col=='Control', 2] #20.38731 - control
reaeration <- redone_shoots_log_reshaped[redone_shoots_log_reshaped$factor_col=='Re-aeration', 2] #19.09779 - re-aeration
mean_OD_df <- data.frame(Density = c(anoxia, control, reaeration), 
                         Condition=c(rep('anoxia',222),rep('control',148),rep('reaeration',148)))

t_test_df=data.frame(comparision=c("Anoxia_vs_Control","Anoxia_vs_Re-aeration", 'Control_vs_Re-aeration'), p.value=1)
t_test_df[1, 2] <- t.test(anoxia, control)$p.value
t_test_df[2, 2] <- t.test(anoxia, reaeration)$p.value
t_test_df[3, 2] <- t.test(control, reaeration)$p.value
t_test_df$p.adj <- p.adjust(t_test_df$p.value, method = "fdr")
write.xlsx(t_test_df,'t_test_OD_shoots.xlsx')


#### exclude rubisco

anoxia_no_rubisco = redone_shoots_log_no_rubisco_reshaped[redone_shoots_log_no_rubisco_reshaped$factor_col=='Anoxia', 2]  #19.11415 - anoxia
control_no_rubisco = redone_shoots_log_no_rubisco_reshaped[redone_shoots_log_no_rubisco_reshaped$factor_col=='Control', 2] #20.38731 - control
reaeration_no_rubisco <- redone_shoots_log_no_rubisco_reshaped[redone_shoots_log_no_rubisco_reshaped$factor_col=='Re-aeration', 2] #19.09779 - re-aeration
mean_OD_df_no_rubisco <- data.frame(Density = c(anoxia_no_rubisco, control_no_rubisco, reaeration_no_rubisco), 
                         Condition=c(rep('anoxia',210),rep('control',140),rep('reaeration',140)))
length(anoxia_no_rubisco) #210
length(control_no_rubisco) #140
length(reaeration_no_rubisco) #140

mean(anoxia_no_rubisco) #18.93872
mean(control_no_rubisco) #20.28655
mean(reaeration_no_rubisco) #19.03892


t_test_df_no_rubisco=data.frame(comparision=c("Anoxia_vs_Control","Anoxia_vs_Re-aeration", 'Control_vs_Re-aeration'), p.value=1)
t_test_df_no_rubisco[1, 2] <- t.test(anoxia_no_rubisco, control_no_rubisco)$p.value
t_test_df_no_rubisco[2, 2] <- t.test(anoxia_no_rubisco, reaeration_no_rubisco)$p.value
t_test_df_no_rubisco[3, 2] <- t.test(control_no_rubisco, reaeration_no_rubisco)$p.value
t_test_df_no_rubisco$p.adj <- p.adjust(t_test_df_no_rubisco$p.value, method = "fdr")
write.xlsx(t_test_df_no_rubisco,'t_test_OD_shoots_no_rubisco.xlsx')



ggplot(redone_shoots_log_no_rubisco_reshaped, aes(y=measurement, x=factor_col, fill=factor_col))+
  geom_violin(alpha=0.9, color='black')+theme_bw()+
  scale_fill_manual(values = c("#678667","#375DB6","#A5432B"))+
  geom_hline(yintercept=18.93872, size=2, linetype='dashed', color='#375DB6', alpha=0.6)+
  geom_hline(yintercept=20.28655, size=2, linetype='dashed', color='#678667')+
  geom_hline(yintercept=19.03892, size=2, linetype='dashed', color='#A5432B', alpha=0.6)+
  xlab('Condition') + ylab('Optical density')+
  theme(axis.text.y = element_text(color='black', 
                                   size=16),
        axis.title.y=element_text(color="black", 
                                  size=18),
        panel.background = element_blank(), 
        axis.line = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(color='black', 
                                   size=16),
        axis.title.x = element_text( color="black", 
                                     size=18),
        legend.title=element_text(size=18), 
        legend.text=element_text(size=16)
  )+ 
  guides(fill = guide_legend(override.aes = list(size = 1),title="Condition"))



#normalize data
redone_shoots_normalized <- cbind(redone_shoots_log[1],  normalizeQuantiles(redone_shoots_log[2:8]))
redone_shoots_normalized_reshaped <- gather(redone_shoots_normalized[2:8], condition, measurement, A1:R2, factor_key=TRUE)
redone_shoots_normalized_reshaped$factor_col <- 'condition'
redone_shoots_normalized_reshaped[grepl("A", redone_shoots_normalized_reshaped$condition), 3] <- 'Anoxia'
redone_shoots_normalized_reshaped[grepl("C", redone_shoots_normalized_reshaped$condition), 3] <- 'Control'
redone_shoots_normalized_reshaped[grepl("R", redone_shoots_normalized_reshaped$condition), 3] <- 'Re-aeration'


#No Rubisco

redone_shoots_normalized_no_rubisco <- cbind(redone_shoots_log_no_rubisco[1],  normalizeQuantiles(redone_shoots_log_no_rubisco[2:8]))
redone_shoots_normalized_reshaped_no_rubisco <- gather(redone_shoots_normalized_no_rubisco[2:8], condition, measurement, A1:R2, factor_key=TRUE)
redone_shoots_normalized_reshaped_no_rubisco$factor_col <- 'condition'
redone_shoots_normalized_reshaped_no_rubisco[grepl("A", redone_shoots_normalized_reshaped_no_rubisco$condition), 3] <- 'Anoxia'
redone_shoots_normalized_reshaped_no_rubisco[grepl("C", redone_shoots_normalized_reshaped_no_rubisco$condition), 3] <- 'Control'
redone_shoots_normalized_reshaped_no_rubisco[grepl("R", redone_shoots_normalized_reshaped_no_rubisco$condition), 3] <- 'Re-aeration'



redone_shoots_normalized_reshaped_no_rubisco$factor_col <- factor(redone_shoots_normalized_reshaped_no_rubisco$factor_col, 
                                                           levels=c('Control', 'Anoxia','Re-aeration'))

redone_shoots_normalized_reshaped_no_rubisco$condition <- factor(redone_shoots_normalized_reshaped_no_rubisco$condition ,
                                                          levels= c('C1','C2','C3','A1','A2','A3','R1','R2'))

#normalized box-plots
ggplot(redone_shoots_normalized_reshaped_no_rubisco, aes(y=measurement, x=condition, fill=factor_col))+
  geom_boxplot(alpha=0.9)+theme_bw()+
  scale_fill_manual(values = c("#678667","#375DB6","#A5432B"))+
  xlab('Sample') + ylab('Optical density')+
  theme(axis.text.y = element_text(color='black', 
                                   size=16),
        axis.title.y=element_text(color="black", 
                                  size=18),
        panel.background = element_blank(), 
        axis.line = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(color='black', 
                                   size=16),
        axis.title.x = element_text( color="black", 
                                     size=18),
        legend.title=element_text(size=18), 
        legend.text=element_text(size=16)
  )+ 
  guides(fill = guide_legend(override.aes = list(size = 1),title="Condition")) 



#annotate points
annotation_shoots_df_no_rubisco <-  read.csv('PQquest_data/PD_quest_shoot_annot_no_rubisco.csv', sep=";", na.strings=c("","NA"), header = TRUE)
annotation_shoots_df_spots_no_rubisco <-annotation_shoots_df_no_rubisco[annotation_shoots_df_no_rubisco$Accession!='-',1]

annotation_shoots_no_rubisco <- redone_shoots_normalized_no_rubisco[redone_shoots_normalized_no_rubisco$SSP %in% annotation_shoots_df_spots_no_rubisco,-1]
annotation_shoots_annot_no_rubisco <- redone_shoots_normalized_no_rubisco[redone_shoots_normalized_no_rubisco$SSP %in% annotation_shoots_df_spots_no_rubisco,]


redone_shoots_normalized_no_rubisco$spot_num <- '0'
redone_shoots_normalized_no_rubisco$spot_acc <- '0'
redone_shoots_normalized_no_rubisco$spot_annot <- '0'

annotation_shoots_df <-  read.csv('PQquest_data/PD_quest_shoot_annot.csv', sep=";", na.strings=c("","NA"), header = TRUE)
annotation_shoots_df_spots <-annotation_shoots_df[annotation_shoots_df$Accession!='-',1]
annotation_shoots <- redone_shoots_normalized[redone_shoots_normalized$SSP %in% annotation_shoots_df_spots,-1]


redone_shoots_normalized$spot_num <- '0'
redone_shoots_normalized$spot_acc <- '0'
redone_shoots_normalized$spot_annot <- '0'


for(ind in 1:nrow(annotation_shoots_df)){
  spot=annotation_shoots_df[ind, 1]
  redone_shoots_normalized[redone_shoots_normalized$SSP==spot, 9] <- annotation_shoots_df[ind,3]
  redone_shoots_normalized[redone_shoots_normalized$SSP==spot, 10] <- annotation_shoots_df[ind,2]
  redone_shoots_normalized[redone_shoots_normalized$SSP==spot, 11] <- annotation_shoots_df[ind,4]
}
rownames(redone_shoots_normalized) <- redone_shoots_normalized$spot_num

write.xlsx(redone_shoots_normalized,'shoots_with_new_log.xlsx')


annotation_shoots_annot_with_accessions <- redone_shoots_normalized[redone_shoots_normalized$SSP %in% annotation_shoots_df_spots,-c(1,9,11)]
annotation_shoots_annot_with_accessions$fixed_name <- 'name'

for (ind in 1:nrow(annotation_shoots_annot_with_accessions)){
  acc_str=annotation_shoots_annot_with_accessions[ind,8]
  acc_split=strsplit(acc_str,split='[', fixed=TRUE)[[1]][1]
  acc_split=strsplit(acc_split,split=', chloroplastic', fixed=TRUE)[[1]][1]
  annotation_shoots_annot_with_accessions[ind,9] <- acc_split
}

annotation_shoots_annot_with_accessions_means <- data.frame(control=1:41, anoxia=1:41, reaeration=1:41)
annotation_shoots_annot_with_accessions_means$anoxia <- rowMeans(annotation_shoots_annot_with_accessions[,c(1,2,3)])
annotation_shoots_annot_with_accessions_means$control <- rowMeans(annotation_shoots_annot_with_accessions[,c(4,5)])
annotation_shoots_annot_with_accessions_means$reaeration <- rowMeans(annotation_shoots_annot_with_accessions[,c(6,7)])
annotation_shoots_annot_with_accessions_means$accessions <- annotation_shoots_annot_with_accessions$fixed_name

annotation_shoots_annot_with_accessions_means <- annotation_shoots_annot_with_accessions_means %>% 
  group_by(accessions) %>% 
  summarize(OD_control=mean(control), OD_anoxia=mean(anoxia), OD_reaeration =mean(reaeration)) %>% as.data.frame()
annotation_shoots_annot_with_accessions_means <- annotation_shoots_annot_with_accessions_means[-c(3,5,6,8,11,14,15),]


redone_roots_normalized$spot_num <- '0'
redone_roots_normalized$spot_acc <- '0'
redone_roots_normalized$spot_annot <- '0'

for(ind in 1:nrow(annotation_roots_df)){
  spot=annotation_roots_df[ind, 1]
  redone_roots_normalized[redone_roots_normalized$SSP==spot, 11] <- annotation_roots_df[ind,3]
  redone_roots_normalized[redone_roots_normalized$SSP==spot, 12] <- annotation_roots_df[ind,2]
  redone_roots_normalized[redone_roots_normalized$SSP==spot, 13] <- annotation_roots_df[ind,4]
}
rownames(redone_roots_normalized) <- redone_roots_normalized$spot_num



#All clustering Conditions
#elbow method optimization
sc_shoots<-data.frame(t(redone_shoots_normalized[2:8]))


wss <- (nrow(sc_shoots)-1)*sum(apply(sc_shoots,2,var))
for (i in 1:4) wss[i] <- sum(kmeans(sc_shoots,
                                    centers=i)$withinss,nstart=25,iter.max=1000)

wss_gf=data.frame(clust=c(1:4), wss=wss)

ggplot(wss_gf, aes(clust,wss))+
  geom_line()+theme_bw()+ xlab('Number of clusters') + ylab('WSS')+
  geom_vline(xintercept =3, linetype="dashed")+
  theme(axis.text.y = element_text(color='black', 
                                   size=12),
        axis.title.y=element_text(color="black", 
                                  size=16),
        panel.background = element_blank(), 
        axis.line = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(color='black', 
                                   size=12),
        axis.title.x = element_text( color="black", 
                                     size=16)
  ) 


#All clustering Conditions no rubisco
#elbow method optimization
sc_shoots<-data.frame(t(redone_shoots_normalized_no_rubisco[2:8]))


wss <- (nrow(sc_shoots)-1)*sum(apply(sc_shoots,2,var))
for (i in 1:4) wss[i] <- sum(kmeans(sc_shoots,
                                    centers=i)$withinss,nstart=25,iter.max=1000)

wss_gf=data.frame(clust=c(1:4), wss=wss)

ggplot(wss_gf, aes(clust,wss))+
  geom_line()+theme_bw()+ xlab('Number of clusters') + ylab('WSS')+
  geom_vline(xintercept =3, linetype="dashed")+
  theme(axis.text.y = element_text(color='black', 
                                   size=12),
        axis.title.y=element_text(color="black", 
                                  size=16),
        panel.background = element_blank(), 
        axis.line = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(color='black', 
                                   size=12),
        axis.title.x = element_text( color="black", 
                                     size=16)
  ) 


#k-means clusterization
shoots_kmeans <- kmeans(sc_shoots,centers=3)
shoots_kmeans$cluster
sc_shoots$name=rownames(sc_shoots)

autoplot(shoots_kmeans, 
         data=sc_shoots, label = TRUE,frame = TRUE,
         frame.type = 'norm',label.size = 6,label.col='black',alpha=0)+
  scale_colour_manual(values = c('black','black','black'))+
  scale_fill_manual(values = c("#A9A9A9","#A60B0B",'red'))+
  theme_bw()+
  theme(axis.text.y = element_text(color='black', 
                                   size=24),
        axis.title.y=element_text(color="black", 
                                  size=28),
        panel.background = element_blank(), 
        axis.line = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(color='black', 
                                   size=24),
        axis.title.x = element_text(color="black", 
                                    size=28),
        legend.position = 'none'
  ) + 
  scale_x_continuous(expand = c(.1, .1)) +
scale_y_continuous(expand = c(.1, .1))
write.xlsx(redone_roots,'shoots_with_new.xlsx')


#redo functions
k_p_val <- function(x) { vals_df <- x %>% gather(condition, measurement, A1:R2, factor_key=TRUE)
vals_df$factor_col <- 'factor'
vals_df[grepl("A", vals_df$condition), 3] <- 'Anoxia'
vals_df[grepl("C", vals_df$condition), 3] <- 'Control'
vals_df[grepl("R", vals_df$condition), 3] <- 'Re-aeration'
tryCatch(kruskal.test(measurement ~ factor_col, data = vals_df)$p.value,
         error = function(e) NA)
}

#mixed linear model function
lme_p_val <- function(x) { vals_df <- x %>% gather(condition, measurement, A1:R2, factor_key=TRUE)
vals_df$factor_col <- 'factor'
vals_df[grepl("A", vals_df$condition), 3] <- 'Anoxia'
vals_df[grepl("C", vals_df$condition), 3] <- 'Control'
vals_df[grepl("R", vals_df$condition), 3] <- 'Re-aeration'
fit.lme <- lme(measurement ~ factor_col, random=~1|condition, data=vals_df)
ret_vec <- c(as.data.frame(anova(fit.lme))$`p-value`[2])
fit_sum <- as.data.frame(summary(fit.lme)$tTable)
ret_vec <- c(ret_vec, fit_sum$`p-value`[c(2,3)])
tryCatch(ret_vec,
         error = function(e) NA)
}

pvals_shoot_kraskell <- c()
pvals_shoot_lme <- data.frame(all_pval=c(1), cont_pval = c(1), re_pval = c(1))


#find p-vals for tests
for (i in 1:nrow(redone_shoots_normalized)){
  pvals_shoot_kraskell <- c(pvals_shoot_kraskell, k_p_val(redone_shoots_normalized[i,2:8]))
  pvals_shoot_lme <- rbind(pvals_shoot_lme, lme_p_val(redone_shoots_normalized[i,2:8]))
}
pvals_shoot_lme <- pvals_shoot_lme[-1,]

#adjust p-values
pval_df_shoots=data_frame(spot=redone_shoots_normalized[,1], kraskell_raw=pvals_shoot_kraskell, 
                   kraskell_adj = p.adjust(pvals_shoot_kraskell, method = "fdr"),
                   lme_all_raw = pvals_shoot_lme[,1],
                   lme_all_adj = p.adjust(pvals_shoot_lme[,1], method = "fdr"),
                   lme_cont_raw = pvals_shoot_lme[,2],
                   lme_cont_adj = p.adjust(pvals_shoot_lme[,2], method = "fdr"),
                   lme_re_raw = pvals_shoot_lme[,3],
                   lme_re_adj = p.adjust(pvals_shoot_lme[,3], method = "fdr"))

#adjustment results
#FDR
sum(pval_df_shoots$kraskell_raw <= 0.05, na.rm = TRUE) #0 proteins with significant 
sum(pval_df_shoots$lme_all_raw <= 0.05, na.rm = TRUE) #17 proteins with significant
sum(pval_df_shoots$lme_cont_raw <= 0.05, na.rm = TRUE) #12 proteins with significant 
sum(pval_df_shoots$lme_re_raw <= 0.05, na.rm = TRUE) #15 proteins with significant
sum(pval_df_shoots$kraskell_adj <= 0.05, na.rm = TRUE) #0 proteins with significant 
sum(pval_df_shoots$lme_all_adj <= 0.05, na.rm = TRUE) #5 proteins with significant
sum(pval_df_shoots$lme_cont_adj <= 0.05, na.rm = TRUE) #2 proteins with significant 
sum(pval_df_shoots$lme_re_adj <= 0.05, na.rm = TRUE) #3 proteins with significant

#holm
sum(pval_df$kraskell_adj <= 0.05, na.rm = TRUE) #0 proteins with significant 
sum(pval_df$lme_all_adj <= 0.05, na.rm = TRUE) #0 proteins with significant
sum(pval_df$lme_cont_adj <= 0.05, na.rm = TRUE) #0 proteins with significant 
sum(pval_df$lme_re_adj <= 0.05, na.rm = TRUE) #0 proteins with significant

#BH
sum(pval_df$kraskell_adj <= 0.05, na.rm = TRUE) #0 proteins with significant 
sum(pval_df$lme_all_adj <= 0.05, na.rm = TRUE) #0 proteins with significant
sum(pval_df$lme_cont_adj <= 0.05, na.rm = TRUE) #0 proteins with significant 
sum(pval_df$lme_re_adj <= 0.05, na.rm = TRUE) #0 proteins with significant

#Modded_linear_model
shoot.fac <- data.frame(Condition=c('anoxia','anoxia','anoxia','control','control', 're-aeration', 're-aeration'))
rownames(shoot.fac) <- colnames( redone_shoots_normalized[,-1])
rownames(redone_shoots_normalized) <- gsub(" ", "_", redone_shoots_normalized[,1])

X <- model.matrix(~ Condition, data = shoot.fac)
fit <- lmFit(redone_shoots_normalized[,-1], design = X, method = "robust", maxit = 10000)
efit <- eBayes(fit)

#modelled_p-values adjustment
num_spots <- nrow(redone_shoots_normalized)
pval_df_shoots$eBayes_cont_raw <- topTable(efit, coef = 2, number = num_spots,
                                    sort.by = "none")$P.Value
pval_df_shoots$eBayes_cont_adj <- topTable(efit, coef = 2, number = num_spots,
                                    sort.by = "none")$adj.P.Val
pval_df_shoots$eBayes_re_raw <- topTable(efit, coef = 3, number = num_spots,
                                  sort.by = "none")$P.Value
pval_df_shoots$eBayes_re_adj <- topTable(efit, coef = 3, number = num_spots,
                                  sort.by = "none")$adj.P.Val

#the number of signifficant proteins
sum(pval_df_shoots$eBayes_cont_adj<=0.05) #14
sum(pval_df_shoots$eBayes_cont_raw<=0.05) #22
sum(pval_df_shoots$eBayes_re_adj<=0.05) #21
sum(pval_df_shoots$eBayes_re_raw<=0.05) #26

write.xlsx(pval_df,'dots_pvals_shoots.xlsx')

#extract significant spots
signif_pval_shoots <- redone_shoots_normalized[pval_df_shoots$eBayes_cont_adj<=0.05 | pval_df_shoots$eBayes_re_adj<=0.05,]
rownames(signif_pval_shoots) <- redone_shoots_normalized[pval_df$eBayes_cont_adj<=0.05 | pval_df$eBayes_re_adj<=0.05,9][-c(26,27)]

write.xlsx(signif_pval_shoots,'rice_shoots_signif_spots_annot.xlsx')

#K-means Clusters for significant points
sc_shoots<-data.frame(t(signif_pval_shoots[2:8]))


#elbow method visualization
wss <- (nrow(sc_shoots)-1)*sum(apply(sc_shoots,2,var))
for (i in 1:4) wss[i] <- sum(kmeans(sc_shoots,
                                    centers=i)$withinss,nstart=25,iter.max=1000)

wss_gf=data.frame(clust=c(1:4), wss=wss)

ggplot(wss_gf, aes(clust,wss))+
  geom_line()+theme_bw()+ xlab('Number of clusters') + ylab('WSS')+
  geom_vline(xintercept =3, linetype="dashed")+
  theme(axis.text.y = element_text(color='black', 
                                   size=12),
        axis.title.y=element_text(color="black", 
                                  size=16),
        panel.background = element_blank(), 
        axis.line = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(color='black', 
                                   size=12),
        axis.title.x = element_text( color="black", 
                                     size=16)
  ) 

#k-means clustering
shoots_kmeans <- kmeans(sc_shoots,centers=3)
shoots_kmeans$cluster
sc_shoots$name=rownames(sc_shoots)

#k-means visualization
autoplot(shoots_kmeans, 
         data=sc_shoots, label = TRUE,frame = TRUE,
         frame.type = 'norm',label.size = 6,label.col='black',alpha=0)+
  scale_colour_manual(values = c('black','black','black'))+
  scale_fill_manual(values = c("#A9A9A9","#A60B0B",'red'))+
  theme_bw()+
  theme(axis.text.y = element_text(color='black', 
                                   size=24),
        axis.title.y=element_text(color="black", 
                                  size=28),
        panel.background = element_blank(), 
        axis.line = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(color='black', 
                                   size=24),
        axis.title.x = element_text(color="black", 
                                    size=28),
        legend.position = 'none'
  ) + 
  scale_x_continuous(expand = c(.1, .1)) 



#annotated k-means

#Clusters for significant points
sc_shoots<-data.frame(t(annotation_shoots))


#elbow method visualization
wss <- (nrow(sc_shoots)-1)*sum(apply(sc_shoots,2,var))
for (i in 1:4) wss[i] <- sum(kmeans(sc_shoots,
                                    centers=i)$withinss,nstart=25,iter.max=1000)

wss_gf=data.frame(clust=c(1:4), wss=wss)

ggplot(wss_gf, aes(clust,wss))+
  geom_line()+theme_bw()+ xlab('Number of clusters') + ylab('WSS')+
  geom_vline(xintercept =3, linetype="dashed")+
  theme(axis.text.y = element_text(color='black', 
                                   size=12),
        axis.title.y=element_text(color="black", 
                                  size=16),
        panel.background = element_blank(), 
        axis.line = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(color='black', 
                                   size=12),
        axis.title.x = element_text( color="black", 
                                     size=16)
  ) 

#k-means clustering
shoots_kmeans <- kmeans(sc_shoots,centers=3)
shoots_kmeans$cluster
sc_shoots$name=rownames(sc_shoots)

#k-means visualization
autoplot(shoots_kmeans, 
         data=sc_shoots, label = TRUE,frame = TRUE,
         frame.type = 'norm',label.size = 6,label.col='black',alpha=0)+
  scale_colour_manual(values = c('black','black','black'))+
  scale_fill_manual(values = c("#A9A9A9","#A60B0B",'red'))+
  theme_bw()+
  theme(axis.text.y = element_text(color='black', 
                                   size=24),
        axis.title.y=element_text(color="black", 
                                  size=28),
        panel.background = element_blank(), 
        axis.line = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(color='black', 
                                   size=24),
        axis.title.x = element_text(color="black", 
                                    size=28),
        legend.position = 'none'
  ) + 
  scale_x_continuous(expand = c(.1, .1)) 


#annotated + significant k-means
annotated_sigif_shoots <- signif_pval_shoots[signif_pval_shoots$SSP %in% annotation_shoots_df_spots,]

#Clusters for significant points
sc_shoots<-data.frame(t(annotated_sigif_shoots[2:8]))


#elbow method visualization
wss <- (nrow(sc_shoots)-1)*sum(apply(sc_shoots,2,var))
for (i in 1:4) wss[i] <- sum(kmeans(sc_shoots,
                                    centers=i)$withinss,nstart=25,iter.max=1000)

wss_gf=data.frame(clust=c(1:4), wss=wss)

ggplot(wss_gf, aes(clust,wss))+
  geom_line()+theme_bw()+ xlab('Number of clusters') + ylab('WSS')+
  geom_vline(xintercept =3, linetype="dashed")+
  theme(axis.text.y = element_text(color='black', 
                                   size=12),
        axis.title.y=element_text(color="black", 
                                  size=16),
        panel.background = element_blank(), 
        axis.line = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(color='black', 
                                   size=12),
        axis.title.x = element_text( color="black", 
                                     size=16)
  ) 

#k-means clustering
shoots_kmeans <- kmeans(sc_shoots,centers=3)
shoots_kmeans$cluster
sc_shoots$name=rownames(sc_shoots)

#k-means visualization
autoplot(shoots_kmeans, 
         data=sc_shoots, label = TRUE,frame = TRUE,
         frame.type = 'norm',label.size = 6,label.col='black',alpha=0)+
  scale_colour_manual(values = c('black','black','black'))+
  scale_fill_manual(values = c("#A9A9A9","#A60B0B",'red'))+
  theme_bw()+
  theme(axis.text.y = element_text(color='black', 
                                   size=24),
        axis.title.y=element_text(color="black", 
                                  size=28),
        panel.background = element_blank(), 
        axis.line = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(color='black', 
                                   size=24),
        axis.title.x = element_text(color="black", 
                                    size=28),
        legend.position = 'none'
  ) + 
  scale_x_continuous(expand = c(.1, .1)) 



## New heatmap

rownames(redone_shoots_normalized) <- redone_shoots_normalized$spot_num
row_dend = as.dendrogram(hclust(dist((redone_shoots_normalized[,2:8]))))
cluster_dend <- hclust(dist((redone_shoots_normalized[,2:8])))
plot(cluster_dend)
clusterCut <- cutree(cluster_dend, k = 5)

names_frame <- data.frame(keyName=names(clusterCut), value=clusterCut, row.names=NULL) %>% arrange(value)
colnames(names_frame) <- c('Spot', 'Cluster')
names_frame$Control <- 0
names_frame$Anoxia <- 0
names_frame$Reaeration <- 0
names_frame$Annotation <- 0
names_frame$Accession <- 0

for (ind in 1:nrow(names_frame)){
  normalized_row <- c()
  normalized_row <- redone_shoots_normalized[redone_shoots_normalized$spot_num==names_frame[ind, 1], -1]
  names_frame[ind, 3] <- mean(as.numeric(normalized_row[4:5]))
  names_frame[ind, 4] <- mean(as.numeric(normalized_row[1:3]))
  names_frame[ind, 5] <- mean(as.numeric(normalized_row[6:7]))
  names_frame[ind, 6] <- normalized_row[9]
  names_frame[ind, 7] <- normalized_row[10]
}

write.xlsx(names_frame,'h_clusters_shoot_all.xlsx')


col_fun = colorRamp2(c(min(redone_shoots_normalized[,2:8]), 
                       max(redone_shoots_normalized[,2:8])-min(redone_shoots_normalized[,2:8])+4,
                       max(redone_shoots_normalized[,2:8])), c("#746A20", "#F8F5EF", "#31860E"))
col_fun(seq(min(redone_shoots_normalized[,2:8]),  max(redone_shoots_normalized[,2:8])))


Heatmap(redone_shoots_normalized[,2:8], name = "mat",
        heatmap_width  = unit(19, "cm"), heatmap_height = unit(27, "cm"),
        cluster_rows = row_dend, 
        row_names_gp = gpar(fontsize = 9),
        column_names_gp = gpar(fontsize = 18),
        col = col_fun)



##Mean_acc heatmap

row.names(annotation_shoots_annot_with_accessions_means) <- annotation_shoots_annot_with_accessions_means$accessions

row_dend = as.dendrogram(hclust(dist((annotation_shoots_annot_with_accessions_means[2:4]))))

col_fun = colorRamp2(c(min(annotation_shoots_annot_with_accessions_means[2:4])+1, 
                       max(annotation_shoots_annot_with_accessions_means[2:4])-min(annotation_shoots_annot_with_accessions_means[2:4])+7,
                       max(annotation_shoots_annot_with_accessions_means[2:4])), c( "#F8F5EF","#746A20", "#31860E"))
col_fun(seq(min(annotation_shoots_annot_with_accessions_means[2:4]),  max(annotation_shoots_annot_with_accessions_means[2:4])))



Heatmap(annotation_shoots_annot_with_accessions_means[2:4], name = "mat",
        heatmap_width  = unit(15, "cm"), heatmap_height = unit(15, "cm"),
        cluster_rows = row_dend, 
        row_names_gp = gpar(fontsize = 12),
        column_names_gp = gpar(fontsize = 18),
        col = col_fun,
        cluster_columns = FALSE)


##Signif heatmap

rownames(signif_pval_shoots) <- signif_pval_shoots$spot_num
row_dend = as.dendrogram(hclust(dist((signif_pval_shoots[,2:8]))))
cluster_dend <- hclust(dist((signif_pval_shoots[,2:8])))
plot(cluster_dend)
clusterCut <- cutree(cluster_dend, k = 4)

names_frame <- data.frame(keyName=names(clusterCut), value=clusterCut, row.names=NULL) %>% arrange(value)
colnames(names_frame) <- c('Spot', 'Cluster')
names_frame$Control <- 0
names_frame$Anoxia <- 0
names_frame$Reaeration <- 0
names_frame$Annotation <- 0
names_frame$Accession <- 0

for (ind in 1:nrow(names_frame)){
  normalized_row <- c()
  normalized_row <- signif_pval_shoots[signif_pval_shoots$spot_num==names_frame[ind, 1], -1]
  names_frame[ind, 3] <- mean(as.numeric(normalized_row[4:5]))
  names_frame[ind, 4] <- mean(as.numeric(normalized_row[1:3]))
  names_frame[ind, 5] <- mean(as.numeric(normalized_row[6:7]))
  names_frame[ind, 6] <- normalized_row[9]
  names_frame[ind, 7] <- normalized_row[10]
}

write.xlsx(names_frame,'h_clusters_shoot_signif.xlsx')



col_fun = colorRamp2(c(min(signif_pval_shoots[,2:8]), 
                       max(signif_pval_shoots[,2:8])-min(signif_pval_shoots[,2:8])+4,
                       max(signif_pval_shoots[,2:8])), c("#746A20", "#F8F5EF", "#31860E"))

col_fun(seq(min(signif_pval_shoots[,2:8]),  max(signif_pval_shoots[,2:8])))

Heatmap(signif_pval_shoots[,2:8], name = "mat",
        heatmap_width  = unit(15, "cm"), heatmap_height = unit(15, "cm"),
        cluster_rows = row_dend, 
        row_names_gp = gpar(fontsize = 12),
        column_names_gp = gpar(fontsize = 18),
        col = col_fun)

##Annot heatmap
annotation_shoots_annot <- redone_shoots_normalized[redone_shoots_normalized$SSP %in% annotation_shoots_df_spots,]
rownames(annotation_shoots_annot) <- annotation_shoots_annot$spot_num
row_dend = as.dendrogram(hclust(dist((annotation_shoots_annot[,2:8]))))
cluster_dend <- hclust(dist((annotation_shoots_annot[,2:8])))
plot(cluster_dend)
clusterCut <- cutree(cluster_dend, h=14)

names_frame <- data.frame(keyName=names(clusterCut), value=clusterCut, row.names=NULL) %>% arrange(value)
colnames(names_frame) <- c('Spot', 'Cluster')
names_frame$Control <- 0
names_frame$Anoxia <- 0
names_frame$Reaeration <- 0
names_frame$Annotation <- 0
names_frame$Accession <- 0

for (ind in 1:nrow(names_frame)){
  normalized_row <- c()
  normalized_row <- annotation_shoots_annot[annotation_shoots_annot$spot_num==names_frame[ind, 1], -1]
  names_frame[ind, 3] <- mean(as.numeric(normalized_row[4:5]))
  names_frame[ind, 4] <- mean(as.numeric(normalized_row[1:3]))
  names_frame[ind, 5] <- mean(as.numeric(normalized_row[6:7]))
  names_frame[ind, 6] <- normalized_row[9]
  names_frame[ind, 7] <- normalized_row[10]
}

write.xlsx(names_frame,'h_clusters_shoot_annot.xlsx')



col_fun = colorRamp2(c(min(annotation_shoots_annot[,2:8]), 
                       max(annotation_shoots_annot[,2:8])-min(annotation_shoots_annot[,2:8])+4,
                       max(annotation_shoots_annot[,2:8])), c("#746A20", "#F8F5EF", "#31860E"))

col_fun(seq(min(annotation_shoots_annot[,2:8]),  max(annotation_shoots_annot[,2:8])))

Heatmap(annotation_shoots_annot[,2:8], name = "mat",
        heatmap_width  = unit(18, "cm"), heatmap_height = unit(18, "cm"),
        cluster_rows = row_dend, 
        row_names_gp = gpar(fontsize = 12),
        column_names_gp = gpar(fontsize = 18),
        col = col_fun)

##Annot signif heatmap
rownames(annotated_sigif_shoots) <- annotated_sigif_shoots$spot_num
row_dend = as.dendrogram(hclust(dist((annotated_sigif_shoots[,2:8]))))
cluster_dend <- hclust(dist((annotated_sigif_shoots[,2:8])))
plot(cluster_dend)
clusterCut <- cutree(cluster_dend, h=10)

names_frame <- data.frame(keyName=names(clusterCut), value=clusterCut, row.names=NULL) %>% arrange(value)
colnames(names_frame) <- c('Spot', 'Cluster')
names_frame$Control <- 0
names_frame$Anoxia <- 0
names_frame$Reaeration <- 0
names_frame$Annotation <- 0
names_frame$Accession <- 0

for (ind in 1:nrow(names_frame)){
  normalized_row <- c()
  normalized_row <- annotated_sigif_shoots[annotated_sigif_shoots$spot_num==names_frame[ind, 1], -1]
  names_frame[ind, 3] <- mean(as.numeric(normalized_row[4:5]))
  names_frame[ind, 4] <- mean(as.numeric(normalized_row[1:3]))
  names_frame[ind, 5] <- mean(as.numeric(normalized_row[6:7]))
  names_frame[ind, 6] <- normalized_row[9]
  names_frame[ind, 7] <- normalized_row[10]
}

write.xlsx(names_frame,'h_clusters_shoot_signif_annot.xlsx')


col_fun = colorRamp2(c(min(annotated_sigif_shoots[,2:8]), 
                       max(annotated_sigif_shoots[,2:8])-min(annotated_sigif_shoots[,2:8])+5,
                       max(annotated_sigif_shoots[,2:8])), c("#746A20", "#F8F5EF", "#31860E"))

col_fun(seq(min(annotated_sigif_shoots[,2:8]),  max(annotated_sigif_shoots[,2:8])))

Heatmap(annotated_sigif_shoots[,2:8], name = "mat",
        heatmap_width  = unit(18, "cm"), heatmap_height = unit(18, "cm"),
        cluster_rows = row_dend, 
        row_names_gp = gpar(fontsize = 12),
        column_names_gp = gpar(fontsize = 18),
        col = col_fun)

#K-means All spots 
sc_shoots<-data.frame(redone_shoots_normalized[2:8])


#elbow method visualization
wss <- (nrow(sc_shoots)-1)*sum(apply(sc_shoots,2,var))
for (i in 1:4) wss[i] <- sum(kmeans(sc_shoots,
                                    centers=i)$withinss,nstart=25,iter.max=1000)

wss_gf=data.frame(clust=c(1:4), wss=wss)

ggplot(wss_gf, aes(clust,wss))+
  geom_line()+theme_bw()+ xlab('Number of clusters') + ylab('WSS')+
  geom_vline(xintercept =2, linetype="dashed")+
  theme(axis.text.y = element_text(color='black', 
                                   size=12),
        axis.title.y=element_text(color="black", 
                                  size=16),
        panel.background = element_blank(), 
        axis.line = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(color='black', 
                                   size=12),
        axis.title.x = element_text( color="black", 
                                     size=16)
  ) 

#k-means clustering
shoots_kmeans <- kmeans(sc_shoots,centers=2)
shoots_kmeans$cluster
sc_shoots$name=rownames(sc_shoots)

#k-means visualization
autoplot(shoots_kmeans, 
         data=sc_shoots, label = TRUE,frame = TRUE,
         frame.type = 'norm',label.size = 6,label.col='black',alpha=0)+
  scale_colour_manual(values = c('black','black','black'))+
  scale_fill_manual(values = c("#A9A9A9","#A60B0B",'red'))+
  theme_bw()+
  theme(axis.text.y = element_text(color='black', 
                                   size=24),
        axis.title.y=element_text(color="black", 
                                  size=28),
        panel.background = element_blank(), 
        axis.line = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(color='black', 
                                   size=24),
        axis.title.x = element_text(color="black", 
                                    size=28),
        legend.position = 'none'
  ) + 
  scale_x_continuous(expand = c(.1, .1)) 



clusterCut <-shoots_kmeans$cluster
names_frame <- data.frame(keyName=names(clusterCut), value=clusterCut, row.names=NULL) %>% arrange(value)
colnames(names_frame) <- c('Spot', 'Cluster')
names_frame$Control <- 0
names_frame$Anoxia <- 0
names_frame$Reaeration <- 0
names_frame$Annotation <- 0
names_frame$Accession <- 0

for (ind in 1:nrow(names_frame)){
  normalized_row <- c()
  normalized_row <- redone_shoots_normalized[redone_shoots_normalized$spot_num==names_frame[ind, 1], -1]
  names_frame[ind, 3] <- mean(as.numeric(normalized_row[4:5]))
  names_frame[ind, 4] <- mean(as.numeric(normalized_row[1:3]))
  names_frame[ind, 5] <- mean(as.numeric(normalized_row[6:7]))
  names_frame[ind, 6] <- normalized_row[9]
  names_frame[ind, 7] <- normalized_row[10]
}

write.xlsx(names_frame,'kmeans_clusters_shoot_all.xlsx')


#K-means signif spots
sc_shoots<-data.frame(signif_pval_shoots[2:8])


#elbow method visualization
wss <- (nrow(sc_shoots)-1)*sum(apply(sc_shoots,2,var))
for (i in 1:4) wss[i] <- sum(kmeans(sc_shoots,
                                    centers=i)$withinss,nstart=25,iter.max=1000)

wss_gf=data.frame(clust=c(1:4), wss=wss)

ggplot(wss_gf, aes(clust,wss))+
  geom_line()+theme_bw()+ xlab('Number of clusters') + ylab('WSS')+
  geom_vline(xintercept =2, linetype="dashed")+
  theme(axis.text.y = element_text(color='black', 
                                   size=12),
        axis.title.y=element_text(color="black", 
                                  size=16),
        panel.background = element_blank(), 
        axis.line = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(color='black', 
                                   size=12),
        axis.title.x = element_text( color="black", 
                                     size=16)
  ) 

#k-means clustering
shoots_kmeans <- kmeans(sc_shoots,centers=2)
shoots_kmeans$cluster
sc_shoots$name=rownames(sc_shoots)

#k-means visualization
autoplot(shoots_kmeans, 
         data=sc_shoots, label = TRUE,frame = TRUE,
         frame.type = 'norm',label.size = 6,label.col='black',alpha=0)+
  scale_colour_manual(values = c('black','black','black'))+
  scale_fill_manual(values = c("#A9A9A9","#A60B0B",'red'))+
  theme_bw()+
  theme(axis.text.y = element_text(color='black', 
                                   size=24),
        axis.title.y=element_text(color="black", 
                                  size=28),
        panel.background = element_blank(), 
        axis.line = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(color='black', 
                                   size=24),
        axis.title.x = element_text(color="black", 
                                    size=28),
        legend.position = 'none'
  ) + 
  scale_x_continuous(expand = c(.1, .1)) 



clusterCut <- shoots_kmeans$cluster
names_frame <- data.frame(keyName=names(clusterCut), value=clusterCut, row.names=NULL) %>% arrange(value)
colnames(names_frame) <- c('Spot', 'Cluster')
names_frame$Control <- 0
names_frame$Anoxia <- 0
names_frame$Reaeration <- 0
names_frame$Annotation <- 0
names_frame$Accession <- 0

for (ind in 1:nrow(names_frame)){
  normalized_row <- c()
  normalized_row <- signif_pval_shoots[signif_pval_shoots$spot_num==names_frame[ind, 1], -1]
  names_frame[ind, 3] <- mean(as.numeric(normalized_row[4:5]))
  names_frame[ind, 4] <- mean(as.numeric(normalized_row[1:3]))
  names_frame[ind, 5] <- mean(as.numeric(normalized_row[6:7]))
  names_frame[ind, 6] <- normalized_row[9]
  names_frame[ind, 7] <- normalized_row[10]
}

write.xlsx(names_frame,'kmeans_clusters_shoot_signif.xlsx')


#K-means annot spots
sc_shoots<-data.frame(annotation_shoots_annot[2:8])


#elbow method visualization
wss <- (nrow(sc_shoots)-1)*sum(apply(sc_shoots,2,var))
for (i in 1:4) wss[i] <- sum(kmeans(sc_shoots,
                                    centers=i)$withinss,nstart=25,iter.max=1000)

wss_gf=data.frame(clust=c(1:4), wss=wss)

ggplot(wss_gf, aes(clust,wss))+
  geom_line()+theme_bw()+ xlab('Number of clusters') + ylab('WSS')+
  geom_vline(xintercept =2, linetype="dashed")+
  theme(axis.text.y = element_text(color='black', 
                                   size=12),
        axis.title.y=element_text(color="black", 
                                  size=16),
        panel.background = element_blank(), 
        axis.line = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(color='black', 
                                   size=12),
        axis.title.x = element_text( color="black", 
                                     size=16)
  ) 

#k-means clustering
shoots_kmeans <- kmeans(sc_shoots,centers=2)
shoots_kmeans$cluster
sc_shoots$name=rownames(sc_shoots)

#k-means visualization
autoplot(shoots_kmeans, 
         data=sc_shoots, label = TRUE,frame = TRUE,
         frame.type = 'norm',label.size = 6,label.col='black',alpha=0)+
  scale_colour_manual(values = c('black','black','black'))+
  scale_fill_manual(values = c("#A9A9A9","#A60B0B",'red'))+
  theme_bw()+
  theme(axis.text.y = element_text(color='black', 
                                   size=24),
        axis.title.y=element_text(color="black", 
                                  size=28),
        panel.background = element_blank(), 
        axis.line = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(color='black', 
                                   size=24),
        axis.title.x = element_text(color="black", 
                                    size=28),
        legend.position = 'none'
  ) + 
  scale_x_continuous(expand = c(.1, .1)) 


clusterCut <- shoots_kmeans$cluster
names_frame <- data.frame(keyName=names(clusterCut), value=clusterCut, row.names=NULL) %>% arrange(value)
colnames(names_frame) <- c('Spot', 'Cluster')
names_frame$Control <- 0
names_frame$Anoxia <- 0
names_frame$Reaeration <- 0
names_frame$Annotation <- 0
names_frame$Accession <- 0

for (ind in 1:nrow(names_frame)){
  normalized_row <- c()
  normalized_row <- annotation_shoots_annot[annotation_shoots_annot$spot_num==names_frame[ind, 1], -1]
  names_frame[ind, 3] <- mean(as.numeric(normalized_row[4:5]))
  names_frame[ind, 4] <- mean(as.numeric(normalized_row[1:3]))
  names_frame[ind, 5] <- mean(as.numeric(normalized_row[6:7]))
  names_frame[ind, 6] <- normalized_row[9]
  names_frame[ind, 7] <- normalized_row[10]
}

write.xlsx(names_frame,'kmeans_clusters_shoot_annot.xlsx')


#K-means signif annot spots
sc_shoots<-data.frame(annotated_sigif_shoots[2:8])


#elbow method visualization
wss <- (nrow(sc_shoots)-1)*sum(apply(sc_shoots,2,var))
for (i in 1:4) wss[i] <- sum(kmeans(sc_shoots,
                                    centers=i)$withinss,nstart=25,iter.max=1000)

wss_gf=data.frame(clust=c(1:4), wss=wss)

ggplot(wss_gf, aes(clust,wss))+
  geom_line()+theme_bw()+ xlab('Number of clusters') + ylab('WSS')+
  geom_vline(xintercept =2, linetype="dashed")+
  theme(axis.text.y = element_text(color='black', 
                                   size=12),
        axis.title.y=element_text(color="black", 
                                  size=16),
        panel.background = element_blank(), 
        axis.line = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(color='black', 
                                   size=12),
        axis.title.x = element_text( color="black", 
                                     size=16)
  ) 

#k-means clustering
shoots_kmeans <- kmeans(sc_shoots,centers=2)
shoots_kmeans$cluster
sc_shoots$name=rownames(sc_shoots)

#k-means visualization
autoplot(shoots_kmeans, 
         data=sc_shoots, label = TRUE,frame = TRUE,
         frame.type = 'norm',label.size = 6,label.col='black',alpha=0)+
  scale_colour_manual(values = c('black','black','black'))+
  scale_fill_manual(values = c("#A9A9A9","#A60B0B",'red'))+
  theme_bw()+
  theme(axis.text.y = element_text(color='black', 
                                   size=24),
        axis.title.y=element_text(color="black", 
                                  size=28),
        panel.background = element_blank(), 
        axis.line = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(color='black', 
                                   size=24),
        axis.title.x = element_text(color="black", 
                                    size=28),
        legend.position = 'none'
  ) + 
  scale_x_continuous(expand = c(.1, .1)) 


clusterCut <- shoots_kmeans$cluster
names_frame <- data.frame(keyName=names(clusterCut), value=clusterCut, row.names=NULL) %>% arrange(value)
colnames(names_frame) <- c('Spot', 'Cluster')
names_frame$Control <- 0
names_frame$Anoxia <- 0
names_frame$Reaeration <- 0
names_frame$Annotation <- 0
names_frame$Accession <- 0

for (ind in 1:nrow(names_frame)){
  normalized_row <- c()
  normalized_row <- annotated_sigif_shoots[annotated_sigif_shoots$spot_num==names_frame[ind, 1], -1]
  names_frame[ind, 3] <- mean(as.numeric(normalized_row[4:5]))
  names_frame[ind, 4] <- mean(as.numeric(normalized_row[1:3]))
  names_frame[ind, 5] <- mean(as.numeric(normalized_row[6:7]))
  names_frame[ind, 6] <- normalized_row[9]
  names_frame[ind, 7] <- normalized_row[10]
}

write.xlsx(names_frame,'kmeans_clusters_shoot_signif_annot.xlsx')



###Dots types

root_sumary = data.frame(Spots=c("Annotated",'Signifficant',
                                 'Annotated and Signifficant', 'Other'), 
                         num=c(31,17,10,35))

shoot_sumary = data.frame(Spots=c("Annotated",'Signifficant',
                                  'Annotated and Signifficant', 'Other'), 
                          num=c(26,10,15,23))

ggplot(root_sumary, aes(x='',y=num,fill=Spots)) +
  geom_bar(stat="identity", width=1, alpha=0.55, col='black')+
  coord_polar(theta = "y")+
  theme_bw()+
  geom_text(aes(label = num),
            position = position_stack(vjust = 0.5), size=10)+
  geom_text(aes(label = num),
            position = position_stack(vjust = 0.5), size=10)+
  geom_text(aes(label = num),
            position = position_stack(vjust = 0.5), size=10)+
  scale_fill_manual(values = c("#a60b0b", '#2980b9', '#9A9EAB','#bdbd00'))+
  scale_size(guide = 'none')+
  theme( axis.text.x = element_blank(),
         axis.title.x=element_blank(),
         panel.background = element_blank(), 
         axis.line = element_line(colour = "black"),
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         axis.text.y = element_blank(),
         axis.title.y = element_blank(),
         legend.title=element_text(face="bold",size=14), 
         legend.text=element_text(size=12),
         axis.ticks = element_blank())+
  guides(fill= guide_legend(title="Spot type"))

### ALL t.tests

anoxia_roots = redone_roots_log_reshaped[redone_roots_log_reshaped$factor_col=='Anoxia', 2] #19.8538 - anoxia
control_roots = redone_roots_log_reshaped[redone_roots_log_reshaped$factor_col=='Control', 2] #20.76237 - control
reaeration_roots <- redone_roots_log_reshaped[redone_roots_log_reshaped$factor_col=='Re-aeration', 2] #21.6658 - re-aeration
mean_OD_df_roots <- data.frame(Density = c(anoxia_roots, control_roots, reaeration_roots), 
                         Condition=c(rep('anoxia_shoots',279),rep('control_shoots',279),rep('reaeration_shoots',279)))

anoxia_shoots = redone_shoots_log_reshaped[redone_shoots_log_reshaped$factor_col=='Anoxia', 2]  #19.11415 - anoxia
control_shoots = redone_shoots_log_reshaped[redone_shoots_log_reshaped$factor_col=='Control', 2] #20.38731 - control
reaeration_shoots <- redone_shoots_log_reshaped[redone_shoots_log_reshaped$factor_col=='Re-aeration', 2] #19.09779 - re-aeration
mean_OD_df_shoots <- data.frame(Density = c(anoxia_shoots, control_shoots, reaeration_shoots), 
                         Condition=c(rep('anoxia_shoots',222),rep('control_shoots',148),rep('reaeration_shoots',148)))



t_test_df_all=data.frame(comparision=c("Anoxia_roots_vs_Control_roots","Anoxia_roots_vs_Re-aeration_roots", 'Control_roots_vs_Re-aeration_roots',
                                       "Anoxia_shoots_vs_Control_shoots","Anoxia_shoots_vs_Re-aeration_shoots", 'Control_shots_vs_Re-aeration_shoots',
                                       "Anoxia_roots_vs_Anoxia_shoots","Control_roots_vs_Control_shoots","Re-aeration_roots_vs_Re-aeration_shoots",
                                       "Anoxia_roots_vs_Control_shoots","Anoxia_roots_vs_Re-aeration_shoots","Control_roots_vs_Re-aeration_shoots",
                                       "Anoxia_shoots_vs_Control_roots","Anoxia_shoots_vs_Re-aeration_roots","Control_shoots_vs_Re-aeration_roots"
                                       ), p.value=1)

t_test_df_all[1, 2] <- t.test(anoxia_roots, control_roots)$p.value
t_test_df_all[2, 2] <- t.test(anoxia_roots, reaeration_roots)$p.value
t_test_df_all[3, 2] <- t.test(control_roots, reaeration_roots)$p.value

t_test_df_all[4, 2] <- t.test(anoxia_shoots, control_shoots)$p.value
t_test_df_all[5, 2] <- t.test(anoxia_shoots, reaeration_shoots)$p.value
t_test_df_all[6, 2] <- t.test(control_shoots, reaeration_shoots)$p.value

t_test_df_all[7, 2] <- t.test(anoxia_roots, anoxia_shoots)$p.value
t_test_df_all[8, 2] <- t.test(control_roots, control_shoots)$p.value
t_test_df_all[9, 2] <- t.test(reaeration_roots, reaeration_shoots)$p.value

t_test_df_all[10, 2] <- t.test(anoxia_roots, control_shoots)$p.value
t_test_df_all[11, 2] <- t.test(anoxia_roots, reaeration_shoots)$p.value
t_test_df_all[12, 2] <- t.test(control_roots, reaeration_shoots)$p.value

t_test_df_all[13, 2] <- t.test(anoxia_shoots, control_roots)$p.value
t_test_df_all[14, 2] <- t.test(anoxia_shoots, reaeration_shoots)$p.value
t_test_df_all[15, 2] <- t.test(control_shoots, reaeration_roots)$p.value


t_test_df_all$p.adj <- p.adjust(t_test_df_all$p.value, method = "fdr")
write.xlsx(t_test_df_all,'t_test_shoots_vs_roots.xlsx')



t.test(mean_OD_df_roots$Density, mean_OD_df_shoots$Density) #p-value < 2.2e-16

OD_organs_df <- data.frame(Density=c(mean_OD_df_roots$Density, mean_OD_df_shoots$Density),
                           factor_col=c(rep('roots',837), rep('shoots',518)))

#Mean OD violins for organs
ggplot(OD_organs_df, aes(y=Density, x=factor_col, fill=factor_col))+
  geom_violin(alpha=0.9, color='black')+theme_bw()+
  xlab('Condition') + ylab('Optical density')+
  theme(axis.text.y = element_text(color='black', 
                                   size=16),
        axis.title.y=element_text(color="black", 
                                  size=18),
        panel.background = element_blank(), 
        axis.line = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(color='black', 
                                   size=16),
        axis.title.x = element_text( color="black", 
                                     size=18),
        legend.title=element_text(size=18), 
        legend.text=element_text(size=16)
  )+ 
  guides(fill = guide_legend(override.aes = list(size = 1),title="Condition"))





#______________________ signifficant shoots no rubicso


#redo functions
k_p_val <- function(x) { vals_df <- x %>% gather(condition, measurement, A1:R2, factor_key=TRUE)
vals_df$factor_col <- 'factor'
vals_df[grepl("A", vals_df$condition), 3] <- 'Anoxia'
vals_df[grepl("C", vals_df$condition), 3] <- 'Control'
vals_df[grepl("R", vals_df$condition), 3] <- 'Re-aeration'
tryCatch(kruskal.test(measurement ~ factor_col, data = vals_df)$p.value,
         error = function(e) NA)
}

#mixed linear model function
lme_p_val <- function(x) { vals_df <- x %>% gather(condition, measurement, A1:R2, factor_key=TRUE)
vals_df$factor_col <- 'factor'
vals_df[grepl("A", vals_df$condition), 3] <- 'Anoxia'
vals_df[grepl("C", vals_df$condition), 3] <- 'Control'
vals_df[grepl("R", vals_df$condition), 3] <- 'Re-aeration'
fit.lme <- lme(measurement ~ factor_col, random=~1|condition, data=vals_df)
ret_vec <- c(as.data.frame(anova(fit.lme))$`p-value`[2])
fit_sum <- as.data.frame(summary(fit.lme)$tTable)
ret_vec <- c(ret_vec, fit_sum$`p-value`[c(2,3)])
tryCatch(ret_vec,
         error = function(e) NA)
}

pvals_shoot_kraskell <- c()
pvals_shoot_lme <- data.frame(all_pval=c(1), cont_pval = c(1), re_pval = c(1))


#find p-vals for tests
for (i in 1:nrow(redone_shoots_normalized_no_rubisco)){
  pvals_shoot_kraskell <- c(pvals_shoot_kraskell, k_p_val(redone_shoots_normalized_no_rubisco[i,2:8]))
  pvals_shoot_lme <- rbind(pvals_shoot_lme, lme_p_val(redone_shoots_normalized_no_rubisco[i,2:8]))
}
pvals_shoot_lme <- pvals_shoot_lme[-1,]

#adjust p-values
pval_df_shoots=data_frame(spot=redone_shoots_normalized_no_rubisco[,1], kraskell_raw=pvals_shoot_kraskell, 
                          kraskell_adj = p.adjust(pvals_shoot_kraskell, method = "fdr"),
                          lme_all_raw = pvals_shoot_lme[,1],
                          lme_all_adj = p.adjust(pvals_shoot_lme[,1], method = "fdr"),
                          lme_cont_raw = pvals_shoot_lme[,2],
                          lme_cont_adj = p.adjust(pvals_shoot_lme[,2], method = "fdr"),
                          lme_re_raw = pvals_shoot_lme[,3],
                          lme_re_adj = p.adjust(pvals_shoot_lme[,3], method = "fdr"))

#adjustment results
#FDR
sum(pval_df_shoots$kraskell_raw <= 0.05, na.rm = TRUE) #0 proteins with significant 
sum(pval_df_shoots$lme_all_raw <= 0.05, na.rm = TRUE) #17 proteins with significant
sum(pval_df_shoots$lme_cont_raw <= 0.05, na.rm = TRUE) #12 proteins with significant 
sum(pval_df_shoots$lme_re_raw <= 0.05, na.rm = TRUE) #15 proteins with significant
sum(pval_df_shoots$kraskell_adj <= 0.05, na.rm = TRUE) #0 proteins with significant 
sum(pval_df_shoots$lme_all_adj <= 0.05, na.rm = TRUE) #5 proteins with significant
sum(pval_df_shoots$lme_cont_adj <= 0.05, na.rm = TRUE) #2 proteins with significant 
sum(pval_df_shoots$lme_re_adj <= 0.05, na.rm = TRUE) #3 proteins with significant

#holm
sum(pval_df_shoots$kraskell_adj <= 0.05, na.rm = TRUE) #0 proteins with significant 
sum(pval_df_shoots$lme_all_adj <= 0.05, na.rm = TRUE) #0 proteins with significant
sum(pval_df_shoots$lme_cont_adj <= 0.05, na.rm = TRUE) #0 proteins with significant 
sum(pval_df_shoots$lme_re_adj <= 0.05, na.rm = TRUE) #0 proteins with significant

#BH
sum(pval_df_shoots$kraskell_adj <= 0.05, na.rm = TRUE) #0 proteins with significant 
sum(pval_df_shoots$lme_all_adj <= 0.05, na.rm = TRUE) #0 proteins with significant
sum(pval_df_shoots$lme_cont_adj <= 0.05, na.rm = TRUE) #0 proteins with significant 
sum(pval_df_shoots$lme_re_adj <= 0.05, na.rm = TRUE) #0 proteins with significant

#Modded_linear_model
shoot.fac <- data.frame(Condition=c('anoxia','anoxia','anoxia','control','control', 're-aeration', 're-aeration'))
rownames(shoot.fac) <- colnames( redone_shoots_normalized_no_rubisco[,-1])
rownames(redone_shoots_normalized_no_rubisco) <- gsub(" ", "_", redone_shoots_normalized_no_rubisco[,1])

X <- model.matrix(~ Condition, data = shoot.fac)
fit <- lmFit(redone_shoots_normalized_no_rubisco[,2:8], design = X, method = "robust", maxit = 10000)
efit <- eBayes(fit)

#modelled_p-values adjustment
num_spots <- nrow(redone_shoots_normalized_no_rubisco)
pval_df_shoots$eBayes_cont_raw <- topTable(efit, coef = 2, number = num_spots,
                                           sort.by = "none")$P.Value
pval_df_shoots$eBayes_cont_adj <- topTable(efit, coef = 2, number = num_spots,
                                           sort.by = "none")$adj.P.Val
pval_df_shoots$eBayes_re_raw <- topTable(efit, coef = 3, number = num_spots,
                                         sort.by = "none")$P.Value
pval_df_shoots$eBayes_re_adj <- topTable(efit, coef = 3, number = num_spots,
                                         sort.by = "none")$adj.P.Val

#the number of signifficant proteins
sum(pval_df_shoots$eBayes_cont_adj<=0.05) #14
sum(pval_df_shoots$eBayes_cont_raw<=0.05) #22
sum(pval_df_shoots$eBayes_re_adj<=0.05) #21
sum(pval_df_shoots$eBayes_re_raw<=0.05) #26

write.xlsx(pval_df,'dots_pvals_shoots.xlsx')

#extract significant spots
signif_pval_shoots <- redone_shoots_normalized_no_rubisco[pval_df_shoots$eBayes_cont_adj<=0.05 | pval_df_shoots$eBayes_re_adj<=0.05,]
rownames(signif_pval_shoots) <- redone_shoots_normalized_no_rubisco[pval_df_shoots$eBayes_cont_adj<=0.05 | pval_df_shoots$eBayes_re_adj<=0.05,9][-c(26,27)]

write.xlsx(signif_pval_shoots,'rice_shoots_signif_spots_annot.xlsx')

#K-means Clusters for significant points
sc_shoots<-data.frame(t(signif_pval_shoots[2:8]))


#elbow method visualization
wss <- (nrow(sc_shoots)-1)*sum(apply(sc_shoots,2,var))
for (i in 1:4) wss[i] <- sum(kmeans(sc_shoots,
                                    centers=i)$withinss,nstart=25,iter.max=1000)

wss_gf=data.frame(clust=c(1:4), wss=wss)

ggplot(wss_gf, aes(clust,wss))+
  geom_line()+theme_bw()+ xlab('Number of clusters') + ylab('WSS')+
  geom_vline(xintercept =3, linetype="dashed")+
  theme(axis.text.y = element_text(color='black', 
                                   size=12),
        axis.title.y=element_text(color="black", 
                                  size=16),
        panel.background = element_blank(), 
        axis.line = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(color='black', 
                                   size=12),
        axis.title.x = element_text( color="black", 
                                     size=16)
  ) 

#k-means clustering
shoots_kmeans <- kmeans(sc_shoots,centers=3)
shoots_kmeans$cluster
sc_shoots$name=rownames(sc_shoots)

#k-means visualization
autoplot(shoots_kmeans, 
         data=sc_shoots, label = TRUE,frame = TRUE,
         frame.type = 'norm',label.size = 6,label.col='black',alpha=0)+
  scale_colour_manual(values = c('black','black','black'))+
  scale_fill_manual(values = c("#A9A9A9","#A60B0B",'red'))+
  theme_bw()+
  theme(axis.text.y = element_text(color='black', 
                                   size=24),
        axis.title.y=element_text(color="black", 
                                  size=28),
        panel.background = element_blank(), 
        axis.line = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(color='black', 
                                   size=24),
        axis.title.x = element_text(color="black", 
                                    size=28),
        legend.position = 'none'
  ) + 
  scale_x_continuous(expand = c(.1, .1)) 




rownames(signif_pval_shoots) <- signif_pval_shoots$spot_num
row_dend = as.dendrogram(hclust(dist((signif_pval_shoots[,2:8]))))
cluster_dend <- hclust(dist((signif_pval_shoots[,2:8])))
plot(cluster_dend)
clusterCut <- cutree(cluster_dend, k = 4)

names_frame <- data.frame(keyName=names(clusterCut), value=clusterCut, row.names=NULL) %>% arrange(value)
colnames(names_frame) <- c('Spot', 'Cluster')
names_frame$Control <- 0
names_frame$Anoxia <- 0
names_frame$Reaeration <- 0
names_frame$Annotation <- 0
names_frame$Accession <- 0

for (ind in 1:nrow(names_frame)){
  normalized_row <- c()
  normalized_row <- signif_pval_shoots[signif_pval_shoots$spot_num==names_frame[ind, 1], -1]
  names_frame[ind, 3] <- mean(as.numeric(normalized_row[4:5]))
  names_frame[ind, 4] <- mean(as.numeric(normalized_row[1:3]))
  names_frame[ind, 5] <- mean(as.numeric(normalized_row[6:7]))
  names_frame[ind, 6] <- normalized_row[9]
  names_frame[ind, 7] <- normalized_row[10]
}



col_fun = colorRamp2(c(min(signif_pval_shoots[,2:8]), 
                       max(signif_pval_shoots[,2:8])-min(signif_pval_shoots[,2:8])+4,
                       max(signif_pval_shoots[,2:8])), c("#746A20", "#F8F5EF", "#31860E"))

col_fun(seq(min(signif_pval_shoots[,2:8]),  max(signif_pval_shoots[,2:8])))

Heatmap(signif_pval_shoots[,2:8], name = "mat",
        heatmap_width  = unit(15, "cm"), heatmap_height = unit(15, "cm"),
        cluster_rows = row_dend, 
        row_names_gp = gpar(fontsize = 12),
        column_names_gp = gpar(fontsize = 18),
        col = col_fun)
