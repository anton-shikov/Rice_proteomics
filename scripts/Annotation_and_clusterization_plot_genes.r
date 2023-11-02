library(dplyr)
library(ggplot2)
library(BiocManager)
library(ALL)
library(hgu95av2.db)
library(topGO)
library(reshape2)
library(xlsx)
library(wesanderson)
require(grid)


setwd("../data")

manhattan_file = read.csv('For_processing_scripts/Coordinates.csv', header = T, sep = "\t") 
manhattan_new = data.frame(chrom=c(), coord=c(), intens=c(), organ=c(), type=c())

for (ind in 1:nrow(manhattan_file)){
  anox <- as.data.frame(c(manhattan_file[ind, 1], (manhattan_file[ind, 3]+manhattan_file[ind, 2])/2, manhattan_file[ind, c(6,5)], 'Anoxia'))
  cont <- as.data.frame(c(manhattan_file[ind, 1], (manhattan_file[ind, 3]+manhattan_file[ind, 2])/2, manhattan_file[ind, c(7,5)], 'Control'))
  rear <- as.data.frame(c(manhattan_file[ind, 1], (manhattan_file[ind, 3]+manhattan_file[ind, 2])/2, manhattan_file[ind, c(8,5)], 'Re-aeration'))
  
  colnames(anox)<-c('chrom', 'coord', 'intens', 'organ', 'type')
  colnames(cont)<-c('chrom', 'coord', 'intens', 'organ', 'type')
  colnames(rear)<-c('chrom', 'coord', 'intens', 'organ', 'type')
  manhattan_new <- rbind(manhattan_new, anox)
  manhattan_new <- rbind(manhattan_new, cont)
  manhattan_new <- rbind(manhattan_new, rear)
}


# #75522B - roots, #1F732D - shoots
colnames(manhattan_new) = c('CHR', "BP", "INT",'ORGAN', 'TYPE')

don <- manhattan_new %>% 
  group_by(CHR) %>% 
  summarise(chr_len=max(BP)) %>% 
  mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
  dplyr::select(-chr_len) %>%
  left_join(manhattan_new, ., by=c("CHR"="CHR")) %>%
  arrange(CHR, BP) %>%
  mutate( BPcum=BP+tot)
axisdf = don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

don <- don[don$ORGAN!='ref',]

don$ORGAN <- factor(don$ORGAN, levels = c('shoot','root'))
levels(don$ORGAN) <-c('Shoot','Root')
don$TYPE <- factor(don$TYPE, levels = c('Control','Anoxia','Re-aeration'))

ggplot(don, aes(x=BPcum, y=INT, shape=ORGAN)) +
  geom_point( aes(color=as.factor(TYPE)), alpha=0.9, size=2.7) +
  scale_color_manual(values = rep(c("#678667","#375DB6","#A5432B"), 22)) +
  scale_x_continuous(label = axisdf$CHR, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 0),limits = c(10, 27)) +
  coord_fixed(1 * max(don$BPcum) / (max(don$INT)+30))+
  xlab('Genomic coordinate') + ylab('Optical density')+
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size=16,angle = 45),
    axis.text.y = element_text(size=16),
    axis.title.x =  element_text(size=18),
    axis.title.y= element_text(size=18),
    axis.line = element_line(colour = "black")
  )+ guides(color=guide_legend(title="Condition"),shape=guide_legend(title="Plant organ"))

manhattan_stat <-as.data.frame( manhattan_file[manhattan_file$Strand!='ref',c(1,5)] %>%table())
manhattan_stat_filt <- manhattan_stat[manhattan_stat$Type!='ref',]


COG_stat_raw <- read.table('Functional_annotation/all_annots.tsv',  
                                    header=F, sep='\t', stringsAsFactors = F)
colnames(COG_stat_raw) <- c('Type','Group')
COG_stat_ref <- data.frame(Type=c(),Group=c())

for (ind in 1:nrow(COG_stat_raw)){
  if (nchar(COG_stat_raw[ind,]$Type)==1){
    COG_stat_ref <- rbind(COG_stat_ref, COG_stat_raw[ind,])
  } else{
    Group <- COG_stat_raw[ind,2]
    COG_split <- strsplit(COG_stat_raw[ind,1], "")[[1]]
    for (cog in COG_split){
      COG_stat_ref <- rbind(COG_stat_ref, data.frame(Type=cog, Group=Group))
    }
  }
}


COG_stat_fin <- COG_stat_ref %>% group_by(Group) %>% table() %>% as.data.frame() 

ggplot(COG_stat_fin[COG_stat_fin$Freq!=0,], aes(x=reorder(Type,-Freq),y=Freq,fill=Group)) +scale_size(guide = 'none')+
  geom_bar(stat="identity", width=1, alpha=0.75, col='black')+
  theme_bw()+scale_fill_manual(values = c("#75522B", '#1F732D'))+
  theme(axis.title.x=element_text(face="bold", color="black", 
                                  size=20),
        panel.background =  element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(color='black', 
                                   size=16),
        axis.title.y = element_text(face="bold", color="black", 
                                    size=18),
        legend.title=element_text(face="bold",size=14), 
        legend.text=element_text(size=12),
        axis.text.x = element_text(color='black', 
                                   angle = 60, vjust = 1, 
                                   size = 20, hjust = 1))+
  guides(fill= guide_legend(title="Plant organ"))+
  ylab('Sum of hits')+
  xlab("COG type")+
  facet_wrap(~Group, scales='free')



GO_roots <- read.table('Functional_annotation/roots.emapper.annotations.tsv',  
                           header=F, sep='\t', stringsAsFactors = F)
GO_shoots <- read.table('Functional_annotation/shoots.emapper.annotations.tsv',  
                        header=F, sep='\t', stringsAsFactors = F)

GO_roots <- GO_roots[,c(1,10)]
GO_shoots <- GO_shoots[,c(1,10)]
GO_all <-rbind (GO_roots, GO_shoots)

GO_root_samples <- as.vector(rep(FALSE,33))
GO_shoot_samples <- as.vector(rep(FALSE,33))

GO_roots_list <- list()
for (ind in 1:nrow(GO_roots)){
  GO_vec=c()
  GO_split <- strsplit(GO_roots[ind,2], ",")[[1]]
  for (GO in GO_split){
    if (GO=='-'){
      GO_vec <- c(GO_vec, '')} else {
        GO_root_samples[ind] <- TRUE
        GO_vec <- c(GO_vec, GO)
      }
    
  }
  GO_roots_list[[GO_roots[ind,1]]] <- GO_vec
  
}

GO_shoots_list <- list()
for (ind in 1:nrow(GO_shoots)){
  GO_vec=c()
  GO_split <- strsplit(GO_shoots[ind,2], ",")[[1]]
  for (GO in GO_split){
    if (GO=='-'){
      GO_vec <- c(GO_vec, '')} else {
        GO_vec <- c(GO_vec, GO)
        GO_shoot_samples[ind+15] <- TRUE
      }
    
  }
  GO_shoots_list[[GO_shoots[ind,1]]] <- GO_vec
}

All_GO_annot <- c(GO_roots_list, GO_shoots_list) 
geneNames_roots <- names(All_GO_annot)

myInterestingGenes_roots  <- geneNames_roots[GO_root_samples]
geneList_roots  <- factor(as.integer(geneNames_roots  %in% myInterestingGenes_roots ))
names(geneList_roots ) <- geneNames_roots 
str(geneList_roots )

#without cycle:
#geneList_roots <- factor(as.integer(geneNames_roots  %in% GO_roots$V1 ))
#names(geneList_roots ) <- geneNames_roots 

GOdata_roots  <- new("topGOdata", ontology = "MF", allGenes = geneList_roots ,
              annot = annFUN.gene2GO, gene2GO = All_GO_annot)

resultFisher_roots  <- runTest(GOdata_roots, algorithm = "classic", statistic = "fisher")

allRes_roots  <- GenTable(GOdata_roots , classicFisher = resultFisher_roots )

write.table(allRes_roots, 'GO_all_roots.tsv', sep="\t", row.names = FALSE, quote=FALSE)


geneNames_shoots <- names(All_GO_annot)

myInterestingGenes_shoots   <- geneNames_roots[GO_shoot_samples]
geneList_shoots  <- factor(as.integer(geneNames_shoots  %in% myInterestingGenes_shoots))
names(geneList_shoots) <- geneNames_shoots
str(geneList_shoots)
GOdata_shoots <- new("topGOdata", ontology = "MF", allGenes = geneList_shoots ,
                     annot = annFUN.gene2GO, gene2GO = All_GO_annot)

resultFisher_shoots  <- runTest(GOdata_shoots, algorithm = "classic", statistic = "fisher")

allRes_shoots  <- GenTable(GOdata_shoots, classicFisher = resultFisher_shoots)
write.table(allRes_shoots, 'GO_all_shoots.tsv', sep="\t", row.names = FALSE, quote=FALSE)


#Signifficant GO
#XP_015615023.1, XP_015647554.1, AAG44566.1 , XP_015635326.1, 
#XP_015635326.1, XP_015647554.1, XP_015647554.1,XP_015617521.1, XP_015613809.1, XP_015643864.1

signif_roots_names <- c('XP_015615023.1', 'XP_015647554.1', 'AAG44566.1' , 'XP_015635326.1',
                        'XP_015635326.1', 'XP_015647554.1', 'XP_015647554.1','XP_015617521.1', 'XP_015613809.1', 'XP_015643864.1')

myInterestingGenes_roots  <- geneNames_roots[geneNames_roots %in% signif_roots_names]
geneList_roots  <- factor(as.integer(geneNames_roots  %in% myInterestingGenes_roots ))
names(geneList_roots ) <- geneNames_roots 
str(geneList_roots )
GOdata_roots  <- new("topGOdata", ontology = "MF", allGenes = geneList_roots ,
                     annot = annFUN.gene2GO, gene2GO = All_GO_annot)

resultFisher_roots  <- runTest(GOdata_roots, algorithm = "classic", statistic = "fisher")

allRes_roots  <- GenTable(GOdata_roots , classicFisher = resultFisher_roots )

write.table(allRes_roots, 'GO_signif_roots.tsv', sep="\t", row.names = FALSE, quote=FALSE)

signif_shoots_names <- c('XP_015627405.1', 'XP_015619408', 'XP_015639965.1', 'XP_015643023.1', 'XP_015646728.1', '2002393A', 'XP_015643207.1', 
                         'QFR04205.1', 'XP_015619408', 'XP_015640756.1', 'XP_015625719.1', 'XP_015633459.1', 'XP_015616895.1', 'XP_015616895.1')

myInterestingGenes_shoots   <- geneNames_roots[geneNames_roots %in% signif_shoots_names]
geneList_shoots  <- factor(as.integer(geneNames_shoots  %in% myInterestingGenes_shoots))
names(geneList_shoots) <- geneNames_shoots
str(geneList_shoots)
GOdata_shoots <- new("topGOdata", ontology = "MF", allGenes = geneList_shoots ,
                     annot = annFUN.gene2GO, gene2GO = All_GO_annot)

resultFisher_shoots  <- runTest(GOdata_shoots, algorithm = "classic", statistic = "fisher")

allRes_shoots  <- GenTable(GOdata_shoots, classicFisher = resultFisher_shoots)

write.table(allRes_shoots, 'GO_signif_shoots.tsv', sep="\t", row.names = FALSE, quote=FALSE)


AT_content <- read.table('For_processing_scripts/AT_content_sum.csv',  
                        header=T, sep='\t', stringsAsFactors = F)

t.test(AT_content[AT_content$Organ=='Root',2],AT_content[AT_content$Organ=='Shoot',2]) #0.4124 - insignifficant

t.test(AT_content[AT_content$Organ=='Root' & AT_content$Gene.type=='Genes',2],
       AT_content[AT_content$Organ=='Root' & AT_content$Gene.type=='Promoters',2]) #p-value = 4.156e-06

t.test(AT_content[AT_content$Organ=='Shoot' & AT_content$Gene.type=='Genes',2],
       AT_content[AT_content$Organ=='Shoot' & AT_content$Gene.type=='Promoters',2]) #p-value = 0.0004747

AT_content_mean = data.frame(AT_content=c(mean(AT_content[AT_content$Organ=='Root' & AT_content$Gene.type=='Genes',2]),
                                          mean(AT_content[AT_content$Organ=='Root' & AT_content$Gene.type=='Promoters',2]),
                                          mean(AT_content[AT_content$Organ=='Shoot' & AT_content$Gene.type=='Genes',2]),
                                          mean(AT_content[AT_content$Organ=='Shoot' & AT_content$Gene.type=='Promoters',2])),
                             Organ=c('Root','Root','Shoot', 'Shoot'),
                             Type= c('Genes', 'Promoters','Genes', 'Promoters')
                             )

AT_content_mean$Organ <- factor(AT_content_mean$Organ, levels = c('Shoot','Root'))

ggplot(AT_content_mean, aes(x= Organ, y=AT_content, fill=Type)) +
  geom_bar(position="dodge", stat="identity",width=0.6, alpha=0.9, col='black')+
  theme_bw()+scale_fill_manual(values = c("#CAC583", '#95539F'))+
  theme(axis.title.x=element_text(face="bold", color="black", 
                                  size=20),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(color='black', 
                                   size=16),
        axis.title.y = element_text(face="bold", color="black", 
                                    size=18),
        legend.title=element_text(face="bold",size=18), 
        legend.text=element_text(size=16),
        axis.text.x =  element_text(color='black', 
                                    size=16))+
  ylab("AT content, %")+
  xlab("Plant organ")


Spots_clusters <- read.table('Clusterization/Clusters_intensity.tsv',  
                         header=T, sep='\t', stringsAsFactors = F)
Spots_clusters$Cluster <- as.factor(Spots_clusters$Cluster)

hclust_clusters_shoots <- Spots_clusters[Spots_clusters$Method=='h_clust' & Spots_clusters$Organ=='Shoot',]
hclust_clusters_shoots$Condition <- factor(hclust_clusters_shoots$Condition, levels = c('Control','Anoxia','Reaeration')) 

kmeans_clusters_shoots <- Spots_clusters[Spots_clusters$Method=='K-means' & Spots_clusters$Organ=='Shoot',]
kmeans_clusters_shoots$Condition <- factor(kmeans_clusters_shoots$Condition, levels = c('Control','Anoxia','Reaeration')) 

hclust_clusters_roots <- Spots_clusters[Spots_clusters$Method=='h_clust' & Spots_clusters$Organ=='Root',]
hclust_clusters_roots$Condition <- factor(hclust_clusters_roots$Condition, levels = c('Control','Anoxia','Reaeration')) 

kmeans_clusters_roots <- Spots_clusters[Spots_clusters$Method=='K-means' & Spots_clusters$Organ=='Root',]
kmeans_clusters_roots$Condition <- factor(kmeans_clusters_roots$Condition, levels = c('Control','Anoxia','Reaeration')) 


#Shoot_hclust_intensity_supp

ggplot(hclust_clusters_roots[hclust_clusters_roots$Group %in% c('All', 'Significant'),], aes(x = Cluster, y = Intensity, fill = Condition))+
  geom_col(position = position_dodge())+
  scale_fill_manual(values = c("#678667","#375DB6","#A5432B"))+
  facet_wrap(~Group, strip.position = 'top', scales='free_x')+
  theme(strip.placement = "outside")+
  theme_bw()+
  theme(axis.title.x=element_text(face="bold", color="black", 
                                  size=20),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(color='black', 
                                   size=16),
        axis.title.y = element_text(face="bold", color="black", 
                                    size=18),
        legend.title=element_text(face="bold",size=18), 
        legend.text=element_text(size=16),
        axis.text.x =  element_text(color='black', 
                                    size=16),
        strip.text.x = element_text(size = 18, colour = "black"))



Shoot_Kmeans_heatmap = read.table('Clusterization/heatmaps/Shoot_k-means_clusters_similarity.tsv',  
                                  row.names=1,header=T, sep='\t', stringsAsFactors = F)
Shoot_h_clust_heatmap = read.table('Clusterization/heatmaps/Shoot_hclust_clusters_similarity.tsv',  
                                   row.names=1,header=T, sep='\t', stringsAsFactors = F)
Root_Kmeans_heatmap = read.table('Clusterization/heatmaps/Root_k-means_clusters_similarity.tsv',  
                                 row.names=1,header=T, sep='\t', stringsAsFactors = F)
Root_h_clust_heatmap = read.table('Clusterization/heatmaps/Root_hclust_clusters_similarity.tsv',  
                                  row.names=1,header=T, sep='\t', stringsAsFactors = F)


melted_dist <- melt(as.matrix(Shoot_Kmeans_heatmap), na.rm = TRUE)
melted_dist$value <- round(melted_dist$value, 2)

melted_dist <- melt(as.matrix(Shoot_h_clust_heatmap), na.rm = TRUE)
melted_dist$value <- round(melted_dist$value, 2)

melted_dist <- melt(as.matrix(Root_Kmeans_heatmap), na.rm = TRUE)
melted_dist$value <- round(melted_dist$value, 2)

melted_dist <- melt(as.matrix(Root_h_clust_heatmap), na.rm = TRUE)
melted_dist$value <- round(melted_dist$value, 2)

ggplot(data = melted_dist, aes(Var2, Var1, fill = value))+
  geom_tile(color = 'black')+ 
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 4.5)+
  scale_fill_gradient(low = "#DCDCDC", high = "#A60B0B", limit = c(0.3,1), space = "Lab", 
                      name="Tree simillarity") +
  theme_bw() + xlab('Group') + ylab('Group') +
  theme( axis.text.y = element_text(color='black', 
                                    size=12),
         axis.title.y=element_text(color="black", 
                                   size=14),
         panel.background = element_rect(fill = "#DCDCDC"),
         panel.grid.minor = element_blank(), 
         panel.grid.major = element_line(size = 0.1, linetype = 'dashed',
                                         colour = "black"),
         legend.position = "none",
         axis.text.x = element_text(color='black', 
                                    angle = 60, vjust = 1, 
                                    size = 12, hjust = 1),
         axis.title.x = element_text( color="black", 
                                      size=14))+
  coord_fixed()


inter_clust_similarity = read.table('Clusterization/Inter_clusters_similarity.tsv',  
                                    header=F, sep='\t', stringsAsFactors = F)

ggplot(inter_clust_similarity, aes(x=V2,y=V3,fill=V2)) +
  geom_bar(stat="identity", width=1, alpha=0.75, col='black')+
  theme_bw()+scale_fill_manual(values = c("#a60b0b", '#2980b9', '#9A9EAB','#bdbd00'))+
  theme(axis.title.x=element_text(face="bold", color="black", 
                                  size=20),
        panel.background =  element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(color='black', 
                                   size=16),
        axis.title.y = element_text(face="bold", color="black", 
                                    size=18),
        legend.title=element_text(face="bold",size=14), 
        legend.text=element_text(size=12),
        axis.text.x = element_text(color='black', 
                                   angle = 60, vjust = 1, 
                                   size = 12, hjust = 1))+
  guides(fill= guide_legend(title="Spots group"))+
  ylab('Mean Jaccard similarity')+
  xlab("Spots group")+
  facet_wrap(~V1)


protein_props = read.table('For_processing_scripts/proteins_mass_pI.csv',  
                                    header=T, sep='\t', stringsAsFactors = F)
protein_props$Organ <- factor(protein_props$Organ, levels = c('Shoot','Root'))

ggplot(protein_props, aes(x=pI,y=real.pI, fill=Organ, color=Organ)) +
         geom_point(shape=21, col='black', size=6, alpha=0.7)+
  geom_smooth( method='lm', se=F, size=2)+
  scale_fill_manual(values = c('#1F732D',"#75522B"))+
  scale_color_manual(values = c('#1F732D',"#75522B"))+
  theme_bw()+
  theme(axis.title.x=element_text(face="bold", color="black", 
                                  size=20),
        panel.background =  element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(color='black', 
                                   size=14),
        axis.title.y = element_text(face="bold", color="black", 
                                    size=20),
        legend.title=element_text(face="bold",size=14), 
        legend.text=element_text(size=12),
        axis.text.x = element_text(color='black', 
                                   angle = 60, vjust = 1, 
                                   size = 14, hjust = 1))+
  guides( fill=guide_legend(show = FALSE))+
  ylab('Expected pI')+
  xlab("Observed pI")

cor.test(protein_props$pI, protein_props$real.pI, method = c("pearson", "kendall", "spearman")) # p-value = 1.697e-07
cor.test(protein_props$Mass, protein_props$real.mass, method = c("pearson", "kendall", "spearman")) #p-value = 6.777e-11


ggplot(protein_props, aes(x=Mass,y=real.mass, fill=Organ, color=Organ)) +
  geom_point(shape=21, col='black', size=6, alpha=0.7)+
  geom_smooth( method='lm', se=F, size=2)+
  scale_fill_manual(values = c('#1F732D',"#75522B"))+
  scale_color_manual(values = c('#1F732D',"#75522B"))+
  theme_bw()+
  theme(axis.title.x=element_text(face="bold", color="black", 
                                  size=20),
        panel.background =  element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(color='black', 
                                   size=14),
        axis.title.y = element_text(face="bold", color="black", 
                                    size=20),
        legend.title=element_text(face="bold",size=14), 
        legend.text=element_text(size=12),
        axis.text.x = element_text(color='black', 
                                   angle = 60, vjust = 1, 
                                   size = 14, hjust = 1))+
  guides( fill=guide_legend(show = FALSE))+
  ylab('Expected mass, Da')+
  xlab("Observed mass, Da")


TF_sites_rice_all = read.table('TF_predictions/Rice/Coords/Coords_rice_promotor_sites_coords_all.tsv',  
                           header=T, sep='\t', stringsAsFactors = F)
Signif_TF_sites_coords = read.table('TF_predictions/Rice/Coords/Coords_rice_promotor_sites_coords_signif.tsv',  
                                          header=T, sep='\t', stringsAsFactors = F)

#TF_sites_rice_roots <- TF_sites_rice_all[TF_sites_rice_all$Root_genes!=0,]
#TF_sites_rice_shoots <- TF_sites_rice_all[TF_sites_rice_all$Shoot_genes!=0,]

TF_sites_rice_roots <- Signif_TF_sites_coords[Signif_TF_sites_coords$Root_genes!=0,]
TF_sites_rice_shoots <- Signif_TF_sites_coords[Signif_TF_sites_coords$Shoot_genes!=0,]


ggplot(TF_sites_rice_shoots, aes(x=Site,y=Shoot_genes, size=Root_calls, fill=TF)) +
  geom_point(shape=21, col='black', alpha=0.9)+
  theme_bw()+
  theme(axis.title.x=element_text(face="bold", color="black", 
                                  size=20),
        panel.background =  element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(color='black', 
                                   size=14),
        axis.title.y = element_text(face="bold", color="black", 
                                    size=20),
        legend.title=element_text(face="bold",size=14), 
        legend.text=element_text(size=12),
        axis.text.x = element_text(color='black', 
                                   angle = 60, vjust = 1, 
                                   size = 14, hjust = 1))+
  guides( size="none", fill=guide_legend(override.aes = list(size=5)))+
  ylab('Number of signals')+
  xlab("Upstream coordinate, bp")+
  scale_fill_manual(values = wes_palette("Darjeeling2", 34, type = "continuous"))
  

ggplot(TF_sites_rice_shoots, aes(x=Site,y=Shoot_genes, size=Shoot_calls, fill=TF)) +
  geom_point(shape=21, col='black', alpha=0.7)+
  theme_bw()+
  theme(axis.title.x=element_text(face="bold", color="black", 
                                  size=20),
        panel.background =  element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(color='black', 
                                   size=14),
        axis.title.y = element_text(face="bold", color="black", 
                                    size=20),
        legend.title=element_text(face="bold",size=14), 
        legend.text=element_text(size=12),
        axis.text.x = element_text(color='black', 
                                   angle = 60, vjust = 1, 
                                   size = 14, hjust = 1))+
  guides( size="none", fill=guide_legend(override.aes = list(size=5)))+
  ylab('Number of signals')+
  xlab("Upstream coordinate, bp")


Signif_TF_sites_coords$Calls_all <- Signif_TF_sites_coords$Root_calls+ Signif_TF_sites_coords$Shoot_calls
Signif_TF_sites_coords$Genes_all <- Signif_TF_sites_coords$Root_genes+ Signif_TF_sites_coords$Shoot_genes


ggplot(Signif_TF_sites_coords, aes(x=Site,y=Genes_all, size=Calls_all, fill=TF)) +
  geom_point(shape=21, col='black', alpha=0.7)+
  theme_bw()+
  theme(axis.title.x=element_text(face="bold", color="black", 
                                  size=20),
        panel.background =  element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(color='black', 
                                   size=14),
        axis.title.y = element_text(face="bold", color="black", 
                                    size=20),
        legend.title=element_text(face="bold",size=14), 
        legend.text=element_text(size=12),
        axis.text.x = element_text(color='black', 
                                   angle = 60, vjust = 1, 
                                   size = 14, hjust = 1))+
  guides( size="none", fill=guide_legend(override.aes = list(size=5)))+
  ylab('Number of signals')+
  xlab("Upstream coordinate, bp")

#TF_genes_total_num <- TF_sites_rice_all %>% group_by(Site) %>% summarize(Root = sum (Root_genes),
#                                                                         Shoot = sum(Shoot_genes),
#                                                                         All = sum(Genes_all))
#TF_sites_total_num <- TF_sites_rice_all %>% group_by(Site) %>% summarize(Root = sum (Root_calls),
#                                                                         Shoot = sum(Shoot_calls),
#                                                                         All = sum(Calls_all))

TF_genes_total_num <- Signif_TF_sites_coords %>% group_by(Site) %>% summarize(Root = sum (Root_genes),
                                                                         Shoot = sum(Shoot_genes),
                                                                         All = sum(Genes_all))
TF_sites_total_num <- Signif_TF_sites_coords %>% group_by(Site) %>% summarize(Root = sum (Root_calls),
                                                                         Shoot = sum(Shoot_calls),
                                                                         All = sum(Calls_all))

ggplot(TF_genes_total_num, aes(x=Site,y=Root)) +
  geom_bar(stat='identity',fill='darkgrey' ,alpha=0.9)+ #alpha=0.7)
  theme_bw()+
  theme(axis.title.x=element_text(face="bold", color="black", 
                                  size=20),
        panel.background =  element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(color='black', 
                                   size=14),
        axis.title.y = element_text(face="bold", color="black", 
                                    size=20),
        legend.title=element_text(face="bold",size=14), 
        legend.text=element_text(size=12),
        axis.text.x = element_text(color='black', 
                                   angle = 60, vjust = 1, 
                                   size = 14, hjust = 1))+
  guides( size="none", fill=guide_legend(override.aes = list(size=5)))+
  ylab('Number of signals')+
  xlab("Upstream coordinate, bp")


TF_sites_sum <- TF_sites_rice_all%>% group_by(TF) %>%  summarize(root_sum=sum(Root_calls),
                                                                               shoot_sum=sum(Shoot_calls))
TF_sites_sum$All_TFs <- TF_sites_sum$root_sum+TF_sites_sum$shoot_sum

TF_sites_sum_melt <- melt(TF_sites_sum, id.vars = c('TF'))



TF_sites_rice_summary = read.table('TF_predictions/Rice/rice_promotor_sites_numbers_all.tsv',  
                               header=T, sep='\t', stringsAsFactors = F)

TF_sites_rice_summary_melt_sites <- melt(TF_sites_rice_summary[,c(1,3,5,7)], id.vars = c('TF'))
TF_sites_rice_summary_melt_genes<- melt(TF_sites_rice_summary[,c(1,2,4,6)], id.vars = c('TF'))
Signif_TF_sites_rice_summary = read.table('TF_predictions/Rice/rice_promotor_sites_annotated_numbers_signif.tsv',  
                                   header=T, sep='\t', stringsAsFactors = F)

Signif_TF_sites_rice_summary_sum <- Signif_TF_sites_rice_summary[, c(1,2,3,4,8,9,10)]%>% group_by(TF) %>% summarise(Root_genes=sum(Root_genes),
                                                                                                                    Shoot_genes=sum(Shoot_genes),
                                                                                                                    All_genes=sum(All_genes),
                                                                                                                    Shoot_sites=sum(Shoot_sites),
                                                                                                                    Root_sites=sum(Root_sites),
                                                                                                                    All_sites=sum(All_sites))
Signif_TF_sites_rice_summary_melt_sites <- melt(Signif_TF_sites_rice_summary_sum[,c(1,5,6,7)], id.vars = c('TF'))
Signif_TF_sites_rice_summary_melt_genes <- melt(Signif_TF_sites_rice_summary_sum[,c(1,2,3,4)], id.vars = c('TF'))

TF_genes_roots_and_shoots <- TF_sites_rice_summary_melt_genes[TF_sites_rice_summary_melt_genes$variable %in% c('Root_genes','Shoot_genes'),]
levels(TF_genes_roots_and_shoots$variable) <- c("Root", "Shoot",'NA')


TF_genes_roots_and_shoots$percent = 0
TF_genes_roots_and_shoots[TF_genes_roots_and_shoots$variable=='Root',4] <- TF_genes_roots_and_shoots[TF_genes_roots_and_shoots$variable=='Root',3]/17
TF_genes_roots_and_shoots[TF_genes_roots_and_shoots$variable=='Shoot',4] <- TF_genes_roots_and_shoots[TF_genes_roots_and_shoots$variable=='Shoot',3]/19


TF_sites_roots_and_shoots <- TF_sites_rice_summary_melt_sites[TF_sites_rice_summary_melt_sites$variable %in% c('Root_sites','Shoot_sites'),]
levels(TF_sites_roots_and_shoots$variable) <- c("Root", "Shoot",'NA')



TF_sites_roots_and_shoots$percent = TF_sites_roots_and_shoots$value/TF_genes_roots_and_shoots$value

cor.test(TF_sites_roots_and_shoots[TF_sites_roots_and_shoots$variable=='Root', 3], 
         TF_sites_roots_and_shoots[TF_sites_roots_and_shoots$variable=='Shoot', 3], method='pearson') #p-value = 1.178e-14

cor.test(TF_sites_roots_and_shoots[TF_sites_roots_and_shoots$variable=='Root', 4], 
         TF_sites_roots_and_shoots[TF_sites_roots_and_shoots$variable=='Shoot', 4], method='pearson')


TF_genes_roots_and_shoots$variable <- factor(TF_genes_roots_and_shoots$variable, levels = c('Shoot','Root'))
TF_genes_roots_and_shoots$variable <- factor(TF_genes_roots_and_shoots$variable, levels = c('Shoot','Root'))

ggplot(TF_genes_roots_and_shoots, aes(x=reorder(TF,-percent),y=percent*100, fill=variable))+
  geom_bar(position="dodge",stat="identity", width=1, alpha=0.65, col='black')+ #fill=#75522B, #1F732D, #D8A466
  theme_bw()+
  scale_fill_manual(values = c('#1F732D',"#75522B"))+
  theme(axis.title.x=element_text(face="bold", color="black", 
                                  size=20),
        panel.background =  element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(color='black', 
                                   size=16),
        axis.title.y = element_text(face="bold", color="black", 
                                    size=18),
        legend.title=element_text(face="bold",size=14), 
        legend.text=element_text(size=12),
        axis.text.x = element_text(color='black', 
                                   angle = 60, vjust = 1, 
                                   size = 11, hjust = 1))+
  guides(fill= guide_legend(title="Organ"))+
  ylab('Percent of genes')+ #Percentage of genes    Number of sites
  xlab("TF") #Num_sites_signif_all

ggplot(Signif_TF_sites_rice_summary_melt_genes[Signif_TF_sites_rice_summary_melt_genes$variable=='Root_genes',], aes(x=reorder(TF,-value),y=value))+
  geom_bar(stat="identity", width=1, alpha=0.65, col='black', fill='#75522B')+ #fill=#75522B, #1F732D, #D8A466
  theme_bw()+
  theme(axis.title.x=element_text(face="bold", color="black", 
                                  size=20),
        panel.background =  element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(color='black', 
                                   size=16),
        axis.title.y = element_text(face="bold", color="black", 
                                    size=18),
        legend.title=element_text(face="bold",size=14), 
        legend.text=element_text(size=12),
        axis.text.x = element_text(color='black', 
                                   angle = 60, vjust = 1, 
                                   size = 11, hjust = 1))+
  guides(fill= guide_legend(title="Group"))+
  ylab('Number of genes')+
  xlab("TF")


All_sites<-with(TF_sites_rice_summary_melt_sites[TF_sites_rice_summary_melt_sites$variable=='All_sites', c(1, 3)],  
        TF_sites_rice_summary_melt_sites[TF_sites_rice_summary_melt_sites$variable=='All_sites', c(1, 3)][order( TF) , ])

Sig_sites<-with(Signif_TF_sites_rice_summary_melt_sites[Signif_TF_sites_rice_summary_melt_sites$variable=='All_sites', c(1, 3)],  
                Signif_TF_sites_rice_summary_melt_sites[Signif_TF_sites_rice_summary_melt_sites$variable=='All_sites', c(1, 3)][order( TF) , ])

All_genes<-with(TF_sites_rice_summary_melt_genes[TF_sites_rice_summary_melt_genes$variable=='All_genes', c(1, 3)],  
                TF_sites_rice_summary_melt_genes[TF_sites_rice_summary_melt_genes$variable=='All_genes', c(1, 3)][order( TF) , ])

Sig_genes<-with(Signif_TF_sites_rice_summary_melt_genes[Signif_TF_sites_rice_summary_melt_genes$variable=='All_genes', c(1, 3)],  
                Signif_TF_sites_rice_summary_melt_genes[Signif_TF_sites_rice_summary_melt_genes$variable=='All_genes', c(1, 3)][order( TF) , ])                               



Root_genes<-with(TF_sites_rice_summary_melt_genes[TF_sites_rice_summary_melt_genes$variable=='Root_genes', c(1, 3)],  
                TF_sites_rice_summary_melt_genes[TF_sites_rice_summary_melt_genes$variable=='Root_genes', c(1, 3)][order( TF) , ])

Shoot_genes<-with(TF_sites_rice_summary_melt_genes[TF_sites_rice_summary_melt_genes$variable=='Shoot_genes', c(1, 3)],  
                  TF_sites_rice_summary_melt_genes[TF_sites_rice_summary_melt_genes$variable=='Shoot_genes', c(1, 3)][order( TF) , ])                               

Root_sites<-with(TF_sites_rice_summary_melt_sites[TF_sites_rice_summary_melt_sites$variable=='Root_sites', c(1, 3)],  
                 TF_sites_rice_summary_melt_sites[TF_sites_rice_summary_melt_sites$variable=='Root_sites', c(1, 3)][order( TF) , ])

Shoot_sites<-with(TF_sites_rice_summary_melt_sites[TF_sites_rice_summary_melt_sites$variable=='Shoot_sites', c(1, 3)],  
                  TF_sites_rice_summary_melt_sites[TF_sites_rice_summary_melt_sites$variable=='Shoot_sites', c(1, 3)][order( TF) , ])                               

Root_shoot_comparision <- cbind(Root_genes, Root_sites, Shoot_genes, Shoot_sites)
Root_shoot_comparision <- Root_shoot_comparision[,c(1,2,4,6,8)]
colnames(Root_shoot_comparision) <- c('TF', 'Root_genes','Root_sites',
                                      'Shoot_genes','Shoot_sites')


Signif_all_comparision <- cbind(All_sites, Sig_sites, All_genes, Sig_genes)
Signif_all_comparision <- Signif_all_comparision[,c(1,2,4,6,8)]
colnames(Signif_all_comparision) <- c('TF', 'All_sites','Sig_sites',
                                      'All_genes','Sig_genes')

cor.test(Signif_all_comparision$All_sites, Signif_all_comparision$Sig_sites, method='pearson') #p-value < 2.2e-16
cor.test(Signif_all_comparision$All_genes, Signif_all_comparision$Sig_genes, method='pearson') #p-value < 2.2e-16

cor.test(Signif_all_comparision$All_genes, Signif_all_comparision$All_sites, method='pearson') #p-value = 4.769e-06
cor.test(Signif_all_comparision$Sig_genes, Signif_all_comparision$Sig_sites, method='pearson')

ggplot(Root_shoot_comparision, aes(x=Root_genes,y=Shoot_genes)) +
  geom_point(shape=21, alpha=0.75, col='black', size = 5, fill='#09ADAE')+
  geom_line(stat="smooth", method='lm', se=F, size=2, col='grey', alpha=0.8)+
  theme_bw()+
  theme(axis.title.x=element_text(face="bold", color="black", 
                                  size=18),
        panel.background =  element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(color='black', 
                                   size=16),
        axis.title.y = element_text(face="bold", color="black", 
                                    size=18),
        legend.title=element_text(face="bold",size=14), 
        legend.text=element_text(size=12),
        axis.text.x = element_text(color='black', 
                                   angle = 60, vjust = 1, 
                                   size = 16, hjust = 1),
        legend.position = 'none')+
  ylab('Shoot genes')+ # Significant sites Root genes Shoot genes
  xlab("Root genes") #genes


Signif_conditions_melt_genes <- melt(Signif_TF_sites_rice_summary[,c(1:5)], id.vars = c('TF','Type'))
Signif_conditions_melt_sites <- melt(Signif_TF_sites_rice_summary[,c(1,5,8:10)], id.vars = c('TF','Type'))


ggplot(Signif_conditions_melt_sites[Signif_conditions_melt_sites$variable=='Shoot_sites' & Signif_conditions_melt_sites$Type=='Re-aeration',], 
       aes(x=reorder(TF,-value),y=value))+
  geom_bar(stat="identity", width=1, alpha=0.65, col='black', fill='#A5432B')+ #c("#375DB6","#678667","#A5432B")
  theme_bw()+
  theme(axis.title.x=element_text(face="bold", color="black", 
                                  size=20),
        panel.background =  element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(color='black', 
                                   size=16),
        axis.title.y = element_text(face="bold", color="black", 
                                    size=18),
        legend.title=element_text(face="bold",size=14), 
        legend.text=element_text(size=12),
        axis.text.x = element_text(color='black', 
                                   angle = 60, vjust = 1, 
                                   size = 11, hjust = 1))+
  guides(fill= guide_legend(title="Group"))+
  ylab('Number of genes')+
  xlab("TF")



AT_signif_all <- read.table('TF_predictions/AT/signif/All_signif_mearged.csv',  
                            header=T, sep='\t', stringsAsFactors = F)
AT_signif_roots <- read.table('TF_predictions/AT/signif/Root_signif_mearged.csv',  
                            header=T, sep='\t', stringsAsFactors = F)
AT_signif_shoots  <- read.table('TF_predictions/AT/signif/Shoot_signif_mearged.csv',  
                            header=T, sep='\t', stringsAsFactors = F)

AT_signif_anoxia_roots <- read.table('TF_predictions/AT/Conditions/Root_anoxia.csv',  
                            header=T, sep='\t', stringsAsFactors = F)
AT_signif_anoxia_shoots <- read.table('TF_predictions/AT/Conditions/Shoot_anoxia.csv',  
                                     header=T, sep='\t', stringsAsFactors = F)
AT_signif_anoxia_all <- read.table('TF_predictions/AT/Conditions/All_anoxia.csv',  
                                     header=T, sep='\t', stringsAsFactors = F)

AT_signif_control_roots <- read.table('TF_predictions/AT/Conditions/Root_control.csv',  
                                     header=T, sep='\t', stringsAsFactors = F)
AT_signif_control_shoots <- read.table('TF_predictions/AT/Conditions/Shoot_control.csv',  
                                      header=T, sep='\t', stringsAsFactors = F)
AT_signif_control_all <- read.table('TF_predictions/AT/Conditions/All_control.csv',  
                                   header=T, sep='\t', stringsAsFactors = F)

AT_signif_reaeration_shoots <- read.table('TF_predictions/AT/Conditions/Shoot_reaeration.csv',  
                                    header=T, sep='\t', stringsAsFactors = F)

AT_res <- AT_signif_all %>% group_by(Family) %>%  summarize(gene_num=max(Number.of.genes ))
AT_res <- AT_signif_roots %>% group_by(Family) %>%  summarize(gene_num=max(Number.of.genes ))
AT_res <- AT_signif_shoots %>% group_by(Family) %>%  summarize(gene_num=max(Number.of.genes ))

AT_res <- AT_signif_anoxia_roots %>% group_by(Family) %>%  summarize(gene_num=max(Number.of.genes ))
AT_res <- AT_signif_anoxia_shoots %>% group_by(Family) %>%  summarize(gene_num=max(Number.of.genes ))
AT_res <- AT_signif_anoxia_all %>% group_by(Family) %>%  summarize(gene_num=max(Number.of.genes ))

AT_res <- AT_signif_control_roots %>% group_by(Family) %>%  summarize(gene_num=max(Number.of.genes ))
AT_res <- AT_signif_control_shoots %>% group_by(Family) %>%  summarize(gene_num=max(Number.of.genes ))
AT_res <- AT_signif_control_all %>% group_by(Family) %>%  summarize(gene_num=max(Number.of.genes ))
AT_res <- AT_signif_reaeration_shoots %>% group_by(Family) %>%  summarize(gene_num=max(Number.of.genes ))


AT_res <- AT_signif_reaeration_shoots %>% group_by(Family) %>%  summarize(gene_num=max(Number.of.genes ))


#AT_root_anoxia
ggplot(AT_res[AT_res$Family!='- other',], aes(x=reorder(Family,-gene_num),y=gene_num))+
  geom_bar(stat="identity", width=1, alpha=0.65, col='black', fill='#A5432B')+ # "#375DB6","#678667","#A5432B"
  theme_bw()+
  theme(axis.title.x=element_text(face="bold", color="black", 
                                  size=20),
        panel.background =  element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(color='black', 
                                   size=16),
        axis.title.y = element_text(face="bold", color="black", 
                                    size=18),
        legend.title=element_text(face="bold",size=14), 
        legend.text=element_text(size=12),
        axis.text.x = element_text(color='black', 
                                   angle = 60, vjust = 1, 
                                   size = 11, hjust = 1))+
  guides(fill= guide_legend(title="Group"))+
  ylab('Number of genes')+
  xlab("TF")


ggplot(AT_res[AT_res$Family!='- other',], 
       aes(x=reorder(Family,-gene_num),y=gene_num))+
  geom_bar(stat="identity", width=1, alpha=0.65, col='black', fill='#D8A466')+ #fill=#75522B, #1F732D, #D8A466
  theme_bw()+
  theme(axis.title.x=element_text(face="bold", color="black", 
                                  size=20),
        panel.background =  element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(color='black', 
                                   size=16),
        axis.title.y = element_text(face="bold", color="black", 
                                    size=18),
        legend.title=element_text(face="bold",size=14), 
        legend.text=element_text(size=12),
        axis.text.x = element_text(color='black', 
                                   angle = 60, vjust = 1, 
                                   size = 11, hjust = 1))+
  guides(fill= guide_legend(title="Group"))+
  ylab('Number of genes')+
  xlab("TF")


AT_dist_root <- read.table('For_processing_scripts/root_AT_content_per_site',  
                           header=T, sep='\t', stringsAsFactors = F)
AT_dist_shoot <- read.table('For_processing_scripts/shoot_AT_content_per_site',  
                           header=T, sep='\t', stringsAsFactors = F)


ggplot(AT_dist_shoot, 
       aes(y=AT_content, x=Site))+
  geom_line(stat="identity", alpha=0.65, size=0.5, col='#1F732D')+
  geom_line(stat="smooth", method='lm', se=F, size=2, col='black', alpha=0.4)+
  theme_bw()+
  theme(axis.title.x=element_text(face="bold", color="black", 
                                  size=20),
        panel.background =  element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(color='black', 
                                   size=16),
        axis.title.y = element_text(face="bold", color="black", 
                                    size=18),
        legend.title=element_text(face="bold",size=14), 
        legend.text=element_text(size=12),
        axis.text.x = element_text(color='black', 
                                   angle = 60, vjust = 1, 
                                   size = 11, hjust = 1))+
  guides(fill= guide_legend(title="Group"))+
  ylab('AT content, %')+
  xlab("Upstream coordinate, bp")


AA_content <- read.table('For_processing_scripts/AA_percent_all_dist.csv',  
                         header=T, sep='\t', stringsAsFactors = F)

AA_content$Organ <- factor(AA_content$Organ, levels =c('shoot','root'))
levels(AA_content$Organ) <- c('Shoot','Root')

ggplot(AA_content, aes(x= reorder(Amino_acid,-Cum_freq), y=Cum_freq, fill=Organ)) +
  geom_bar(position="dodge", stat="identity",width=0.9, alpha=0.7, col='black')+
  theme_bw()+
  scale_fill_manual(values = c( '#1F732D', "#75522B"))+
  theme(axis.title.x=element_text(face="bold", color="black", 
                                  size=20),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(color='black', 
                                   size=16),
        axis.title.y = element_text(face="bold", color="black", 
                                    size=18),
        legend.title=element_text(face="bold",size=18), 
        legend.text=element_text(size=16),
        axis.text.x =  element_text(color='black', 
                                    size=16))+
  ylab("Amino acid content, %")+
  xlab("Plant organ")


extr_AA <-read.table('For_processing_scripts/extr_amino_acid_content.csv',  
                     header=T, sep='\t', stringsAsFactors = F)

extr_AA_df <- data.frame(group = c('K', 'A', 'G', 'E'), 
                         p.value = c(0,0,0,0), 
                         p.adj = c(0,0,0,0))

extr_AA_df[1,2] <- t.test(extr_AA[extr_AA$AA=='K' & extr_AA$organ=='shoot',1], extr_AA[extr_AA$AA=='K' & extr_AA$organ=='root',1])$p.value
extr_AA_df[2,2] <-t.test(extr_AA[extr_AA$AA=='A' & extr_AA$organ=='shoot',1], extr_AA[extr_AA$AA=='A' & extr_AA$organ=='root',1])$p.value
extr_AA_df[3,2] <-t.test(extr_AA[extr_AA$AA=='G' & extr_AA$organ=='shoot',1], extr_AA[extr_AA$AA=='G' & extr_AA$organ=='root',1])$p.value
extr_AA_df[4,2] <-t.test(extr_AA[extr_AA$AA=='E' & extr_AA$organ=='shoot',1], extr_AA[extr_AA$AA=='E' & extr_AA$organ=='root',1])$p.value

extr_AA_df$p.adj <- p.adjust(extr_AA_df$p.value, method='fdr')

write.table(extr_AA_df, 'AA_test_pval.tsv', sep="\t", row.names = FALSE, quote=FALSE)


#TF by condition All
Condition_TF_sites_rice_summary = read.table('TF_predictions/Rice/rice_promotor_sites_numbers_all.tsv',  
                                          header=T, sep='\t', stringsAsFactors = F)

#Condition_TF_sites_rice_summary = read.table('/home/anton/fbi/TF_processing_results/rice_promotor_sites_annotated_numbers_all_adjusted_roots.tsv',  
#                                                                                       header=T, sep='\t', stringsAsFactors = F)

Condition_TF_sites_rice_summary_sum <- Condition_TF_sites_rice_summary[, c(1,2,3,4,5,8,9,10)]%>% group_by(TF, Type) %>% summarise(Root_genes=sum(Root_genes),
                                                                                                                    Shoot_genes=sum(Shoot_genes),
                                                                                                                    All_genes=sum(All_genes),
                                                                                                                    Shoot_sites=sum(Shoot_sites),
                                                                                                                    Root_sites=sum(Root_sites),
                                                                                                                    All_sites=sum(All_sites))
#num genes:
#Root: anoxia - 3(1), control - 8(2), re-aeration - 6(5)
#Shoot: anoxia - 3(2), control - 8(4), re-aeration - 8(7)


Condition_TF_sites_rice_summary_melt_genes<- melt(Condition_TF_sites_rice_summary_sum[,c(1,2,3,4)], id.vars = c('TF', 'Type'))

Condition_TF_sites_rice_summary_melt_genes$percent = 0


##### Root
Condition_TF_sites_rice_summary_melt_genes[Condition_TF_sites_rice_summary_melt_genes$Type=='Anoxia' 
                                           & Condition_TF_sites_rice_summary_melt_genes$variable=='Root_genes' 
                                           ,5] <- Condition_TF_sites_rice_summary_melt_genes[Condition_TF_sites_rice_summary_melt_genes$Type=='Anoxia' 
                                                                                             & Condition_TF_sites_rice_summary_melt_genes$variable=='Root_genes' 
                                                                                             ,4]/1

Condition_TF_sites_rice_summary_melt_genes[Condition_TF_sites_rice_summary_melt_genes$Type=='Control' 
                                           & Condition_TF_sites_rice_summary_melt_genes$variable=='Root_genes' 
                                           ,5] <- Condition_TF_sites_rice_summary_melt_genes[Condition_TF_sites_rice_summary_melt_genes$Type=='Control' 
                                                                                             & Condition_TF_sites_rice_summary_melt_genes$variable=='Root_genes' 
                                                                                             ,4]/2


Condition_TF_sites_rice_summary_melt_genes[Condition_TF_sites_rice_summary_melt_genes$Type=='Re-aeration' 
                                           & Condition_TF_sites_rice_summary_melt_genes$variable=='Root_genes' 
                                           ,5] <- Condition_TF_sites_rice_summary_melt_genes[Condition_TF_sites_rice_summary_melt_genes$Type=='Re-aeration' 
                                                                                             & Condition_TF_sites_rice_summary_melt_genes$variable=='Root_genes' 
                                                                                             ,4]/5

##### Shoot
Condition_TF_sites_rice_summary_melt_genes[Condition_TF_sites_rice_summary_melt_genes$Type=='Anoxia' 
                                           & Condition_TF_sites_rice_summary_melt_genes$variable=='Shoot_genes' 
                                           ,5] <- Condition_TF_sites_rice_summary_melt_genes[Condition_TF_sites_rice_summary_melt_genes$Type=='Anoxia' 
                                                                                             & Condition_TF_sites_rice_summary_melt_genes$variable=='Shoot_genes' 
                                                                                             ,4]/2

Condition_TF_sites_rice_summary_melt_genes[Condition_TF_sites_rice_summary_melt_genes$Type=='Control' 
                                           & Condition_TF_sites_rice_summary_melt_genes$variable=='Shoot_genes' 
                                           ,5] <- Condition_TF_sites_rice_summary_melt_genes[Condition_TF_sites_rice_summary_melt_genes$Type=='Control' 
                                                                                             & Condition_TF_sites_rice_summary_melt_genes$variable=='Shoot_genes' 
                                                                                             ,4]/4


Condition_TF_sites_rice_summary_melt_genes[Condition_TF_sites_rice_summary_melt_genes$Type=='Re-aeration' 
                                           & Condition_TF_sites_rice_summary_melt_genes$variable=='Shoot_genes' 
                                           ,5] <- Condition_TF_sites_rice_summary_melt_genes[Condition_TF_sites_rice_summary_melt_genes$Type=='Re-aeration' 
                                                                                             & Condition_TF_sites_rice_summary_melt_genes$variable=='Shoot_genes' 
                                                                                             ,4]/7


Condition_TF_sites_rice_summary_melt_sites <- melt(Condition_TF_sites_rice_summary_sum[,c(1,2,6,7)], id.vars = c('TF', 'Type'))

#### Sites (190 rows - all) 1:95 176 - signif (1:88)
Condition_TF_sites_rice_summary_melt_sites$percent = Condition_TF_sites_rice_summary_melt_sites$value/rbind(Condition_TF_sites_rice_summary_melt_genes[c(96:190 ),],Condition_TF_sites_rice_summary_melt_genes[c(1:95),])[,4]
Condition_TF_sites_rice_summary_melt_sites[is.infinite(Condition_TF_sites_rice_summary_melt_sites$percent),5] <- 0
Condition_TF_sites_rice_summary_melt_sites[is.nan(Condition_TF_sites_rice_summary_melt_sites$percent),5] <- 0

Condition_TF_sites_rice_summary_melt_sites$Type <- factor(Condition_TF_sites_rice_summary_melt_sites$Type, levels = c("Control", "Anoxia", "Re-aeration"))
Condition_TF_sites_rice_summary_melt_genes$Type <- factor(Condition_TF_sites_rice_summary_melt_genes$Type, levels = c("Control", "Anoxia", "Re-aeration"))


write.xlsx(Condition_TF_sites_rice_summary_melt_genes, "Rice_TF_genes_fixed.xlsx", sheetName="Sheet1",
           col.names=TRUE, row.names=FALSE, append=FALSE)

write.xlsx(Condition_TF_sites_rice_summary_melt_sites, "Rice_TF_sites_fixed.xlsx", sheetName="Sheet1",
           col.names=TRUE, row.names=FALSE, append=FALSE)


#Add zeros
organs_vec <- c('Root_genes','Shoot_genes')
Cond_vec <- c("Control", "Anoxia", "Re-aeration")

for (TF in Condition_TF_sites_rice_summary_melt_genes$TF){
  for (check_org in organs_vec){
    for (chek_cond in Cond_vec){
      sub_df <- Condition_TF_sites_rice_summary_melt_genes[Condition_TF_sites_rice_summary_melt_genes$TF==TF & 
                                                             Condition_TF_sites_rice_summary_melt_genes$variable==check_org & 
                                                             Condition_TF_sites_rice_summary_melt_genes$Type==chek_cond,] 
      if (nrow(sub_df)<1){
        Condition_TF_sites_rice_summary_melt_genes <- rbind(Condition_TF_sites_rice_summary_melt_genes, 
                                                            data.frame(TF=TF, Type = chek_cond,
                                                                       variable=check_org,
                                                                       value=0,percent=0))
      }
    }
  }
  
}


organs_vec_sites <- c('Root_sites','Shoot_sites')
for (TF in Condition_TF_sites_rice_summary_melt_sites$TF){
  for (check_org in organs_vec_sites){
    for (chek_cond in Cond_vec){
      sub_df <- Condition_TF_sites_rice_summary_melt_sites[Condition_TF_sites_rice_summary_melt_sites$TF==TF & 
                                                             Condition_TF_sites_rice_summary_melt_sites$variable==check_org & 
                                                             Condition_TF_sites_rice_summary_melt_sites$Type==chek_cond,] 
      if (nrow(sub_df)<1){
        Condition_TF_sites_rice_summary_melt_sites <- rbind(Condition_TF_sites_rice_summary_melt_sites, 
                                                            data.frame(TF=TF, Type = chek_cond,
                                                                       variable=check_org,
                                                                       value=0, percent=0))
      }
    }
  }
  
}

#20 width for PDF
ggplot(Condition_TF_sites_rice_summary_melt_sites, aes(x= reorder(TF,-percent), y = percent, fill = Type))+
  geom_bar(position="dodge", stat="identity",width=0.9, alpha=0.7, col='black')+
  scale_fill_manual(values = c("#678667","#375DB6", "#A5432B"))+
  facet_wrap(~variable, strip.position = 'top')+
  theme(strip.placement = "outside")+
  theme_bw()+
  theme(axis.title.x=element_text(face="bold", color="black", 
                                  size=20),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(color='black', 
                                   size=16),
        axis.title.y = element_text(face="bold", color="black", 
                                    size=18),
        legend.title=element_text(face="bold",size=18), 
        legend.text=element_text(size=16),
        axis.text.x =  element_text(color='black', 
                                    size=14, angle = 45,
                                    vjust = 1, hjust=1),
        strip.text.x = element_text(size = 18, colour = "black"))+
  ylab("Number of genes")+
  xlab("TF")

#AT all organs



AT_signif_roots <- read.table('TF_predictions/AT/All/roots_common_sites.csv',  
                              header=T, sep='\t', stringsAsFactors = F)
AT_signif_shoots  <- read.table('TF_predictions/AT/All/shoots_common_sites.csv',  
                                header=T, sep='\t', stringsAsFactors = F)

AT_signif_roots$Type <- 'Root'
AT_signif_shoots$Type <- 'Shoot'
AT_mearged <- rbind(AT_signif_roots, AT_signif_shoots)
AT_mearged <- AT_mearged[AT_mearged$Family!='- other',]

AT_res <- AT_mearged %>% group_by(Family, Type) %>%  summarize(gene_num=max(Number.of.genes),
                                                                               Sites_sum=sum(Sum.of.TFBSs.in.total))

AT_res$genes_percent <- 0
AT_res[AT_res$Type=='Root',5] <- AT_res[AT_res$Type=='Root',3]/17
AT_res[AT_res$Type=='Shoot',5] <- AT_res[AT_res$Type=='Shoot',3]/19

AT_res$sites_percent <- 0
AT_res[,6] <- AT_res$Sites_sum/AT_res$gene_num

signle_families <- AT_res %>%  group_by(Family) %>% summarize(num=n()) %>% as.data.frame() %>% 
  filter(num==1) %>% dplyr::select(Family)

AT_res <- AT_res %>% as.data.frame()

for (TF in signle_families$Family){
  pres_organ <- AT_res[AT_res$Family==TF, 2] %>% as.vector()
  if (!'Root' %in% pres_organ) {
    AT_res <- rbind(AT_res, data.frame(Family=TF, Type='Root',gene_num=0,Sites_sum=0,
                                       genes_percent=0, sites_percent=0))
  } else {
    AT_res <- rbind(AT_res, data.frame(Family=TF, Type='Shoot',gene_num=0,Sites_sum=0,
                                       genes_percent=0, sites_percent=0))
  }
}

# 
AT_res$Type <- factor(AT_res$Type, levels=c('Shoot','Root'))

ggplot(AT_res, aes(x=reorder(Family,-genes_percent),y=genes_percent*100, fill=Type))+
  geom_bar(position = position_dodge2(width = 1, preserve = "single"), 
           stat="identity", alpha=0.65, col='black')+ #fill=#75522B, #1F732D, #D8A466
  theme_bw()+
  scale_fill_manual(values = c( '#1F732D',"#75522B"))+
  theme(axis.title.x=element_text(face="bold", color="black", 
                                  size=20),
        panel.background =  element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(color='black', 
                                   size=16),
        axis.title.y = element_text(face="bold", color="black", 
                                    size=18),
        legend.title=element_text(face="bold",size=14), 
        legend.text=element_text(size=12),
        axis.text.x = element_text(color='black', 
                                   angle = 60, vjust = 1, 
                                   size = 11, hjust = 1))+
  guides(fill= guide_legend(title="Organ"))+
  ylab('Number of sites')+
  xlab("TF") #Num_sites_signif_all



AT_anoxia_roots <- read.table('TF_predictions/AT/Conditions/AT_root_anoxia.csv',  
                                     header=T, sep='\t', stringsAsFactors = F)
AT_anoxia_shoots <- read.table('TF_predictions/AT/Conditions/AT_shoot_anoxia.csv',  
                                      header=T, sep='\t', stringsAsFactors = F)

AT_control_roots <- read.table('TF_predictions/AT/Conditions/AT_root_control.csv',  
                                      header=T, sep='\t', stringsAsFactors = F)
AT_control_shoots <- read.table('TF_predictions/AT/Conditions/AT_shoot_control.csv',  
                                       header=T, sep='\t', stringsAsFactors = F)

AT_reaeration_roots <- read.table('TF_predictions/AT/Conditions/AT_root_re-aeration.csv',  
                                          header=T, sep='\t', stringsAsFactors = F)
AT_reaeration_shoots <- read.table('TF_predictions/AT/Conditions/AT_shoot_re-aeration.csv',  
                                          header=T, sep='\t', stringsAsFactors = F)

AT_anoxia_roots$organ <- 'Root'
AT_anoxia_roots$condition <- 'Anoxia'

AT_anoxia_shoots$organ <- 'Shoot'
AT_anoxia_shoots$condition <- 'Anoxia'

AT_control_roots$organ <- 'Root'
AT_control_roots$condition <- 'Control'

AT_control_shoots$organ <- 'Shoot'
AT_control_shoots$condition <- 'Control'

AT_reaeration_roots$organ <- 'Root'
AT_reaeration_roots$condition <- 'Re-aeration'

AT_reaeration_shoots$organ <- 'Shoot'
AT_reaeration_shoots$condition <- 'Re-aeration'


AT_mearged <- rbind(AT_anoxia_roots,AT_anoxia_shoots, AT_control_roots, AT_control_shoots, AT_reaeration_roots, AT_reaeration_shoots)

AT_summary <- AT_mearged %>% group_by(Family, organ, condition) %>%  summarize(gene_num=max(Number.of.genes),
                                                              site_num=sum(Sum.of.TFBSs.in.total)/max(Number.of.genes)) %>% as.data.frame()


AT_summary[AT_summary$organ=='Root'& AT_summary$condition=='Anoxia',4] <- AT_summary[AT_summary$organ=='Root'& AT_summary$condition=='Anoxia',4]/5
AT_summary[AT_summary$organ=='Root'& AT_summary$condition=='Control',4] <- AT_summary[AT_summary$organ=='Root'& AT_summary$condition=='Control',4]/9
AT_summary[AT_summary$organ=='Root'& AT_summary$condition=='Re-aeration',4] <- AT_summary[AT_summary$organ=='Root'& AT_summary$condition=='Re-aeration',4]/3


AT_summary[AT_summary$organ=='Shoot'& AT_summary$condition=='Anoxia',4] <- AT_summary[AT_summary$organ=='Shoot'& AT_summary$condition=='Anoxia',4]/3
AT_summary[AT_summary$organ=='Shoot'& AT_summary$condition=='Control',4] <- AT_summary[AT_summary$organ=='Shoot'& AT_summary$condition=='Control',4]/8
AT_summary[AT_summary$organ=='Shoot'& AT_summary$condition=='Re-aeration',4] <- AT_summary[AT_summary$organ=='Shoot'& AT_summary$condition=='Re-aeration',4]/8


AT_summary$condition <- factor(AT_summary$condition, levels = c("Control", "Anoxia", "Re-aeration"))


write.xlsx(AT_summary, "AT_TF_genes.xlsx", sheetName="Sheet1",
           col.names=TRUE, row.names=FALSE, append=FALSE)

organs_vec <- c('Root','Shoot')
Cond_vec <- c("Control", "Anoxia", "Re-aeration")

single_TF_AT <- AT_summary %>%  group_by(Family) %>% summarize(num=n()) %>% as.data.frame() %>% 
  filter(num!=6) %>% dplyr::select(Family)

for (TF in AT_summary$Family){
      #AT_summary <- rbind(AT_summary, data.frame(Family=TF, organ=m_org,gene_num=0,Sites_sum=0,
      #                                           condition = m_cond))
  for (check_org in organs_vec){
    for (chek_cond in Cond_vec){
      sub_df <- AT_summary[AT_summary$Family==TF & AT_summary$organ==check_org & AT_summary$condition==chek_cond,] 
      if (nrow(sub_df)<1){
        AT_summary <- rbind(AT_summary, data.frame(Family=TF, organ=check_org,
                                                   condition = chek_cond,gene_num=0,site_num=0))
      }
    }
  }

}

#20 width for PDF
ggplot(AT_summary, aes(x= reorder(Family,-gene_num), y = gene_num*100, fill = condition))+
  geom_bar(position="dodge", stat="identity",width=0.9, alpha=0.7, col='black')+
  scale_fill_manual(values = c("#678667","#375DB6","#A5432B"))+
  facet_wrap(~organ, strip.position = 'top')+
  theme(strip.placement = "outside")+
  theme_bw()+
  theme(axis.title.x=element_text(face="bold", color="black", 
                                  size=20),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(color='black', 
                                   size=16),
        axis.title.y = element_text(face="bold", color="black", 
                                    size=18),
        legend.title=element_text(face="bold",size=18), 
        legend.text=element_text(size=16),
        axis.text.x =  element_text(color='black', 
                                    size=14, angle = 45,
                                    vjust = 1, hjust=1),
        strip.text.x = element_text(size = 18, colour = "black"))+
  ylab("Number of sites")+
  xlab("TF")


Lysine_per_prots_roots <- read.table('For_processing_scripts/Lysin_coords_root.csv',  
                              header=F, sep='\t', stringsAsFactors = F)

Lysine_per_prots_shoots <- read.table('For_processing_scripts/Lysin_coords_shoot.csv',  
                                     header=F, sep='\t', stringsAsFactors = F)


lys_root_percent <- Lysine_per_prots_roots %>% group_by(V3) %>% summarize(prot_percent=V1/V4) %>% as.data.frame()
lys_root_percent$organ <- 'Root'

lys_shoot_percent <- Lysine_per_prots_shoots %>% group_by(V3) %>% summarize(prot_percent=V1/V4) %>% as.data.frame()
lys_shoot_percent$organ <- 'Shoot'

Lysine_all_prots <- rbind(lys_root_percent, lys_shoot_percent)
Lysine_all_prots$V3

Lysine_all_prots <- Lysine_all_prots[Lysine_all_prots$V3 %in% c('fructose-bisphosphate aldolase 3, cytoplasmic ',
                                                                '20 kDa chaperonin',
                                                                'heat shock cognate 70 kDa protein 2 ', 
                                                                'oxygen-evolving complex protein 1',
                                                                'oxygen-evolving enhancer protein 2',
                                                                'stromal 70 kDa heat shock-related protein'),]


ggplot(Lysine_all_prots, aes(x=prot_percent,y=V3, fill=organ)) +
  geom_point(shape=21, col='black', alpha=0.6, size=3)+
  theme_bw()+
  scale_fill_manual(values = c("#75522B", '#1F732D'))+
  theme(axis.title.x=element_text(face="bold", color="black", 
                                  size=20),
        panel.background =  element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(color='black', 
                                   size=14),
        axis.title.y = element_text(face="bold", color="black", 
                                    size=20),
        legend.title=element_text(face="bold",size=14), 
        legend.text=element_text(size=12),
        axis.text.x = element_text(color='black', 
                                   angle = 60, vjust = 1, 
                                   size = 14, hjust = 1))+
  guides( size="none", fill=guide_legend(override.aes = list(size=5)))+
  ylab('Protein')+
  xlab("Coordinate")

##OD for gene groups

OD_genes_with_errors <- read.table('PQquest_data/err_mean_ODS_genes.csv',  
                                     header=T, sep='\t', stringsAsFactors = F)

OD_genes_with_errors <- OD_genes_with_errors[OD_genes_with_errors$Signif=='s',-c(9,11)]

colnames(OD_genes_with_errors)[4] <- 'Re-aeration'
colnames(OD_genes_with_errors)[8] <- "err_Re-aeration"


OD_genes_with_errors_melt <- melt(OD_genes_with_errors[,c(1,2,3,4,9)], id.vars = c('Accession','Organ')) %>% as.data.frame()

for (ind in 1:nrow(OD_genes_with_errors)){
  acc_str=OD_genes_with_errors[ind,5]
  acc_split=strsplit(acc_str,split='[', fixed=TRUE)[[1]][1]
  OD_genes_with_errors[ind,5] <- acc_split
}

colnames(OD_genes_with_errors_melt)[3] <- 'Condition'
colnames(OD_genes_with_errors_melt)[4] <- 'OD'

OD_genes_with_errors_melt$error <- 0
OD_genes_with_errors_melt$annotation <-'-'



for (ind in 1:nrow(OD_genes_with_errors_melt)){
  acc_str <- OD_genes_with_errors_melt[ind,1]
  condition <- OD_genes_with_errors_melt[ind,3]
  
  annot_name <- OD_genes_with_errors[OD_genes_with_errors$Accession==acc_str,5]
  OD_genes_with_errors_melt[ind,6] <- annot_name
  
  error_anoxia <- OD_genes_with_errors[OD_genes_with_errors$Accession==acc_str,6]
  error_control <- OD_genes_with_errors[OD_genes_with_errors$Accession==acc_str,7]
  error_reaeration <- OD_genes_with_errors[OD_genes_with_errors$Accession==acc_str,8]

  if (condition=='Anoxia'){
    OD_genes_with_errors_melt[ind,5] <- error_anoxia
  } else if (condition=='Control'){
    OD_genes_with_errors_melt[ind,5] <- error_control    
  } else if (condition=='Re-aeration'){
    OD_genes_with_errors_melt[ind,5] <- error_reaeration
  }
  
}

OD_genes_with_errors_melt$Condition <- factor(OD_genes_with_errors_melt$Condition, levels = c('Control','Anoxia','Re-aeration')) 


group1 <-c('XP_015616895.1','XP_015619408','XP_015625719.1', 'XP_015633459.1', 'XP_015643864.1','XP_015613809.1')
group2 <-c('QFR04205.1','XP_015616495.1','AAQ19031.1')
group3 <-c('2002393A','XP_015639965.1','XP_015640756.1','XP_015643023.1','XP_015643207.1',
           'XP_015646728.1','XP_015615023.1','XP_015647554.1',
           'AAG44566.1','XP_015635326.1')
group4 <-c('XP_015627405.1','XP_015617521.1')

OD_genes_with_errors_melt1<-OD_genes_with_errors_melt[OD_genes_with_errors_melt$Accession %in% group1,]
OD_genes_with_errors_melt2<-OD_genes_with_errors_melt[OD_genes_with_errors_melt$Accession %in% group2,]
OD_genes_with_errors_melt3<-OD_genes_with_errors_melt[OD_genes_with_errors_melt$Accession %in% group3,]
OD_genes_with_errors_melt4<-OD_genes_with_errors_melt[OD_genes_with_errors_melt$Accession %in% group4,]

#group1_ODs
  
#17/11
#10/11
#21/11
#7/11
ggplot(OD_genes_with_errors_melt, aes(x=reorder(annotation,-OD) , y = OD, fill = Condition))+
  geom_bar(position="dodge", stat="identity",width=0.9, alpha=0.7, col='black')+
  scale_fill_manual(values = c("#678667","#375DB6","#A5432B"))+
  geom_errorbar(aes(ymin=OD-error, ymax=OD+error), width=.2,
                position=position_dodge(.9)) +
  theme_bw()+
  theme(axis.title.x=element_text(face="bold", color="black", 
                                  size=20),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(color='black', 
                                   size=16),
        axis.title.y = element_text(face="bold", color="black", 
                                    size=18),
        legend.title=element_text(face="bold",size=18), 
        legend.text=element_text(size=16),
        axis.text.x =  element_text(color='black', 
                                    size=14, angle = 65,
                                    vjust = 1, hjust=1),
        strip.text.x = element_text(size = 18, colour = "black"))+
  ylab("Mean OD")+
  xlab("Protein") +
  coord_cartesian(ylim = c(14, 23))+
  facet_grid(. ~ Organ, scales = "free", space = "free")


