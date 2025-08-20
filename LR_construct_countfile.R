setwd("/scratch/wsspaces/rsellers-GL96_phase2/trans_analysis_HP/")
.libPaths("/data/compbio/resources/R4_3libs/")

col24<-c("#72da56",
         "#783cca",
         "#cfd749",
         "#d553c8",
         "#5e9441",
         "#726bd4",
         "#b89844",
         "#3b246d",
         "#67d3a8",
         "#d94980",
         "#c7d19d",
         "#8c3d89",
         "#d97c37",
         "#647fba",
         "#d2433c",
         "#89c7d8",
         "#89482d",
         "#c99dd0",
         "#4a582e",
         "#8c3f5c",
         "#577d79",
         "orange",
         "#d39992",
         "slateblue2",
         "springgreen1")

library(ggplot2)
library(ggrepel)
library(DESeq2)
library(edgeR)
library(viridis)
library(stringi)
library(VennDiagram)
library(ComplexHeatmap)
library(sva)

set.seed(777)
proj<-"phase2_trans"

files<-list.files("trans_process/","*_abundance.tsv",full.names=T,recursive=T)
ldat<-lapply(files,function(x){
  read.table(x,sep="\t",header=T)
})
cf_rn<-unique(do.call(c,lapply(ldat,function(x){
  unname(unlist(x$transcript_name))
})))
dat<-do.call(cbind,lapply(ldat,function(x){
  hold<-data.frame(count=rep(0,length(cf_rn)),row.names = cf_rn)
  hold$count[match(x$transcript_name,rownames(hold))]<-x$est_count
  return(hold)
}))
colnames(dat)<-make.names(sub("_abundance.tsv","",basename(files)),unique=T)

dat[]<-lapply(dat,function(x){as.numeric(x)})

fanot<-read.table("annotations.csv",sep=",",header=T)
fbcds<-read.table("barcoding.csv",sep=",",header=F,skip=1)
fbcds$UniqID<-fbcds$V3
fbcds$barcode<-fbcds$V17
idx<-match(fanot$UniqID,fbcds$UniqID)
fanot$barcode<-fbcds$barcode[idx]
fanot<-fanot[!is.na(fanot$barcode),]
fanot$MatchID<-gsub("-",".",paste0("SQK-NBD114-96_",sub("NB","barcode",fanot$barcode),
                      dplyr::case_when(fanot$RUN=="run1" ~ "",
                                       fanot$RUN=="run2" ~ ".1",
                                       fanot$RUN=="run3" ~ ".2")))
dat<-dat[,which(!is.na(match(colnames(dat),fanot$MatchID)))]
colnames(dat)[match(fanot$MatchID,colnames(dat))]<-fanot$UniqID
fanot<-fanot[order(fanot$parent.cluster),]
dat<-dat[,match(fanot$UniqID,colnames(dat))]

summary(colSums(dat))
anot<-data.frame(ID=colnames(dat),
                 cond=fanot$parent.cluster,
                 run=fanot$RUN)
anot$LS<-apply(dat,2,sum)
anot$nGenes_1<-apply(dat,2,function(x){length(which(x>=1))})
anot$nGenes_5<-apply(dat,2,function(x){length(which(x>=5))})
anot$nGenes_10<-apply(dat,2,function(x){length(which(x>=10))})
anot$nGenes_50<-apply(dat,2,function(x){length(which(x>=50))})
anot$nGenes_100<-apply(dat,2,function(x){length(which(x>=100))})
anot$top100pct<-apply(dat,2,function(x){sum(x[order(x,decreasing = T)[1:100]])/sum(x)*100})
write.table(anot,paste(proj,"stats_table.txt",sep="_"),sep="\t",row.names=F)

write.table(data.frame(geneID=rownames(dat),dat),paste(proj,"Estcount.txt",sep="_"),sep="\t",row.names=F)

ndat<-do.call(cbind,lapply(ldat,function(x){
  hold<-data.frame(count=rep(0,length(cf_rn)),row.names = cf_rn)
  hold$count[match(x$transcript_name,rownames(hold))]<-x$tpm
  return(hold)
}))
colnames(ndat)<-make.names(sub("_abundance.tsv","",basename(files)),unique=T)

ndat[]<-lapply(ndat,function(x){as.numeric(x)})

ndat<-ndat[,which(!is.na(match(colnames(ndat),fanot$MatchID)))]
colnames(ndat)[match(fanot$MatchID,colnames(ndat))]<-fanot$UniqID
fanot<-fanot[order(fanot$parent.cluster),]
ndat<-ndat[,match(fanot$UniqID,colnames(ndat))]

write.table(data.frame(geneID=rownames(ndat),ndat),paste(proj,"tpm.txt",sep="_"),sep="\t",row.names=F)

##
anot<-anot[!scuttle::isOutlier(anot$LS,type = "lower",nmads = 1),]
dat<-dat[,match(anot$ID,colnames(dat))]
dat[]<-lapply(dat,function(x){round(as.numeric(x))})
##

dds<-DESeqDataSetFromMatrix(dat,anot,~run+cond)
des<-DESeq(dds)
Rdat<-vst(des,blind = T)
rdat<-assay(Rdat)
ndat<-counts(des,normalized=T)
colnames(ndat)<-paste("DESeq-normalised",colnames(ndat),sep="_")
ldat<-log2(dat+1)
colnames(ldat)<-paste("Log2-normalised",colnames(ldat),sep="_")
cdat<-edgeR::cpm(dat)
colnames(cdat)<-paste("CPM-normalised",colnames(cdat),sep="_")

colv<-setNames(col24[1:length(unique(anot$cond))],unique(anot$cond))

pca_l<-prcomp(t(rdat))

png(paste0(proj,"_QC_PCA_cond.png"),res=300,height=3000,width=3000)
ggplot(as.data.frame(pca_l$x)) + geom_point(aes(PC1,PC2,col=anot$cond)) +
  theme_classic() + geom_text_repel(aes(PC1,PC2),label = rownames(pca_l$x), size=3) + 
  xlab(paste0("PC1: ",round(summary(pca_l)$importance[2,1]*100,2),"%")) +
  ylab(paste0("PC2: ",round(summary(pca_l)$importance[2,2]*100,2),"%")) + theme(legend.position = "bottom") +
  scale_color_manual(values=colv) + labs(col="cond")
dev.off()

png(paste0(proj,"_QC_PCA_15K_cond.png"),res=300,height=3000,width=3000)
plotPCA(Rdat,intgroup="cond",ntop=15000) + geom_text_repel(label=anot$ID) + scale_color_manual(values=colv) + theme_classic()
dev.off()

png(paste0(proj,"_QC_PCA_10K_cond.png"),res=300,height=3000,width=3000)
plotPCA(Rdat,intgroup="cond",ntop=10000) + geom_text_repel(label=anot$ID) + scale_color_manual(values=colv) + theme_classic()
dev.off()

png(paste0(proj,"_QC_PCA_5K_cond.png"),res=300,height=3000,width=3000)
plotPCA(Rdat,intgroup="cond",ntop=5000) + geom_text_repel(label=anot$ID) + scale_color_manual(values=colv) + theme_classic()
dev.off()

png(paste0(proj,"_QC_PCA_1K_cond.png"),res=300,height=3000,width=3000)
plotPCA(Rdat,intgroup="cond",ntop=1000) + geom_text_repel(label=anot$ID) + scale_color_manual(values=colv) + theme_classic()
dev.off()

png(paste0(proj,"_QC_PCA_500_cond.png"),res=300,height=3000,width=3000)
plotPCA(Rdat,intgroup="cond",ntop=500) + geom_text_repel(label=anot$ID) + scale_color_manual(values=colv) + theme_classic()
dev.off()

png(paste0(proj,"_QC_PCA_500_cond_nL.png"),res=300,height=3000,width=3000)
plotPCA(Rdat,intgroup="cond",ntop=500) + scale_color_manual(values=colv) + theme_classic()
dev.off()

png(paste0(proj,"_QC_PCA_500_run_nL.png"),res=300,height=3000,width=3000)
plotPCA(Rdat,intgroup="run",ntop=500) + scale_color_manual(values=rev(col24)) + theme_classic()
dev.off()

png(paste0(proj,"_QC_PCA_500_LS_nL.png"),res=300,height=3000,width=3000)
plotPCA(Rdat,intgroup="LS",ntop=500) + scale_color_viridis(option = 'E') + theme_classic()
dev.off()

write.table(data.frame(geneID=rownames(ndat),ndat),paste(proj,"DESeq2_count.txt",sep="_"),sep="\t",row.names=F)
write.table(data.frame(geneID=rownames(rdat),rdat),paste(proj,"rlog_count.txt",sep="_"),sep="\t",row.names=F)
write.table(data.frame(geneID=rownames(ldat),ldat),paste(proj,"log2_count.txt",sep="_"),sep="\t",row.names=F)
