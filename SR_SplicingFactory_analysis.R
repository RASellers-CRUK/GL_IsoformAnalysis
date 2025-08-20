library(SplicingFactory)
library(SummarizedExperiment)
library(org.Mm.eg.db)
library(tximport)
library(rtracklayer)

gtf_path <- "/scratch/wsspaces/rsellers-GL96_phase2/resources/gencode.vM37.annotation.gtf.gz"
txdb <- makeTxDbFromGFF(file=gtf_path)
k <- AnnotationDbi::keys(txdb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")
write.table(tx2gene,"/scratch/wsspaces/rsellers-GL96_phase2/resources/gencode.vM37.annotation.tx2gene.txt",sep="\t",row.names=F)

gtf<-import(gtf_path)

anot<-read.table("../../HP_analysis/HP_cells.txt", sep = "\t", header = T, stringsAsFactors = F)
anot<-anot[order(anot$condition),]
files<-list.files("../salmon","quant.sf",full.names = T,recursive = T)
names(files)<-basename(dirname(files))

files<-files[match(anot$sampleID,names(files))]
coldata<-data.frame(files,names=names(files),condition=anot$condition[match(names(files),anot$sampleID)])

salmonQuant<-tximport(files, type = "salmon", tx2gene = tx2gene,txOut = TRUE)

genes <- rownames(salmonQuant$abundance)

readcounts<-salmonQuant$counts
tokeep <- rowSums(readcounts > 5) > 5

readcounts <- readcounts[tokeep, ]
genes <- genes[tokeep]

abund<-salmonQuant$counts
abund<-abund[tokeep,]

laplace_entropy <- calculate_diversity(readcounts, tx2gene$GENEID[match(genes,tx2gene$TXNAME)],  method = "laplace",
                                       norm = TRUE, verbose = TRUE,)

colData(laplace_entropy) <- cbind(colData(laplace_entropy),
                                  sample_type = anot$condition[match(colData(laplace_entropy)$samples,anot$sampleID)])
conds<-unique(colData(laplace_entropy)$sample_type)
comps<-expand.grid(conds,conds)
comps<-comps[which(comps$Var1!=comps$Var2),]
comps<-as.data.frame(t(apply(comps,1,sort)))
comps<-dplyr::distinct(comps)
colnames(comps)<-c("cond1","cond2")

calc_diffs<-lapply(1:nrow(comps),function(x){
  cond1<-comps$cond1[x]
  cond2<-comps$cond2[x]
  hold<-laplace_entropy[,which(laplace_entropy$sample_type==cond1 | laplace_entropy$sample_type==cond2)]
  entropy_significance <- calculate_difference(x = hold, samples = "sample_type",
                                               control = cond1,
                                               method = "mean", test = "wilcoxon",
                                               verbose = TRUE)
  entropy_significance$symbol<-gtf$gene_name[match(entropy_significance$genes,gtf$gene_id)]
  return(entropy_significance)
})
names(calc_diffs)<-apply(comps,1,function(x){paste(x[1],"vs",x[2],sep="_")})
calc_diffs<-lapply(calc_diffs,function(x){data.frame(x,sig=x$adjusted_p_values<=0.05)})
saveRDS(calc_diffs,"SplicingFactory_Laplace_diffs.rds")

dir.create("SplicingFactory_QC_plots",showWarnings = F)
dir.create("SplicingFactory_res_tables",showWarnings = F)
for(x in 1:length(calc_diffs)) {
  png(paste("SplicingFactory_QC_plots",paste0(names(calc_diffs)[x],"_QC_plot_entropy_hist.png"),sep="/"),res=300,width=1800,height=1600)
  print(ggplot(calc_diffs[[x]]) + geom_histogram(aes(log2_fold_change,fill=sig),bins = 100,alpha=0.7,col="white") + xlim(c(-1,1)) +
    theme_classic() + coord_cartesian(expand=F) + theme(legend.position="top") + scale_fill_manual(values=c("grey88","grey9")) + ggtitle(gsub("_"," ",names(calc_diffs)[x])))
  dev.off()
  write.table(calc_diffs[[x]],paste("SplicingFactory_res_tables",paste0(names(calc_diffs)[x],"_res.txt"),sep="/"),sep="\t",row.names=F)
  calc_diffs[[x]]$mean<-apply(calc_diffs[[x]][,c(2,3)],1,function(y){mean(y)})
  
  png(paste("SplicingFactory_QC_plots",paste0(names(calc_diffs)[x],"_QC_plot_entropy.png"),sep="/"),res=300,width=1800,height=1600)
  print(ggplot(calc_diffs[[x]]) +
    geom_point(aes(x = mean, y = mean_difference,col=sig), size = 1,alpha=0.7) +
    theme_minimal() +
    labs(title = paste(gsub("_"," ",names(calc_diffs)[x]),"Normalized Laplace entropy"),
         subtitle = "Wilcoxon signed rank test",
         x = "Mean entropy",
         y = "Mean difference") +
    scale_color_manual(values=c("grey88","grey9")))
  dev.off()
}

table(calc_diffs$AGM_vs_Ysb$adjusted_p_values<0.05)
table(calc_diffs$AGM_vs_Ysm$adjusted_p_values<0.05)
table(calc_diffs$Ysb_vs_Ysm$adjusted_p_values<0.05)

cons<-sort(intersect(unique(calc_diffs$AGM_vs_Ysb$genes[calc_diffs$AGM_vs_Ysb$adjusted_p_values<=0.05]),
                     unique(calc_diffs$AGM_vs_Ysm$genes[calc_diffs$AGM_vs_Ysm$adjusted_p_values<=0.05])))

Ysm_cons<-calc_diffs$AGM_vs_Ysm[which(!is.na(match(calc_diffs$AGM_vs_Ysm$genes,cons))),]
Ysb_cons<-calc_diffs$AGM_vs_Ysb[which(!is.na(match(calc_diffs$AGM_vs_Ysb$genes,cons))),]

comp<-data.frame(ensID=Ysm_cons$genes,
                 symbol=gtf$gene_name[match(Ysm_cons$genes,gtf$gene_id)],
                 Ysm_p=Ysm_cons$adjusted_p_values,
                 Ysb_p=Ysb_cons$adjusted_p_values,
                 Ysm_prnk=rank(-log10(Ysm_cons$adjusted_p_values)),
                 Ysb_prnk=rank(-log10(Ysb_cons$adjusted_p_values)),
                 prnksm=rank(-log10(Ysm_cons$adjusted_p_values)) + rank(-log10(Ysb_cons$adjusted_p_values)),
                 prnksm_rnk=rank(rank(-log10(Ysm_cons$adjusted_p_values)) + rank(-log10(Ysb_cons$adjusted_p_values))))
write.table(comp,"HMP_LMP_vs_AGM_cons_SplicingFactory.txt",sep="\t",row.names=F)
