setwd("/scratch/wsspaces/rsellers-GL96_phase2/rsellers-GL96_phase2-1753600156/trans_analysis_HP/")
.libPaths("/data/compbio/resources/R4_3libs/")

library(SplicingFactory)
library(SummarizedExperiment)
library(org.Mm.eg.db)
library(tximport)
library(rtracklayer)
library(GenomicFeatures)
library(viridis)

anot<-read.table("GL96_phase2_trans_stats_table.txt", sep = "\t", header = T, stringsAsFactors = F)
counts<-read.table("GL96_phase2_trans_Estcount.txt", sep = "\t", header = T, stringsAsFactors = F,row.names = 1)
anot$ID<-gsub("-",".",anot$ID)

anot<-anot[!scuttle::isOutlier(anot$LS,type = "lower",nmads = 1),]
counts<-counts[,match(anot$ID,colnames(counts))]

anot$cond_MLaL<-dplyr::case_when(grepl("AGM_3|AGM_4",anot$cond,perl=T) ~ "HE_agm",
                                 grepl("YS_m",anot$cond,perl=T) ~ "HE_lmp",
                                 grepl("YS_b",anot$cond,perl=T) ~ "HE_emp",
                                 grepl("YS_p",anot$cond,perl=T) ~ "progenitors",
                                 grepl("AGM_5",anot$cond,perl=T) ~ "agm_hb")

genes<-do.call(rbind,strsplit(rownames(counts),"[|]"))[,6]

tokeep <- rowSums(counts > 5) > 5

counts <- counts[tokeep,]
genes <- genes[tokeep]

laplace_entropy <- calculate_diversity(counts, genes,  method = "laplace",
                                       norm = TRUE, verbose = TRUE)

colData(laplace_entropy) <- cbind(colData(laplace_entropy),
                                  sampleID = anot$ID,
                                  sample_type = anot$cond_MLaL)
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
  
  png(paste("SplicingFactory_QC_plots",paste0(names(calc_diffs)[x],"_QC_plot_entropy_cont.png"),sep="/"),res=300,width=1800,height=1600)
  print(ggplot(calc_diffs[[x]]) +
          geom_point(aes(x = mean, y = mean_difference,col=raw_p_values), size = 1,alpha=0.7) +
          theme_minimal() +
          labs(title = paste(gsub("_"," ",names(calc_diffs)[x]),"Normalized Laplace entropy"),
               subtitle = "Wilcoxon signed rank test",
               x = "Mean entropy",
               y = "Mean difference") +
          scale_color_viridis(option = 'F',direction = -1))
  dev.off()
}

table(calc_diffs$HE_agm_vs_HE_emp$adjusted_p_values<0.05)
table(calc_diffs$HE_agm_vs_HE_lmp$adjusted_p_values<0.05)
table(calc_diffs$HE_emp_vs_HE_lmp$adjusted_p_values<0.05)

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
