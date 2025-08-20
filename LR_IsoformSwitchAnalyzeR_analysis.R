setwd("/scratch/wsspaces/rsellers-GL96_phase2/rsellers-GL96_phase2-1753600156/trans_analysis_HP/")
.libPaths("/data/compbio/resources/R4_3libs/")

library(IsoformSwitchAnalyzeR)
library(rtracklayer)

col24<-c("#72da56",
         "#783cca",
         "#cfd749",
         "#d553c8",
         "#5e9441",
         "#726bd4",
         "#b89844",
         "#3b246d",
         "#67d3a8",
         "#c7d19d",
         "#8c3d89",
         "#d97c37",
         "#647fba",
         "#d2433c",
         "#89c7d8",
         "#89482d",
         "#d94980",
         "#c99dd0",
         "#4a582e",
         "#8c3f5c",
         "#577d79",
         "orange",
         "#d39992",
         "slateblue2",
         "springgreen1")

anot<-read.table("GL96_phase2_trans_stats_table.txt", sep = "\t", header = T, stringsAsFactors = F)
anot$ID<-gsub("-",".",anot$ID)
counts<-read.table("GL96_phase2_trans_Estcount.txt", sep = "\t", header = T, stringsAsFactors = F,row.names = 1)

ndat<-read.table("GL96_phase2_trans_DESeq2_count.txt", sep = "\t", header = T, stringsAsFactors = F,row.names = 1)
colnames(ndat)<-sub("DESeq.normalised_","",colnames(ndat))

anot<-anot[!scuttle::isOutlier(anot$LS,type = "lower",nmads = 1),]
counts<-counts[,match(anot$ID,colnames(counts))]
ndat<-ndat[,match(anot$ID,colnames(ndat))]

anot$cond_MLaL<-dplyr::case_when(grepl("AGM_3|AGM_4",anot$cond,perl=T) ~ "HE_agm",
                                 grepl("YS_m",anot$cond,perl=T) ~ "HE_lmp",
                                 grepl("YS_b",anot$cond,perl=T) ~ "HE_emp",
                                 grepl("YS_p",anot$cond,perl=T) ~ "progenitors",
                                grepl("AGM_5",anot$cond,perl=T) ~ "agm_hb")

myDesign <- data.frame(
  sampleID = anot$ID,
  condition = anot$cond_MLaL
)

gtf<-import("resources/gencode.vM37.annotation.gtf.gz")
gtf<-gtf[which(!is.na(match(gtf$transcript_id,stringi::stri_split_fixed(rownames(counts),"|",2,simplify=T)[,1]))),]
export(gtf,"resources/gencode.vM37.annotation.sub_HP.gtf.gz")

counts<-counts[which(!is.na(match(stringi::stri_split_fixed(rownames(counts),"|",2,simplify=T)[,1],gtf$transcript_id))),]

aSwitchList <- importRdata(
  isoformCountMatrix   = counts,
  isoformRepExpression = ndat,
  designMatrix         = myDesign,
  isoformExonAnnoation = "resources/gencode.vM37.annotation.sub_HP.gtf.gz",
  isoformNtFasta       = "resources/gencode.vM37.transcripts.fa.gz",
  showProgress = TRUE,
  removeNonConvensionalChr = TRUE
)

SwitchListFiltered <- preFilter(
  switchAnalyzeRlist = aSwitchList,
  geneExpressionCutoff = 0,
  isoformExpressionCutoff = 0,
  removeSingleIsoformGenes = TRUE
)
 
SwitchListAnalyzed <- isoformSwitchTestSatuRn(
  switchAnalyzeRlist = SwitchListFiltered,
  reduceToSwitchingGenes=TRUE,
  reduceFurtherToGenesWithConsequencePotential = FALSE,
  alpha = 0.05,
  dIFcutoff = 0.1,
  onlySigIsoforms = FALSE
)

SwitchListFiltered <- analyzeORF(SwitchListFiltered,
                                 orfMethod='longest')


extractSwitchSummary(SwitchListAnalyzed)

SwitchListAnalyzed$isoformSwitchAnalysis$symbol<-SwitchListAnalyzed$isoformFeatures$gene_id[match(SwitchListAnalyzed$isoformSwitchAnalysis$isoform_id,SwitchListAnalyzed$isoformFeatures$isoform_id)]
SwitchListAnalyzed$isoformSwitchAnalysis$dIF<-SwitchListAnalyzed$isoformFeatures$dIF[match(SwitchListAnalyzed$isoformSwitchAnalysis$iso_ref,SwitchListAnalyzed$isoformFeatures$iso_ref)]

conds<-unique(SwitchListAnalyzed$conditions$condition)
comps<-expand.grid(conds,conds)
comps<-comps[which(comps$Var1!=comps$Var2),]
comps<-as.data.frame(t(apply(comps,1,sort)))
comps<-dplyr::distinct(comps)
colnames(comps)<-c("cond1","cond2")

calc_diffs<-lapply(1:nrow(comps),function(x){
  cond1<-comps$cond1[x]
  cond2<-comps$cond2[x]
  hold<-SwitchListAnalyzed$isoformSwitchAnalysis[which((SwitchListAnalyzed$isoformSwitchAnalysis$condition_1==cond1 & SwitchListAnalyzed$isoformSwitchAnalysis$condition_2==cond2) | (SwitchListAnalyzed$isoformSwitchAnalysis$condition_1==cond2 & SwitchListAnalyzed$isoformSwitchAnalysis$condition_2==cond1)),]
  hold$sig<-hold$padj<=0.05
  return(hold)
})
names(calc_diffs)<-apply(comps,1,function(x){paste(x[1],"vs",x[2],sep="_")})
saveRDS(calc_diffs,"IsoformSwitch_diffs.rds")

dir.create("IsoformSwitch_QC_plots",showWarnings = F)
dir.create("IsoformSwitch_res_tables",showWarnings = F)
for(x in 1:length(calc_diffs)) {
  png(paste("IsoformSwitch_QC_plots",paste0(names(calc_diffs)[x],"_QC_plot_hist.png"),sep="/"),res=300,width=1800,height=1600)
  print(ggplot(calc_diffs[[x]]) + geom_histogram(aes(estimates,fill=sig),bins = 100,alpha=0.7,col="white") + scale_y_log10() +
          xlim(c(-max(abs(SwitchListAnalyzed$isoformSwitchAnalysis$estimates)),max(abs(SwitchListAnalyzed$isoformSwitchAnalysis$estimates)))) +
          theme_classic() + coord_cartesian(expand=F) + theme(legend.position="top") + scale_fill_manual(values=c("grey88","grey9")) + ggtitle(gsub("_"," ",names(calc_diffs)[x])))
  dev.off()
  
  png(paste("IsoformSwitch_QC_plots",paste0(names(calc_diffs)[x],"_QC_plot_dIF_hist.png"),sep="/"),res=300,width=1800,height=1600)
  print(ggplot(calc_diffs[[x]]) + geom_histogram(aes(dIF,fill=sig),bins = 100,alpha=0.7,col="white") + scale_y_log10() +
          theme_classic() + coord_cartesian(expand=F) + theme(legend.position="top") + scale_fill_manual(values=c("grey88","grey9")) + ggtitle(gsub("_"," ",names(calc_diffs)[x])))
  dev.off()
  
  write.table(calc_diffs[[x]],paste("IsoformSwitch_res_tables",paste0(names(calc_diffs)[x],"_res.txt"),sep="/"),sep="\t",row.names=F)

  png(paste("IsoformSwitch_QC_plots",paste0(names(calc_diffs)[x],"_QC_plot_dIF.png"),sep="/"),res=300,width=1800,height=1600)
  print(ggplot(calc_diffs[[x]]) + geom_point(aes(dIF,estimates,col=sig),alpha=0.7) +
          ylim(c(-max(abs(SwitchListAnalyzed$isoformSwitchAnalysis$estimates)),max(abs(SwitchListAnalyzed$isoformSwitchAnalysis$estimates)))) +
          theme_classic() + coord_cartesian(expand=F) + theme(legend.position="top") + scale_color_manual(values=c("grey88","grey9")) + ggtitle(gsub("_"," ",names(calc_diffs)[x])))
  dev.off()
  
  png(paste("IsoformSwitch_QC_plots",paste0(names(calc_diffs)[x],"_QC_plot_dIF_padj.png"),sep="/"),res=300,width=1800,height=1600)
  print(ggplot(calc_diffs[[x]]) + geom_point(aes(dIF,-log10(padj),col=sig),alpha=0.7) +
          theme_classic() + coord_cartesian(expand=F) + theme(legend.position="top") + scale_color_manual(values=c("grey88","grey9")) + ggtitle(gsub("_"," ",names(calc_diffs)[x])))
  dev.off()
  
  calc_diffs[[x]]$sig<-calc_diffs[[x]]$padj<=0.05 & abs(calc_diffs[[x]]$dIF)>0.1
  png(paste("IsoformSwitch_QC_plots",paste0(names(calc_diffs)[x],"_QC_plot_dIF_padj_guides.png"),sep="/"),res=300,width=1800,height=1600)
  print(ggplot(calc_diffs[[x]]) + geom_point(aes(dIF,-log10(padj),col=sig),alpha=0.7) +
          geom_hline(yintercept = -log10(0.05),lty="dashed",col="red") +
          geom_vline(xintercept = -0.1,lty="dashed",col="red") + geom_vline(xintercept = 0.1,lty="dashed",col="red") +
          theme_classic() + coord_cartesian(expand=F) + theme(legend.position="top") + scale_color_manual(values=c("grey88","grey9")) + ggtitle(gsub("_"," ",names(calc_diffs)[x])))
  dev.off()
}

extractSwitchSummary(SwitchListAnalyzed)
write.csv(SwitchListAnalyzed$isoformFeatures, "isoformfeatures_2.2_est_sva_MLAL.tr.csv")
saveRDS(SwitchListAnalyzed,"isoformfeatures_2.2_est_sva_MLAL.tr.rds")
save.image("IsoformSwitchAnalyzeR_analysis_est_MLAL.RData")

pdf(file = 'Sox4_switchPlot.pdf', onefile = FALSE, height=6, width = 9)
switchPlot(SwitchListAnalyzed,
           gene='Sox4',
           condition1 = 'HE_agm',
           condition2 = 'HE_emp',
           localTheme = theme_bw(base_size = 13))
dev.off()
