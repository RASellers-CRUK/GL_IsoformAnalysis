library(IsoformSwitchAnalyzeR)
library(rtracklayer)

anot<-read.table("../../HP_analysis/HP_cells.txt", sep = "\t", header = T, stringsAsFactors = F)
anot<-anot[order(anot$condition),]
files<-list.files("../salmon","quant.sf",full.names = T,recursive = T)
names(files)<-basename(dirname(files))

myDesign <- data.frame(
  sampleID = anot$sampleID,
  condition = anot$condition
)

files<-files[match(anot$sampleID,names(files))]

salmonQuant <- importIsoformExpression(
  sampleVector = files
)

gtf<-import("resources/gencode.vM37.annotation.gtf.gz")
gtf<-gtf[which(!is.na(match(gtf$transcript_id,salmonQuant$abundance$isoform_id))),]
export(gtf,"resources/gencode.vM37.annotation.sub.gtf.gz")

salmonQuant_sub<-list(abundance=salmonQuant$abundance[which(!is.na(match(salmonQuant$abundance$isoform_id,gtf$transcript_id))),],
                  counts=salmonQuant$counts[which(!is.na(match(salmonQuant$counts$isoform_id,gtf$transcript_id))),],
                  length=salmonQuant$length[which(!is.na(match(salmonQuant$length$isoform_id,gtf$transcript_id))),],
                  importOptions=salmonQuant$importOptions)

aSwitchList <- importRdata(
  isoformCountMatrix   = salmonQuant_sub$counts,
  isoformRepExpression = salmonQuant_sub$abundance,
  designMatrix         = myDesign,
  isoformExonAnnoation = "resources/gencode.vM37.annotation.sub.gtf.gz",
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

SwitchListFiltered <- analyzeORF(SwitchListFiltered,
                                 orfMethod='longest')

SwitchListAnalyzed <- isoformSwitchTestSatuRn(
  switchAnalyzeRlist = SwitchListFiltered,
  reduceToSwitchingGenes=FALSE,
  reduceFurtherToGenesWithConsequencePotential = FALSE,
  alpha = 0.05,
  dIFcutoff = 0.1,
  onlySigIsoforms = FALSE
)
SwitchListAnalyzed$isoformSwitchAnalysis$symbol<-SwitchListAnalyzed$isoformFeatures$gene_id[match(SwitchListAnalyzed$isoformSwitchAnalysis$isoform_id,SwitchListAnalyzed$isoformFeatures$isoform_id)]
SwitchListAnalyzed$isoformSwitchAnalysis$dIF<-SwitchListAnalyzed$isoformFeatures$dIF[match(SwitchListAnalyzed$isoformSwitchAnalysis$iso_ref,SwitchListAnalyzed$isoformFeatures$iso_ref)]

SwitchListAnalyzed_RS <- isoformSwitchTestSatuRn(
  switchAnalyzeRlist = SwitchListFiltered,
  reduceToSwitchingGenes=TRUE,
  reduceFurtherToGenesWithConsequencePotential = TRUE,
  alpha = 0.05,
  dIFcutoff = 0.1,
  onlySigIsoforms = FALSE
)
SwitchListAnalyzed_RS$isoformSwitchAnalysis$symbol<-SwitchListAnalyzed_RS$isoformFeatures$gene_id[match(SwitchListAnalyzed_RS$isoformSwitchAnalysis$isoform_id,SwitchListAnalyzed_RS$isoformFeatures$isoform_id)]
SwitchListAnalyzed_RS$isoformSwitchAnalysis$dIF<-SwitchListAnalyzed_RS$isoformFeatures$dIF[match(SwitchListAnalyzed_RS$isoformSwitchAnalysis$iso_ref,SwitchListAnalyzed_RS$isoformFeatures$iso_ref)]

extractSwitchSummary(SwitchListAnalyzed)
extractSwitchSummary(SwitchListAnalyzed_RS)

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
          theme_classic() + coord_cartesian(expand=F) + theme(legend.position="top") + scale_fill_manual(values=c("grey9","purple")) + ggtitle(gsub("_"," ",names(calc_diffs)[x])))
  dev.off()
  
  png(paste("IsoformSwitch_QC_plots",paste0(names(calc_diffs)[x],"_QC_plot_dIF_hist.png"),sep="/"),res=300,width=1800,height=1600)
  print(ggplot(calc_diffs[[x]]) + geom_histogram(aes(dIF,fill=sig),bins = 100,alpha=0.7,col="white") + scale_y_log10() +
          theme_classic() + coord_cartesian(expand=F) + theme(legend.position="top") + scale_fill_manual(values=c("grey88","grey9")) + ggtitle(gsub("_"," ",names(calc_diffs)[x])))
  dev.off()

  write.table(calc_diffs[[x]],paste("IsoformSwitch_res_tables",paste0(names(calc_diffs)[x],"_res.txt"),sep="/"),sep="\t",row.names=F)

  png(paste("IsoformSwitch_QC_plots",paste0(names(calc_diffs)[x],"_QC_plot_dIF.png"),sep="/"),res=300,width=1800,height=1600)
  print(ggplot(calc_diffs[[x]]) + geom_point(aes(dIF,estimates,col=sig),alpha=0.7) +
          ylim(c(-max(abs(SwitchListAnalyzed$isoformSwitchAnalysis$estimates)),max(abs(SwitchListAnalyzed$isoformSwitchAnalysis$estimates)))) +
          theme_classic() + coord_cartesian(expand=F) + theme(legend.position="top") + scale_color_manual(values=c("grey9","purple")) + ggtitle(gsub("_"," ",names(calc_diffs)[x])))
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

write.csv(SwitchListAnalyzed$isoformFeatures, "isoformfeatures_2.2_est_sva_MLAL_HP.tr.csv")
write.csv(SwitchListAnalyzed_RS$isoformFeatures, "isoformfeatures_2.2_est_sva_MLAL_HP_RS.tr.csv")
saveRDS(SwitchListAnalyzed,"isoformfeatures_2.2_est_sva_MLAL_HP.tr.rds")
saveRDS(SwitchListAnalyzed_RS,"isoformfeatures_2.2_est_sva_MLAL_HP_RS.tr.rds")
save.image("IsoformSwitchAnalyzeR_analysis_est_sva_MLAL_HP.RData")

Ysm<-SwitchListAnalyzed$isoformSwitchAnalysis[which(SwitchListAnalyzed$isoformSwitchAnalysis$condition_1=="AGM" & SwitchListAnalyzed$isoformSwitchAnalysis$condition_2=="Ysm"),]
Ysm<-Ysm[which(Ysm$padj<=0.05),]
sort(unique(Ysm$symbol))

Ysb<-SwitchListAnalyzed$isoformSwitchAnalysis[which(SwitchListAnalyzed$isoformSwitchAnalysis$condition_1=="AGM" & SwitchListAnalyzed$isoformSwitchAnalysis$condition_2=="Ysb"),]
Ysb<-Ysb[which(Ysb$padj<=0.05),]
sort(unique(Ysb$symbol))

Ysm_sym<-do.call(rbind,lapply(split(Ysm,Ysm$symbol),function(x){x[which.min(x$padj),]}))
Ysb_sym<-do.call(rbind,lapply(split(Ysb,Ysb$symbol),function(x){x[which.min(x$padj),]}))

cons<-sort(intersect(unique(Ysb$symbol),unique(Ysm$symbol)))

Ysm_cons<-Ysm_sym[which(!is.na(match(Ysm_sym$symbol,cons))),]
Ysb_cons<-Ysb_sym[which(!is.na(match(Ysb_sym$symbol,cons))),]

comp<-data.frame(sym=Ysm_cons$symbol,
                 Ysm_p=Ysm_cons$padj,
                 Ysb_p=Ysb_cons$padj,
                 Ysm_prnk=rank(-log10(Ysm_cons$padj)),
                 Ysb_prnk=rank(-log10(Ysb_cons$padj)),
                 prnksm=rank(-log10(Ysm_cons$padj)) + rank(-log10(Ysb_cons$padj)),
                 prnksm_rnk=rank(rank(-log10(Ysm_cons$padj)) + rank(-log10(Ysb_cons$padj))))

cons_iso<-unique(sort(intersect(unique(Ysb$isoform_id),unique(Ysm$isoform_id))))

Ysm_cons_iso<-Ysm[which(!is.na(match(Ysm$isoform_id,cons_iso))),]
Ysb_cons_iso<-Ysb[which(!is.na(match(Ysb$isoform_id,cons_iso))),]

comp_iso<-data.frame(sym=SwitchListAnalyzed$isoformFeatures$gene_name[match(Ysm_cons_iso$isoform_id,SwitchListAnalyzed$isoformFeatures$isoform_id)],
                     iso=Ysm_cons_iso$isoform_id,
                 Ysm_p=Ysm_cons_iso$padj,
                 Ysb_p=Ysb_cons_iso$padj,
                 Ysm_fc=Ysm_cons_iso$estimates,
                 Ysb_fc=Ysb_cons_iso$estimates,
                 Ysm_d=ifelse(Ysm_cons_iso$estimates>0,"pos","neg"),
                 Ysb_d=ifelse(Ysb_cons_iso$estimates>0,"pos","neg"),
                 Ysm_prnk=rank(-log10(Ysm_cons_iso$padj)),
                 Ysb_prnk=rank(-log10(Ysb_cons_iso$padj)),
                 prnksm=rank(-log10(Ysm_cons_iso$padj)) + rank(-log10(Ysb_cons_iso$padj)),
                 prnksm_rnk=rank(rank(-log10(Ysm_cons_iso$padj)) + rank(-log10(Ysb_cons_iso$padj))))

## Plots
pdf(file = 'Runx1_switchPlot_sva_HP_AGM_Ysm.pdf', onefile = FALSE, height=6, width = 9)
switchPlot(SwitchListAnalyzed,
           gene='Runx1',
           condition1 = 'AGM',
           condition2 = 'Ysb',
           localTheme = theme_bw(base_size = 13))
dev.off()

pdf(file = 'Runx1_switchPlot_sva_HP_AGM_Ysb.pdf', onefile = FALSE, height=6, width = 9)
switchPlot(SwitchListAnalyzed,
           gene='Runx1',
           condition1 = 'AGM',
           condition2 = 'Ysm',
           localTheme = theme_bw(base_size = 13))
dev.off()

pdf(file = 'Mecom_switchPlot_sva_HP_AGM_Ysm.pdf', onefile = FALSE, height=6, width = 9)
switchPlot(SwitchListAnalyzed,
           gene='Mecom',
           condition1 = 'AGM',
           condition2 = 'Ysb',
           localTheme = theme_bw(base_size = 13))
dev.off()

pdf(file = 'Mecom_switchPlot_sva_HP_AGM_Ysb.pdf', onefile = FALSE, height=6, width = 9)
switchPlot(SwitchListAnalyzed,
           gene='Mecom',
           condition1 = 'AGM',
           condition2 = 'Ysm',
           localTheme = theme_bw(base_size = 13))
dev.off()

dir.create("cons_switchplots_sva",showWarnings = F,recursive = T)
for(x in order(comp$prnksm_rnk,decreasing=T)[1:100]){
  gene<-comp$sym[x]
  
  pdf(file = paste0("cons_switchplots_sva/","rank_",comp$prnksm_rnk[comp$sym %in% gene],"_",gene,'_switchPlot_HP_AGM_Ysm.pdf'),
      onefile = FALSE, height=6, width = 9)
  switchPlot(SwitchListAnalyzed,
             gene=gene,
             condition1 = 'AGM',
             condition2 = 'Ysm',
             localTheme = theme_bw(base_size = 13))
  dev.off()
  
  pdf(file = paste0("cons_switchplots_sva/","rank_",comp$prnksm_rnk[comp$sym %in% gene],"_",gene,'_switchPlot_HP_AGM_Ysb.pdf'),
      onefile = FALSE, height=6, width = 9)
  switchPlot(SwitchListAnalyzed,
             gene=gene,
             condition1 = 'AGM',
             condition2 = 'Ysb',
             localTheme = theme_bw(base_size = 13))
  dev.off()
}

save.image("IsoformSwitchAnalyzeR_analysis.RData")
