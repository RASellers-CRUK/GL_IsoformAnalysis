#!/bin/bash -l
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --mem=72G
#SBATCH -t 48:00:00
#SBATCH --mail-type=END
#SBATCH --mail-user=robert.sellers@manchester.ac.uk

ml apps/fastqc
ml apps/minimap2
ml apps/bbmap
ml apps/samtools
ml apps/qualimap
ml dbdata/gatk/hg38

FASTQ=$(sed -n "${SLURM_ARRAY_TASK_ID} p" ${SLURM_SUBMIT_DIR}/${1})
REF=${2}
SAMP=$(basename ${FASTQ} | cut -f1 -d.)
ODIR=${SLURM_SUBMIT_DIR}/${3}/${SAMP}
mkdir -p ${ODIR}/${SAMP}_bamqc

bbduk.sh literal=bbduk.sh trimq=5 in=${FASTQ} out=${ODIR}/${SAMP}_trimmed.fq literal=AAGCAGTGGTATCAACGCAGAGTACTTTTTTTTTT stats=${ODIR}/${SAMP}_trim.stats ksplit=f tossjunk=t trimpolya=0 minlength=0

conda activate yacrd
minimap2 -t ${SLURM_NTASKS} -x ava-ont -g 500 ${ODIR}/${SAMP}_trimmed.fq ${ODIR}/${SAMP}_trimmed.fq > ${ODIR}/${SAMP}.overlap.paf
yacrd -i ${ODIR}/${SAMP}.overlap.paf -o ${ODIR}/${SAMP}.report.yacrd -c 4 -n 0.4 filter -i ${ODIR}/${SAMP}_trimmed.fq -o ${ODIR}/${SAMP}_trimmed_filtered.fq
rm -rf ${ODIR}/${SAMP}.overlap.paf

porechop-runner.py -i ${ODIR}/${SAMP}_trimmed_filtered.fq -o ${ODIR}/${SAMP}_trimmed_filtered_porechop.fq -v 2 -t ${SLURM_NTASKS} --discard_middle
fastqc -o ${ODIR} ${ODIR}/${SAMP}_trimmed_filtered_porechop.fq --memory 10000

minimap2 -t ${SLURM_NTASKS} -ax map-ont -p 0 -N 10 ${REF} ${ODIR}/${SAMP}_trimmed_filtered_porechop.fq | samtools view -@ ${SLURM_NTASKS} -O BAM -o ${ODIR}/${SAMP}.bam -
samtools sort -@ ${SLURM_NTASKS} -O BAM -o ${ODIR}/${SAMP}.sort.bam ${ODIR}/${SAMP}.bam
qualimap bamqc -bam ${ODIR}/${SAMP}.sort.bam --java-mem-size=${SLURM_MEM_PER_NODE}M -outdir ${ODIR}/${SAMP}_bamqc

conda activate NanoCount
NanoCount -i ${ODIR}/${SAMP}.sort.bam -o ${ODIR}/${SAMP}_abundancePlus.tsv -p align_score -a --extra_tx_info -n

exit 0
