#!/bin/bash
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --mem=32G
#SBATCH -t 48:00:00
#SBATCH --mail-type=END
#SBATCH --mail-user=robert.sellers@manchester.ac.uk

ml apps/salmon/1.10.2

SAMP=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${SLURM_SUBMIT_DIR}/${1})
FQLIST=${SLURM_SUBMIT_DIR}/${2}
R1=$(grep ${SAMP} ${FQLIST} | grep _R1_001)
R2=$(grep ${SAMP} ${FQLIST} | grep _R2_001)
LANE=$(basename ${R1} | cut -f1 -d. | sed "s/_R[12]_001$//")
IND=${3}
ODIR=${SLURM_SUBMIT_DIR}/${4}/salmon/${SAMP}
mkdir -p ${ODIR}

echo ${R1} | tr " " "\n" | xargs -I % zcat % | pigz -c > ${ODIR}/${SAMP}_R1.fq.gz
echo ${R2} | tr " " "\n" | xargs -I % zcat % | pigz -c > ${ODIR}/${SAMP}_R2.fq.gz

salmon quant -l OU -i ${IND} -1 ${ODIR}/${SAMP}_R1.fq.gz -2 ${ODIR}/${SAMP}_R2.fq.gz \
	-o ${ODIR} -p ${SLURM_NTASKS} --validateMappings

exit 0
