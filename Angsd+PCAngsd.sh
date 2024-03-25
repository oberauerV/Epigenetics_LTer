#This assumes that the mapped files are in a folder called "bams".

#Create list of mapped bam files (normally when in the parent directory of "bams"):

ls ./bams/*.bam > list

#use a text editor to clean the list to the part of the name needed for each accession by using search and replace. Save as "indlist" for plotting later.


SCAFF_LIST=($(</lisc/scratch/botany/fabian/Lumbricus/PCA/LumTer.scaffoldlist))
SCAFF=${SCAFF_LIST[${SLURM_ARRAY_TASK_ID}]}

angsd -bam indlist -r ${SCAFF} \
  -GL 1 -doMajorMinor 1 -doMaf 1 -SNP_pval 1e-6 -minMapQ 20 -minQ 20 -minInd 13 -doGlf 2\
  -out ${SCAFF}_nominmaf -P 1

ml pcangsd
pcangsd -beagle -b ../LumTer_nominmaf.beagle.gz -e 3 -o LumTer


