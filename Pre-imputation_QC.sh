#!/bin/bash
set -o nounset #exit if encounter an uninitialized variable
set -o errexit #exit if error

#######
#Input parameters
#######

#Directory containing makefiles
GXAPPI=/GWD/appbase/projects/GXapp/GWASImputationPipeline

#RTP directory containing the plink binary data retrieved from the bank
RTP=/GWD/appbase/projects/RD-MDD-GX_SCRATCH/ABC123456

#Study ID
STUDY_ID=ABC123456

#e-mail address to receive notification when complete
EMAIL=first.last@PAREXEL.com

#root of plink keep format file in phenotypes directory containing the subjects from the largest ancestry group
MAJOR_POP=EUR



#######
#Execution
#######
#Create sub-directories
mkdir -p genotypes phenotypes doc

#Copy plink binary data from RTP
for EXT in bed bim fam; do rsync -az us2us00013:$RTP/$STUDY_ID.$EXT genotypes/; done
#Copy log from RTP
rsync -az us2us00013:$RTP/RETRIEVE.log doc/

#Simple QC using single largest ancestry group for HW check
STUDY=genotypes/$STUDY_ID HGVERSION=19 make SUBGROUPS=$MAJOR_POP -e -j --file=$GXAPPI/gwas-qc-Makefile

#E-mail notification
echo $STUDY" pre-imputation QC is complete, please review results" | mail -s $STUDY" pre-imputation QC is complete" $EMAIL
