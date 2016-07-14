#!/bin/bash

cd /GWD/bioinfo/projects/statgen/GXapp/GWASImputationPipelineUpdates/HMaxiom_1KGP3

########################################################################
#Stage 1. Starting a README.sh and making a working directories
########################################################################
set -o nounset
set -o errexit
mkdir -p genotypes phenotypes doc
GXAPP=/GWD/bioinfo/projects/statgen/GXapp/imputation/GWASImputationPipeline-master
ALIGHN=$GXAPP/imputation/imputation-Makefile
HG19_RESOLUTION=/GWD/appbase/projects/statgen/RD-MDD-GX_PUBLIC/Affymetrix_Product_Files/Axiom_GSKBB2/hg19_b37_na34/hg19_resolution
QC=$GXAPP/imputation/gwas-qc-Makefile
PHASE=$GXAPP/imputation/phase-dispatch.py
HAPI2MACH=$GXAPP/imputation/convert-hapiur-mach.py
IMPUTE=$GXAPP/imputation/chunk-dispatch.py
PLINK=/GWD/appbase/projects/statgen/GXapp/plink/plink-1.07/plink
PLINK2DOSE=$GXAPP/imputation/plink2dose.r

SERVERS=dl580@us1us049,dl580@us1us050,dl580@us1us0166,dl580@us1us0167,dl580@us1us0168,dl580@us1us0169,dl580@us1us0170

########################################################################  
#Stage 5. phase (need direction for phasing of separate chromosomes and/or regions)
#Confirm success/completion of the phasing process. If the .err files are not empty or the .done files are missing, review the chr*.out files in the phased/jobs directory.
#option --chr=22 for specific chromosome
########################################################################
echo -e "\nPhasing......"
$PHASE --genotypes=aligned-hg19/alignment3-clean3 --hapi2mach-path=$HAPI2MACH --submitmem=0.2 --submitopt=--mail_me --dispatch --algorithm=auto

for BASE in  "par" "F" "M" 
do
   if [ -f aligned-hg19/alignment3-clean3-chr23-${BASE}.bim  ]; then 
       $PHASE --genotypes=aligned-hg19/alignment3-clean3-chr23-${BASE} --hapi2mach-path=$HAPI2MACH --chr=23 --name=chr23-${BASE} --submitmem=0.2 --submitopt=--mail_me --dispatch --algorithm=auto
   fi
done

