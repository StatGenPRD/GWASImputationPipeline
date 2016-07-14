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
#Stage 6.  Impute
#option: --submitopt=--memory=4g for large sample sets
#option: --chr=22 for specific chromosome
#option: --minimacopt="--rounds 5 --states 200 --probs" for imputation prob files
########################################################################  
echo -e "\nImputing......"

if [ -d phased-hapiur-hg19 ]; then 
  $IMPUTE --submitqueue="${SERVERS}" --phased=phased-hapiur-hg19  --submitopt="--memory=4g" --dispatch --sleep
fi
if [ -d phased-mach-hg19 ]; then 
  $IMPUTE --submitqueue="${SERVERS}" --phased=phased-mach-hg19  --submitopt="--memory=4g" --dispatch --sleep
fi
for BASE in  "par" "F" "M" 
do
   if [ -d phased-hapiur-hg19 ] && [ -f phased-hapiur-hg19/chr23-${BASE}.snps.gz ]; then 
      $IMPUTE --submitqueue="${SERVERS}" --phased=phased-hapiur-hg19 --chr=23 --name=chr23-${BASE} --submitopt="--memory=4g" --dispatch --sleep
   fi
   if [ -d phased-mach-hg19 ] && [ -f phased-mach-hg19/chr23-${BASE}.snps.gz ]; then 
      $IMPUTE --submitqueue="${SERVERS}" --phased=phased-mach-hg19 --chr=23 --name=chr23-${BASE} --submitopt="--memory=4g" --dispatch --sleep
   fi
done

