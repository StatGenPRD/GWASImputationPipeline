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
#Stage 7.  Genotyped variant not in imputation
########################################################################  
echo -e "Obtaining genotyped vairants not in imputation......"

awk '{print $2 " " $1 ":" $4}' genotypes/hg19plus.bim >genotypes/update-map-name
$PLINK --noweb  --bfile genotypes/hg19plus --update-map genotypes/update-map-name --update-name --make-bed --out  genotypes/hg19plus2 >/dev/null
STUDY=genotypes/hg19plus2 HGVERSION=19 SUBGROUPS=EUR make -e -j --file=$QC
awk '{print $2}' aligned-hg19/alignment3-clean3.bim aligned-hg19/alignment3-clean3-chr23.bim |awk 'NF>0' >genotypes/hg19plus2.exclude
$PLINK --noweb  --bfile genotypes/hg19plus2-clean3 --exclude genotypes/hg19plus2.exclude --make-bed --out  genotypes/no-imputation >/dev/null
$PLINK --noweb  --bfile genotypes/no-imputation --recodeA --out  genotypes/no-imputation >/dev/null
/GWD/bioinfo/tools/bin/R --vanilla --slave --args genotypes/no-imputation imputed-20130502 < $PLINK2DOSE 
