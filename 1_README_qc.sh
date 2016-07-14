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
#Stage 2.  rsync or copy the genotype files, plus other files your analysis will depend on
########################################################################
echo "Copying the genotype/phenotype files......"
for EXT in bed bim fam
do 
    cp /GWD/bioinfo/projects/statgen2/HMaxiom/genotypes/HM_imp.$EXT genotypes/ 
done

cp -r /GWD/bioinfo/projects/statgen2/HMaxiom/phenotypes .

########################################################################
#Stage 3. Convert strand and check alignment (of marker names, positions, and strand) with 1000G reference
########################################################################
echo -e "\nChecking strand and alignment......"
#GSKBB2
#GXAPP=$GXAPP  PLINK=$PLINK TOPSOURCE=genotypes/PGx7656 REF=$HG19_RESOLUTION PLINKFILTERS="--maf 0.000001"  make -e -j --file=$ALIGHN aligned-hg19

#the following 3 lines are specific to GSKBB1
cp -r /GWD/appbase/projects/statgen/RD-MDD-GX_PUBLIC/Affymetrix_Product_Files/Axiom_GSKBB1/hg19_b37_na34/hg19_resolution .
cp hg19_resolution/PLINK_extract.list.orig.no.multi-allelic hg19_resolution/PLINK_extract.list
GXAPP=$GXAPP PLINK=$PLINK TOPSOURCE=genotypes/HM_imp REF=hg19_resolution PLINKFILTERS="--maf 0.000001"  make -e -j --file=$ALIGHN aligned-hg19

########################################################################
#Stage 4.  Check GWAS genotype data quality
#Need to review QC-makefile results- kinship, PCA
#autosome data clean3
#chr x data: non-PAR region (clean3-chr23-F clean3-chr23-M) and PAR region (clean3-chr23-par, empty for GSKBB2)
########################################################################
echo -e "\nChecking GWAS genotype data quality......"
STUDY=aligned-hg19/alignment3 HGVERSION=19 SUBGROUPS=EUR make -e -j --file=$QC

