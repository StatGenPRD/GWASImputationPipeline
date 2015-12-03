#!/bin/bash
set -o nounset #exit if encounter an uninitialized variable
set -o errexit #exit if error

#######
#Input parameters
#######

#Directory containing makefiles
GXAPPI=/GWD/appbase/projects/GXapp/GWASImputationPipeline

#Directory containing strand resolution files
STRAND=/GWD/appbase/projects/RD-MDD-GX_PUBLIC/Affymetrix_Product_Files/Axiom_GSKBB2/hg19_b37_na34/hg19_resolution

#root of plink files, usually the clinical study or project ID
STUDY_ID=ABC123456
#e-mail address to receive notification when complete
EMAIL=first.last@PAREXEL.com


#######
#Execution
#######

#Remove any related subjects
plink --bfile genotypes/$STUDY_ID-clean3 --remove genotypes/remove_related.sub --make-bed --out genotypes/$STUDY_ID-clean4

#"Align" data to 1000 Genomes reference haplotypes
TOPSOURCE=genotypes/$STUDY_ID-clean4 REF=$STRAND make -e -j --file=$GXAPPI/imputation-Makefile aligned-hg19

#Submit phasing jobs to grid with e-mail notification
$GXAPPI/phase-dispatch.py --genotypes=aligned-hg19/alignment3 --dispatch --submitopt=--mail_to=$EMAIL --submitqueue=dl580@us1us049,dl580@us1us050,dl580@us1us0166,dl580@us1us0167,dl580@us1us0168,dl580@us1us0169,dl580@us1us0170 

#Sleep until phasing complete then submit imputation jobs to grid
$GXAPPI/chunk-dispatch.py --sleep --phased=phased-mach-hg19 --dispatch --submitopt=--memory=4g --submitqueue=dl580@us1us049,dl580@us1us050,dl580@us1us0166,dl580@us1us0167,dl580@us1us0168,dl580@us1us0169,dl580@us1us0170

