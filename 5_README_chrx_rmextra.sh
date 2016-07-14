#!/bin/bash

IMPDIR=/GWD/bioinfo/projects/statgen/GXapp/GWASImputationPipelineUpdates/ChrX/HMaxiom_1KGP3

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
#Stage 7.  Remove extra data in chrX imputation files
#chrX dosage files have an extra 3rd column for all chunks except chunk 1
#chrX info and dosage files contains in flanking region
#these extra columns need to be removed manually
########################################################################  
echo -e "\nRemoving extra data in chrX imputation files......"

##-----------------------------------------------------------------------------
## move all chr23 file to subfolder dose_orig
##-----------------------------------------------------------------------------
cd $IMPDIR/imputed-20130502
if [ ! -d chr23_orig ]; then mkdir chr23_orig; fi
mv chr23*chunk[[:digit:]]* chr23_orig/

##-----------------------------------------------------------------------------
## tmp folder to hold dosage files with extra 3rd column removed
## add dose file for chunk 1 at the end
##-----------------------------------------------------------------------------
if [ ! -d chr23_tmp ]; then mkdir chr23_tmp; fi
cp chr23_orig/*.info.gz chr23_tmp/

##-----------------------------------------------------------------------------
## remove column 3 from chr23 chunks except chunk 1
##-----------------------------------------------------------------------------
for file in chr23_orig/*.dose.gz
do
   file=$(basename $file)
   if [[ $file == *chunk1.dose.gz ]]
   then
       echo "coping ${file}"
       cp chr23_orig/${file} chr23_tmp/
   else 
       echo "removing column 3 for ${file}"
       zcat chr23_orig/${file} | cut -f1-2,4- |  gzip > chr23_tmp/${file} &
   fi
done

##-----------------------------------------------------------------------------
## retain target region for each chunk file
## create dose.gz and info.gz under imputation directory
##-----------------------------------------------------------------------------   
for file in chr23_tmp/*.info.gz
do
      file=$(basename $file)
      file="${file%.*.*}"            
      #for each file, get region in jobs/*.sh
      start=$(cat jobs/${file}.sh |grep -i vcfstart| cut -d\   -f6)
      end=$(cat jobs/${file}.sh |grep -i vcfend | cut -d\   -f8)
      echo "python $GXAPP/imputation/extract-dose.py --input=chr23_tmp/${file}  --output=${file} --range --begin=${start} --end=${end} &"
      python $GXAPP/imputation/extract-dose.py --input=chr23_tmp/${file}  --output=${file} --range --begin=${start} --end=${end} 
done

##-----------------------------------------------------------------------------
##remove dose_tmp folder
##-----------------------------------------------------------------------------
rm -r chr23_tmp

