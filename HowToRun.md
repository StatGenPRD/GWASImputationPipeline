Note, this guidance assumes a standardized extraction of Axiom GSKBB2 data from GSK's data bank. Although this workflow is flexible and can support data from any source, for the sake of simplicity, we are only describing the most likely scenario. Please dig into the documentation on each step of the process for details on how to adjust the workflow to handle other scenarios.

**Table of Contents**  *generated with [DocToc](http://doctoc.herokuapp.com/)*

- [Create workspace](#workspace)
- [Create shell scripts](#scripts)
- [Refer to data bank retrieval](#retrieve)
- [Determine ancestry](#ancestry)
- [Execute the pre-imputation Q](#execQC)
- [Review the QC results](#review)
- [Dispatch phase & impute jobs](#dispatch)
- [Monitoring progress](#monitor)
- [Check results](#check)
- [Re-running](#rerun)
	- [Failed jobs](#fail)
	- [Corrected data](#corrected)



## <a name="workspace">Create workspace</a>
You will need a directory on a file share in the GDC with sufficient space for the imputation output (estimate ~75 MB per sample) with a structure like this:
```
> tree /path/to/share/StudyID
`-- AnalysisReadyData
    |-- phenotypes
 ```
where ```StudyID``` is the clinical study ID for a routine instream analysis or an OPT ID for a bespoke or meta-analysis.


## <a name="scripts">Create shell scripts</a>
Although there are only a few commands required, it is best to run as a shell script to maintain a record of what was done. Specifically, it is best to have two scripts as the process naturally breaks into two steps - the pre-imputation QC and the dispatching of the phasing & imputation jobs. You can copy the template scripts here to the ```AnalysisReadyData``` directory to modify.
* [Pre-imputation QC](https://raw.githubusercontent.com/StatGenPRD/GWASImputationPipeline/master/Pre-imputation_QC.sh)
* [Phase & Impute](https://raw.githubusercontent.com/StatGenPRD/GWASImputationPipeline/master/Phase_Impute.sh)


## <a name="retrieve">Refer to data bank retrieval</a>
In the ```Pre-imputation_QC.sh``` file, set the ```RTP```, ```STUDY_ID```, and ```EMAIL``` variables like:
```
#RTP directory containing the plink binary data retrieved from the bank
RTP=/GWD/appbase/projects/RD-MDD-GX_SCRATCH/ABC123456
#root of plink files, usually the clinical study or project ID
STUDY_ID=ABC123456
#e-mail address to receive notification when complete
EMAIL=first.last@PAREXEL.com
```


## <a name="ancestry">Determine ancestry</a>
Use the self-reported ancestry from dmenv (race dataset for race and demo dataset for ethnicity) only considering subjects in the genetic dataset (not the whole study) to determine the largest homogeneous ancestry group for use in Hardy-Weinberg disequilibrium assay removal. This will most likely be non-hispanic European whites. Create a plink keep format file using the USUBJIDs as in the fam file and name as ```group.subjects``` where ```group``` is something like EUR or ASN. Place this file in the ```phenotypes``` directory and then set the ```MAJOR_POP``` variable in the ```Pre-imputation_QC.sh``` file to match the ```group``` like:
```
#root of plink keep format file in phenotypes directory containing the subjects from the largest ancestry group
MAJOR_POP=EUR
```


## <a name="execQC">Execute the pre-imputation QC</a>
Leave the ```Execution``` section of the ```Pre-imputation_QC.sh``` file intact and run using nohup as it will run a while (you should receive an e-mail notice when done unless something went wrong):
```
nohup /path/to/share/StudyID/AnalysisReadyData/Pre-imputation_QC.sh >/path/to/share/StudyID/AnalysisReadyData/Pre-imputation_QC.out 2>/path/to/share/StudyID/AnalysisReadyData/Pre-imputation_QC.err
```
**Do not** add an ```&``` to the end of the command to send it to the background as you will be prompted for your password 4 times. After the fourth time, you can do ctrl-Z to stop the job and then ```bg``` to resume it in the background.

Note, the script assumes data was retrieved from the bank and thus a ```RETRIEVE.log``` file will be in the same ```RTP``` directory. If this isn't the case and you delete/comment out the command that copies this file to the ```doc``` directory, then you should consider placing some other file describing the data provenance in the ```doc``` directory.




## <a name="review">Review the QC results</a>
After execution has finished, check the ```Pre-imputation_QC.out``` and ```.err``` files to confirm everything executed as expected. The final files produced in the ```genotypes``` directory should be (where ```STUDY_ID``` matches what was specified in the script):
* a plink data set ```STUDY_ID-clean3*```
* a matrix of eigenvectors - the first two of which can be plotted in Excel to confirm cluster separate as expected (optionally, you can color these with the self-reported ancestry) ```STUDY_ID-pca.eigenvec```
* pair-wise kinship values from KING ```STUDY_ID-kinship.kin0.gz```.  Since a similar check is done as part of the pre-bank lab QC, you should not expect many (or any) related samples and can check with a command like this:
```
zcat HZC113108-kinship.kin0.gz | awk '{if ($8 > 0.0884) print $0}'
```
Looking at the 8th column (kinship) for values > 0.0884  [see the wiki page on KING for more details](https://connect.gsk.com/sites/genetics/GeneticsWIKI/Wiki%20Pages/Software%20-%20KING.aspx). Remove all duplicates as potentially mis-identified subjects. Remove the lowest call rate of 1st or 2nd degree relative pairs per the ```STUDY_ID-clean2-missing.imiss``` file. Create a plink format remove file named ```remove_related.sub``` and place in the ```genotypes``` directory using USUBJID as in the fam file. If no subjects to remove, create an empty file. If  subjects to remove, create a tab-delimited ```KIN_EXCLUDE.log``` file in the ```docs``` directory like:
```
#USUBJID	Reason
#ABC123456.0002227	Lower call rate in 1st or 2nd degree relative pair
#ABC123456.0002223	Unexpected duplicate
#ABC123456.0002229	Unexpected duplicate
```


## <a name="dispatch">Dispatch phase & impute jobs</a>
In the ```Phase_Impute.sh``` file, set the ```STUDY_ID``` and ```EMAIL``` variables to match that above:
```
#root of plink files, usually the clinical study or project ID
STUDY_ID=ABC123456
#e-mail address to receive notification when complete
EMAIL=first.last@PAREXEL.com
```
Then run using nohup:
```
nohup /path/to/share/StudyID/AnalysisReadyData/Phase_Impute.sh >/path/to/share/StudyID/AnalysisReadyData/Phase_Impute.out 2>/path/to/share/StudyID/AnalysisReadyData/Phase_Impute.err &
```
The e-mail notifications will come at the end of each phasing job (one per autosome) and will include exit status (any non-zero value indicates an error). When you start to receive these, you can check the ```Phase_Impute.out``` and ```.err``` files to ensure everything ran ok. You should expect the following files in the ```aligned-hg19``` directory:
* a plink dataset containing all genotypes going into phasing ```alignment3*```
* a strand/allele alignment report ```report-match.csv```. In this report, fail-filter should be **very** low, non-ACGT and no-match should be < 10k and the vast majority should be in match-plus.
* an allele frequency comparison plot ```report-freq.png``` which should show a strong diagonal.
If anything is wrong with these results, see the [Corrected data](#corrected) section betlow for how to kill the phasing and restart after correcting the data.


## <a name="monitor">Monitoring progress</a>
You will know the imputation is done when there are 406 ```*.done``` files in the ```imputed-20120314``` directory which you can check as follows
```
ls *.done | wc -l
```
Until this happens, it is a good idea to check the SGE queue periodically to confirm that your jobs are running at the rate expected (ideally 100 jobs at a time). Use the following command:
```
qstat
```
State of r indicates running, q or T indicates queued. E indicates something is wrong, check with an administrator.

If it appears too many of your jobs are queued and not enough running, check to see how busy the queue is:
```
qstat -u "*"
```
If there aren't a lot of running jobs from other users, then the queues may be in error state and require resetting, check
```
> qhost -q | grep 'dl580'
   dl580                BIP   0/0/16
   dl580                BIP   0/0/16
   dl580                BIP   0/0/16
   dl580                BIP   0/0/18
   dl580                BIP   0/0/18
   dl580                BIP   0/0/18
   dl580                BIP   0/0/18
   dl580                BIP   0/0/18
   dl580                BIP   0/0/12
   dl580                BIP   0/0/12
```
Where the lack of E afer the N/N/N indicates the queues are fine. If several are in E state, e-mail R&D_IT_Infra_Services@gsk.com and ask to reset. Be sure to include a screen cap of what you have observed and note that you are on the us1us0168 server.


## <a name="check">Check results</a>
To confirm all the phasing and imputation jobs ran as expected, you can check all the ```*.err``` files in the ```phased-*-hg19/jobs``` and ```imputed-20120314/jobs``` directories which should be empty
```
cat phased-*-hg19/jobs/*.err imputed-20120314/jobs/*.err
```
And you can check the ```*.out``` files in the same directories:
```
grep '^Analysis took' phased-*-hg19/jobs/*.out | wc -l
grep '^Run completed' imputed-20120314/jobs/*.out | wc -l
```
where the first command should return 22 and the second 406.


## <a name="rerun">Re-running</a>
Below are some tips for re-running under different scenarios:

### <a name="fail">Failed jobs</a>
If any of the phasing or imputation jobs failed to complete for any reason (e.g. server goes down), you can resume the analysis where it left off by...BE SURE TO ADDRESS LOGGING

### <a name="corrected">Corrected data</a>
If you discover an issue with the assayed genotype data after the phasing and imputation have been initiated and those jobs are still running, you can kill with this command
```
qdel -u mudID
```
where ```mudID``` is your user ID. **Note**, this will terminate all SGE jobs you have running including those not associated with this dataset (e.g. analysis pipeline or the phasing/imputation of another dataset). Consult an administrator if you need to target kill the jobs associated with a specific run of the pipeline. The jobs should die relatively quickly. Once the data has been corrected, you can re-run by...BE SURE TO ADDRESS FORCING OVERWRITE

