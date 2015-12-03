The files here are used to run an automated workflow in GSK's global data center (GDC) making use of the wrappered SGE (abq) to:

1. [QC and reformat the directly assayed array genotypes producing reports for review.](https://github.com/StatGenPRD/GWASImputationPipeline/blob/master/pre-imputation_QC.md)
2. [Phase results from (1) using either MACH if < 1,000 samples or hapi-ur if > 1,000 samples.](https://github.com/StatGenPRD/GWASImputationPipeline/blob/master/phasing.md)
3. [Impute using the 1000 Genomes Phase 1 reference haplotypes into the results from (2) with minimac.](https://github.com/StatGenPRD/GWASImputationPipeline/blob/master/imputation.md) 

Click on each of the above for details or click here for guidance on [how to run](https://github.com/StatGenPRD/GWASImputationPipeline/blob/master/HowToRun.md).
