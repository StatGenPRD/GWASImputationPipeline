#!/usr/bin/python
from optparse import OptionParser
## deprecated, but argparse requires python >= 2.7 does not work on GSK servers
import logging
import sys
import os
import stat
import subprocess
import string
import math
import time

###
### parse command line arguments and check
###

parser = OptionParser(description = 'usage: %prog OPTIONS')
parser.add_option('-d', '--dispatch', help = 'dispatch jobs (default is NOT to dispatch)',
                  action = 'store_true', dest = 'dispatch',
                  default = False)
parser.add_option('-r', '--redo', help = 'redo jobs already completed',
                  action = 'store_true', dest = 'redo',
                  default = False)
parser.add_option('--chr', help = 'impute specified chromosome(s) only', metavar = 'CHR',
                  action = 'append', type = 'string', dest = 'chr', default = [])
parser.add_option('-p', '--phased', help = 'directory for phased haplotypes (default phased-hapiur-hg19)', metavar = 'DIR',
                  type = 'string', dest = 'phased',
                  default = 'phased-hapiur-hg19')
parser.add_option('-i', '--imputed', help = 'directory for imputed chunks (default imputed-20130502)', metavar = 'DIR',
                  type = 'string', dest = 'imputed',
                  default = 'imputed-20130502')
parser.add_option('-n', '--name', help = 'job name (default chr[chr] or JOBNAME for chr 23)', metavar = 'JOBNAME',
                  type = 'string', dest = 'jobname',
                  default = '')
parser.add_option('--chunkMb', help = 'chunk size in Mb (default 4)', metavar = 'MB',
                  type = 'float', dest = 'chunkMb',
                  default = 4.)
parser.add_option('--windowMb', help = 'flanking window size in Mb (default 0.25)', metavar = 'MB',
                  type = 'float', dest = 'windowMb',
                  default = 0.25)
parser.add_option('--minimacopt', help = 'options for minimac (default \'--rounds 5 --states 200 --probs\')', metavar = 'OPTIONS',
                  type = 'string', dest = 'minimacopt',
                  default = '--rounds 5 --states 200 --probs')
parser.add_option('--minimac-path', help = 'full path to minimac executable', metavar = 'MINIMAC',
                  type = 'string', dest = 'minimac',
                  default = '/GWD/appbase/projects/statgen/GXapp/minimac/minimac')
parser.add_option('--sleep', help = 'sleep until phased haplotypes are done',
                  action = 'store_true', dest = 'sleep',
                  default = False)
parser.add_option('--submit', help = 'command to submit job script', metavar = 'COMMAND',
                  type = 'string', dest = 'submit',
                  default = '/GWD/bioinfo/common/scripts/abq-sub')
parser.add_option('--submitqueue', help = '--queue option argument (default dl580)', metavar = 'QUEUE',
                  type = 'string', dest = 'submitqueue',
                  default = 'dl580')
parser.add_option('--submitopt', help = 'append other option(s) for submit command', metavar = 'OPTION',
                  action = 'append', type = 'string', dest = 'submitopt', default = [])
#recommend ['--memory=4g']
#but not hard coded as default optins because then there would be no way to override 
(options, args) = parser.parse_args()

logging.basicConfig(stream=sys.stdout, format='%(levelname)s:%(message)s', level=logging.DEBUG)
logging.info('minimac chunk dispatcher v0.3 Toby.x.Johnson@gsk.com')

if len(options.chr) == 0:
    logging.info('Target imputation of chromosomes 1-22 (default)')
    for chridx in range(22): # C style range 0..21
        options.chr.append(str(chridx + 1))
else:
    logging.info('Using --chr arguments, target impution of chromosomes [ ' + ','.join(options.chr) + ' ]')

if options.chunkMb <= 0:
    logging.error('Invalid --chunkMb argument [ ' + options.chunkMb + ' ], should be positive')

try:
    minimac_chunk_size = int(options.chunkMb * 1e6)
except ValueError:
    logging.error('Could not parse --chunkMb argument [ ' + options.chunkMb + ' ]')
    sys.exit(1)

if options.windowMb <= 0:
    logging.error('Invalid --window Mb argument [ ' + options.windowMb + ' ], should be positive')
    
try:
    minimac_window = str(int(options.windowMb * 1e6)) # string
except ValueError:
    logging.error('Could not parse --windowMb argument [ ' + options.windowMb + ' ]')
    sys.exit(1)

if not os.path.isfile(options.minimac):
    logging.error('Could not find minimac executable [ ' + options.minimac + ' ]')
    sys.exit(1)

if not os.path.isfile(options.submit):
    logging.error('Could not find submit executable [ ' + options.submit + ' ]')
    sys.exit(1)

if not options.submitqueue in ['', 'any', 'default']:
    options.submitopt.append('--queue=' + options.submitqueue)
logging.info('Submit command is [ ' + options.submit + ' ' + ' '.join(options.submitopt) + ' ]')

###
### Hard coded options for nosingletons 1000G full panel
### Hard coded paths for us1us
###

vcfpath = '/GWD/appbase/projects/statgen/RD-MDD-GX_PUBLIC/1KG/share.sph.umich.edu/1000genomes/fullProject/2013.05.02'
vcfpre = 'reduced.ALL.chr'
vcfpost = '.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz'

###
### Get chromosome begin and end positions from data on last variant in VCF files
###

logging.info('Reading chromosome begin and end positions')
        
try:
    chrom_begin_fh = open(os.path.join(vcfpath, 'chrbegin.data'), 'r')
except IOError:
    logging.error('Cannot open [ ' + os.path.join(vcfpath, 'chrbegin.data') + ' ]')
    sys.exit(2)

chrom_begin = dict()
for chrom_read in chrom_begin_fh.readlines():
    chrom_read = chrom_read.strip()
    if len(chrom_read) == 0 or chrom_read[0] == '#':
        continue
    chrom_data = chrom_read.split(None, 2)
    if len(chrom_data) >= 2:
        chrom_begin[chrom_data[0]] = int(chrom_data[1])
    else:
        logging.error('Error reading in [ ' + os.path.join(vcfpath, 'chrbegin.data') + ' ]')

chrom_begin_fh.close()
logging.info('Read begin positions for chromosomes [ ' + ','.join(sorted(chrom_begin.keys(), key=int)) + ' ]')

try:
    chrom_end_fh = open(os.path.join(vcfpath, 'chrend.data'), 'r')
except IOError:
    logging.error('Cannot open [ ' + os.path.join(vcfpath, 'chrend.data') + ' ]')
    sys.exit(2)

chrom_end = dict()
for chrom_read in chrom_end_fh.readlines():
    chrom_read = chrom_read.strip()
    if len(chrom_read) == 0 or chrom_read[0] == '#':
        continue
    chrom_data = chrom_read.split(None, 2)
    if len(chrom_data) >= 2:
        chrom_end[chrom_data[0]] = int(chrom_data[1])
    else:
        logging.error('Error reading in [ ' + os.path.join(vcfpath, 'chrend.data') + ' ]')

chrom_end_fh.close()
logging.info('Read end positions for chromosomes [ ' + ','.join(sorted(chrom_end.keys(), key=int)) + ' ]')

###
### Remove chromosomes without begin/end positions from options.chr
###
chrom_set_tmp = list(options.chr) # temporary copy list, as cannot remove elements from a list we are iterating over
for chrom in chrom_set_tmp:
    if not chrom in chrom_begin:
        options.chr.remove(chrom)
        logging.info('removed')
    if not chrom in chrom_end:
        options.chr.remove(chrom)
        logging.info('removed')

## make more informative logging... XXX
logging.info('Target imputation of chromosomes, with begin and end positions [ ' + ','.join(options.chr) + ' ]') 

###
### Check directory for phased haplotypes
###

if os.path.isabs(options.phased):
    phasedir = options.phased
else:
    phasedir = os.path.abspath(options.phased)
    logging.info('Inferred absolute path to phased haplotypes [ ' + phasedir + ' ]')
if not os.path.isdir(phasedir):
    logging.error('Directory [ ' + phasedir + ' ] does not exist')
    sys.exit(2)

###
### Check existence/make directories for imputation and job scripts
###

if os.path.isabs(options.imputed):
    imputedir = options.imputed
else:
    imputedir = os.path.abspath(options.imputed)
    logging.info('Inferred absolute path for output imputation [ ' + imputedir + ' ]')

jobdir = os.path.join(imputedir, 'jobs')
for thisdir in [imputedir, jobdir]:
    if not os.path.isdir(thisdir):
        try:
            os.mkdir(thisdir)
        except OSError:
            logging.error('Directory [ ' + thisdir + ' ] does not exist and cannot be created')
            sys.exit(2)

logging.info('Setting up jobs to output results in [ ' + imputedir + ' ]')
logging.info('Writing job scripts in [ ' + jobdir + ' ]')


## sleeping loop
chrom_dispatched = list()
jobs_done = list()
jobs_wrote = list()
jobs_disp = list()
while True:
    ## list of chromosomes to sleep to wait for
    chrom_todispatch = list(options.chr)
    for chrom in chrom_dispatched:
        chrom_todispatch.remove(chrom)
    logging.info('Chromosomes already dispatched [ ' + ','.join(chrom_dispatched) + ' ]')
    logging.info('Chromosomes to dispatch [ ' + ','.join(chrom_todispatch) + ' ]')
    chrom_waiting = list()
    for chrom in chrom_todispatch:    	
    	## chrX: modified by Li Li
    	if chrom == '23':
            xregion = options.jobname.split('-')[-1]
            phased = 'chr' + chrom + '-' + xregion
            ref_file = os.path.join(vcfpath, vcfpre + 'X.no.auto' + vcfpost)
            if 'par' in xregion:
            	ref_file = os.path.join(vcfpath, vcfpre + 'X.auto' + vcfpost)
        else:
        	phased = 'chr' + chrom 
        	ref_file = os.path.join(vcfpath, vcfpre + chrom + vcfpost)
        	
        ### Check whether phased haplotypes exist
        if not os.path.isfile(os.path.join(phasedir, phased + '.done')):
            chrom_waiting.append(chrom)
            continue # for chrom

        if (not os.path.isfile(os.path.join(phasedir, phased + '.gz'))
            or not os.path.isfile(os.path.join(phasedir, phased + '.snps.gz'))):
            logging.error('.done files exists but phased chromosome .gz and .snps.gz do not')
            ## print file names XXX
            sys.exit(2)

        logging.info('Found phased haplotypes for chr' + chrom + ' in [ ' + phasedir + ' ]')

        ###
        ### Calculate chunk start and ends, write job scripts and dispatch
        ###
        num_chunk = int(math.ceil(float(chrom_end[chrom] - chrom_begin[chrom] + 1) / minimac_chunk_size))
        logging.info(phased + ':' + str(chrom_begin[chrom]) + '-' + str(chrom_end[chrom]) + ' to be split into ' + str(num_chunk) + ' chunks')
        for chunk in range(num_chunk):
            chunk_id = str(chunk + 1)
            chunk_start = str(chunk*minimac_chunk_size + chrom_begin[chrom])
            chunk_end = str((chunk + 1)*minimac_chunk_size + chrom_begin[chrom] - 1)            
            job_name = phased + 'chunk' + chunk_id
           
            if not options.redo and os.path.isfile(os.path.join(imputedir, job_name + '.done')):
                jobs_done.append(job_name)
            else:
                ## XXX should we overwrite existing job .sh files?
                try:
                    job_fh = open(os.path.join(jobdir, job_name + '.sh'), 'w')
                except IOError:
                    logging.error('Could not write job script [ ' + os.path.join(jobdir, job_name + '.sh') + ' ]')
                    sys.exit(2)

                job_fh.write('#!/bin/bash\n')
                job_fh.write('/bin/uname -a\n')
                job_fh.write('/bin/echo \'Chunk imputation starting\' $(/bin/date)\n')
                job_fh.write('/bin/rm -f ' + os.path.join(imputedir, job_name + '.done') + '\n')
                job_fh.write(options.minimac + ' --vcfReference --refHaps ' + ref_file \
                             + ' --vcfstart ' + chunk_start + ' --vcfend ' + chunk_end + ' --vcfwindow ' + minimac_window \
                             + ' --haps ' + os.path.join(phasedir, phased + '.gz') \
                             + ' --snps ' + os.path.join(phasedir, phased + '.snps.gz') \
                             + ' ' + options.minimacopt \
                             + ' --prefix ' + os.path.join(imputedir, job_name) \
                             + ' --gzip' + ' \\\n')
                job_fh.write('&& /bin/touch ' + os.path.join(imputedir, job_name + '.done') + '\n')
                job_fh.write('/bin/echo \'Chunk imputation finished\' $(/bin/date)\n')
                job_fh.close()
                os.chmod(os.path.join(jobdir, job_name + '.sh'), stat.S_IRUSR | stat.S_IRGRP | stat.S_IROTH | stat.S_IWUSR | stat.S_IXUSR )
                jobs_wrote.append(job_name)
            
                if options.dispatch:
                    subprocess.call([options.submit] + options.submitopt + \
                                    ['--outfile=' + os.path.join(jobdir, job_name + '.out'), \
                                     '--errfile=' + os.path.join(jobdir, job_name + '.err'), \
                                     os.path.join(jobdir, job_name + '.sh')])
                    jobs_disp.append(job_name)

        # for chunk end
        chrom_dispatched.append(chrom)
        logging.info('All chunks dispatched for ' + phased)
    # for chrom end

    if not options.sleep:
        logging.info('Chromosomes [ ' + ','.join(chrom_waiting) + ' ] not yet phased, nothing done')
        break
    if len(chrom_waiting) > 0:
        logging.info('Chromosomes [ ' + ','.join(chrom_waiting) + ' ] not yet phased, ' + time.strftime('%H:%M %d-%m-%Y') + ' now; sleeping for 1 hour')
        time.sleep(3600)
    else:
        break # break while True loop

    # Note, will loop forever at this point unless break above
# while True end
logging.info('All chromosomes targetted have been dispatched [ ' + ','.join(chrom_dispatched) + ' ]')

if not options.redo and len(jobs_done) > 0:
    logging.info('There are ' + str(len(jobs_done)) + ' jobs already done, will only redo with --redo option')

if len(jobs_wrote) > 0:
    if len(jobs_disp) > 0:
        logging.info('Wrote ' + str(len(jobs_wrote)) + ' job descriptions and dispatched ' + str(len(jobs_disp)))
    else:
        logging.info('Wrote ' + str(len(jobs_wrote)) + ' job descriptions, will only dispatch with --dispatch option')

sys.exit(0)



##

## warn/error if no phased haplotypes at all
#        logging.error('Directory [ ' + phasedir + ' ] does not contain phased haplotypes')
#    sys.exit(2)
