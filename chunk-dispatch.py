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
parser.add_option('-i', '--imputed', help = 'directory for imputed chunks (default imputed-20120314)', metavar = 'DIR',
                  type = 'string', dest = 'imputed',
                  default = 'imputed-20120314')
parser.add_option('--chunkMb', help = 'chunk size in Mb (default 7)', metavar = 'MB',
                  type = 'float', dest = 'chunkMb',
                  default = 7.)
parser.add_option('--windowMb', help = 'flanking window size in Mb (default 0.25)', metavar = 'MB',
                  type = 'float', dest = 'windowMb',
                  default = 0.25)
parser.add_option('--minimacopt', help = 'options for minimac (default \'--rounds 5 --states 200 --probs\')', metavar = 'OPTIONS',
                  type = 'string', dest = 'minimacopt',
                  default = '--rounds 5 --states 200 --probs')
parser.add_option('--minimac-path', help = 'full path to minimac binary', metavar = 'MINIMAC',
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
logging.info('minimac chunk dispatcher v0.2 Toby.x.Johnson@gsk.com')

if len(options.chr) == 0:
    logging.info('Imputing chromosomes 1-22 (default)')
    for chridx in range(22): # C style range 0..21
        options.chr.append(str(chridx + 1))
else:
    logging.info('Using --chr arguments, imputing only chromosomes [ ' + ','.join(options.chr) + ' ]')

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
    logging.error('Could not find minimac binary [ ' + options.minimac + ' ]')
    sys.exit(1)

if not os.path.isfile(options.submit):
    logging.error('Could not find submit binary [ ' + options.submit + ' ]')
    sys.exit(1)

if not options.submitqueue in ['', 'any', 'default']:
    options.submitopt.append('--queue=' + options.submitqueue)
logging.info('Submit command is [ ' + options.submit + ' ' + ' '.join(options.submitopt) + ' ]')

###
### Hard coded options for nosingletons 1000G full panel
### Hard coded paths for us1us
###

vcfpath = '/GWD/appbase/projects/statgen/RD-MDD-GX_PUBLIC/1KG/share.sph.umich.edu/1000genomes/fullProject/2012.03.14/nosingletons'
vcfpre = 'chr'
vcfpost = '.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.ALL.vcf.gz'


###
### Check phased haplotypes exist
###

if os.path.isabs(options.phased):
    phasedir = options.phased
else:
    phasedir = os.path.abspath(options.phased)
    logging.info('Inferred absolute path to phased haplotypes [ ' + phasedir + ' ]')

if not os.path.isdir(phasedir):
    logging.error('Directory [ ' + phasedir + ' ] does not exist')
    sys.exit(2)

if options.sleep:
    while True:
        chrom_waiting = list()
        for chrom in options.chr:
            if not os.path.isfile(os.path.join(phasedir, 'chr' + chrom + '.done')):
                chrom_waiting.append(chrom)
        if len(chrom_waiting) > 0:
            logging.info('Chromosomes [ ' + ','.join(chrom_waiting) + ' ] not yet phased, ' + time.strftime('%H:%M %d-%m-%Y') + ' now; sleeping for 1 hour')
            time.sleep(3600)
        else:
            break # break while True loop
    
chrom_set = list()
for chridx in range(22): # C style range 0..21
    chrom = str(chridx + 1)
    if os.path.isfile(os.path.join(phasedir, 'chr' + chrom + '.gz')) and \
       os.path.isfile(os.path.join(phasedir, 'chr' + chrom + '.snps.gz')):
        chrom_set.append(chrom)

if len(chrom_set) == 0:
    logging.error('Directory [ ' + phasedir + ' ] does not contain phased haplotypes')
    sys.exit(2)
    
logging.info('Directory [ ' + phasedir + ' ] contains phased haplotypes')
logging.info('Haplotypes for chromosomes [ ' + ','.join(chrom_set) + ' ]')

###
### Get chromosome begin and end positions from data on last variant in VCF files
###

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
        if chrom_set.count(chrom_data[0]):
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
        if chrom_set.count(chrom_data[0]):
            chrom_end[chrom_data[0]] = int(chrom_data[1])
    else:
        logging.error('Error reading in [ ' + os.path.join(vcfpath, 'chrend.data') + ' ]')

chrom_end_fh.close()
logging.info('Read end positions for chromosomes [ ' + ','.join(sorted(chrom_end.keys(), key=int)) + ' ]')

# make temporary list as cannot remove elements from list we are iterating over
chrom_set_tmp = list()
for chrom in chrom_set:
    chrom_set_tmp.append(chrom)

for chrom in chrom_set_tmp:
    if not chrom in options.chr:
        chrom_set.remove(chrom)
    if not chrom in chrom_begin:
        chrom_set.remove(chrom)
    if not chrom in chrom_end:
        chrom_set.remove(chrom)

logging.info('Imputing for chromosomes specified and with haplotypes, begin and end positions [ ' + ','.join(chrom_set) + ' ]') 

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

###
### Calculate chunk start and ends, write job scripts and dispatch
###

jobs_done = list()
jobs_wrote = list()
jobs_disp = list()
for chrom in chrom_set:
    num_chunk = int(math.ceil(float(chrom_end[chrom] - chrom_begin[chrom] + 1) / minimac_chunk_size))
    logging.info('chr' + chrom + ':' + str(chrom_begin[chrom]) + '-' + str(chrom_end[chrom]) + ' to be split into ' + str(num_chunk) + ' chunks')
    for chunk in range(num_chunk):
        chunk_id = str(chunk + 1)
        chunk_start = str(chunk*minimac_chunk_size + chrom_begin[chrom])
        chunk_end = str((chunk + 1)*minimac_chunk_size + chrom_begin[chrom] - 1)
        job_name = 'chr' + chrom + 'chunk' + chunk_id

        if not options.redo and os.path.isfile(os.path.join(imputedir, job_name + '.done')):
            jobs_done.append(job_name)
        else:
            try:
                job_fh = open(os.path.join(jobdir, job_name + '.sh'), 'w')
            except IOError:
                logging.error('Could not write job script [ ' + os.path.join(jobdir, job_name + '.sh') + ' ]')
                sys.exit(2)

            job_fh.write('#!/bin/bash\n')
            job_fh.write('/bin/uname -a\n')
            job_fh.write('/bin/echo \'Chunk imputation starting\' $(/bin/date)\n')
            job_fh.write('/bin/rm -f ' + os.path.join(imputedir, job_name + '.done') + '\n')
            job_fh.write(options.minimac + ' --vcfReference --refHaps ' + os.path.join(vcfpath, vcfpre + chrom + vcfpost) \
                         + ' --vcfstart ' + chunk_start + ' --vcfend ' + chunk_end + ' --vcfwindow ' + minimac_window \
                         + ' --haps ' + os.path.join(phasedir, 'chr' + chrom + '.gz') \
                         + ' --snps ' + os.path.join(phasedir, 'chr' + chrom + '.snps.gz') \
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

if not options.redo and len(jobs_done) > 0:
    logging.info('There are ' + str(len(jobs_done)) + ' jobs already done, will only redo with --redo option')

if len(jobs_wrote) > 0:
    if len(jobs_disp) > 0:
        logging.info('Wrote ' + str(len(jobs_wrote)) + ' job descriptions and dispatched ' + str(len(jobs_disp)))
    else:
        logging.info('Wrote ' + str(len(jobs_wrote)) + ' job descriptions, will only dispatch with --dispatch option')

sys.exit(0)
