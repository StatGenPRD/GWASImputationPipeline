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
parser.add_option('--chr', help = 'phase specified chromosome(s) only', metavar = 'CHR',
                  action = 'append', type = 'string', dest = 'chr', default = [])
parser.add_option('-g', '--genotypes', help = 'basename of PLINK fileset of measured genotypes (default plink)', metavar = 'BASENAME',
                  type = 'string', dest = 'genotypes',
                  default = 'plink')
parser.add_option('-a', '--algorithm', help = 'phasing algorithm (auto, hapiur, or mach, default auto)', metavar = 'NAME',
                  type = 'string', dest = 'algorithm',
                  default = 'auto')
parser.add_option('-p', '--phased', help = 'directory for phased haplotypes (default phased-ALGORITHM-hg19)', metavar = 'DIR',
                  type = 'string', dest = 'phased',
                  default = '')
parser.add_option('-n', '--name', help = 'job name (default chr[chr] or JOBNAME for chr 23)', metavar = 'JOBNAME',
                  type = 'string', dest = 'jobname',
                  default = '')
parser.add_option('--plinkopt', help = 'options for PLINK', metavar = 'OPTIONS',
                  type = 'string', dest = 'plinkopt', default = '')
parser.add_option('--hapiur-win', help = 'HAPI-UR window size (default 0 = automatic)', metavar = 'INTEGER',
                  type = 'int', dest = 'hapiurwin',
                  default = 0)
parser.add_option('--machopt', help = 'options for mach (default \'--rounds 20 --states 200\')', metavar = 'OPTIONS',
                  type = 'string', dest = 'machopt',
                  default = '--rounds 20 --states 200')
parser.add_option('--plink-path', help = 'full path to plink binary', metavar = 'PLINK',
                  type = 'string', dest = 'plink',
                  default = '/GWD/bioinfo/projects/RD-TSci-Software/CB/linuxbrew/bin/plink')
parser.add_option('--hapiur-path', help = 'full path to hapiur binary', metavar = 'HAPIUR',
                  type = 'string', dest = 'hapiur',
                  default = '/GWD/appbase/projects/statgen/GXapp/hapiur/hapi-ur')
parser.add_option('--hapi2mach-path', help = 'full path to convert-hapiur-mach.py script', metavar = 'HAPI2MACH',
                  type = 'string', dest = 'hapi2mach',
                  default = '/GWD/appbase/projects/statgen/GXapp/imputation/convert-hapiur-mach.py')
parser.add_option('--mach-path', help = 'full path to mach binary', metavar = 'MACH',
                  type = 'string', dest = 'mach',
                  default = '/GWD/appbase/projects/statgen/GXapp/mach1/mach1')
parser.add_option('--submit', help = 'command to submit job script', metavar = 'COMMAND',
                  type = 'string', dest = 'submit',
                  default = '/GWD/bioinfo/projects/lsf/SGE/6.2u5/bin/lx24-amd64/qsub')
parser.add_option('--submitqueue', help = '--queue option argument (default rhel7)', metavar = 'QUEUE',
                  type = 'string', dest = 'submitqueue',
                  default = 'rhel7')
parser.add_option('--submitmem', help = 'job memory (Mbyte/marker, default 0.1)', metavar = 'VALUE',
                  type = 'float', dest = 'submitmem',
                  default = 0.1)
parser.add_option('--submitopt', help = 'append other option(s) for submit command', metavar = 'OPTION',
                  action = 'append', type = 'string', dest = 'submitopt', default = [])
(options, args) = parser.parse_args()

logging.basicConfig(stream=sys.stdout, format='%(levelname)s:%(message)s', level=logging.DEBUG)
logging.info('haplotype phasing dispatcher v0.1 Toby.x.Johnson@gsk.com')

## check what chromosomes are to be phased
##
if len(options.chr) == 0:
    logging.info('Imputing chromosomes 1-22 (default)')
    for chridx in range(22): # C style range 0..21
        options.chr.append(str(chridx + 1))
else:
    logging.info('Using --chr arguments, imputing only chromosomes [ ' + ','.join(options.chr) + ' ]')

## check binaries and scripts exist
##
if not os.path.isfile(options.plink):
    logging.error('Could not find plink binary [ ' + options.plink + ' ]')
    sys.exit(1)

if options.algorithm in ['auto', 'hapiur', 'mach']:
    if not options.algorithm == 'mach':
        if not os.path.isfile(options.hapiur):
            logging.error('Could not find hapiur binary [ ' + options.hapiur + ' ]')
            sys.exit(1)
        if not os.path.isfile(options.hapi2mach):
            logging.error('Could not find hapi2mach script [ ' + options.hapi2mach + ' ]')
            sys.exit(1)
    if not options.algorithm == 'hapiur':
        if not os.path.isfile(options.mach):
            logging.error('Could not find mach binary [ ' + options.mach + ' ]')
            sys.exit(1)
else:
    logging.error('Algorithm [ ' + options.algorithm + ' ] not recognised; use --algorithm=auto, hapiur, or mach')
    sys.exit(1)

if not os.path.isfile(options.submit):
    logging.error('Could not find submit binary [ ' + options.submit + ' ]')
    sys.exit(1)

if not options.submitqueue in ['', 'any', 'default']:
    options.submitopt.extend(['-q', options.submitqueue])
if options.submitmem > 0:
    logging.info('Submit command is [ ' + options.submit + ' ' + ' '.join(options.submitopt) + ' -l h_data=xxG ]')
    logging.info(' where xxG is >= ' + str(options.submitmem) + ' Mbyte/marker')
else:
    logging.info('submit command is [ ' + options.submit +  ' '.join(options.submitopt) + ' ]')
    
## check input genotype files exist
##
if not os.path.isabs(options.genotypes):
    options.genotypes = os.path.abspath(options.genotypes)
    logging.info('Inferred absolute path to measured genotypes [ ' + options.genotypes + ' ]')

logging.info('Checking measured genotypes with basename [ ' + options.genotypes + ' ]')
for ext in ['.bed', '.bim', '.fam']:
    if not os.path.isfile(options.genotypes + ext):
        logging.error('File [ ' + options.genotypes + ext + ' ] does not exist')
        sys.exit(2)

## if automatic algorithm selection, count subjects
##
if options.algorithm == 'auto':
    logging.info('Using automatic selection of phasing algorithm')
    logging.info('Counting subjects in [ ' + options.genotypes + '.fam ]')
    try:
        fam_fh = open(options.genotypes + '.fam', 'r')
    except IOError:
        logging.error('Cannot open [ ' + options.genotypes + '.fam ]')
        sys.exit(2)
    subjects = len(fam_fh.readlines()) # is this memory expensive?
    fam_fh.close()
    logging.info('Counted [ ' + str(subjects) + ' ] subjects in [ ' + options.genotypes + '.fam ]')
    if subjects >= 1000:
        logging.info('Using [ hapiur ] algorithm for phasing >=1000 subjects')
        options.algorithm = 'hapiur'
    elif subjects >0:
        logging.info('Using [ mach ] algorithm for phasing <1000 subjects')
        options.algorithm = 'mach'
    else:
        logging.error('Zero subjects in [ ' + options.genotypes + '.fam ]')
        sys.exit(2)
else:
    logging.info('Using [ ' + options.algorithm +' ] algorithm for phasing')
        
## read .bim file and count markers by autosome
##
logging.info('Counting markers in [ ' + options.genotypes + '.bim ]')
chrom_counts = dict()
for chridx in range(23): # C style range 0..21 #add 22 for chrX by Li Li
    chrom = str(chridx + 1)
    chrom_counts[chrom] = 0

try:
    bim_fh = open(options.genotypes + '.bim', 'r')
except IOError:
    logging.error('Cannot open [ ' + options.genotypes + '.bim ]')
    sys.exit(2)

for snp_read in bim_fh:
    snp_data = snp_read.strip().split(None, 1)
    if len(snp_data) >= 2:
        if snp_data[0] in chrom_counts:
            chrom_counts[snp_data[0]] += 1

bim_fh.close()
logging.info('Counted markers in [ ' + options.genotypes + '.bim ]')

chrom_set = list()
marker_total = 0
for chrom in sorted(chrom_counts.keys(), key=int):
    if chrom_counts[chrom] > 0:
        chrom_set.append(chrom)
        marker_total += chrom_counts[chrom]
logging.info('There are ' + str(marker_total) + ' markers on ' + str(len(chrom_set)) + ' chromosomes') #changed autosomes to chromosomes

## Williams et al. recommend extrapolating linearly using their marker densities and window sizes
## so that for A autosomal markers, window size would be (90 - 64)/(755008 - 386353) * (A - 386353) + 64
## Hence --win 100 should be sufficient for 1M and OmniExpressExome arrays, but --win 300 might be needed
## for 5M data with ~4M markers passing QC
## By private email, Amy Williams suggested --win 200 might be enough for 5M data
if options.algorithm == 'hapiur':
    if len(chrom_set) == 22:
        auto_win = int(float(90 - 64)/float(755008 - 386353) * float(marker_total - 386353)) + 64
        if auto_win < 100:
            auto_win = 100
        logging.info('Extrapolated HAPI-UR recommended window size = max(' + str(auto_win) + ', 200)')
        if auto_win > 200:
            auto_win = 200
        if options.hapiurwin > 0:
            logging.info('Using command line specified window size [ ' + str(options.hapiurwin) + ' ]')
        else:
            options.hapiurwin = auto_win
            logging.info('Using automatic window size [ ' + str(options.hapiurwin) + ' ]')
    else:
        logging.info('Cannot recommend window size without marker counts for all 22 autosomes')
        if options.hapiurwin > 0:
            logging.info('Using command line specified window size [ ' + str(options.hapiurwin) + ' ]')
        else:
            options.hapiurwin = 100
            logging.info('Using default window size [ ' + str(options.hapiurwin) + ' ]')

if len(chrom_set) == 0:
    logging.error('Measured genotypes [ ' + options.genotypes + '.bim ] do not contain any autosomal markers')
    sys.exit(2)

# make temporary list as cannot remove elements from list we are iterating over
chrom_set_tmp = list()
for chrom in chrom_set:
    chrom_set_tmp.append(chrom)

for chrom in chrom_set_tmp:
    if not chrom in options.chr:
        chrom_set.remove(chrom)

###
### Check existence/make directories for phasing and job scripts
###

if options.phased == '':
    options.phased = 'phased-' + options.algorithm + '-hg19'
    logging.info('Set automatic output phasing directory [ ' + options.phased + ' ]')

if os.path.isabs(options.phased):
    phasedir = options.phased
else:
    phasedir = os.path.abspath(options.phased)
    logging.info('Inferred absolute path for output phasing [ ' + phasedir + ' ]')

jobdir = os.path.join(phasedir, 'jobs')
for thisdir in [phasedir, jobdir]:
    if not os.path.isdir(thisdir):
        try:
            os.mkdir(thisdir)
        except OSError:
            logging.error('Directory [ ' + thisdir + ' ] does not exist and cannot be created')
            sys.exit(2)

logging.info('Setting up jobs to output results in [ ' + phasedir + ' ]')
logging.info('Writing job scripts in [ ' + jobdir + ' ]')

###
### Write job scripts and dispatch
###

jobs_done = list()
jobs_wrote = list()
jobs_disp = list()
for chrom in chrom_set:
    job_name = 'chr' + chrom
    if options.jobname != '':
    	job_name = options.jobname
    if not options.redo and os.path.isfile(os.path.join(phasedir, job_name + '.done')):
        jobs_done.append(job_name)
    else:
        try:
            job_fh = open(os.path.join(jobdir, job_name + '.sh'), 'w')
        except IOError:
            logging.error('Could not write job script [ ' + os.path.join(jobdir, job_name + '.sh') + ' ]')
            sys.exit(2)
            
        job_fh.write('#!/bin/bash\n')
        job_fh.write('set -o nounset\n')
        job_fh.write('set -o errexit\n')
        job_fh.write('/bin/uname -a\n')
        if options.algorithm == 'mach':
            job_fh.write('/bin/echo \'MACH phasing starting\' $(/bin/date)\n')
            job_fh.write('/bin/rm -f ' + os.path.join(phasedir, job_name + '.done') + '\n')
            job_fh.write(options.plink + ' ' + options.plinkopt \
                         + ' --bfile ' + options.genotypes + ' --chr ' + chrom \
                         + ' --recode --out ' + os.path.join(phasedir, job_name) \
                         + '\n')
            job_fh.write('/bin/gzip -f ' + os.path.join(phasedir, job_name + '.ped') + ' && \\\n')
            job_fh.write('(/bin/echo "T PHENOTYPE";'
                         + ' /bin/awk \'{print "M " $1 ":" $4}\' ' + os.path.join(phasedir, job_name + '.map') + ')' \
                         + ' | /bin/gzip - >' + os.path.join(phasedir, job_name + '.dat.gz') \
                         + '\n')
            job_fh.write(options.mach \
                         + ' --datfile ' + os.path.join(phasedir, job_name + '.dat.gz') \
                         + ' --pedfile ' + os.path.join(phasedir, job_name + '.ped.gz') \
                         + ' --phase ' + options.machopt \
                         + ' --prefix ' + os.path.join(phasedir, job_name) \
## no log file because overwrites PLINK --recode logfile; merge back into main Makefile
                         + '\n')
            job_fh.write('/bin/zcat ' + os.path.join(phasedir, job_name + '.dat.gz') \
                         + ' | /bin/awk \'$1 == "M" {gsub("23:", "X:"); print $2}\'' \
                         + ' | /bin/gzip - >' + os.path.join(phasedir, job_name + '.snps.gz') \
                         + '\n')
            job_fh.write('/bin/grep \'Analysis took\' ' \
                         + os.path.join(phasedir, job_name + '.log') \
                         + ' >' + os.path.join(phasedir, job_name + '.done') + '\n')
        elif options.algorithm == 'hapiur':
            job_fh.write('/bin/echo \'HAPI-UR phasing starting\' $(/bin/date)\n')
            job_fh.write('/bin/rm -f ' + os.path.join(phasedir, job_name + '.done') + '\n')
            job_fh.write(options.plink + ' ' + options.plinkopt \
                         + ' --bfile ' + options.genotypes + ' --chr ' + chrom \
                         + ' --make-bed --out ' + os.path.join(phasedir, job_name) \
                         + '\n')
            ## $(R) --vanilla --args $(subst .gz,.bim,$@) <$(ADDGENETICMAP) >/dev/null 2>/dev/null
            job_fh.write(options.hapiur + ' --win ' + str(options.hapiurwin) \
                         + ' --plink ' + os.path.join(phasedir, job_name) \
                         + ' --out ' + os.path.join(phasedir, job_name + '-hapiur') \
                         + ' >' + os.path.join(phasedir, job_name + '-hapiur.progress') \
                         + '\n')
            job_fh.write('/bin/gzip -f ' + os.path.join(phasedir, job_name + '-hapiur.phgeno') + '\n')
            job_fh.write(options.hapi2mach \
                         + ' ' + os.path.join(phasedir, job_name + '-hapiur') \
                         + ' ' + os.path.join(phasedir, job_name) \
                         + '\n')
            job_fh.write('/bin/grep \'Total wall clock time\' ' \
                         + os.path.join(phasedir, job_name + '-hapiur.progress') \
                         + ' >' + os.path.join(phasedir, job_name + '.done') + '\n')
        else:
            logging.error('Algorithm [ ' + options.algorithm + ' ] not implemented; must be --algorithm=mach or --algorithm=hapiur')
            ## this should never happen because options.algorithm already checked
            sys.exit(1)

        job_fh.write('/bin/touch ' + os.path.join(phasedir, job_name + '.done') + '\n')
        if options.algorithm == 'mach':
            job_fh.write('/bin/echo \'MACH phasing finished\' $(/bin/date)\n')
        elif options.algorithm == 'hapiur':
            job_fh.write('/bin/echo \'HAPI-UR phasing finished\' $(/bin/date)\n')
        else:
            logging.error('Algorithm [ ' + options.algorithm + ' ] not implemented; must be --algorithm=mach or --algorithm=hapiur')
            ## this should never happen because options.algorithm already checked
            sys.exit(1)
        job_fh.close()
        os.chmod(os.path.join(jobdir, job_name + '.sh'), stat.S_IRUSR | stat.S_IRGRP | stat.S_IROTH | stat.S_IWUSR | stat.S_IXUSR )
        jobs_wrote.append(job_name)

        thisopt = list(options.submitopt) # must not modify options for other chromosomes
        if options.submitmem > 0:
            submitmem = ['-l', 'h_data=' + str(max(int(float(chrom_counts[chrom])*options.submitmem*1e-3) + 1, 2)) + 'G']
            logging.info(job_name + ' appending option for submit command [ ' + ' '.join(submitmem) + ' ]')
            thisopt.extend(submitmem)
        
        if options.dispatch:
            subprocess.call([options.submit] + thisopt + \
                            ['-o', os.path.join(jobdir, job_name + '.out'), \
                             '-e', os.path.join(jobdir, job_name + '.err'), \
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
