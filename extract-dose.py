#!/usr/bin/python
from optparse import OptionParser
## deprecated, but argparse requires python >= 2.7 does not work on GSK servers
import logging
import gzip
import sys
import os
import stat
import subprocess
import string
import math

###
### parse command line arguments and check
###

parser = OptionParser(description = 'usage: %prog OPTIONS')
parser.add_option('-i', '--input', help = 'input file basename (without .dose.gz and .info.gz extensions)', metavar = 'BASENAME',
                  type = 'string', dest = 'input',
                  default = '')
parser.add_option('-o', '--output', help = 'output file basename (without .dose.gz and .info.gz extensions)', metavar = 'BASENAME',
                  type = 'string', dest = 'output',
                  default = '')
parser.add_option('-p', '--pos', help = 'extract variant by individual position', metavar = 'POS',
                  action = 'append', type = 'string', dest = 'pos', default = [])
parser.add_option('-r', '--range', help = 'extract variants by position range', 
                  action = 'store_true', dest = 'range',
                  default = False)
parser.add_option('-b', '--begin', help = 'set range begin position', metavar = 'BEGIN',
                  type = 'int', dest = 'begin',
                  default = '0')
parser.add_option('-e', '--end', help = 'set range end position', metavar = 'END',
                  type = 'int', dest = 'end',
                  default = '0')
(options, args) = parser.parse_args()

logging.basicConfig(stream=sys.stdout, format='%(levelname)s:%(message)s', level=logging.DEBUG)
logging.info('dose extractor v0.1 Toby.x.Johnson@gsk.com')

## check end >= begin, warn if begin/end used without range 

###
###
###

posdict = dict()
for pos in options.pos:
    posdict[pos] = 1

if len(posdict) > 0:
    logging.info('Extracting dosage for variants at positions [ ' + ','.join(posdict.keys()) + ' ]')

if options.range:
    logging.info('Extracting dosage for variants at position range [ ' + str(options.begin) + ' ] to [ ' + str(options.end) + ' ]')

try:
    info_in_fh = gzip.open(options.input + '.info.gz', 'r')
except IOError:
    logging.error('Cannot open [ ' + options.input + '.info.gz ]')
    sys.exit(2)
logging.info('Reading info from [ ' + options.input + '.info.gz ]')

try:
    info_out_fh = gzip.open(options.output + '.info.gz', 'w')
except IOError:
    logging.error('Cannot open [ ' + options.output + '.info.gz ] for writing')
    sys.exit(2)
logging.info('Writing info to [ ' + options.output + '.info.gz ]')
    
linesel = list()
snp_line = 0 # count lines EXCLUDING HEADER
first_line = True
range_flag = 0 # 0=before 1=in 2=after
for info_in in info_in_fh:
    info_in_data = info_in.strip().split(None, 2)
    if first_line:
        if info_in_data[0] != 'SNP':
            logging.error('First line [ ' + options.input + '.info.gz ] does not begin \'SNP\'')
        first_line = False
        info_out_fh.write(info_in)
    else:
        snp_line = snp_line + 1
        snp_data = info_in_data[0].split(":", 3)
        if len(snp_data) >= 2:
            if (snp_data[1] in posdict):
                linesel.append(snp_line)
                info_out_fh.write(info_in)
            if options.range:
                if (range_flag == 0 and int(snp_data[1]) >= options.begin):
                    range_flag = 1
                if (range_flag == 1 and int(snp_data[1]) > options.end):
                    range_flag = 2
                if (range_flag == 1):
                    linesel.append(snp_line)
                    info_out_fh.write(info_in)

info_out_fh.close()
info_in_fh.close()

logging.info('Read info for [ ' + str(snp_line) + ' ] SNPs from [ ' + options.input + '.info.gz ]')
logging.info('Wrote info for [ ' + str(len(linesel)) + ' ] SNPs to [ ' + options.output + '.info.gz ]')

###
###
###

try:
    dose_in_fh = gzip.open(options.input + '.dose.gz', 'r')
except IOError:
    logging.error('Cannot open [ ' + options.input + '.dose.gz ]')
    sys.exit(2)
logging.info('Reading dosages from [ ' + options.input + '.dose.gz ]')

try:
    dose_out_fh = gzip.open(options.output + '.dose.gz', 'w')
except IOError:
    logging.error('Cannot open [ ' + options.output + '.dose.gz ] for writing')
    sys.exit(2)
logging.info('Writing dosages to [ ' + options.output + '.dose.gz ]')

dose_line = 0
for dose_in in dose_in_fh:
    dose_in_data = dose_in.strip().split()
    dose_line = dose_line + 1
    dose_out_fh.write(dose_in_data[0] + '\t' + dose_in_data[1])
    for sel in linesel:
        dose_out_fh.write('\t' + dose_in_data[sel + 1])  # sel-th SNP is (sel+2)-th column (1-begin) is (sel+1)-th element (0-begin)
    dose_out_fh.write('\n')

dose_out_fh.close()
dose_in_fh.close()

logging.info('Read dosages for [ ' + str(dose_line) + ' ] subjects from [ ' + options.input + '.dose.gz ]')
logging.info('Wrote dosages for [ ' + str(len(linesel)) + ' ] SNPs and [ ' + str(dose_line) + ' ] subjects to [ ' + options.output + '.dose.gz ]')

sys.exit(0)
