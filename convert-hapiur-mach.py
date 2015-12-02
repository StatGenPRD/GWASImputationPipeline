#!/usr/bin/python
import gzip
import logging
import string
import sys

logging.basicConfig(stream=sys.stdout, format='%(levelname)s:%(message)s', level=logging.DEBUG)
logging.info('HAPI-UR to MACH haplotype format conversion tool v0.2 Toby.x.Johnson@gsk.com')

### HAPI-UR to MACH haplotype format conversion tool v0.1
### Toby.x.Johnson@gsk.com 04 Feb 2013
###
### Requires two command line arguments:
### First argument INROOT is the root name of the input HAPI-UR format haplotypes
### Second argument OUTROOT is the root name of the output MACH format haplotypes
###
### This script converts fileset INROOT.phsnp, INROOT.phind and INROOT.phgeno.gz to
### fileset OUTROOT.snps.gz and OUTROOT.gz
###
### Note that HAPI-UR v1.01 writes a .phgeno file that has to be gzipped before
### running this script
###
### v0.2 allows underscores in subject IDs

###
### parse command line arguments
###

if len(sys.argv) == 3:
    inroot = sys.argv[1]
    outroot = sys.argv[2]
else:
    logging.error('Need two arguments: INROOT, OUTROOT (to which .phsnp, .phind, .phgeno.gz, .snps.gz, .gz extensions will be added')
    sys.exit(2)

try:
    mapfh = open(inroot + '.phsnp', 'r')
except IOError:
    logging.error('Cannot open file [ %s.phsnp ]', inroot)
    sys.exit(1)

logging.info('Reading map information from [ %s.phsnp ] ...', inroot)
snpname = list()
snpalleles = list()
for mapline in mapfh:
    mapinfo = mapline.strip().split()
    if len(mapinfo) != 6:
        logging.error('Expecting 6 items but found [ %i ] in map file [ %s.phsnp ] line [ %i ]', len(mapinfo), inroot, len(snpname) + 1)
        sys.exit(1)
    snpname.append(mapinfo[0])
    snpalleles.append([mapinfo[5], mapinfo[4]]) # columns 5,4 are PLINK alleles A2,A1 coded 0,1 respectively

logging.info('Read map information for [ %i ] SNPs from [ %s.phsnp ]', len(snpname), inroot)
mapfh.close()

try:
    subjfh = open(inroot + '.phind', 'r')
except IOError:
    logging.error('Cannot open file [ %s.phind ]', inroot)
    sys.exit(1)

logging.info('Reading subject information from [ %s.phind ] ...', inroot)
subjid = list()
for subjline in subjfh:
    subjinfo = subjline.strip().split()
    if len(subjinfo) != 3:
        logging.error('Expecting 3 items but found [ %i ] in subject file [ %s.phind ] line [ %i ]', len(subjinfo), inroot, len(subjid) + 1)
        sys.exit(1)

    # parse HAPI-UR format '1:HZA106829.0010840_A' to MACH format '1->HZA106829.0010840 HAPLO1'
    if subjinfo[0].count(':') != 1:
        logging.error('Subject ID [ %s ] does not contain exactly one colon (:), in subject file [ %s.phind ] line [ %i ]', subjinfo[0], inroot, len(subjid) + 1)
        sys.exit(1)
#    if subjinfo[0].count('_') != 1:
#        logging.error('Subject ID [ %s ] does not contain exactly one underscore (_), in subject file [ %s.phind ] line [ %i ]', subjinfo[0], inroot, len(subjid) + 1)
#        sys.exit(1)
    if subjinfo[0].endswith('_A') == subjinfo[0].endswith('_B'): # logical !XOR
        logging.error('Subject ID [ %s ] does not end with \'_A\' or \'_B\', in subject file [ %s.phind ] line [ %i ]', subjinfo[0], inroot, len(subjid) + 1)
        sys.exit(1)
# previously
# subjid.append(subjinfo[0].replace('_B', ' HAPLO2', 1).replace('_A', ' HAPLO1', 1).replace(':', '->', 1))
# only worked if single underscore
    if subjinfo[0].endswith('_A'):
        subjid.append((subjinfo[0][:subjinfo[0].rfind('_A')]).replace(':', '->', 1) + ' HAPLO1')
    elif subjinfo[0].endswith('_B'):    
        subjid.append((subjinfo[0][:subjinfo[0].rfind('_B')]).replace(':', '->', 1) + ' HAPLO2')
    else:
        # should never happen
        logging.error('Internal error: Subject ID [ %s ] does not end with \'_A\' or \'_B\', in subject file [ %s.phind ] line [ %i ]', subjinfo[0], inroot, len(subjid) + 1)
        sys.exit(1)

logging.info('Read subject information for [ %i ] haplotypes from [ %s.phind ]', len(subjid), inroot)
subjfh.close()

try:
    hapfh = gzip.open(inroot + '.phgeno.gz', 'r')
except IOError:
    logging.error('Cannot open file [ %s.phgeno.gz ]', inroot)
    sys.exit(1)

logging.info('Reading haplotypes from [ %s.phgeno.gz ] ...', inroot)
hap = list()
for hapline in hapfh:
    thishap = hapline.strip()
##    logging.info('line [%i] length [%i]', len(hap) + 1, len(thishap))
    if len(thishap) != len(subjid):
        logging.error('Expecting line length [ %i ] but found [ %i ] in haplotype file [ %s.phgeno.gz ] line [ %i ]', len(subjid), len(thishap), inroot, len(hap) + 1)
        sys.exit(1)
    hap.append(thishap)

if len(hap) != len(snpname):
    logging.error('Expecting [ %i ] lines but found [ %i ] in haplotype file [ %s.phgeno.gz ]', len(snpname), len(hap), inroot)
    sys.exit(1)

logging.info('Read haplotypes from [ %s.phgeno.gz ]', inroot)
hapfh.close()

logging.info('Writing SNP names to [ %s.snps.gz ] ...', outroot)
try:
    osnpfh = gzip.open(outroot + '.snps.gz', 'w')
except IOError:
    logging.error('Cannot write to file [ %s.snps.gz ]', outroot)
    sys.exit(1)

for j in range(len(snpname)):
    osnpfh.write(snpname[j] + '\n')

osnpfh.close()
logging.info('Wrote SNP names for [ %i ] SNPs to [ %s.snps.gz ]', len(snpname), outroot)

logging.info('Writing haplotypes to [ %s.gz ] ...', outroot)
try:
    ohapfh = gzip.open(outroot + '.gz', 'w')
except IOError:
    logging.error('Cannot write to file [ %s.gz ]', outroot)
    sys.exit(1)

for i in range(len(subjid)):
    ohapfh.write(subjid[i] + ' ')
    for j in range(len(snpname)):
        ohapfh.write(snpalleles[j][int(hap[j][i])])
    ohapfh.write('\n')

ohapfh.close()
logging.info('Wrote haplotypes for [ %i ] haplotypes and [ %i ] SNPs to [ %s.gz ]', len(subjid), len(snpname), outroot)

sys.exit(0)
