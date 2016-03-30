#!/usr/bin/python
import gzip
import logging
import string
import sys

logging.basicConfig(format='%(levelname)s:%(message)s') # can add level=logging.DEBUG
## logging to stderr
logging.info('1KG VCF alignment (lookup/liftover) tool v0.2 Toby.x.Johnson@gsk.com') # should these go into output?

### notes:
### MATCHBY=NAME matches on lower case of SNP ID
### matches on BOTH alleles, hence SOFT and HARD liftovers both drop monomorphic SNPs that PLINK calls 0/A etc


###
### parse command line arguments
###

if len(sys.argv) == 5:
    matchby = sys.argv[1]
    mode = sys.argv[2]
    vcffile = sys.argv[3]
    markerfile = sys.argv[4]
else:
    logging.error('Need four arguments: MATCHBY [name|chrpos], MODE [hard|soft], VCFFILE [reference], and MARKERFILE [target .bim or .frq file]')
    sys.exit(2)

matchby = matchby.strip().upper()
if matchby in ['NAME', 'CHRPOS'] :
    logging.info('Match mode is %s', matchby)
    if matchby == 'NAME' :
        logging.info('Target markers will be matched to reference markers by name (case insensitive)')
    elif matchby == 'CHRPOS' :
        if markerfile.endswith('.frq'):
            logging.error('Cannot match by position using target .frq file')
            sys.exit(2)
        logging.info('Target markers will be matched to reference markers by chromosome and position (build sensitive)')
else:
    logging.error('First argument (MATCHBY) must be either \'NAME\' or \'CHRPOS\'')
    sys.exit(2)
    
mode = mode.strip().upper()
if mode in ['HARD', 'SOFT', 'TEST'] :
    logging.info('Liftover mode is %s', mode)
    if mode == 'HARD' :
        logging.info('Only A/C, A/G, C/T and G/T SNPs will be kept; strand will be flipped to match reference')
        safealleles = ['A/C', 'A/G', 'C/A', 'C/T', 'G/A', 'G/T', 'T/C', 'T/G']
    elif mode == 'SOFT' :
        logging.info('All (A/C, A/G, A/T, C/G, C/T, G/T) SNPs will be kept, if alleles match reference WITHOUT strand flipping')
        safealleles = ['A/C', 'A/G', 'A/T', 'C/A', 'C/G', 'C/T', 'G/A', 'G/C', 'G/T', 'T/A', 'T/C', 'T/G']
    elif mode == 'TEST' :
        logging.info('All SNPs that match reference by name will be kept; USE ONLY FOR PRE-LIFTOVER TESTING')
else:
    logging.error('Second argument (MODE) must be either \'HARD\', \'SOFT\', or \'TEST\'')
    sys.exit(2)

###
### build dictionary of SNPs to lookup
### python guarantees(?) lookups average complexity O(1)    
###   

try:
    markerfh = open(markerfile, 'r')
except IOError:
    logging.error('Cannot open file [ %s ]', markerfile)
    sys.exit(1)

logging.info('Reading SNP IDs from [ %s ] ...', markerfile)
snpdict = dict() # keys are ID to match by (lowercase of SNP IDs markerfile, or chr:pos); values are SNP IDs as found in markerfile
alleledict = dict()
alleleflipdict = dict()
flipper = string.maketrans('ACGT', 'TGCA')

if markerfile.endswith('.bim'):
    for markerline in markerfh:
        markerinfo = markerline.strip().split()
        if matchby == 'NAME' :
            thisname = markerinfo[1].lower()
            if thisname in snpdict :
                logging.warning('Excluding repeated instance of marker [ %s ] with non-unique name in [ %s ]', thisname, markerfile)
            else :
                snpdict[thisname] = markerinfo[1]
                alleledict[thisname] = markerinfo[4].upper() + '/' + markerinfo[5].upper()
                alleleflipdict[thisname] = markerinfo[4].upper().translate(flipper) + '/' + markerinfo[5].upper().translate(flipper)
        elif matchby == 'CHRPOS' and markerinfo[0] != '0' and markerinfo[3] != '0' :
            thischrpos = markerinfo[0] + ':' + markerinfo[3]
            if thischrpos in snpdict :
                logging.warning('Excluding marker [ %s ] because at same position [ %s ] as marker [ %s ] in [ %s ]', markerinfo[1], thischrpos, snpdict[thischrpos], markerfile)
            else :
                snpdict[thischrpos] = markerinfo[1]
                alleledict[thischrpos] = markerinfo[4].upper() + '/' + markerinfo[5].upper()
                alleleflipdict[thischrpos] = markerinfo[4].upper().translate(flipper) + '/' + markerinfo[5].upper().translate(flipper)
elif markerfile.endswith('.frq'):
    for markerline in markerfh:
        markerinfo = markerline.strip().split()
        if markerinfo[1].upper() == 'SNP': continue
        thisname = markerinfo[1].lower()
        snpdict[thisname] = markerinfo[1]
        ## .frq file does not contain position information
        alleledict[thisname] = markerinfo[2].upper() + '/' + markerinfo[3].upper()
        alleleflipdict[thisname] = markerinfo[2].upper().translate(flipper) + '/' + markerinfo[3].upper().translate(flipper)

markerfh.close()
logging.info('Successfully read [ %i ] SNP IDs from [ %s ]', len(snpdict), markerfile)

###
### read VCF file and output data matching SNPs to lookup
###

try:
    vcffh = gzip.open(vcffile, 'r')
except IOError:
    logging.error('Cannot open file [ %s ]', vcffile)
    sys.exit(1)

logging.info('Lookup in [ %s ] ...', vcffile)
gotfielddefs = 0
numread = 0;

for vcfline in vcffh:
    vcfline = vcfline.lstrip() # lstrip() faster than strip() when we don't care about right end

    if vcfline.startswith('##'):
        # we ignore meta-information lines
        continue

    if vcfline.startswith('#'):
        # field definition line
        if gotfielddefs:
                logging.error('File [ %s ] has multiple field definition lines', vcffile)
                sys.exit(1)

        fielddefs = vcfline.lstrip('#').split()
        if len(fielddefs) < 8:
                logging.error('File [ %s ] must have >=8 field definitions', vcffile)
                sys.exit(1)
        if fielddefs[0] != 'CHROM':
                logging.error('File [ %s ] field definition 1 must be CHROM', vcffile)
                sys.exit(1)
        if fielddefs[1] != 'POS':
                logging.error('File [ %s ] field definition 2 must be POS', vcffile)
                sys.exit(1)
        if fielddefs[2] != 'ID':
                logging.error('File [ %s ] field definition 3 must be ID', vcffile)
                sys.exit(1)
        if fielddefs[3] != 'REF':
                logging.error('File [ %s ] field definition 4 must be REF', vcffile)
                sys.exit(1)
        if fielddefs[4] != 'ALT':
                logging.error('File [ %s ] field definition 5 must be ALT', vcffile)
                sys.exit(1)
        if fielddefs[5] != 'QUAL':
                logging.error('File [ %s ] field definition 6 must be QUAL', vcffile)
                sys.exit(1)
        if fielddefs[6] != 'FILTER':
                logging.error('File [ %s ] field definition 7 must be FILTER', vcffile)
                sys.exit(1)
        if fielddefs[7] != 'INFO':
                logging.error('File [ %s ] field definition 8 must be INFO', vcffile)
                sys.exit(1)

        print '#SNP ALLELES STRAND ID CHROM POS REF ALT AC AN'
        gotfielddefs = 1
        continue

    if gotfielddefs == 0 :
        continue
    
    ## counter for number of variant (data) lines
    numread += 1

    fields = vcfline.split(None, 8) # 'None' means split on whitespace
    if matchby == 'NAME' :
        thisid = fields[2].lower()
    elif matchby == 'CHRPOS' :
    	if fields[0] == 'X': #modified by Li Li for chromosome X
    		fields[0] = '23'  #modified by Li Li for chromosome X
        thisid = fields[0] + ':' + fields[1]
    else :
        logging.error('Internal error matchby not NAME or CHRPOS')
        sys.exit(1)

    if thisid in snpdict:

        if mode != 'TEST' and alleledict[thisid] not in safealleles :
            ## okay that safealleles undefined when mode == 'TEST'
            ## no need to check alleleflipdict, assuming safealleles contains flipped pairs
            continue

        strandmatch = '.'
        if (alleledict[thisid] == fields[3] + '/' + fields[4]) or (alleledict[thisid] == fields[4] + '/' + fields[3]) :
            strandmatch = '+'
        elif (alleleflipdict[thisid] == fields[3] + '/' + fields[4]) or (alleleflipdict[thisid] == fields[4] + '/' + fields[3]) :
            strandmatch = '-'
            if mode == 'SOFT' :
                continue
        elif mode != 'TEST' :
            continue
        
        infos = fields[7].split(';')
        infoAC = '.'
        infoAN = '.'
        for info in infos:
            if info.startswith('AC='):
                infoAC = info[3:]
            elif info.startswith('AN='):
                infoAN = info[3:]

        print snpdict[thisid], alleledict[thisid], strandmatch, fields[2], fields[0], fields[1], fields[3], fields[4], infoAC, infoAN

   
vcffh.close()

sys.exit(0)
