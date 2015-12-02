#!/usr/bin/python
import gzip
import logging
import time
import string
import sys

logging.basicConfig(stream=sys.stdout, format='%(levelname)s:%(message)s', level=logging.DEBUG)
logging.info('ped/map consensus-call merge tool v0.3 Toby.x.Johnson@gsk.com')

###
### parse command line arguments
###

if len(sys.argv) == 3:
    inroot = sys.argv[1]
    outroot = sys.argv[2]
else:
    logging.error('Need two arguments: INROOT, OUTROOT (to which .map and .ped extensions will be added')
    sys.exit(2)

try:
    mapfh = open(inroot + '.map', 'r')
except IOError:
    logging.error('Cannot open file [ %s.map ]', inroot)
    sys.exit(1)

logging.info('Analysis started: %s', time.asctime())

logging.info('Reading map information from [ %s.map ] ...', inroot)
snpchr = list()
snpname = list()
snpbp = list()
snppos = list()
for mapline in mapfh:
    mapinfo = mapline.strip().split()
    if len(mapinfo) != 4:
        logging.error('Expecting 4 items but found [ %i ] in map file [ %s.map ] line [ %i ]', len(mapinfo), inroot, len(snppos) + 1)
        sys.exit(1)
    snpchr.append(mapinfo[0])
    snpname.append(mapinfo[1])
    snpbp.append(mapinfo[3])
    snppos.append(mapinfo[0] + ":" + mapinfo[3])
logging.info('Read map information for [ %i ] SNPs from [ %s.map ]', len(snppos), inroot)

mapfh.close()

logging.info('Writing new map information to [ %s.map ] ...', outroot)
try:
    omapfh = open(outroot + '.map', 'w')
except IOError:
    logging.error('Cannot write to file [ %s.map ]', outroot)
    sys.exit(1)

snpdupcount = list() # for each SNP, 0 means merge to later SNP, >0 means merge with this number of previous SNP(s)
snpcname = list() # if snpdupcount>0, the name for the merged group
snpgroup = list() # comma separated names of all merged SNPs
thisname = ''
thisgroup = ''
dupcount = 0
for i in range(len(snppos)):
    if thisgroup == '':
        thisgroup = snpname[i]
    else:
        thisgroup += ','+snpname[i]
    if thisname == '' or (snpname[i][0:2] == 'rs' and thisname[0:2] != 'rs'):
        thisname = snpname[i]
    dupcount += 1
    if i == len(snppos)-1 or snppos[i+1] != snppos[i]:
        snpdupcount.append(dupcount) # guaranteed dupcount >= 1
        snpcname.append(thisname)
        snpgroup.append(thisgroup)
        omapfh.write(snpchr[i] + ' ' + thisname + ' 0 ' + snpbp[i] + '\n') ## output all info...
        thisname = ''
        thisgroup = ''
        dupcount = 0
    else:
        snpdupcount.append(0)
        snpcname.append("<NA>") # this should never be printed
        snpgroup.append("<NA>") # this should never be printed

omapfh.close()
logging.info('Wrote new map information to [ %s.map ]', outroot)

numdup = dict()
for i in range(len(snppos)):
    if snpdupcount[i]:
        if snpdupcount[i] in numdup:
            numdup[snpdupcount[i]] += 1
        else:
            numdup[snpdupcount[i]] = 1
for nd, i in numdup.iteritems():
    if nd > 1:
        addstring = 's to be consensus-merged'
    else:
        addstring = ''
    logging.info('There are [ %i ] unique markers with %i genotype%s', i, nd, addstring)

countccall = [0]*len(snppos)
countcall = [0]*len(snppos)

logging.info('Reading genotypes from [ %s.ped ] ...', inroot)
try:
    pedfh = open(inroot + '.ped', 'r')
except IOError:
    logging.error('Cannot open file [ %s.ped ]', inroot)
    sys.exit(1)

logging.info('Writing genotypes to [ %s.ped ] ...', outroot)
try:
    opedfh = open(outroot + '.ped', 'w')
except IOError:
    logging.error('Cannot write to file [ %s.ped ]', outroot)
    sys.exit(1)

numindiv = 0
for pedline in pedfh:
    pedinfo = pedline.strip().split()
    if len(pedinfo) != 6+2*len(snppos):
        logging.error('Expecting [ %i = 6+2*%i ] items but found [ %i ] in ped file [ %s.ped ] line [ %i ]', 6+2*len(snppos), len(snppos), len(pedinfo), inroot, numindiv + 1)
        sys.exit(1)
    numindiv += 1

    opedfh.write(pedinfo[0] + ' ' + pedinfo[1] + ' ' + pedinfo[2] + ' ' + pedinfo[3] + ' ' + pedinfo[4] + ' ' + pedinfo[5])

    numcall = 0
    callset = dict()
    for i in range(len(snppos)):
        if snpdupcount[i] == 1: # handle this quickly
            opedfh.write(' ' + pedinfo[2*i+6] + ' ' + pedinfo[2*i+7])
        else:
            a1 = pedinfo[2*i+6]
            a2 = pedinfo[2*i+7]
            if a1 == '0' and a2 == '0':
                pass # do nothing
            else:
                if a1 > a2:
                    call = a2 + ' ' + a1
                else:
                    call = a1 + ' ' + a2
                numcall += 1
                if call in callset:
                    callset[call] += 1
                else:
                    callset[call] = 1
            if snpdupcount[i] == 0:
                pass # do nothing
            else:
                ccall = '0 0'
                for call, numcons in callset.iteritems():
                    if numcons == numcall: # require all noncalls to be the same
                        ccall = call
                        countccall[i] += numcons
                        # could add break statement here
                countcall[i] += numcall
                opedfh.write(' '  + ccall)
                numcall = 0
                callset = dict()
            
    opedfh.write('\n')
            
logging.info('Read genotypes for [ %i ] individuals', numindiv)
logging.info('Wrote genotypes for [ %i ] individuals', numindiv)

pedfh.close()
opedfh.close()

logging.info('Writing consensus information to [ %s.info ] ...', outroot)
totcall = 0
totccall = 0
try:
    oinfofh = open(outroot + '.info', 'w')
except IOError:
    logging.error('Cannot write to file [ %s.info ]', outroot)
    sys.exit(1)

oinfofh.write('SNP CHR BP GROUP CALLS CONSENSUS F_MISMATCH\n')
for i in range(len(snppos)):
    if snpdupcount[i] > 1:
        if countcall[i] > 0:
            fmismatch = str(round(float(countcall[i] - countccall[i])/float(countcall[i]), 6))
        else:
            fmismatch = 'NA'
        totcall += countcall[i]
        totccall += countccall[i]
        oinfofh.write(snpcname[i] + ' ' + snpchr[i] + ' ' + snpbp[i] + ' ' + snpgroup[i] + ' ' + str(countcall[i]) + ' ' + str(countccall[i]) + ' ' + fmismatch + '\n')

oinfofh.close()
logging.info('Wrote consensus information to [ %s.info ]', outroot)
if totcall > 0:
    fmismatch = str(round(float(totcall - totccall)/float(totcall), 6))
else:
    fmismatch = 'NA'
logging.info('Overall mismatch rate [ %s ]', fmismatch)

logging.info('Analysis finished: %s', time.asctime())

sys.exit(0)
