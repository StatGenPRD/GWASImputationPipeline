#!/usr/bin/python
import gzip
import logging
import string
import sys

logging.basicConfig(format='%(levelname)s:%(message)s') # can add level=logging.DEBUG

## pipe tool to filter based on one column matching big dictionary, v0.2 Toby.x.Johnson@gsk.com
## v0.2: input lines with <COLUMN whitespace separated tokens are skipped

###
### parse command line arguments
###

if len(sys.argv) == 3:
    columnnum = sys.argv[1]
    includefile = sys.argv[2]
else:
    logging.error('Need two arguments: COLUMN [number >=1] and INCLUDEFILE [filename]')
    sys.exit(1)

try:
    columnnum = int(columnnum) - 1
except ValueError:
    logging.error('First argument (COLUMN) must be integer')
    sys.exit(1)

if columnnum<0 :    
    logging.error('First argument (COLUMN) must be positive integer')
    sys.exit(1)

###
### build dictionary
### python guarantees(?) lookups average complexity O(1)    
###   

try:
    includefh = open(includefile, 'r')
except IOError:
    logging.error('Cannot read file [ %s ]', includefile)
    sys.exit(1)
   
includedict = dict()
for includeline in includefh:
    includedict[includeline.strip()] = 1

includefh.close()

###
### scan stdin for matches and print
###

for inputline in sys.stdin:
    inputitems = inputline.strip().split(None, columnnum + 1)
    if len(inputitems) > columnnum and inputitems[columnnum] in includedict:
        print inputline.strip()

###
### done
###

sys.exit(0)
