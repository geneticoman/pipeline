#!/usr/bin/env python3

import sys
import os
import glob
import re

misoanno = sys.argv[1]
misosummary = sys.argv[2]
misoout = misosummary + ".ugenes"

# read miso annoation file and make a dictionary of events with gene symbol
misof = open(misoanno, "r")
altevents = misof.readlines()
misof.close()
annotations = {}
for event in altevents:
    eventfields = event.split("\t")
    if eventfields[2] == "gene":
        data = eventfields[-1]
        m = re.search('gsymbol=(.+?);', data)
        gsymbol = m.group(1)
        m = re.search('gid=(.+?);', data)
        gid = m.group(1)
        annotations[gid] = gsymbol

# read miso summary and loop over events extracting event name and printing corresponding gene symbol
summaryf = open(misosummary, "r")
summaryevents = summaryf.readlines()
summaryf.close()
summaryevents = summaryevents[1:]
misogenes = []
for sevent in summaryevents:
    seventfields = sevent.split("\t")
    gid = seventfields[0]
    gsymbol = annotations[gid]
    misogenes.append(gsymbol)

# sort found genes and find unique
geneset = set(misogenes)
uniquegenes = list(geneset)
uniquegenes.sort()

# output unique genes
outfile = open(misoout, "w")
for gene in uniquegenes:
    outfile.write("%s\n" %gene)
outfile.close()
