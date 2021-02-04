#!/usr/bin/env python3

import sys
import os
import glob
import re
import numpy as np

summaryfile = sys.argv[1]
m = re.search('(.*)-miso', summaryfile)
samplepath = m.group(1)
dirlist = samplepath.split("/")
samplename = dirlist[-1]

summaryin = open(summaryfile, "r")
summarydata = summaryin.readlines()
summarydata = summarydata[1:]

thresholdlist = np.arange(0.0, 1.0, 0.1)
thresholdevents = []
for threshold in thresholdlist:
    threshold = round(threshold, 1)
    numevents = 0
    for line in summarydata:
        line = line.strip()
        eventdata = line.split("\t")
        psi = float(eventdata[1])
        if psi > threshold:
            numevents = numevents + 1
    thresholdevents.append(numevents)

print(samplename, "\t".join(map(str, thresholdevents)), sep="\t")
