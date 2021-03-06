#!/usr/bin/env python3

import sys
import os
import glob
import re

fastqdir1=sys.argv[1]
fastqdir2=sys.argv[2]
outputdir=sys.argv[3]

fastqlist1=glob.glob(fastqdir1 + "/" + "*.fastq.gz")
fastqlist2=glob.glob(fastqdir2 + "/" + "*.fastq.gz")
fastqlist1.sort()
fastqlist2.sort()

i = 0
for i in range(0,len(fastqlist1)):
    dirlist1 = fastqlist1[i].split("/")
    fastq1 = dirlist1[-1]
    m = re.search('(.*)_001',fastq1)
    fastq = ""
    if m:
        fastq = m.group(1)
    jobname = fastq + "_merge"
    fname = fastq + "_merge" + ".sh"
    outname = fastq + "_merge" + ".fastq.gz"
    script = open(fname,"w")
    script.write("#!/bin/bash\n")
    script.write("#PBS -N %s\n" %jobname)
    script.write("#PBS -l walltime=24:00:00\n")
    script.write("#PBS -l nodes=1\n")
    script.write("#PBS -d %s\n" % (outputdir))
    script.write("\nmodule load gcc/6.2.0\n")
    script.write("gunzip -c %s %s | gzip > %s\n" % (fastqlist1[i],fastqlist2[i],outname))
