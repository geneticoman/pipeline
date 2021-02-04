#!/usr/bin/env python3

import sys
import os
import glob
import re

bamdir=sys.argv[1]

bamlist=glob.glob(bamdir + "/" + "*.bam")

for bam in bamlist:
    bamdirs = bam.split("/")
    bamname = bamdirs[-1]
    #m = re.search('(.+?)Aligned.sortedByCoord.out.bam',bamname)
    m = re.search('(.+?)Aligned.out.bam',bamname)
    samplename = m.group(1)
    jobname = samplename + "_" + "sort_name"
    fname = jobname + ".sh"
    outname = bamdir + "/" + samplename + ".sorted.names.bam"
    script = open(fname,"w")
    script.write("#!/bin/bash\n")
    script.write("#PBS -N %s\n" %jobname)
    script.write("#PBS -l walltime=24:00:00\n")
    script.write("#PBS -l nodes=1\n")
    script.write("#PBS -d %s\n" % bamdir)
    script.write("\nmodule load gcc/6.2.0\n")
    script.write("module load samtools\n")   
    script.write("samtools sort -n -T %s -o %s %s\n" % (samplename, outname, bam))
