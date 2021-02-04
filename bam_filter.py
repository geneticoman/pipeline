#!/usr/bin/env python3

import sys
import os
import glob
import re

bamdir = sys.argv[1]
length = int(sys.argv[2])

bamlist = glob.glob(bamdir + "/" + "*Aligned.sortedByCoord.out.bam")

for bam in bamlist:
    bamdirs = bam.split("/")
    bamname = bamdirs[-1]
    m = re.search('(.+?)Aligned.sortedByCoord.out.bam', bamname)
    samplename = m.group(1)
    jobname = samplename + "_" + "filterbam"
    fname = jobname + ".sh"
    outname = bamdir + "/" + samplename + "Aligned.sortedByCoord.out.filter.bam"
    script = open(fname,"w")
    script.write("#!/bin/bash\n")
    script.write("#PBS -N %s\n" %jobname)
    script.write("#PBS -l walltime=24:00:00\n")
    script.write("#PBS -l nodes=1\n")
    script.write("#PBS -d %s\n" % bamdir)
    script.write("\nmodule load gcc/6.2.0\n")
    script.write("module load samtools\n")
    script.write("samtools view -h %s | awk 'length($10)==%d || $1 ~ /^@/' | samtools view -bS - > %s\n" % (bam, length, outname))
    script.write("samtools index %s\n" % (outname))
