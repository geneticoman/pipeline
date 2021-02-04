#!/usr/bin/env python3

import sys
import os
import glob
import re

bamdir = sys.argv[1]
bamout = sys.argv[2]
bedfile = sys.argv[3]

bamlist = glob.glob(bamdir + "/" + "*Aligned.sortedByCoord.out.bam")

f = open(bedfile, "r")
lines = f.readlines()
f.close()

for bam in bamlist:
    bamdirs = bam.split("/")
    bamname = bamdirs[-1]
    m = re.search('(.+?)Aligned.sortedByCoord.out.bam', bamname)
    samplename = m.group(1)
    jobname = samplename + "_" + "extractbam"
    fname = jobname + ".sh"
    outprefix = bamout + "/" + samplename
    script = open(fname,"w")
    script.write("#!/bin/bash\n")
    script.write("#PBS -N %s\n" %jobname)
    script.write("#PBS -l walltime=24:00:00\n")
    script.write("#PBS -l nodes=1\n")
    script.write("#PBS -d %s\n" % bamdir)
    script.write("\nmodule load gcc/6.2.0\n")
    script.write("module load samtools\n")
    for line in lines:
        line = line.strip()
        data = line.split("\t")
        gene = data[3]
        loc = data[0] + ":" + data[1] + "-" + data[2]
        outname = outprefix + "_" + gene + ".bam"
        script.write("samtools view -b %s \"%s\" > %s\n" % (bam, loc, outname))
