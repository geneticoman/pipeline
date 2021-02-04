#!/usr/bin/env python3

import sys
import os
import glob
import re

bamdir=sys.argv[1]
anno = sys.argv[2]

bamlist=glob.glob(bamdir + "/" + "*.bam")
for bam in bamlist:
    bamdirs = bam.split("/")
    bamname = bamdirs[-1]
    m = re.search('(.+?).bam',bamname)
    samplename = m.group(1)
    jobname = samplename + "_" + "htseq"
    fname = jobname + ".sh"
    outname = bamdir + "/" + samplename + ".htseq.counts"
    script = open(fname,"w")
    script.write("#!/bin/bash\n")
    script.write("#PBS -N %s\n" %jobname)
    script.write("#PBS -l walltime=24:00:00\n")
    script.write("#PBS -l nodes=1\n")
    script.write("#PBS -d %s\n" % bamdir)
    script.write("\nmodule load gcc/6.2.0\n")
    script.write("module load miniconda3\n")
    script.write("source activate htseq\n")
    script.write("htseq-count -f bam -q -m intersection-nonempty -s reverse -t exon -i gene_id %s %s > %s\n" % (bam, anno, outname))
