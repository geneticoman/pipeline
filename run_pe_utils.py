#!/usr/bin/env python3

import sys
import os
import glob
import re

bamdir=sys.argv[1]

bamlist=glob.glob(bamdir + "/" + "*Aligned.sortedByCoord.out.bam")

for bam in bamlist:
    bamdirs = bam.split("/")
    bamname = bamdirs[-1]
    m = re.search('(.+?)Aligned.sortedByCoord.out.bam',bamname)
    samplename = m.group(1)
    jobname = samplename + "_" + "pe_utils"
    fname = jobname + ".sh"
    outdir = samplename + "-insert-dist"
    script = open(fname,"w")
    script.write("#!/bin/bash\n")
    script.write("#PBS -N %s\n" %jobname)
    script.write("#PBS -l walltime=24:00:00\n")
    script.write("#PBS -l nodes=1\n")
    script.write("#PBS -d %s\n" % bamdir)
    script.write("\nmodule load gcc/6.2.0\n")
    script.write("module load python/2.7.13\n")
    script.write("module load samtools/1.6.0\n")
    script.write("module load bedtools\n")
    script.write("PATH=$PATH:/home/martinezc/.local/bin\n")
    script.write("export PATH\n")
    script.write("mkdir %s\n" % (outdir))
    script.write("pe_utils --compute-insert-len %s /gpfs/data/pitroda-lab/ReferenceData/anno/miso/Homo_sapiens_exons/Homo_sapiens.GRCh38.84.min_1000.const_exons.gff --output-dir %s/\n" %(bam, outdir))
