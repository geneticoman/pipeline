#!/usr/bin/env python3

import sys
import os
import glob
import re

fastqdir=sys.argv[1]
outputdir=sys.argv[2]
RSEMdir =sys.argv[3]

fastqlist=glob.glob(fastqdir + "/" + "*.fastq.gz")
fastqlist.sort()
numfastq = len(fastqlist)

for i in range(0,numfastq,2):
    pathdirs1 = fastqlist[i].split("/")
    pathdirs2 = fastqlist[i+1].split("/")
    fastqname1 = pathdirs1[-1]
    fastqname2 = pathdirs2[-1]
    m = re.search('(.+?)_.*.fastq.gz',fastqname1)
    samplename = m.group(1)
    jobname = samplename + "_" + "RSEM"
    outfnameprefix = outputdir + "/" + samplename
    fname = jobname + ".sh"
    script = open(fname,"w")
    script.write("#!/bin/bash\n")
    script.write("#PBS -N %s\n" %jobname)
    script.write("#PBS -l walltime=24:00:00\n")
    script.write("#PBS -l nodes=1:ppn=10\n")
    script.write("#PBS -l mem=40GB\n")
    script.write("#PBS -d %s\n" % outputdir)
    script.write("\nmodule load gcc/6.2.0\n")
    script.write("module load RSEM\n")
    script.write("module load STAR/2.6.1d\n")
    script.write("module load samtools\n")
    script.write("rsem-calculate-expression --paired-end --strandedness reverse --star --star-gzipped-read-file --output-genome-bam --sort-bam-by-read-name -p 10 %s %s %s %s\n" % (fastqlist[i], fastqlist[i+1], RSEMdir, samplename))
