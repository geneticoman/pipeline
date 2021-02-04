#!/usr/bin/env python3

import sys
import os
import glob
import re

fastqdir=sys.argv[1]
outputdir=sys.argv[2]
indexfile =sys.argv[3]

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
    jobname = samplename + "_" + "kallisto"
    sampleout = outputdir + "/" + samplename
    fname = jobname + ".sh"
    script = open(fname,"w")
    script.write("#!/bin/bash\n")
    script.write("#PBS -N %s\n" %jobname)
    script.write("#PBS -l walltime=24:00:00\n")
    script.write("#PBS -l nodes=1:ppn=10\n")
    script.write("#PBS -l mem=40GB\n")
    script.write("#PBS -d %s\n" % outputdir)
    script.write("\nmodule load gcc/6.2.0\n")
    script.write("module load kallisto\n")
    script.write("cd %s\n" % outputdir)
    script.write("mkdir %s\n" % sampleout)
    script.write("mkdir %s/kallisto\n" % sampleout)
    #script.write("kallisto quant -i %s --fr-stranded -t 10 -b 10 -o %s %s %s\n" % (indexfile, sampleout, fastqlist[i], fastqlist[i+1]))
    script.write("kallisto quant -i %s --rf-stranded -t 10 -b 10 -o %s/kallisto %s %s\n" % (indexfile, sampleout, fastqlist[i], fastqlist[i+1]))
    
