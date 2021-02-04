#!/usr/bin/env python3

import sys
import os
import glob
import re

fastqdir=sys.argv[1]
outputdir=sys.argv[2]
IRFinder_index =sys.argv[3]

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
    jobname = samplename + "_" + "IRFinder"
    outfnameprefix = outputdir + "/" + samplename
    fname = jobname + ".sh"
    script = open(fname,"w")
    script.write("#!/bin/bash\n")
    script.write("#PBS -N %s\n" %jobname)
    script.write("#PBS -l walltime=24:00:00\n")
    script.write("#PBS -l nodes=1:ppn=10\n")
    script.write("#PBS -l mem=50GB\n")
    script.write("#PBS -d %s\n" % outputdir)
    script.write("\nmodule load gcc/6.2.0\n")
    script.write("module load STAR/2.6.1a\n")
    script.write("module load bedtools\n")
    script.write("module load pigz\n")
    script.write("module unload perl\n")
    script.write("cd %s\n" % outputdir)
    script.write("mkdir %s\n" % samplename)
    script.write("/home/martinezc/Github/IRFinder/bin/IRFinder -m FastQ -t 10 -d %s -s NoSharedMemory -r %s %s %s\n" % (outfnameprefix, IRFinder_index, fastqlist[i], fastqlist[i+1]))
