#!/usr/bin/env python3

import sys
import os

fastqdir=sys.argv[1]
outputdir=sys.argv[2]

fastqlist=os.listdir(fastqdir)
fastqlist.sort()
numlist = len(fastqlist)
for i in range(0,numlist,2):
    fqname1 = fastqlist[i]
    fqname2 = fastqlist[i+1]
    fname = fqname1.replace("_R1_merge.fastq.gz","")
    jobname = fname + "_" + "trim"
    fname = fname + "_" + "trim.sh"
    script = open(fname,"w")
    script.write("#!/bin/bash\n")
    script.write("#PBS -N %s\n" %jobname)
    script.write("#PBS -d %s\n" %fastqdir)
    script.write("#PBS -l walltime=12:00:00\n")
    script.write("#PBS -l nodes=1\n")
    script.write("\nmodule load gcc/6.2.0\n")
    script.write("module load trimgalore/0.4.5\n")
    script.write("trim_galore -q 30 --length 20 --paired -o %s %s %s\n" % (outputdir, fqname1, fqname2))
