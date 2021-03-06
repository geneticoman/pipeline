#!/usr/bin/env python3

import sys
import os
import glob
import re

sradir = sys.argv[1]
fastqdir = sys.argv[2]

sralist = glob.glob(sradir + "/" + "*.sra")
sralist.sort()

for sra in sralist:
    pathdirs = sra.split("/")
    sraname = pathdirs[-1]
    m = re.search('(.+?).sra',sraname)
    samplename = m.group(1)
    jobname = samplename + "_" + "fastq-dump"
    fname = jobname + ".sh"
    script = open(fname,"w")
    script.write("#!/bin/bash\n")
    script.write("#PBS -N %s\n" %jobname)
    script.write("#PBS -l walltime=24:00:00\n")    
    script.write("#PBS -l nodes=1\n")
    script.write("#PBS -d %s\n" % fastqdir)
    script.write("\nmodule load gcc/6.2.0\n")
    script.write("module load sra-tools\n")
    script.write("fastq-dump --gzip --skip-technical --readids --dumpbase --split-files --clip %s\n" % samplename)

    
