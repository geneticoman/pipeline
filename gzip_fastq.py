#!/usr/bin/env python3

import sys
import os
import glob
import re

fastqdir=sys.argv[1]
fastqlist=glob.glob(fastqdir + "/" + "*.fastq")

for fastq in fastqlist:
    dirlist = fastq.split("/")
    fastqname = dirlist[-1]
    m = re.search('(.*).fastq',fastqname)
    sample = ""
    if m:
        sample = m.group(1)
    jobname = sample + "_gzip"
    fname = sample + "_gzip" + ".sh"
    script = open(fname,"w")
    script.write("#!/bin/bash\n")
    script.write("#PBS -N %s\n" %jobname)
    script.write("#PBS -l walltime=24:00:00\n")
    script.write("#PBS -l nodes=1\n")
    script.write("#PBS -d %s\n" % (fastqdir))
    script.write("\nmodule load gcc/6.2.0\n")
    script.write("gzip %s\n" % (fastq))
