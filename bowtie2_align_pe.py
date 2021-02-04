#!/usr/bin/env python3

import sys
import os
import glob
import re

fastqdir=sys.argv[1]
outputdir=sys.argv[2]

fastqlist=glob.glob(fastqdir + "/" + "*.fq.gz")
fastqlist.sort()
numfastq = len(fastqlist)

for i in range(0,numfastq,2):
    pathdirs1 = fastqlist[i].split("/")
    pathdirs2 = fastqlist[i+1].split("/")
    fastqname1 = pathdirs1[-1]
    fastqname2 = pathdirs2[-1]
    m = re.search('(.+?)_.*.fq.gz',fastqname1)
    samplename = m.group(1)
    jobname = samplename + "_" + "bowtie2"
    outfname = outputdir + "/" + samplename + ".sam"
    outbam = outputdir + "/" + samplename + ".bam"
    outsorted = outputdir + "/" + samplename + "_sorted.bam"
    fname = jobname + ".sh"
    script = open(fname,"w")
    script.write("#!/bin/bash\n")
    script.write("#PBS -N %s\n" %jobname)
    script.write("#PBS -l walltime=24:00:00\n")
    script.write("#PBS -l nodes=1:ppn=10\n")
    script.write("#PBS -l mem=16GB\n")
    script.write("#PBS -d %s\n" % outputdir)
    script.write("\nmodule load gcc/6.2.0\n")
    script.write("module load bowtie2\n")
    script.write("bowtie2 -x /gpfs/data/pitroda-lab/ReferenceData/Bowtie2/hg38 -1 %s -2 %s -S %s -p 10\n" % (fastqlist[i], fastqlist[i+1], outfname))
    script.write("\nmodule load samtools\n")
    script.write("samtools view -Sb %s > %s\n" % (outfname, outbam))
    script.write("samtools sort %s > %s\n" % (outbam, outsorted))
    script.write("samtools index %s\n" % outsorted)
