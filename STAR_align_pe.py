#!/usr/bin/env python3

import sys
import os
import glob
import re

fastqdir=sys.argv[1]
outputdir=sys.argv[2]
STARdir =sys.argv[3]

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
    jobname = samplename + "_" + "STAR"
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
    script.write("module load STAR/2.6.1d\n")
    #script.write("STAR --runMode alignReads --genomeDir /gpfs/data/pitroda-lab/ReferenceData/STAR_indexes/hg38_gencode.v29_76 --runThreadN 10 --readFilesIn %s %s --readFilesCommand zcat -c --outFileNamePrefix %s --outSAMtype BAM Unsorted --chimSegmentMin 20 --quantMode TranscriptomeSAM --outReadsUnmapped Fastq --outMultimapperOrder Random --outSAMattrRGline ID:%s --outFilterMultimapNmax 20 --outFilterMismatchNmax 2\n" % (fastqlist[i], fastqlist[i+1], outfnameprefix, samplename))
    script.write("STAR --runMode alignReads --genomeDir %s --runThreadN 10 --readFilesIn %s %s --readFilesCommand zcat -c --outFileNamePrefix %s --outSAMtype BAM Unsorted --chimSegmentMin 20 --outReadsUnmapped Fastq --outMultimapperOrder Random --outSAMattrRGline ID:%s --outFilterMultimapNmax 20 --outFilterMismatchNmax 2\n" % (STARdir, fastqlist[i], fastqlist[i+1], outfnameprefix, samplename))
