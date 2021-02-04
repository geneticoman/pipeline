#!/usr/bin/env python3

import sys
import os
import glob
import re

bamdir = sys.argv[1]
settings = sys.argv[2]
readlen = int(sys.argv[3])
#sample_readlen = sys.argv[3]

bamlist = glob.glob(bamdir + "/" + "*Aligned.sortedByCoord.out.bam")
misoindexdir = "/gpfs/data/pitroda-lab/ReferenceData/anno/miso/indexed"

#readlenf = open(sample_readlen,"r")
#readlendat = readlenf.readlines()
#readlengths = {}
#for line in readlendat:
#    line_array = line.split()
#    readlengths[line_array[0]] = int(line_array[-1]) - 1

for bam in bamlist:
    bamdirs = bam.split("/")
    bamname = bamdirs[-1]
    m = re.search('(.+?)Aligned.sortedByCoord.out.bam',bamname)
    samplename = m.group(1)
    jobname = samplename + "_" + "miso"
    fname = jobname + ".sh"
    outdir = samplename + "-miso"
    inlendir = bamdir + "/" + samplename + "-insert-dist/"
    #inlenfn = inlendir + bamname + ".insert_len"
    inlenfn = inlendir + samplename + "Aligned.sortedByCoord.out.bam.insert_len"
    inlenf = open(inlenfn,"r")
    inlendata = inlenf.readlines()
    inlen1 = inlendata[0]
    inlen1 = inlen1[1:]
    inlena = inlen1.split(",")
    m = re.search('mean=(.+)', inlena[0])
    mean = float(m.group(1))
    m = re.search('sdev=(.+)', inlena[1])
    sdev = float(m.group(1))
    inlenf.close()
#    readlen = readlengths[samplename]
    script = open(fname,"w")
    script.write("#!/bin/bash\n")
    script.write("#PBS -N %s\n" %jobname)
    script.write("#PBS -l walltime=168:00:00\n")
    script.write("#PBS -l nodes=1\n")
    script.write("#PBS -l mem=16GB\n")
    script.write("#PBS -d %s\n" % bamdir)
    script.write("\nmodule load gcc/6.2.0\n")
    script.write("module load python/2.7.13\n")
    script.write("module load samtools/1.6.0\n")
    script.write("module load bedtools\n")
    script.write("PATH=$PATH:/home/martinezc/.local/bin\n")
    script.write("export PATH\n")
    script.write("mkdir %s\n" % (outdir))
    script.write("miso --run %s/SE %s --output-dir %s/SE --read-len %d --paired-end %f %f --settings-filename=%s\n" %(misoindexdir, bamname, outdir, readlen, mean, sdev, settings))
    script.write("miso --run %s/RI %s --output-dir %s/RI --read-len %d --paired-end %f %f --settings-filename=%s\n" %(misoindexdir, bamname, outdir, readlen, mean, sdev, settings))
    script.write("miso --run %s/A3SS %s --output-dir %s/A3SS --read-len %d --paired-end %f %f --settings-filename=%s\n" %(misoindexdir, bamname, outdir, readlen, mean, sdev, settings))
    script.write("miso --run %s/A5SS %s --output-dir %s/A5SS --read-len %d --paired-end %f %f --settings-filename=%s\n" %(misoindexdir, bamname, outdir, readlen, mean, sdev, settings))
    script.write("miso --run %s/MXE %s --output-dir %s/MXE --read-len %d --paired-end %f %f --settings-filename=%s\n" %(misoindexdir, bamname, outdir, readlen, mean, sdev, settings))
    script.write("summarize_miso --summarize-samples %s/SE %s/SE\n" %(outdir, outdir))
    script.write("summarize_miso --summarize-samples %s/RI %s/RI\n" %(outdir, outdir))
    script.write("summarize_miso --summarize-samples %s/A3SS %s/A3SS\n" %(outdir, outdir))
    script.write("summarize_miso --summarize-samples %s/A5SS %s/A5SS\n" %(outdir, outdir))
    script.write("summarize_miso --summarize-samples %s/MXE %s/MXE\n" %(outdir, outdir))
    script.close()
