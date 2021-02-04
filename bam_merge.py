#!/usr/bin/env python3

import sys
import os
import glob
import re

bamdir = sys.argv[1]
outdir = sys.argv[2]
bedfile = sys.argv[3]
groupfile = sys.argv[4]

bf = open(bedfile, "r")
blines = bf.readlines()
bf.close()

gf = open(groupfile, "r")
glines = gf.readlines()
gf.close()

for bedline in blines:
    bedline = bedline.strip()
    bdata = bedline.split("\t")
    gene = bdata[3]
    jobname = gene + "_merge"
    outlist1 = outdir + "/" + gene + "_" + "bamlist_g2.txt"
    outlist2 = outdir + "/" + gene + "_" + "bamlist_g1.3.txt"
    outbam1 = outdir + "/" + gene + "_g2.bam"
    outbam2 = outdir + "/" + gene + "_g1.3.bam"
    scriptname = gene + ".sh"
    fout1 = open(outlist1, "w")
    fout2 = open(outlist2, "w")
    script = open(scriptname, "w")
    for groupline in glines:
        groupline = groupline.strip()
        gdata = groupline.split("\t")
        sample = gdata[0]
        group = gdata[2]
        bam = bamdir + "/" + sample + "_" + gene + ".bam"
        if group=="2":
            fout1.write("%s\n" % bam)
        elif group=="1,3":
            fout2.write("%s\n" % bam)
    fout1.close()
    fout2.close()
    script.write("#!/bin/bash\n")
    script.write("#PBS -N %s\n" %jobname)
    script.write("#PBS -l walltime=24:00:00\n")
    script.write("#PBS -l nodes=1\n")
    script.write("#PBS -d %s\n" % bamdir)
    script.write("\nmodule load gcc/6.2.0\n")
    script.write("module load bamtools\n")
    script.write("bamtools merge -list %s -out %s\n" % (outlist1, outbam1))
    script.write("bamtools merge -list %s -out %s\n" % (outlist2, outbam2))
    script.close()

