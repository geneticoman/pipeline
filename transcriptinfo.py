#!/usr/bin/env python3

import sys
import os
import glob
import re

gtfanno = sys.argv[1]

with open(gtfanno, "r") as gtf_file:
    for line in gtf_file:
        if line[0] == "#":
            continue
        line = line.strip()
        gtfdata = line.split("\t")
        if gtfdata[2] == "transcript":
            anno = gtfdata[-1]
            #print(anno)
            m = re.search('gene_id "(.*?)"', anno)
            gene_id = m.group(1)
            m = re.search('transcript_id "(.*?)"', anno)
            transcript_id = m.group(1)
            m = re.search('gene_type "(.*?)"', anno)
            gene_type = m.group(1)
            m = re.search('gene_name "(.*?)"', anno)
            gene_name = m.group(1)
            m = re.search('transcript_type "(.*?)"', anno)
            transcript_type = m.group(1)
            m = re.search('transcript_name "(.*?)"', anno)
            transcript_name = m.group(1)
            m = re.search('level (\d)', anno)
            level = m.group(1)
            m = re.search('protein_id "(.*?)"', anno)
            if not m:
                protein_id = "-"
            else:
                protein_id = m.group(1)
            print(gene_id, transcript_id, gene_type, gene_name, transcript_type, transcript_name, level, protein_id, sep = "\t")
            
