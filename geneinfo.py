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
        if gtfdata[2] == "gene":
            anno = gtfdata[-1]
            m = re.search('gene_id "(.*?)"', anno)
            gene_id = m.group(1)
            m = re.search('gene_type "(.*?)"', anno)
            gene_type = m.group(1)
            m = re.search('gene_name "(.*?)"', anno)
            gene_name = m.group(1)
            m = re.search('level (\d)', anno)
            level = m.group(1)
            m = re.search('havana_gene "(.*?)"', anno)
            if not m:
                havana_gene = "-"
            else:
                havana_gene = m.group(1)
            print(gene_id, gene_type, gene_name, level, havana_gene, sep = "\t")
            
