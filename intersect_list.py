#!/usr/bin/env python3

import sys
import os
import glob
import re
from scipy.stats import hypergeom

filename1 = sys.argv[1]
filename2 = sys.argv[2]
popsize = int(sys.argv[3])

file1 = open(filename1, "r")
file2 = open(filename2, "r")

genelist1 = file1.readlines()
genelist2 = file2.readlines()
    
file1.close()
file2.close()

# make a set of all the genes in genelist2 separating comma-delimeted gene multiplets
geneset2 = set()
for gene in genelist2:
    gene = gene.strip()
    genemult = gene.split(",")
    for element in genemult:
        geneset2.add(element)

# loop through genes in genelist1 and make a set of all genes in genelist2
commonset = set()
for gene in genelist1:
    gene = gene.strip()
    genemult = gene.split(",")
    for element in genemult:
        if element in geneset2:
            commonset.add(element)
            break

commongenes = list(commonset)
commongenes.sort()

numgenes1 = len(genelist1)
numgenes2 = len(genelist2)
numcommon = len(commongenes)

# Calculate p-value of overlap
pval = hypergeom.sf(numcommon-1, popsize, numgenes1, numgenes2)
print("%d\t%d\t%d\t%g" % (numgenes1, numgenes2, numcommon, pval))
