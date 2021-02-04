#!/usr/bin/env python3

import sys
import os
import glob
import re

tarlist = sys.argv[1]
outputdir = sys.argv[2]

listfile = open(tarlist, "r")
tarlines = listfile.readlines()
listfile.close()

for line in tarlines:
    line = line.strip()
    filepath = line.split("/")
    fileID = filepath[-2]
    jobname = fileID + "_untar"
    scriptname = fileID + "_untar.sh"
    script = open(scriptname, "w")
    script.write("#!/bin/bash\n")
    script.write("#PBS -N %s\n" %jobname)
    script.write("#PBS -l walltime=24:00:00\n")
    script.write("#PBS -l nodes=1\n")
    script.write("#PBS -d %s\n" % (outputdir))
    script.write("\nmodule load gcc/6.2.0\n")
    script.write("tar -xvzf %s\n" % (line))
    script.close()
