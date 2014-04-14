
import statePlanar as sp
import itertools
import sys
import numpy as np
import copy
import os


d = 5
error_type = "xyz" # "y" or "xyz"


qend = 0.01
pend = 0.15
nsteps = 20
home = os.environ['HOME']
results_file = "%s/precalculatedDecoder/%doutput_%s_lies3.txt"%(home,d,error_type)


pmultiplier = pend/nsteps

for p in [pmultiplier*i for i in range(20)]:
    

    f = open("fsub2","w")
    f.write("#!/bin/sh\n#PBS -l walltime=1:00:00\n#PBS -l mem=1000mb\n#PBS -l ncpus=1\n\n")
    f.write("python $HOME/precalculatedDecoder/part_compute_lies.py %d %s %f %f %d %s"%(d,error_type,p,qend,nsteps,results_file))
    f.close()

    os.system("qsub fsub2")


