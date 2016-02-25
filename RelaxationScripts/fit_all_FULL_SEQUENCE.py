#!/bin/python
# Do fitting for M samples of K columns from correlation data files (???.txt) 
# For residues in RESIDS, and calculate appropriate R1,R2,R2/R1 for each sample
# using fit3.sh script

##### IMPORTS #####
import numpy.random as nprnd
import glob
import subprocess
import multiprocessing

###### PARAMS ######
M=20 # number of samples
K=50 # number of columns to sample
RESIDS=range(251,376,1) # resids to sample from
GENERATE_COL=1

##### METHODS #####
def call_fit(resid,columns): 
  # call fit3 for resid.txt file, with columns specified as a list of integers
  cmd="bash ./fit7_intervals.sh %d.txt %d %d %d %d" % (resid, columns[0], columns[1],columns[2], columns[3])
  print cmd
  subprocess.call(cmd,shell=True);

##### Main: #####
# Save list of random columns:
if GENERATE_COL:
  C=[nprnd.randint(2,11,size=K) for r in xrange(M)];
  F=open('COLUMNS','w');
  for columns in C:
    for c in columns:
      print >>F, "%3d" % c,
    print >>F
  F.close() 
# Process columns for each specified residue
pool=multiprocessing.Pool(processes=2)
F=open('COLUMNS','r');
for line in F:
  columns=[int(x) for x in line.split()]   # convert to int to verify proper input
  for resid in RESIDS:
    pool.apply_async(call_fit, [resid,columns])
#    pool.apply(call_fit, [resid,columns])
pool.close()
pool.join()
F.close()
