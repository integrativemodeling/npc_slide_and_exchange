#!/bin/python
# Do fitting for M samples of K columns from correlation data files (???.txt) 
# For residues in RESIDS, and calculate appropriate R1,R2,R2/R1 for each sample
# using fit3.sh script, and based on samples in COLUMNS file (see below)
#
# If GENERATE_COL is 1, also generates a COLUMNS file that holds the columns to sample from

##### IMPORTS #####
import numpy.random as nprnd
import glob
import subprocess
import multiprocessing

###### PARAMS ######
M=50 # number of samples
K=4 # number of columns to sample
RESIDS=range(1,20,1) # resids to sample from
GENERATE_COL=1 # if zero, use an old COLUMN file

##### METHODS #####
def call_fit(resid,columns): 
  # call fit3 for resid.txt file, with columns specified as a list of integers
  cmd="bash ./fit7_intervals.sh F%d.txt %d %d %d %d" % (resid, columns[0], columns[1],columns[2], columns[3])
  print cmd
  subprocess.call(cmd,shell=True);

##### Main: #####
if __name__=='__main__':
  # Make F files - averaged over all six repeats in each of the time interaval
  cmd='for k in `seq 1 19`; do paste $((251+$k)).txt $((270+$k)).txt $((289+$k)).txt $((308+$k)).txt $((327+$k)).txt $((346+$k)).txt | awk \'{printf("%.2f ",$1); for(i=2;i<=10;i++){printf "%.4f ",(($(i)+$(i+10)+$(i+20)+$(i+30)+$(i+40)+$(i+50))/6.0)} print ""}\' > F$k.txt; done'
  subprocess.call(cmd,shell=True)

  # Save list of random columns:
  if GENERATE_COL:
    C=[(nprnd.choice(9,K,replace=False)+2) for r in xrange(M)];
    F=open('COLUMNS','w');
    for columns in C:
      for c in columns:
        print >>F, "%3d" % c,
      print >>F
    F.close() 
  # Process columns for each specified residue
  pool=multiprocessing.Pool(processes=1)
  F=open('COLUMNS','r');
  for line in F:
    columns=[int(x) for x in line.split()]   # convert to int to verify proper input
    for resid in RESIDS:
      pool.apply_async(call_fit, [resid,columns])
  pool.close()
  pool.join()
  F.close()
