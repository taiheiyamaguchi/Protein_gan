#!/usr/bin/python3

import os

for i in range(1000):
  fname='final.pdb'
  if(not os.path.exists(fname)):
    break
  fin=open(fname)
  fout=open('final_mirror.pdb', 'w')
  for line in fin:
    z=float(line[46:54])
    fout.write(line[0:46])
    fout.write('%8.3f' % (-z))
    fout.write(line[54:])
  fin.close()
  fout.close()
