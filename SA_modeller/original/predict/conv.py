#!/usr/bin/python3

import numpy as np

cut=22.0

mat=np.load('input.npy')
num=mat.shape[0]
fp=open('rest.dat','w')
for i in range(num):
  for j in range(i+1,num):
    d=mat[i,j]
    if(d < cut):
      fp.write('%d %d %.1f %.1f 1.0 1.0\n' % (i,j,d-0.5,d+0.5))
    else:
      fp.write('%d %d %.1f %.1f 1.0 1.0\n' % (i,j,cut,1000.0))
fp.close()

