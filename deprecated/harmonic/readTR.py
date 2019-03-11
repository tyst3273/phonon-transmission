#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 27 21:12:11 2018

@author: ty
"""
import numpy as np

infile = 'tr.dat'

with open(infile, 'r') as fid:
    num_lines = sum(1 for line in fid)-6 #get number of line in file
    fid.seek(0) #return cursor to start of file
    
    data = np.zeros((num_lines,15)).astype(complex)
    
    for i in range(4):
        fid.readline() #skip comments
    
    columns = fid.readline().strip().split()
    fid.readline() #skip comment
    
    for i in range(num_lines):
        tmp = fid.readline().strip().split()
        for j in range(15):
            data[i,j] = complex(tmp[j])
            
del infile, i, j, num_lines, tmp
    
    
