#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 15 16:39:11 2019

@author: ty
"""
from timeit import default_timer as timer
import numpy as np
import scipy.linalg.blas as blas

a = np.array(np.random.randint(0,100,size=(100,5000)),order='C')
b = np.array(np.random.randint(0,100,size=(5000,100)),order='C')
c = np.array(np.random.randint(0,100,size=(100,100)),order='C')

start = timer()
d = np.dot(b,np.dot(c,a))
stop = timer()
print((stop-start))

start = timer()
f = blas.dgemm(1,b,blas.dgemm(1,c,a))
stop = timer()
print((stop-start))
