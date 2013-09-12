#!/usr/bin/env python

from scipy.io import loadmat
import sys
a = loadmat(sys.argv[1])
params = a['params'][0]
print '\n'.join(str(params).strip().strip('{ }').split(','))

