#!/usr/bin/env python

from scipy.io import loadmat
import sys
filename = sys.argv[1]
a = loadmat(filename)
params = a['params'][0]
print "## YAML script used to generate simulation in %s ##" % filename
print '\n'.join(str(params).strip().strip('{ }').split(','))

