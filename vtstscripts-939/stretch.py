#!/usr/bin/env python
import aselite
from sys import exit, argv
from os.path import isfile
import numpy as np

if '-h' in argv or len(argv) < 5:
    print 'usage: stretch.py FILENAME INDEX1 INDEX2 DISTANCE'
    print '       Stretches the atom INDEX2 by the bond defined by atoms INDEX1'
    print '       and INDEX2 by DISTANCE. Saves the result to FILENAME_stretch.'
    print
    exit(0)

filename = argv[1]
index1 = int(argv[2])
index2 = int(argv[3])
distance = float(argv[4])

if not isfile(filename):
    print 'No such file: %s' % filename
    exit(1)

atoms = aselite.read_any(filename)
r = atoms.get_positions()
bond = r[index2] - r[index1]
bond /= np.linalg.norm(bond)

r[index2] += distance*bond
atoms.set_positions(r)

atoms.write(filename+'_stretch')
