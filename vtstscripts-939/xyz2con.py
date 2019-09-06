#!/usr/bin/env python
import aselite
import numpy as np
from sys import argv, exit

if len(argv) < 3 or '-h' in argv:
    print "usage: xyz2con.py FILENAME BOXSIZE\n"
    exit(0)
    
filename = argv[1]
atoms = aselite.read_xyz(filename)

atoms.positions -= np.min(atoms.positions)
a = float(argv[2])
atoms.set_cell((a,a,a))


aselite.write_con(filename.replace('xyz', 'con'), atoms)
