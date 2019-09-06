#!/usr/bin/env python
import aselite
from sys import argv

if '-h' in argv or len(argv) != 2:
    print 'usage: pos2xyz.py POSCAR'
    print
    exit(1)

atoms = aselite.read_vasp(argv[1])
aselite.write_xyz('%s.xyz' % argv[1], atoms)
