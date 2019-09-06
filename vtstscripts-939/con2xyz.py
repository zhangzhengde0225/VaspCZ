#!/usr/bin/env python

import aselite

from sys import argv, exit

if len(argv) < 2 or '-h' in argv:
    print "usage: con2xyz.py FILENAME\n"
    exit(0)
    
filename = argv[1]
atoms = aselite.read_con(filename)

aselite.write_xyz(filename.replace('con', 'xyz'), atoms)
