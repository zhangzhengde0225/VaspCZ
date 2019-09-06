#!/usr/bin/env python

import aselite

from sys import argv, exit

if len(argv) < 3:
    print 'usage: center.py FILE STDDEV'
    print '       randomly displaces the atoms in FILE by a gaussian with' 
    print '       a standard deviation of STDDEV' 
    print
    exit(0)
    
filename = argv[1]
stddev = float(argv[2])
atoms = aselite.read_any(filename)
atoms.rattle(stddev)

atoms.write(filename)
