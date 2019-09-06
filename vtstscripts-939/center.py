#!/usr/bin/env python
import aselite
from sys import argv, exit

if len(argv) < 2 or len(argv) > 3 or '-h' in argv:
    print 'usage: center.py FILE [DISTANCE]'
    print '       centers the structure in the current box and'
    print '       optionally adds DISTANCE amount of vacuum to FILE'
    print
    exit(0)
    
filename = argv[1]
if len(argv) == 3:
    distance = float(argv[2])
else:
    distance = None

atoms = aselite.read_any(filename)
if distance:
    atoms.center(distance)
else:
    atoms.center()

atoms.write(filename)
