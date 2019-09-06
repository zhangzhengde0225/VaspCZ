#!/usr/bin/env python
import sys
import aselite

if len(sys.argv) < 3:
    print 'usage: 2con.py IN OUT'
    print '       converts file IN of type POSCAR or xyz to con and saves it to OUT'
    print
    sys.exit(0)

traj = aselite.read_any(sys.argv[1])
aselite.write_con(sys.argv[2], traj)
