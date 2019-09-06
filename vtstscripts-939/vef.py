#!/usr/bin/env python
import aselite
from os.path import isfile, dirname, abspath, join
from os import system
from sys import exit, argv

if '-h' in argv:
    print 'usage: vef.py'
    print '       prints the force and energy for each ionic step of a vasp run'
    print
    exit(0)

vtst_path = dirname(abspath(__file__))

filename = 'OUTCAR'
if not isfile(filename):
    print 'No such file: %s' % filename
    exit(1)

traj = aselite.read_vasp_out(filename)
if len(traj) == 0:
    exit(0)

fe = open('fe.dat', 'w')
for i, atoms in enumerate(traj):
    e = atoms.get_potential_energy()
    if i == 0:
        e0 = e
    f = atoms.get_max_atom_force()
    str = '%5i %20.8f %20.6f %20.6g ' % (i,f,e,e-e0)
    print str
    fe.write(str+'\n')
fe.close()

if i > 1:
    system('gnuplot %s' % join(vtst_path, 'vef.gnu'))
