#!/usr/bin/env python
import numpy as np
#import ase
import aselite as ase

### READ DOSCAR ###
def read_dosfile():
    f = open("DOSCAR", 'r')
    lines = f.readlines()
    f.close()
    index = 0
    natoms = int(lines[index].strip().split()[0])
    index = 5
    nedos = int(lines[index].strip().split()[2])
    efermi = float(lines[index].strip().split()[3])
    print natoms, nedos, efermi

    return lines, index, natoms, nedos, efermi

### READ POSCAR or CONTCAR and save pos
def read_posfile():
    from ase.io import read

    try:
        atoms = read('POSCAR')
    except IOError:
        print "[__main__]: Couldn't open input file POSCAR, atomic positions will not be written...\n"
        atoms = []

    return atoms

### WRITE DOS0 CONTAINING TOTAL DOS ###
def write_dos0(lines, index, nedos, efermi):

    fdos = open("DOS0", 'w')
    line = lines[index+1].strip().split()
    ncols = int(len(line))
    # fdos.write('# %d \n' % (ncols)) #GH not sure why this is here

    for n in xrange(nedos):
        index +=1
        e = float(lines[index].strip().split()[0])
        e_f = e-efermi
        fdos.write('%15.8f ' % (e_f))

        for col in xrange(1, ncols):
            dos = float(lines[index].strip().split()[col])
            fdos.write('%15.8f ' % (dos))
        fdos.write('\n')
    return index

### LOOP OVER SETS OF DOS, NATOMS ###
def write_nospin(lines, index, nedos, natoms, ncols, efermi):
    
    atoms = read_posfile()
    if len(atoms) < natoms:
    	pos = np.zeros((natoms, 3))
    else:
        pos = atoms.get_positions()

    for i in xrange(1,natoms+1):
        si = str(i)

    ## OPEN DOSi FOR WRITING ##
        fdos = open("DOS"+si, 'w')
        index += 1
        ia = i-1
        # fdos.write('# %d \n' % (ncols))
        fdos.write('# %15.8f %15.8f %15.8f \n' % (pos[ia,0], pos[ia,1], pos[ia,2]))

    ### LOOP OVER NEDOS ###
        for n in xrange(nedos):
            index += 1
            e = float(lines[index].strip().split()[0])
            e_f = e-efermi
            fdos.write('%15.8f ' % (e_f))

            for col in xrange(1, ncols):
                dos = float(lines[index].strip().split()[col])
                fdos.write('%15.8f ' % (dos))
            fdos.write('\n')
    fdos.close()

def write_spin(lines, index, nedos, natoms, ncols, efermi):
    #pos=[]
    atoms = read_posfile()
    if len(atoms) < natoms:
        pos = np.zeros((natoms, 3))
    else:
        pos = atoms.get_positions()

    nsites = (ncols -1)/2

    for i in xrange(1,natoms+1):
        si = str(i)
    ## OPEN DOSi FOR WRITING ##
        fdos = open("DOS"+si, 'w')
        index += 1
        ia = i-1
        fdos.write('# %d \n' % (ncols))
        fdos.write('# %15.8f %15.8f %15.8f \n' % (pos[ia,0], pos[ia,1], pos[ia,2]))

    ### LOOP OVER NEDOS ###
        for n in xrange(nedos):
            index +=1   
            e = float(lines[index].strip().split()[0])
            e_f = e-efermi
            fdos.write('%15.8f ' % (e_f))

            for site in xrange(nsites):
                dos_up = float(lines[index].strip().split()[site*2+1])
                dos_down = float(lines[index].strip().split()[site*2+2])*-1
                fdos.write('%15.8f %15.8f ' % (dos_up, dos_down))
            fdos.write('\n')
        fdos.close()

#
if __name__ == '__main__':
    import sys
    import os
    import datetime
    import time
    import optparse

    lines, index, natoms, nedos, efermi = read_dosfile()
    index = write_dos0(lines, index, nedos, efermi)
    ## Test if a spin polarized calculation was performed ##
    line = lines[index+2].strip().split()
    ncols = int(len(line)) 
    if ncols==7 or ncols==19 or ncols==9 or ncols==33:
        write_spin(lines, index, nedos, natoms, ncols, efermi)
        is_spin=True
    else: 
        write_nospin(lines, index, nedos, natoms, ncols, efermi)
        is_spin=False
    print "Spin unrestricted calculation: ", is_spin
