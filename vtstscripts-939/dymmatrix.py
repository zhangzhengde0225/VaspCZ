#!/usr/bin/env python
import numpy
numpy.seterr(all='raise')
import os
import glob
import aselite
from sys import exit

def check_file(filename):
    if not os.path.isfile(filename):
        print 'No such file: %s' % filename
        exit(1)


def dymmatrix(displacecars, outcars):

    all_displacements = []
    displacements = None
    for displacecar in displacecars:
        check_file(displacecar)
        print 'Reading %s' % displacecar
        d = numpy.loadtxt(displacecar)[:,:3].flatten()
        all_displacements.extend(d)
        if displacements == None:
            displacements = d
        else:
            displacements += d

    ndisp = numpy.count_nonzero(displacements)
    #indices of the nonzero displacements
    di = numpy.nonzero(all_displacements)[0] % ndisp
    print 'Number of displacements: %i' % ndisp

    traj = []
    for outcar in outcars:
        check_file(outcar)
        print 'Reading %s' % outcar
        images = aselite.read_vasp_out(outcar)
        atoms = images[0]
        traj += images[1:]

    reference_force = atoms.get_forces().flatten()
    masses = traj[0].get_masses()
    displacement_masses = []
    for mass in masses:
        displacement_masses.extend([mass,mass,mass])
    masses = numpy.array(displacement_masses)

    if len(traj) != ndisp:
        print 'ERROR: number of displacements (%i) !=' % ndisp,
        print 'number of calculations (%i)' % len(traj)
        exit(1)

    print 'Building dynamical matrix'
    dymmat = numpy.zeros((ndisp,ndisp))
    hessian = numpy.zeros_like(dymmat)

    for i in range(ndisp):
        f1 = traj[i].get_forces().flatten()[di]
        f0 = reference_force[di]
        dymmat[i] = -(f1-f0)
        dymmat[i] /= displacements[di]
        hessian[i] = dymmat[i]
        for j in range(ndisp):
            dymmat[i,j] /= numpy.sqrt(masses[di][i]*masses[di][j])

    #symmetrize
    dymmat = (dymmat + dymmat.transpose()) / 2.0
    hessian = (hessian + hessian.transpose()) / 2.0

    numpy.savetxt('freq.mat', dymmat, fmt='%16.8f')

    print 'Diagonalizing matrix'

    omegas, ev = numpy.linalg.eigh(dymmat)
    numpy.savetxt('eigs.dat', omegas, fmt='%25.15g')

    f = open('freq.dat', 'w')
    for omega in omegas:
        imag = 0
        if omega < 0:
            imag = 1
        freq = numpy.sqrt(numpy.abs(omega))*521.47
        s = '%12.6f cm^{-1} ... %i ' % (freq, imag)
        f.write(s+'\n')
        print s
    f.close()

    numpy.savetxt('modes.dat', ev, fmt='%16.8f')
    #
    f = open('modes_sqrt_amu.dat', 'w')
    masses_ = atoms.get_masses()
    for i in range(len(atoms)*3):
        evec = ev[:,i].tolist() # eigenvectors are in columns  
        for j in range(len(atoms)):
            dx = evec[3*j]/numpy.sqrt(masses_[j])
            dy = evec[3*j+1]/numpy.sqrt(masses_[j])
            dz = evec[3*j+2]/numpy.sqrt(masses_[j])
            f.write('%10.6f  %10.6f  %10.6f\n' % (dx, dy, dz))
        f.write('\n')
    f.close()
         #

    force_constants, ev = numpy.linalg.eigh(hessian)
    numpy.savetxt('force_constants.dat', force_constants, fmt='%16.12f')
    effective_masses = force_constants/omegas
    numpy.savetxt('effective_masses.dat', effective_masses, fmt='%12.6f')

def usage():
    print 'usage: dymmmatrix.py [DISPLACECAR] [OUTCAR]'
    print '   or  dymmmatrix.py #DISPLACECAR DISPLACECAR1',
    print 'DISPLACECAR2 ...'
    print '                     OUTCAR1 OUTCAR2 OUTCAR3 ...'
    print

if __name__ == '__main__':
    from sys import argv
    if '-h' in argv:
        usage()
        exit(0)
    if len(argv) == 1:
        displacecars = ['DISPLACECAR']
        outcars = ['OUTCAR']
    elif len(argv) == 3:
        displacecars = [argv[1]]
        outcars = [argv[2]]
    elif len(argv) > 3:
        ndisplacecars = int(argv[1])
        displacecars = argv[2:2+ndisplacecars]
        outcars = argv[2+ndisplacecars:]
    else:
        usage()
        exit(1)

    dymmatrix(displacecars, outcars)
