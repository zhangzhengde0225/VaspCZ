#!/usr/bin/env python
#encoding: utf-8

'''aselite is a striped down single file version of ase that retains the
following features: atom and atoms objects, some of ase.io and some of
ase.constraints.'''

# Copyright 2008, 2009 CAMd
# (see accompanying license files for details).

from math import cos, sin
import warnings
import numpy as np
np.seterr(all='raise')
import os

def read_any(filename):
    try:
        return read_vasp(filename)
    except:
        pass
    try:
        return read_xyz(filename)
    except:
        pass
    try:
        return read_con(filename)
    except:
        pass
    raise IOError, "Could not read file %s." % filename

def write_jmol(filename, atoms, eigenvalues, eigenvectors):
    f_xyz = open(filename,'w')
    for i in xrange(len(eigenvectors)):
        mode = eigenvectors[:,i]
        mode.shape = (len(mode)/3,3)

        f_xyz.write("%i\n"%len(atoms))
        f_xyz.write("%f\n"%eigenvalues[i])
        for j,atom in enumerate(atoms):
            f_xyz.write("%s %f %f %f %f %f %f\n" % (atom.symbol, atom.position[0], atom.position[1], atom.position[2], mode[j,0], mode[j,1], mode[j,2]))
    f_xyz.close()


def get_atomtypes(fname):
    """Given a file name, get the atomic symbols. 

    The function can get this information from OUTCAR and POTCAR
    format files.  The files can also be compressed with gzip or
    bzip2.

    """
    atomtypes=[]
    if fname.find('.gz') != -1:
        import gzip
        f = gzip.open(fname)
    elif fname.find('.bz2') != -1:
        import bz2
        f = bz2.BZ2File(fname)
    else:
        f = open(fname)
    for line in f:
        if line.find('TITEL') != -1:
            atomtypes.append(line.split()[3].split('_')[0].split('.')[0])
    return atomtypes

def atomtypes_outpot(posfname, numsyms):
    """Try to retreive chemical symbols from OUTCAR or POTCAR
    
    If getting atomtypes from the first line in POSCAR/CONTCAR fails, it might
    be possible to find the data in OUTCAR or POTCAR, if these files exist.

    posfname -- The filename of the POSCAR/CONTCAR file we're trying to read
    
    numsyms -- The number of symbols we must find

    """
    import os.path as op
    import glob

    # First check files with exactly same name except POTCAR/OUTCAR instead
    # of POSCAR/CONTCAR.
    fnames = [posfname.replace('POSCAR', 'POTCAR').replace('CONTCAR', 
                                                           'POTCAR')]
    fnames.append(posfname.replace('POSCAR', 'OUTCAR').replace('CONTCAR',
                                                               'OUTCAR'))
    # Try the same but with compressed files
    fsc = []
    for fn in fnames:
        fsc.append(fn + '.gz')
        fsc.append(fn + '.bz2')
    for f in fsc:
        fnames.append(f)
    # Finally try anything with POTCAR or OUTCAR in the name
    vaspdir = op.dirname(posfname)
    fs = glob.glob(vaspdir + '*POTCAR*')
    for f in fs:
        fnames.append(f)
    fs = glob.glob(vaspdir + '*OUTCAR*')
    for f in fs:
        fnames.append(f)

    tried = []
    files_in_dir = os.listdir('.')
    for fn in fnames:
        if fn in files_in_dir:
            tried.append(fn)
            at = get_atomtypes(fn)
            if len(at) == numsyms:
                return at

    raise IOError('Could not determine chemical symbols. Tried files ' 
                  + str(tried))


def get_atomtypes_from_formula(formula):
    """Return atom types from chemical formula (optionally prepended
    with and underscore).
    """
    symbols = string2symbols(formula.split('_')[0])
    atomtypes = [symbols[0]]
    for s in symbols[1:]:
        if s != atomtypes[-1]: atomtypes.append(s)
    return atomtypes


def read_vasp(filename='CONTCAR'):
    """Import POSCAR/CONTCAR type file.

    Reads unitcell, atom positions and constraints from the POSCAR/CONTCAR
    file and tries to read atom types from POSCAR/CONTCAR header, if this fails
    the atom types are read from OUTCAR or POTCAR file.
    """
 
    if isinstance(filename, str):
        f = open(filename)
    else: # Assume it's a file-like object
        f = filename

    # First line should contain the atom symbols , eg. "Ag Ge" in
    # the same order
    # as later in the file (and POTCAR for the full vasp run)
    atomtypes = f.readline().split()

    # Sometimes the first line in POSCAR/CONTCAR is of the form
    # "CoP3_In-3.pos". Check for this case and extract atom types
    if len(atomtypes) == 1 and '_' in atomtypes[0]:
        atomtypes = get_atomtypes_from_formula(atomtypes[0])

    lattice_constant = float(f.readline().split()[0])

    # Now the lattice vectors
    a = []
    for ii in range(3):
        s = f.readline().split()
        floatvect = float(s[0]), float(s[1]), float(s[2])
        a.append(floatvect)

    basis_vectors = np.array(a) * lattice_constant

    # Number of atoms. Again this must be in the same order as
    # in the first line
    # or in the POTCAR or OUTCAR file
    atom_symbols = []
    numofatoms = f.readline().split()
    #vasp5.1 has an additional line which gives the atom types
    #the following try statement skips this line
    try:
        int(numofatoms[0])
    except ValueError:
        numofatoms = f.readline().split()

    # check for comments in numofatoms line and get rid of them if necessary
    commentcheck = np.array(['!' in s for s in numofatoms])
    if commentcheck.any():
        # only keep the elements up to the first including a '!':
        numofatoms = numofatoms[:np.arange(len(numofatoms))[commentcheck][0]]
        
    numsyms = len(numofatoms)
    if len(atomtypes) < numsyms:
        # First line in POSCAR/CONTCAR didn't contain enough symbols.
        atomtypes = atomtypes_outpot(f.name, numsyms)
    else:
        try:
            for atype in atomtypes[:numsyms]:
                if not atype in chemical_symbols:
                    raise KeyError
        except KeyError:
            atomtypes = atomtypes_outpot(f.name, numsyms)

    for i, num in enumerate(numofatoms):
        numofatoms[i] = int(num)
        [atom_symbols.append(atomtypes[i]) for na in xrange(numofatoms[i])]

    # Check if Selective dynamics is switched on
    sdyn = f.readline()
    selective_dynamics = sdyn[0].lower() == "s"

    # Check if atom coordinates are cartesian or direct
    if selective_dynamics:
        ac_type = f.readline()
    else:
        ac_type = sdyn
    cartesian = ac_type[0].lower() == "c" or ac_type[0].lower() == "k"
    tot_natoms = sum(numofatoms)
    atoms_pos = np.empty((tot_natoms, 3))
    if selective_dynamics:
        selective_flags = np.empty((tot_natoms, 3), dtype=bool)
    for atom in xrange(tot_natoms):
        ac = f.readline().split()
        atoms_pos[atom] = (float(ac[0]), float(ac[1]), float(ac[2]))
        if selective_dynamics:
            curflag = []
            for flag in ac[3:6]:
                curflag.append(flag == 'F')
            selective_flags[atom] = curflag
    # Done with all reading
    if type(filename) == str:
        f.close()
    if cartesian:
        atoms_pos *= lattice_constant
    atoms = Atoms(symbols = atom_symbols, cell = basis_vectors, pbc = True)
    if cartesian:
        atoms.set_positions(atoms_pos)
    else:
        atoms.set_scaled_positions(atoms_pos)
    if selective_dynamics:
        constraints = []
        indices = []
        for ind, sflags in enumerate(selective_flags):
            if sflags.any() and not sflags.all():
                constraints.append(FixScaled(atoms.get_cell(), ind, sflags))
            elif sflags.all():
                indices.append(ind)
        if indices:
            constraints.append(FixAtoms(indices))
        if constraints:
            atoms.set_constraint(constraints)
    atoms.format = 'vasp'
    return atoms

def read_vasp_out(filename='OUTCAR',index = 'all'):
    """Import OUTCAR type file.

    Reads unitcell, atom positions, energies, and forces from the OUTCAR file
    and attempts to read constraints (if any) from CONTCAR/POSCAR, if present. 
    """
    try:          # try to read constraints, first from CONTCAR, then from POSCAR
        constr = read_vasp('CONTCAR').constraints
    except:
        try:
            constr = read_vasp('POSCAR').constraints
        except:
            constr = None

    if isinstance(filename, str):
        f = open(filename)
    else: # Assume it's a file-like object
        f = filename
    data    = f.readlines()
    natoms  = 0
    images  = []
    atoms   = Atoms(pbc = True, constraint = constr)
    energy  = 0
    species = []
    species_num = []
    symbols = []
    ecount = 0
    poscount = 0
    for n,line in enumerate(data):
        if 'POTCAR:' in line:
            temp = line.split()[2]
            for c in ['.','_','1']:
                if c in temp:
                    temp = temp[0:temp.find(c)]
            species += [temp]
        if 'ions per type' in line:
            species = species[:len(species)/2]
            temp = line.split()
            for ispecies in range(len(species)):
                species_num += [int(temp[ispecies+4])]
                natoms += species_num[-1]
                for iatom in range(species_num[-1]): symbols += [species[ispecies]]
        if 'direct lattice vectors' in line:
            cell = []
            for i in range(3):
                temp = data[n+1+i].split()
                cell += [[float(temp[0]), float(temp[1]), float(temp[2])]]
        if 'energy  without entropy' in line:
            energy = float(data[n].split()[6])
            #energy = float(data[n+2].split()[4])
            if ecount < poscount:
                # reset energy for LAST set of atoms, not current one - VASP 5.11? and up
                images[-1].calc.energy = energy
            ecount += 1
        if 'POSITION          ' in line:
            forces = []
            atoms_symbols = []
            atoms_positions = []
            positions = []
            for iatom in range(natoms):
                temp    = data[n+2+iatom].split()
                atoms_symbols.append(symbols[iatom])
                atoms_positions.append([float(temp[0]),float(temp[1]),float(temp[2])])
                forces += [[float(temp[3]),float(temp[4]),float(temp[5])]]
            atoms = Atoms('H'*natoms, pbc = True, constraint = constr)
            atoms.set_cell(cell)
            atoms.set_chemical_symbols(atoms_symbols)
            atoms.set_positions(atoms_positions)
            atoms.set_calculator(SinglePointCalculator(energy,forces,None,None,atoms))
            images += [atoms]
            poscount += 1

        if 'HIPREC TOTAL-FORCE' in line:
            forces = []
            for line in data[n+2:n+2+natoms]:
                fields = line.split()
                force = []
                for i in range(3):
                    force.append(float(fields[i]))
                forces.append(force)
            images[-1].calc.forces = np.array(forces)

    # return requested images, code borrowed from ase/io/trajectory.py
    if isinstance(index, int):
        return images[index]
    elif index == 'all':
        return images
    else:
        step = index.step or 1
        if step > 0:
            start = index.start or 0
            if start < 0:
                start += len(images)
            stop = index.stop or len(images)
            if stop < 0:
                stop += len(images)
        else:
            if index.start is None:
                start = len(images) - 1
            else:
                start = index.start
                if start < 0:
                    start += len(images)
            if index.stop is None:
                stop = -1
            else:
                stop = index.stop
                if stop < 0:
                    stop += len(images)
        return [images[i] for i in range(start, stop, step)]

def write_vasp(filename, atoms, label='', direct=False, sort=None, symbol_count = None, long_format=True):
    """Method to write VASP position (POSCAR/CONTCAR) files.

    Writes label, scalefactor, unitcell, # of various kinds of atoms,
    positions in cartesian or scaled coordinates (Direct), and constraints
    to file. Cartesian coordiantes is default and default label is the 
    atomic species, e.g. 'C N H Cu'.
    """
    
    if isinstance(filename, str):
        f = open(filename, 'w')
    else: # Assume it's a 'file-like object'
        f = filename
    
    if isinstance(atoms, (list, tuple)):
        if len(atoms) > 1:
            raise RuntimeError("Don't know how to save more than "+
                               "one image to VASP input")
        else:
            atoms = atoms[0]

    # Write atom positions in scaled or cartesian coordinates
    if direct:
        coord = atoms.get_scaled_positions()
    else:
        coord = atoms.get_positions()

    if atoms.constraints:
        sflags = np.zeros((len(atoms), 3), dtype=bool)
        for constr in atoms.constraints:
            if isinstance(constr, FixScaled):
                sflags[constr.a] = constr.mask
            elif isinstance(constr, FixAtoms):
                sflags[constr.index] = [True, True, True]

    if sort:
        ind = np.argsort(atoms.get_chemical_symbols())
        symbols = np.array(atoms.get_chemical_symbols())[ind]
        coord = coord[ind]
        if atoms.constraints:
            sflags = sflags[ind]
    else:
        symbols = atoms.get_chemical_symbols()

    # Create a list sc of (symbol, count) pairs
    if symbol_count:
        sc = symbol_count
    else:
        sc = []
        psym = symbols[0]
        count = 0
        for sym in symbols:
            if sym != psym:
                sc.append((psym, count))
                psym = sym
                count = 1
            else:
                count += 1
        sc.append((psym, count))

    # Create the label
    if label == '':
        for sym, c in sc:
            label += '%2s ' % sym
    f.write(label + '\n')

    # Write unitcell in real coordinates and adapt to VASP convention 
    # for unit cell
    # ase Atoms doesn't store the lattice constant separately, so always
    # write 1.0.
    f.write('%19.16f\n' % 1.0)
    if long_format:
        latt_form = ' %21.16f'
    else:
        latt_form = ' %11.6f'
    for vec in atoms.get_cell():
        f.write(' ')
        for el in vec:
            f.write(latt_form % el)
        f.write('\n')

    # Numbers of each atom
    for sym, count in sc:
        f.write(' %3i' % count)
    f.write('\n')

    if atoms.constraints:
        f.write('Selective dynamics\n')

    if direct:
        f.write('Direct\n')
    else:
        f.write('Cartesian\n')

    if long_format:
        cform = ' %19.16f'
    else:
        cform = ' %9.6f'
    for iatom, atom in enumerate(coord):
        for dcoord in atom:
            f.write(cform % dcoord)
        if atoms.constraints:
            for flag in sflags[iatom]:
                if flag:
                    s = 'F'
                else:
                    s = 'T'
                f.write('%4s' % s)
        f.write('\n')

    if type(filename) == str:
        f.close()
def length_angle_to_box(boxlengths, angles):
    box = np.zeros( (3,3) )
    angles *= np.pi/180.0
    box[0][0] = 1.0
    box[1][0] = np.cos(angles[0])
    box[1][1] = np.sin(angles[0])
    box[2][0] = np.cos(angles[1])
    box[2][1] = (np.cos(angles[2]) - box[1][0] * box[2][0])/box[1][1]
    box[2][2] = np.sqrt(1.0 - box[2][0]**2 - box[2][1]**2)
    box[0,:]*=boxlengths[0]
    box[1,:]*=boxlengths[1]
    box[2,:]*=boxlengths[2]
    return box

def box_to_length_angle(box):
    lengths = np.zeros(3)
    lengths[0] = np.linalg.norm(box[0,:])
    lengths[1] = np.linalg.norm(box[1,:])
    lengths[2] = np.linalg.norm(box[2,:])
    angles = np.zeros(3)
    angles[0] = np.arccos(np.dot(box[0,:]/lengths[0],box[1,:]/lengths[1]))
    angles[1] = np.arccos(np.dot(box[0,:]/lengths[0],box[2,:]/lengths[2]))
    angles[2] = np.arccos(np.dot(box[1,:]/lengths[1],box[2,:]/lengths[2]))
    angles *= 180.0/np.pi
    return lengths, angles
    
    
def read_con(filename):
    f = open(filename, 'r')
    lines = f.readlines()
    f.close()
    trajectory = []
    line_index = 0
    while True:
        try:
            boxlengths = np.array([float(length) for length in lines[line_index+2].split()])
            boxangles = np.array([float(angle) for angle in lines[line_index+3].split()])
            cell = length_angle_to_box(boxlengths, boxangles)
            num_types = int(lines[line_index+6].strip())
            num_each_type = [int(n) for n in lines[line_index+7].split()]
            mass_each_type = [float(n) for n in lines[line_index+8].split()]
            a = Atoms('H'*sum(num_each_type))
            a.format = 'con'
            a.cell = cell
            a.set_pbc((True, True, True))
            frozen = []
            positions = []
            symbols = []
            masses = []
            line_index += 9
            atom_index = 0
            for i in range(num_types):
                symbol = lines[line_index].strip()
                mass = mass_each_type[i]
                line_index += 2
                for j in range(num_each_type[i]):
                    split = lines[line_index].split()
                    positions.append([float(s) for s in split[0:3]])
                    symbols.append(symbol)
                    masses.append(mass)
                    if split[3] != '0':
                        frozen.append(atom_index)
                    atom_index += 1
                    line_index += 1
            a.set_chemical_symbols(symbols)
            a.set_positions(positions)
            a.set_masses(masses)
            a.set_constraint(FixAtoms(frozen))
        except:
            if len(trajectory) == 1:
                return trajectory[0]
            if len(trajectory) == 0:
                raise
            return trajectory
        trajectory.append(a)        
            
            
def write_con(filename, p, w = 'w'):
    con = open(filename, w)
    print >> con, "Generated by vtstscripts"
    print >> con
    lengths, angles = box_to_length_angle(p.cell)
    print >> con, " ".join(['%12.6f' % s for s in lengths])
    print >> con, " ".join(['%12.6f' % s for s in angles])
    print >> con
    print >> con
    atom_count = {}
    name_order = []
    for i in range(len(p)):
        name = p[i].symbol
        if name not in name_order:
            name_order.append(name)
        if name in atom_count:
            atom_count[name] += 1
        else:
            atom_count[name] = 1
    print >> con, len(name_order)
    print >> con, " ".join([str(atom_count[i]) for i in name_order])
    printmasses = []
    index = 0
    for i in range(len(name_order)):
        printmasses.append(p[index].mass)
        index += atom_count[name_order[i]]
    print >> con, " ".join(["%12.6f"% i for i in printmasses])
    index = 0
    for i in range(len(name_order)):
        print >> con, name_order[i]
        print >> con, "Coordinates of Component", i+1
        for j in range(atom_count[name_order[i]]):
            free = 0
            if len(p.constraints) > 0:
                if index in p.constraints[0].index:
                    free = 1
            con.write("%12.6f %12.6f %12.6f %d %d\n" % (p[index].position[0],
                      p[index].position[1], p[index].position[2], free, index))
            index += 1

def read_xdatcar(fileName):
    f = open(fileName, 'r')
    lines = f.readlines()
    f.close()
    lattice_constant = float(lines[1].strip())
    cell = np.array([[float(x) * lattice_constant for x in lines[2].split()], 
                        [float(x) * lattice_constant for x in lines[3].split()], 
                        [float(x) * lattice_constant for x in lines[4].split()]])
    elements = lines[5].split()
    natoms = [int(x) for x in lines[6].split()]
    nframes = (len(lines)-7)/(sum(natoms) + 1)
    trajectory = []
    for i in range(nframes):
        a = Atoms('H'*sum(natoms))
        a.masses = [1.0] * len(a)
        a.set_chemical_symbols(''.join([n*e for (n, e) in zip(natoms, elements)]))
        a.cell = cell.copy()
        j = 0
        for N, e in zip(natoms, elements):
            for k in range(N):
                split = lines[8 + i * (sum(natoms) + 1) + j].split()
                a[j].position = [float(l) for l in split[0:3]]
                j += 1
        a.positions = np.dot(a.positions, cell)
        trajectory.append(a)
    return trajectory

def read_xyz(fileobj, index=-1):
    if isinstance(fileobj, str):
        fileobj = open(fileobj)

    lines = fileobj.readlines()
    L1 = lines[0].split()
    if len(L1) == 1:
        del lines[:2]
        natoms = int(L1[0])
    else:
        natoms = len(lines)
    images = []
    while len(lines) >= natoms:
        positions = []
        symbols = []
        for line in lines[:natoms]:
            symbol, x, y, z = line.split()[:4]
            symbol = symbol.lower().capitalize()
            symbols.append(symbol)
            positions.append([float(x), float(y), float(z)])
        images.append(Atoms(symbols=symbols, positions=positions))
        images[-1].format = 'xyz'
        del lines[:natoms + 2]
    return images[index]

def write_xyz(fileobj, images):
    if isinstance(fileobj, str):
        fileobj = open(fileobj, 'w')

    if not isinstance(images, (list, tuple)):
        images = [images]

    symbols = images[0].get_chemical_symbols()
    natoms = len(symbols)
    for atoms in images:
        fileobj.write('%d\n\n' % natoms)
        for s, (x, y, z) in zip(symbols, atoms.get_positions()):
            fileobj.write('%-2s %22.15f %22.15f %22.15f\n' % (s, x, y, z))

class Units:
    from math import pi, sqrt
    # Constants from Konrad Hinsen's PhysicalQuantities module (1986 CODATA):
    _c = 299792458.              # speed of light, m/s
    _mu0 = 4.e-7 * pi            # permeability of vacuum
    _eps0 = 1 / _mu0 / _c**2     # permittivity of vacuum
    _Grav = 6.67259e-11          # gravitational constant
    _hplanck = 6.6260755e-34     # Planck constant, J s
    _hbar = _hplanck / (2 * pi)  # Planck constant / 2pi, J s
    _e = 1.60217733e-19          # elementary charge
    _me = 9.1093897e-31          # electron mass
    _mp = 1.6726231e-27          # proton mass
    _Nav = 6.0221367e23          # Avogadro number
    _k = 1.380658e-23            # Boltzmann constant, J/K
    _amu = 1.6605402e-27         # atomic mass unit, kg

    Ang = Angstrom = 1.0
    nm = 10.0
    Bohr = 4e10 * pi * _eps0 * _hbar**2 / _me / _e**2  # Bohr radius

    eV = 1.0
    Hartree = _me * _e**3 / 16 / pi**2 / _eps0**2 / _hbar**2
    kJ = 1000.0 / _e
    kcal = 4.184 * kJ
    mol = _Nav
    Rydberg = 0.5 * Hartree
    Ry = Rydberg
    Ha = Hartree

    second = 1e10 * sqrt(_e / _amu)
    fs = 1e-15 * second

    kB = _k / _e                 # Boltzmann constant, eV/K

    Pascal = (1 / _e) / 1e30  # J/m^3
    GPa = 1e9 * Pascal

    Debye = 1e11 *_e * _c
    alpha = _e**2 / (4 * pi * _eps0) / _hbar / _c # fine structure constant

    # Derived atomic units that have no assigned name:
    _aut = _hbar / (alpha**2 * _me * _c**2)      # atomic unit of time, s
    _auv =  _e**2 / _hbar / (4 * pi * _eps0)     # atomic unit of velocity, m/s
    _auf = alpha**3 * _me**2 * _c**3 / _hbar     # atomic unit of force, N
    _aup = alpha**5 * _me**4 * _c**5 / _hbar**3  # atomic unit of pressure, Pa

    AUT = second * _aut
units = Units()


chemical_symbols = ['X',  'H',  'He', 'Li', 'Be',
                    'B',  'C',  'N',  'O',  'F',
                    'Ne', 'Na', 'Mg', 'Al', 'Si',
                    'P',  'S',  'Cl', 'Ar', 'K',
                    'Ca', 'Sc', 'Ti', 'V',  'Cr',
                    'Mn', 'Fe', 'Co', 'Ni', 'Cu',
                    'Zn', 'Ga', 'Ge', 'As', 'Se',
                    'Br', 'Kr', 'Rb', 'Sr', 'Y',
                    'Zr', 'Nb', 'Mo', 'Tc', 'Ru',
                    'Rh', 'Pd', 'Ag', 'Cd', 'In',
                    'Sn', 'Sb', 'Te', 'I',  'Xe',
                    'Cs', 'Ba', 'La', 'Ce', 'Pr',
                    'Nd', 'Pm', 'Sm', 'Eu', 'Gd',
                    'Tb', 'Dy', 'Ho', 'Er', 'Tm',
                    'Yb', 'Lu', 'Hf', 'Ta', 'W',
                    'Re', 'Os', 'Ir', 'Pt', 'Au',
                    'Hg', 'Tl', 'Pb', 'Bi', 'Po',
                    'At', 'Rn', 'Fr', 'Ra', 'Ac',
                    'Th', 'Pa', 'U',  'Np', 'Pu',
                    'Am', 'Cm', 'Bk', 'Cf', 'Es',
                    'Fm', 'Md', 'No', 'Lr']

atomic_numbers = {'A':1}
for Z, symbol in enumerate(chemical_symbols):
    atomic_numbers[symbol] = Z

atomic_names = [
    '', 'Hydrogen', 'Helium', 'Lithium', 'Beryllium', 'Boron',
    'Carbon', 'Nitrogen', 'Oxygen', 'Fluorine', 'Neon', 'Sodium',
    'Magnesium', 'Aluminium', 'Silicon', 'Phosphorus', 'Sulfur',
    'Chlorine', 'Argon', 'Potassium', 'Calcium', 'Scandium',
    'Titanium', 'Vanadium', 'Chromium', 'Manganese', 'Iron',
    'Cobalt', 'Nickel', 'Copper', 'Zinc', 'Gallium', 'Germanium',
    'Arsenic', 'Selenium', 'Bromine', 'Krypton', 'Rubidium',
    'Strontium', 'Yttrium', 'Zirconium', 'Niobium', 'Molybdenum',
    'Technetium', 'Ruthenium', 'Rhodium', 'Palladium', 'Silver',
    'Cadmium', 'Indium', 'Tin', 'Antimony', 'Tellurium',
    'Iodine', 'Xenon', 'Caesium', 'Barium', 'Lanthanum',
    'Cerium', 'Praseodymium', 'Neodymium', 'Promethium',
    'Samarium', 'Europium', 'Gadolinium', 'Terbium',
    'Dysprosium', 'Holmium', 'Erbium', 'Thulium', 'Ytterbium',
    'Lutetium', 'Hafnium', 'Tantalum', 'Tungsten', 'Rhenium',
    'Osmium', 'Iridium', 'Platinum', 'Gold', 'Mercury',
    'Thallium', 'Lead', 'Bismuth', 'Polonium', 'Astatine',
    'Radon', 'Francium', 'Radium', 'Actinium', 'Thorium',
    'Protactinium', 'Uranium', 'Neptunium', 'Plutonium',
    'Americium', 'Curium', 'Berkelium', 'Californium',
    'Einsteinium', 'Fermium', 'Mendelevium', 'Nobelium',
    'Lawrencium', 'Unnilquadium', 'Unnilpentium', 'Unnilhexium']

atomic_masses = np.array([
   0.00000, # X
   1.00794, # H
   4.00260, # He
   6.94100, # Li
   9.01218, # Be
  10.81100, # B
  12.01100, # C
  14.00670, # N
  15.99940, # O
  18.99840, # F
  20.17970, # Ne
  22.98977, # Na
  24.30500, # Mg
  26.98154, # Al
  28.08550, # Si
  30.97376, # P
  32.06600, # S
  35.45270, # Cl
  39.94800, # Ar
  39.09830, # K
  40.07800, # Ca
  44.95590, # Sc
  47.88000, # Ti
  50.94150, # V
  51.99600, # Cr
  54.93800, # Mn
  55.84700, # Fe
  58.93320, # Co
  58.69340, # Ni
  63.54600, # Cu
  65.39000, # Zn
  69.72300, # Ga
  72.61000, # Ge
  74.92160, # As
  78.96000, # Se
  79.90400, # Br
  83.80000, # Kr
  85.46780, # Rb
  87.62000, # Sr
  88.90590, # Y
  91.22400, # Zr
  92.90640, # Nb
  95.94000, # Mo
    np.nan, # Tc
 101.07000, # Ru
 102.90550, # Rh
 106.42000, # Pd
 107.86800, # Ag
 112.41000, # Cd
 114.82000, # In
 118.71000, # Sn
 121.75700, # Sb
 127.60000, # Te
 126.90450, # I
 131.29000, # Xe
 132.90540, # Cs
 137.33000, # Ba
 138.90550, # La
 140.12000, # Ce
 140.90770, # Pr
 144.24000, # Nd
    np.nan, # Pm
 150.36000, # Sm
 151.96500, # Eu
 157.25000, # Gd
 158.92530, # Tb
 162.50000, # Dy
 164.93030, # Ho
 167.26000, # Er
 168.93420, # Tm
 173.04000, # Yb
 174.96700, # Lu
 178.49000, # Hf
 180.94790, # Ta
 183.85000, # W
 186.20700, # Re
 190.20000, # Os
 192.22000, # Ir
 195.08000, # Pt
 196.96650, # Au
 200.59000, # Hg
 204.38300, # Tl
 207.20000, # Pb
 208.98040, # Bi
    np.nan, # Po
    np.nan, # At
    np.nan, # Rn
    np.nan, # Fr
 226.02540, # Ra
    np.nan, # Ac
 232.03810, # Th
 231.03590, # Pa
 238.02900, # U
 237.04820, # Np
    np.nan, # Pu
    np.nan, # Am
    np.nan, # Cm
    np.nan, # Bk
    np.nan, # Cf
    np.nan, # Es
    np.nan, # Fm
    np.nan, # Md
    np.nan, # No
    np.nan])# Lw

# Covalent radii from:
#
#  Covalent radii revisited,
#  Beatriz Cordero, Verónica Gómez, Ana E. Platero-Prats, Marc Revés,
#  Jorge Echeverría, Eduard Cremades, Flavia Barragán and Santiago Alvarez,
#  Dalton Trans., 2008, 2832-2838 DOI:10.1039/B801115J 
missing = 0.2
covalent_radii = np.array([
    missing,  # X
    0.31,  # H
    0.28,  # He
    1.28,  # Li
    0.96,  # Be
    0.84,  # B
    0.76,  # C
    0.71,  # N
    0.66,  # O
    0.57,  # F
    0.58,  # Ne
    1.66,  # Na
    1.41,  # Mg
    1.21,  # Al
    1.11,  # Si
    1.07,  # P
    1.05,  # S
    1.02,  # Cl
    1.06,  # Ar
    2.03,  # K
    1.76,  # Ca
    1.70,  # Sc
    1.60,  # Ti
    1.53,  # V
    1.39,  # Cr
    1.39,  # Mn
    1.32,  # Fe
    1.26,  # Co
    1.24,  # Ni
    1.32,  # Cu
    1.22,  # Zn
    1.22,  # Ga
    1.20,  # Ge
    1.19,  # As
    1.20,  # Se
    1.20,  # Br
    1.16,  # Kr
    2.20,  # Rb
    1.95,  # Sr
    1.90,  # Y
    1.75,  # Zr
    1.64,  # Nb
    1.54,  # Mo
    1.47,  # Tc
    1.46,  # Ru
    1.42,  # Rh
    1.39,  # Pd
    1.45,  # Ag
    1.44,  # Cd
    1.42,  # In
    1.39,  # Sn
    1.39,  # Sb
    1.38,  # Te
    1.39,  # I
    1.40,  # Xe
    2.44,  # Cs
    2.15,  # Ba
    2.07,  # La
    2.04,  # Ce
    2.03,  # Pr
    2.01,  # Nd
    1.99,  # Pm
    1.98,  # Sm
    1.98,  # Eu
    1.96,  # Gd
    1.94,  # Tb
    1.92,  # Dy
    1.92,  # Ho
    1.89,  # Er
    1.90,  # Tm
    1.87,  # Yb
    1.87,  # Lu
    1.75,  # Hf
    1.70,  # Ta
    1.62,  # W
    1.51,  # Re
    1.44,  # Os
    1.41,  # Ir
    1.36,  # Pt
    1.36,  # Au
    1.32,  # Hg
    1.45,  # Tl
    1.46,  # Pb
    1.48,  # Bi
    1.40,  # Po
    1.50,  # At
    1.50,  # Rn
    2.60,  # Fr
    2.21,  # Ra
    2.15,  # Ac
    2.06,  # Th
    2.00,  # Pa
    1.96,  # U
    1.90,  # Np
    1.87,  # Pu
    1.80,  # Am
    1.69,  # Cm
    missing,  # Bk
    missing,  # Cf
    missing,  # Es
    missing,  # Fm
    missing,  # Md
    missing,  # No
    missing,  # Lr
    ])


#         singular,    plural,     default value
names = {'position': ('positions', np.zeros(3)),
         'number':   ('numbers',   0),
         'tag':      ('tags',      0),
         'momentum': ('momenta',   np.zeros(3)),
         'mass':     ('masses',    None),
         'magmom':   ('magmoms',   0.0),
         'charge':   ('charges',   0.0)
         }


def atomproperty(name, doc):
    """Helper function to easily create Atom attribute property."""

    def getter(self):
        return self.get(name)

    def setter(self, value):
        self.set(name, value)

    def deleter(self):
        self.delete(name)

    return property(getter, setter, deleter, doc)


def xyzproperty(index):
    """Helper function to easily create Atom XYZ-property."""

    def getter(self):
        return self.position[index]

    def setter(self, value):
        self.position[index] = value

    return property(getter, setter, doc='XYZ'[index] + '-coordinate')


class Atom(object):
    """Class for representing a single atom.

    Parameters:
    
    symbol: str or int
        Can be a chemical symbol (str) or an atomic number (int).
    position: sequence of 3 floats
        Atomi position.
    tag: int
        Special purpose tag.
    momentum: sequence of 3 floats
        Momentum for atom.
    mass: float
        Atomic mass in atomic units.
    magmom: float or 3 floats
        Magnetic moment.
    charge: float
        Atomic charge.
    """
    __slots__ = ['data', 'atoms', 'index']

    def __init__(self, symbol='X', position=(0, 0, 0),
                 tag=None, momentum=None, mass=None,
                 magmom=None, charge=None,
                 atoms=None, index=None):

        self.data = d = {}

        if atoms is None:
            # This atom is not part of any Atoms object:
            if isinstance(symbol, str):
                d['number'] = atomic_numbers[symbol]
            else:
                d['number'] = symbol
            d['position'] = np.array(position, float)
            d['tag'] = tag
            if momentum is not None:
                momentum = np.array(momentum, float)
            d['momentum'] = momentum
            d['mass'] = mass
            if magmom is not None:
                magmom = np.array(magmom, float)
            d['magmom'] = magmom
            d['charge'] = charge

        self.index = index
        self.atoms = atoms

    def __repr__(self):
        s = "Atom('%s', %s" % (self.symbol, list(self.position))
        for name in ['tag', 'momentum', 'mass', 'magmom', 'charge']:
            value = self.get_raw(name)
            if value is not None:
                if isinstance(value, np.ndarray):
                    value = value.tolist()
                s += ', %s=%s' % (name, value)
        if self.atoms is None:
            s += ')'
        else:
            s += ', index=%d)' % self.index
        return s

    def cut_reference_to_atoms(self):
        """Cut reference to atoms object."""
        for name in names:
            self.data[name] = self.get_raw(name)
        self.index = None
        self.atoms = None
        
    def get_raw(self, name):
        """Get attribute, return None if not explicitely set."""
        if name == 'symbol':
            return chemical_symbols[self.get_raw('number')]

        if self.atoms is None:
            return self.data[name]
        
        plural = names[name][0]
        if plural in self.atoms.arrays:
            return self.atoms.arrays[plural][self.index]
        else:
            return None

    def get(self, name):
        """Get attribute, return default if not explicitely set."""
        value = self.get_raw(name)
        if value is None:
            if name == 'mass':
                value = atomic_masses[self.number]
            else:
                value = names[name][1]
        return value

    def set(self, name, value):
        """Set attribute."""
        if name == 'symbol':
            name = 'number'
            value = atomic_numbers[value]

        if self.atoms is None:
            assert name in names
            self.data[name] = value
        else:
            plural, default = names[name]
            if plural in self.atoms.arrays:
                array = self.atoms.arrays[plural]
                if name == 'magmom' and array.ndim == 2:
                    assert len(value) == 3
                array[self.index] = value
            else:
                if name == 'magmom' and np.asarray(value).ndim == 1:
                    array = np.zeros((len(self.atoms), 3))
                elif name == 'mass':
                    array = self.atoms.get_masses()
                else:
                    default = np.asarray(default)
                    array = np.zeros((len(self.atoms),) + default.shape,
                                     default.dtype)
                array[self.index] = value
                self.atoms.new_array(plural, array)

    def delete(self, name):
        """Delete attribute."""
        assert self.atoms is None
        assert name not in ['number', 'symbol', 'position']
        self.data[name] = None

    symbol = atomproperty('symbol', 'Chemical symbol')
    number = atomproperty('number', 'Atomic number')
    position = atomproperty('position', 'XYZ-coordinates')
    tag = atomproperty('tag', 'Integer tag')
    momentum = atomproperty('momentum', 'XYZ-momentum')
    mass = atomproperty('mass', 'Atomic mass')
    magmom = atomproperty('magmom', 'Initial magnetic moment')
    charge = atomproperty('charge', 'Atomic charge')
    x = xyzproperty(0)
    y = xyzproperty(1)
    z = xyzproperty(2)

    def _get(self, name):
        """Helper function for deprecated get methods."""
        warnings.warn('Use atom.%s' % name, stacklevel=3)
        return getattr(self, name)

    def _set(self, name, value):
        """Helper function for deprecated set methods."""
        warnings.warn('Use atom.%s = ...' % name, stacklevel=3)
        setattr(self, name, value)

    def get_symbol(self): return self._get('symbol')
    def get_atomic_number(self): return self._get('number')
    def get_position(self): return self._get('position')
    def get_tag(self): return self._get('tag')
    def get_momentum(self): return self._get('momentum')
    def get_mass(self): return self._get('mass')
    def get_initial_magnetic_moment(self): return self._get('magmom')
    def get_charge(self): return self._get('charge')

    def set_symbol(self, value): self._set('symbol', value)
    def set_atomic_number(self, value): self._set('number', value)
    def set_position(self, value): self._set('position', value)
    def set_tag(self, value): self._set('tag', value)
    def set_momentum(self, value): self._set('momentum', value)
    def set_mass(self, value): self._set('mass', value)
    def set_initial_magnetic_moment(self, value): self._set('magmom', value)
    def set_charge(self, value): self._set('charge', value)


class Atoms(object):
    """Atoms object.

    The Atoms object can represent an isolated molecule, or a
    periodically repeated structure.  It has a unit cell and
    there may be periodic boundary conditions along any of the three
    unit cell axes.

    Information about the atoms (atomic numbers and position) is
    stored in ndarrays.  Optionally, there can be information about
    tags, momenta, masses, magnetic moments and charges.

    In order to calculate energies, forces and stresses, a calculator
    object has to attached to the atoms object.

    Parameters:

    symbols: str (formula) or list of str
        Can be a string formula, a list of symbols or a list of
        Atom objects.  Examples: 'H2O', 'COPt12', ['H', 'H', 'O'],
        [Atom('Ne', (x, y, z)), ...].
    positions: list of xyz-positions
        Atomic positions.  Anything that can be converted to an
        ndarray of shape (n, 3) will do: [(x1,y1,z1), (x2,y2,z2),
        ...].
    scaled_positions: list of scaled-positions
        Like positions, but given in units of the unit cell.
        Can not be set at the same time as positions.
    numbers: list of int
        Atomic numbers (use only one of symbols/numbers).
    tags: list of int
        Special purpose tags.
    momenta: list of xyz-momenta
        Momenta for all atoms.
    masses: list of float
        Atomic masses in atomic units.
    magmoms: list of float or list of xyz-values
        Magnetic moments.  Can be either a single value for each atom
        for collinear calculations or three numbers for each atom for
        non-collinear calculations.
    charges: list of float
        Atomic charges.
    cell: 3x3 matrix
        Unit cell vectors.  Can also be given as just three
        numbers for orthorhombic cells.  Default value: [1, 1, 1].
    pbc: one or three bool
        Periodic boundary conditions flags.  Examples: True,
        False, 0, 1, (1, 1, 0), (True, False, False).  Default
        value: False.
    constraint: constraint object(s)
        Used for applying one or more constraints during structure
        optimization.
    calculator: calculator object
        Used to attach a calculator for calculating energies and atomic
        forces.
    info: dict of key-value pairs
        Dictionary of key-value pairs with additional information
        about the system.  The following keys may be used by ase:

          - spacegroup: Spacegroup instance
          - unit_cell: 'conventional' | 'primitive' | int | 3 ints
          - adsorbate_info:

        Items in the info attribute survives copy and slicing and can
        be store to and retrieved from trajectory files given that the
        key is a string, the value is picklable and, if the value is a
        user-defined object, its base class is importable.  One should
        not make any assumptions about the existence of keys.

    Examples:

    These three are equivalent:

    >>> d = 1.104  # N2 bondlength
    >>> a = Atoms('N2', [(0, 0, 0), (0, 0, d)])
    >>> a = Atoms(numbers=[7, 7], positions=[(0, 0, 0), (0, 0, d)])
    >>> a = Atoms([Atom('N', (0, 0, 0)), Atom('N', (0, 0, d)])

    FCC gold:

    >>> a = 4.05  # Gold lattice constant
    >>> b = a / 2
    >>> fcc = Atoms('Au',
    ...             cell=[(0, b, b), (b, 0, b), (b, b, 0)],
    ...             pbc=True)

    Hydrogen wire:

    >>> d = 0.9  # H-H distance
    >>> L = 7.0
    >>> h = Atoms('H', positions=[(0, L / 2, L / 2)],
    ...           cell=(d, L, L),
    ...           pbc=(1, 0, 0))
    """

    def __init__(self, symbols=None,
                 positions=None, numbers=None,
                 tags=None, momenta=None, masses=None,
                 magmoms=None, charges=None,
                 scaled_positions=None,
                 cell=None, pbc=None,
                 constraint=None,
                 calculator=None, 
                 info=None):

        atoms = None

        if hasattr(symbols, 'GetUnitCell'):
            from ase.old import OldASEListOfAtomsWrapper
            atoms = OldASEListOfAtomsWrapper(symbols)
            symbols = None
        elif hasattr(symbols, 'get_positions'):
            atoms = symbols
            symbols = None
        elif (isinstance(symbols, (list, tuple)) and
              len(symbols) > 0 and isinstance(symbols[0], Atom)):
            # Get data from a list or tuple of Atom objects:
            data = [[atom.get_raw(name) for atom in symbols]
                    for name in
                    ['position', 'number', 'tag', 'momentum',
                     'mass', 'magmom', 'charge']]
            atoms = self.__class__(None, *data)
            symbols = None

        if atoms is not None:
            # Get data from another Atoms object:
            if scaled_positions is not None:
                raise NotImplementedError
            if symbols is None and numbers is None:
                numbers = atoms.get_atomic_numbers()
            if positions is None:
                positions = atoms.get_positions()
            if tags is None and atoms.has('tags'):
                tags = atoms.get_tags()
            if momenta is None and atoms.has('momenta'):
                momenta = atoms.get_momenta()
            if magmoms is None and atoms.has('magmoms'):
                magmoms = atoms.get_initial_magnetic_moments()
            if masses is None and atoms.has('masses'):
                masses = atoms.get_masses()
            if charges is None and atoms.has('charges'):
                charges = atoms.get_charges()
            if cell is None:
                cell = atoms.get_cell()
            if pbc is None:
                pbc = atoms.get_pbc()
            if constraint is None:
                constraint = [c.copy() for c in atoms.constraints]
            if calculator is None:
                calculator = atoms.get_calculator()

        self.arrays = {}
        
        if symbols is None:
            if numbers is None:
                if positions is not None:
                    natoms = len(positions)
                elif scaled_positions is not None:
                    natoms = len(scaled_positions)
                else:
                    natoms = 0
                numbers = np.zeros(natoms, int)
            self.new_array('numbers', numbers, int)
        else:
            if numbers is not None:
                raise ValueError(
                    'Use only one of "symbols" and "numbers".')
            else:
                self.new_array('numbers', symbols2numbers(symbols), int)

        if cell is None:
            cell = np.eye(3)
        self.set_cell(cell)

        if positions is None:
            if scaled_positions is None:
                positions = np.zeros((len(self.arrays['numbers']), 3))
            else:
                positions = np.dot(scaled_positions, self._cell)
        else:
            if scaled_positions is not None:
                raise RuntimeError('Both scaled and cartesian positions set!')
        self.new_array('positions', positions, float, (3,))

        self.set_constraint(constraint)
        self.set_tags(default(tags, 0))
        self.set_momenta(default(momenta, (0.0, 0.0, 0.0)))
        self.set_masses(default(masses, None))
        self.set_initial_magnetic_moments(default(magmoms, 0.0))
        self.set_charges(default(charges, 0.0))
        if pbc is None:
            pbc = False
        self.set_pbc(pbc)

        if info is None:
            self.info = {}
        else:
            self.info = dict(info)

        self.adsorbate_info = {}

        self.set_calculator(calculator)

    def set_calculator(self, calc=None):
        """Attach calculator object."""
        if hasattr(calc, '_SetListOfAtoms'):
            from ase.old import OldASECalculatorWrapper
            calc = OldASECalculatorWrapper(calc, self)
        if hasattr(calc, 'set_atoms'):
            calc.set_atoms(self)
        self._calc = calc

    def get_calculator(self):
        """Get currently attached calculator object."""
        return self._calc

    def _del_calculator(self):
        self._calc = None

    calc = property(get_calculator, set_calculator, _del_calculator,
                    doc='Calculator object.')
    
    def set_constraint(self, constraint=None):
        """Apply one or more constrains.

        The *constraint* argument must be one constraint object or a
        list of constraint objects."""
        if constraint is None:
            self._constraints = []
        else:
            if isinstance(constraint, (list, tuple)):
                self._constraints = constraint
            else:
                self._constraints = [constraint]

    def _get_constraints(self):
        return self._constraints

    def _del_constraints(self):
        self._constraints = []

    constraints = property(_get_constraints, set_constraint, _del_constraints,
                           'Constraints of the atoms.')
    
    def set_cell(self, cell, scale_atoms=False, fix=None):
        """Set unit cell vectors.

        Parameters:

        cell : 
            Unit cell.  A 3x3 matrix (the three unit cell vectors) or
            just three numbers for an orthorhombic cell.
        scale_atoms : bool
            Fix atomic positions or move atoms with the unit cell?
            Default behavior is to *not* move the atoms (scale_atoms=False).

        Examples:

        Two equivalent ways to define an orthorhombic cell:
        
        >>> a.set_cell([a, b, c])
        >>> a.set_cell([(a, 0, 0), (0, b, 0), (0, 0, c)])

        FCC unit cell:

        >>> a.set_cell([(0, b, b), (b, 0, b), (b, b, 0)])
        """

        if fix is not None:
            raise TypeError('Please use scale_atoms=%s' % (not fix))

        cell = np.array(cell, float)
        if cell.shape == (3,):
            cell = np.diag(cell)
        elif cell.shape != (3, 3):
            raise ValueError('Cell must be length 3 sequence or '
                             '3x3 matrix!')
        if scale_atoms:
            M = np.linalg.solve(self._cell, cell)
            self.arrays['positions'][:] = np.dot(self.arrays['positions'], M)
        self._cell = cell

    def get_cell(self):
        """Get the three unit cell vectors as a 3x3 ndarray."""
        return self._cell.copy()

    def get_reciprocal_cell(self):
        """Get the three reciprocal lattice vectors as a 3x3 ndarray.

        Note that the commonly used factor of 2 pi for Fourier
        transforms is not included here."""
        
        rec_unit_cell = np.linalg.inv(self.get_cell()).transpose()
        return rec_unit_cell

    def set_pbc(self, pbc):
        """Set periodic boundary condition flags."""
        if isinstance(pbc, int):
            pbc = (pbc,) * 3
        self._pbc = np.array(pbc, bool)
        
    def get_pbc(self):
        """Get periodic boundary condition flags."""
        return self._pbc.copy()

    def new_array(self, name, a, dtype=None, shape=None):
        """Add new array.

        If *shape* is not *None*, the shape of *a* will be checked."""
        
        if dtype is not None:
            a = np.array(a, dtype)
        else:
            a = a.copy()
            
        if name in self.arrays:
            raise RuntimeError

        for b in self.arrays.values():
            if len(a) != len(b):
                raise ValueError('Array has wrong length: %d != %d.' %
                                 (len(a), len(b)))
            break

        if shape is not None and a.shape[1:] != shape:
            raise ValueError('Array has wrong shape %s != %s.' %
                             (a.shape, (a.shape[0:1] + shape)))
        
        self.arrays[name] = a

    def get_array(self, name, copy=True):
        """Get an array.

        Returns a copy unless the optional argument copy is false.
        """
        if copy:
            return self.arrays[name].copy()
        else:
            return self.arrays[name]
    
    def set_array(self, name, a, dtype=None, shape=None):
        """Update array.

        If *shape* is not *None*, the shape of *a* will be checked.
        If *a* is *None*, then the array is deleted."""
        
        b = self.arrays.get(name)
        if b is None:
            if a is not None:
                self.new_array(name, a, dtype, shape)
        else:
            if a is None:
                del self.arrays[name]
            else:
                a = np.asarray(a)
                if a.shape != b.shape:
                    raise ValueError('Array has wrong shape %s != %s.' %
                                     (a.shape, b.shape))
                b[:] = a

    def has(self, name):
        """Check for existence of array.

        name must be one of: 'tags', 'momenta', 'masses', 'magmoms',
        'charges'."""
        return name in self.arrays
    
    def set_atomic_numbers(self, numbers):
        """Set atomic numbers."""
        self.set_array('numbers', numbers, int, ())

    def get_atomic_numbers(self):
        """Get integer array of atomic numbers."""
        return self.arrays['numbers'].copy()

    def set_chemical_symbols(self, symbols):
        """Set chemical symbols."""
        self.set_array('numbers', symbols2numbers(symbols), int, ())

    def get_chemical_symbols(self, reduce=False):
        """Get list of chemical symbol strings.

        If reduce is True, a single string is returned, where repeated
        elements have been contracted to a single symbol and a number.
        E.g. instead of ['C', 'O', 'O', 'H'], the string 'CO2H' is returned.
        """
        if not reduce:
            # XXX
            return [chemical_symbols[Z] for Z in self.arrays['numbers']]
        else:
            num = self.get_atomic_numbers()
            N = len(num)
            dis = np.concatenate(([0], np.arange(1, N)[num[1:] != num[:-1]]))
            repeat = np.append(dis[1:], N) - dis
            symbols = ''.join([chemical_symbols[num[d]] + str(r) * (r != 1)
                               for r, d in zip(repeat, dis)])
            return symbols

    def set_tags(self, tags):
        """Set tags for all atoms."""
        self.set_array('tags', tags, int, ())
        
    def get_tags(self):
        """Get integer array of tags."""
        if 'tags' in self.arrays:
            return self.arrays['tags'].copy()
        else:
            return np.zeros(len(self), int)

    def set_momenta(self, momenta):
        """Set momenta."""
        if len(self.constraints) > 0 and momenta is not None:
            momenta = np.array(momenta)  # modify a copy
            for constraint in self.constraints:
                constraint.adjust_forces(self.arrays['positions'], momenta)
        self.set_array('momenta', momenta, float, (3,))

    def set_velocities(self, velocities):
        """Set the momenta by specifying the velocities."""
        self.set_momenta(self.get_masses()[:, np.newaxis] * velocities)
        
    def get_momenta(self):
        """Get array of momenta."""
        if 'momenta' in self.arrays:
            return self.arrays['momenta'].copy()
        else:
            return np.zeros((len(self), 3))
        
    def set_masses(self, masses='defaults'):
        """Set atomic masses.

        The array masses should contain a list of masses.  In case
        the masses argument is not given or for those elements of the
        masses list that are None, standard values are set."""
        
        if masses == 'defaults':
            masses = atomic_masses[self.arrays['numbers']]
        elif isinstance(masses, (list, tuple)):
            newmasses = []
            for m, Z in zip(masses, self.arrays['numbers']):
                if m is None:
                    newmasses.append(atomic_masses[Z])
                else:
                    newmasses.append(m)
            masses = newmasses
        self.set_array('masses', masses, float, ())

    def get_masses(self):
        """Get array of masses."""
        if 'masses' in self.arrays:
            return self.arrays['masses'].copy()
        else:
            return atomic_masses[self.arrays['numbers']]
        
    def set_initial_magnetic_moments(self, magmoms=None):
        """Set the initial magnetic moments.

        Use either one or three numbers for every atom (collinear
        or non-collinear spins)."""
        
        if magmoms is None:
            self.set_array('magmoms', None)
        else:
            magmoms = np.asarray(magmoms)
            self.set_array('magmoms', magmoms, float, magmoms.shape[1:])

    def get_initial_magnetic_moments(self):
        """Get array of initial magnetic moments."""
        if 'magmoms' in self.arrays:
            return self.arrays['magmoms'].copy()
        else:
            return np.zeros(len(self))

    def get_magnetic_moments(self):
        """Get calculated local magnetic moments."""
        if self._calc is None:
            raise RuntimeError('Atoms object has no calculator.')
        if self._calc.get_spin_polarized():
            return self._calc.get_magnetic_moments(self)
        else:
            return np.zeros(len(self))
        
    def get_magnetic_moment(self):
        """Get calculated total magnetic moment."""
        if self._calc is None:
            raise RuntimeError('Atoms object has no calculator.')
        if self._calc.get_spin_polarized():
            return self._calc.get_magnetic_moment(self)
        else:
            return 0.0

    def set_charges(self, charges):
        """Set charges."""
        self.set_array('charges', charges, float, ())

    def get_charges(self):
        """Get array of charges."""
        if 'charges' in self.arrays:
            return self.arrays['charges'].copy()
        else:
            return np.zeros(len(self))

    def set_positions(self, newpositions):
        """Set positions."""
        positions = self.arrays['positions']
        if self.constraints:
            newpositions = np.asarray(newpositions, float)
            for constraint in self.constraints:
                constraint.adjust_positions(positions, newpositions)
                
        self.set_array('positions', newpositions, shape=(3,))

    def get_positions(self):
        """Get array of positions."""
        return self.arrays['positions'].copy()

    def get_calculation_done(self):
        """Let the calculator calculate its thing,
           using the current input.
           """
        if self.calc is None:
            raise RuntimeError('Atoms object has no calculator.')
        self.calc.initialize(self)
        self.calc.calculate(self)

    def get_potential_energy(self):
        """Calculate potential energy."""
        if self._calc is None:
            raise RuntimeError('Atoms object has no calculator.')
        return self._calc.get_potential_energy(self)

    def get_potential_energies(self):
        """Calculate the potential energies of all the atoms.

        Only available with calculators supporting per-atom energies
        (e.g. classical potentials).
        """
        if self._calc is None:
            raise RuntimeError('Atoms object has no calculator.')
        return self._calc.get_potential_energies(self)

    def get_kinetic_energy(self):
        """Get the kinetic energy."""
        momenta = self.arrays.get('momenta')
        if momenta is None:
            return 0.0
        return 0.5 * np.vdot(momenta, self.get_velocities())

    def get_velocities(self):
        """Get array of velocities."""
        momenta = self.arrays.get('momenta')
        if momenta is None:
            return None
        m = self.arrays.get('masses')
        if m is None:
            m = atomic_masses[self.arrays['numbers']]
        return momenta / m.reshape(-1, 1)

    def get_total_energy(self):
        """Get the total energy - potential plus kinetic energy."""
        return self.get_potential_energy() + self.get_kinetic_energy()

    def get_forces(self, apply_constraint=True):
        """Calculate atomic forces.

        Ask the attached calculator to calculate the forces and apply
        constraints.  Use *apply_constraint=False* to get the raw
        forces."""

        if self._calc is None:
            raise RuntimeError('Atoms object has no calculator.')
        forces = self._calc.get_forces(self)
        if apply_constraint:
            for constraint in self.constraints:
                constraint.adjust_forces(self.arrays['positions'], forces)
        return forces

    def get_max_atom_force(self, apply_constraint=True):
        return np.sqrt((self.get_forces()**2).sum(axis=1).max())

    def get_stress(self):
        """Calculate stress tensor.

        Returns an array of the six independent components of the
        symmetric stress tensor, in the traditional order
        (s_xx, s_yy, s_zz, s_yz, s_xz, s_xy).
        """
        if self._calc is None:
            raise RuntimeError('Atoms object has no calculator.')
        stress = self._calc.get_stress(self)
        shape = getattr(stress, 'shape', None)
        if shape == (3, 3):
            return np.array([stress[0, 0], stress[1, 1], stress[2, 2],
                             stress[1, 2], stress[0, 2], stress[0, 1]])
        else:
            # Hopefully a 6-vector, but don't check in case some weird
            # calculator does something else.
            return stress
    
    def get_stresses(self):
        """Calculate the stress-tensor of all the atoms.

        Only available with calculators supporting per-atom energies and
        stresses (e.g. classical potentials).  Even for such calculators
        there is a certain arbitrariness in defining per-atom stresses.
        """
        if self._calc is None:
            raise RuntimeError('Atoms object has no calculator.')
        return self._calc.get_stresses(self)

    def get_dipole_moment(self):
        """Calculate the electric dipole moment for the atoms object.

        Only available for calculators which has a get_dipole_moment()
        method."""
        
        if self._calc is None:
            raise RuntimeError('Atoms object has no calculator.')
        try:
            dipole = self._calc.get_dipole_moment(self)
        except AttributeError:
            raise AttributeError(
                'Calculator object has no get_dipole_moment method.')
        return dipole
    
    def copy(self):
        """Return a copy."""
        import copy
        atoms = self.__class__(cell=self._cell, pbc=self._pbc, info=self.info)

        atoms.arrays = {}
        for name, a in self.arrays.items():
            atoms.arrays[name] = a.copy()
        atoms.constraints = copy.deepcopy(self.constraints)
        atoms.adsorbate_info = copy.deepcopy(self.adsorbate_info)
        return atoms

    def __len__(self):
        return len(self.arrays['positions'])

    def get_number_of_atoms(self):
        """Returns the number of atoms.

        Equivalent to len(atoms) in the standard ASE Atoms class.
        """
        return len(self)
    
    def __repr__(self):
        num = self.get_atomic_numbers()
        N = len(num)
        if N == 0:
            symbols = ''
        elif N <= 60:
            symbols = self.get_chemical_symbols(reduce=True)
        else:
            symbols = ''.join([chemical_symbols[Z] for Z in num[:15]]) + '...'
        s = "%s(symbols='%s', " % (self.__class__.__name__, symbols)
        for name in self.arrays:
            if name == 'numbers':
                continue
            s += '%s=..., ' % name
        if (self._cell - np.diag(self._cell.diagonal())).any():
            s += 'cell=%s, ' % self._cell.tolist()            
        else:
            s += 'cell=%s, ' % self._cell.diagonal().tolist()
        s += 'pbc=%s, ' % self._pbc.tolist()
        if len(self.constraints) == 1:
            s += 'constraint=%s, ' % repr(self.constraints[0])
        if len(self.constraints) > 1:
            s += 'constraint=%s, ' % repr(self.constraints)
        if self._calc is not None:
            s += 'calculator=%s(...), ' % self._calc.__class__.__name__
        return s[:-2] + ')'

    def __add__(self, other):
        atoms = self.copy()
        atoms += other
        return atoms

    def extend(self, other):
        """Extend atoms object by appending atoms from *other*."""
        if isinstance(other, Atom):
            other = self.__class__([other])
            
        n1 = len(self)
        n2 = len(other)
        
        for name, a1 in self.arrays.items():
            a = np.zeros((n1 + n2,) + a1.shape[1:], a1.dtype)
            a[:n1] = a1
            a2 = other.arrays.get(name)
            if a2 is not None:
                a[n1:] = a2
            self.arrays[name] = a

        for name, a2 in other.arrays.items():
            if name in self.arrays:
                continue
            a = np.empty((n1 + n2,) + a2.shape[1:], a2.dtype)
            a[n1:] = a2
            if name == 'masses':
                a[:n1] = self.get_masses()
            else:
                a[:n1] = 0

            self.set_array(name, a)

        return self

    __iadd__ = extend

    def append(self, atom):
        """Append atom to end."""
        self.extend(self.__class__([atom]))

    def __getitem__(self, i):
        """Return a subset of the atoms.

        i -- scalar integer, list of integers, or slice object
        describing which atoms to return.

        If i is a scalar, return an Atom object. If i is a list or a
        slice, return an Atoms object with the same cell, pbc, and
        other associated info as the original Atoms object. The
        indices of the constraints will be shuffled so that they match
        the indexing in the subset returned.

        """
        if isinstance(i, int):
            natoms = len(self)
            if i < -natoms or i >= natoms:
                raise IndexError('Index out of range.')

            return Atom(atoms=self, index=i)
        
        import copy
        
        atoms = self.__class__(cell=self._cell, pbc=self._pbc, info=self.info)
        # TODO: Do we need to shuffle indices in adsorbate_info too?
        atoms.adsorbate_info = self.adsorbate_info
        
        atoms.arrays = {}
        for name, a in self.arrays.items():
            atoms.arrays[name] = a[i].copy()
        
        # Constraints need to be deepcopied, since we need to shuffle
        # the indices
        atoms.constraints = copy.deepcopy(self.constraints)
        condel = []
        for con in atoms.constraints:
            if isinstance(con, FixConstraint):
                try:
                    con.index_shuffle(i)
                except IndexError:
                    condel.append(con)
        for con in condel:
            atoms.constraints.remove(con)
        return atoms

    def __delitem__(self, i):
        check_constraint = np.array([isinstance(c, FixAtoms)
                                     for c in self._constraints])
        if len(self._constraints) > 0 and not check_constraint.all():
            raise RuntimeError('Remove constraint using set_constraint() ' +
                               'before deleting atoms.')
        mask = np.ones(len(self), bool)
        mask[i] = False
        for name, a in self.arrays.items():
            self.arrays[name] = a[mask]
        if len(self._constraints) > 0:
            for n in range(len(self._constraints)):
                self._constraints[n].delete_atom(range(len(mask))[i])

    def pop(self, i=-1):
        """Remove and return atom at index *i* (default last)."""
        atom = self[i]
        atom.cut_reference_to_atoms()
        del self[i]
        return atom
    
    def __imul__(self, m):
        """In-place repeat of atoms."""
        if isinstance(m, int):
            m = (m, m, m)

        M = np.product(m)
        n = len(self)
        
        for name, a in self.arrays.items():
            self.arrays[name] = np.tile(a, (M,) + (1,) * (len(a.shape) - 1))

        positions = self.arrays['positions']
        i0 = 0
        for m0 in range(m[0]):
            for m1 in range(m[1]):
                for m2 in range(m[2]):
                    i1 = i0 + n
                    positions[i0:i1] += np.dot((m0, m1, m2), self._cell)
                    i0 = i1

        if self.constraints is not None:
            self.constraints = [c.repeat(m, n) for c in self.constraints]

        self._cell = np.array([m[c] * self._cell[c] for c in range(3)])

        return self

    def repeat(self, rep):
        """Create new repeated atoms object.

        The *rep* argument should be a sequence of three positive
        integers like *(2,3,1)* or a single integer (*r*) equivalent
        to *(r,r,r)*."""

        atoms = self.copy()
        atoms *= rep
        return atoms

    __mul__ = repeat
    
    def translate(self, displacement):
        """Translate atomic positions.

        The displacement argument can be a float an xyz vector or an
        nx3 array (where n is the number of atoms)."""

        self.arrays['positions'] += np.array(displacement)

    def center(self, vacuum=None, axis=None):
        """Center atoms in unit cell.

        Centers the atoms in the unit cell, so there is the same
        amount of vacuum on all sides.

        Parameters:

        vacuum (default: None): If specified adjust the amount of
        vacuum when centering.  If vacuum=10.0 there will thus be 10
        Angstrom of vacuum on each side.

        axis (default: None): If specified, only act on the specified
        axis.  Default: Act on all axes.
        """
        # Find the orientations of the faces of the unit cell
        c = self.get_cell()
        dirs = np.zeros_like(c)
        for i in range(3):
            dirs[i] = np.cross(c[i - 1], c[i - 2])
            dirs[i] /= np.sqrt(np.dot(dirs[i], dirs[i]))  # normalize
            if np.dot(dirs[i], c[i]) < 0.0:
                dirs[i] *= -1

        # Now, decide how much each basis vector should be made longer
        if axis is None:
            axes = (0, 1, 2)
        else:
            axes = (axis,)
        p = self.arrays['positions']
        longer = np.zeros(3)
        shift = np.zeros(3)
        for i in axes:
            p0 = np.dot(p, dirs[i]).min()
            p1 = np.dot(p, dirs[i]).max()
            height = np.dot(c[i], dirs[i])
            if vacuum is not None:
                lng = (p1 - p0 + 2 * vacuum) - height
            else:
                lng = 0.0  # Do not change unit cell size!
            top = lng + height - p1
            shf = 0.5 * (top - p0)
            cosphi = np.dot(c[i], dirs[i]) / np.sqrt(np.dot(c[i], c[i]))
            longer[i] = lng / cosphi
            shift[i] = shf / cosphi

        # Now, do it!
        translation = np.zeros(3)
        for i in axes:
            nowlen = np.sqrt(np.dot(c[i], c[i]))
            self._cell[i] *= 1 + longer[i] / nowlen
            translation += shift[i] * c[i] / nowlen
        self.arrays['positions'] += translation

    def get_center_of_mass(self, scaled=False):
        """Get the center of mass.

        If scaled=True the center of mass in scaled coordinates
        is returned."""
        m = self.arrays.get('masses')
        if m is None:
            m = atomic_masses[self.arrays['numbers']]
        com = np.dot(m, self.arrays['positions']) / m.sum()
        if scaled:
            return np.linalg.solve(self._cell.T, com)
        else:
            return com

    def get_moments_of_inertia(self, vectors=False):
        """Get the moments of inertia along the principal axes.

        The three principal moments of inertia are computed from the
        eigenvalues of the symmetric inertial tensor. Periodic boundary
        conditions are ignored. Units of the moments of inertia are
        amu*angstrom**2.
        """
        com = self.get_center_of_mass()
        positions = self.get_positions()
        positions -= com  # translate center of mass to origin
        masses = self.get_masses()

        #initialize elements of the inertial tensor
        I11 = I22 = I33 = I12 = I13 = I23 = 0.0
        for i in range(len(self)):
            x, y, z = positions[i]
            m = masses[i]

            I11 += m * (y**2 + z**2)
            I22 += m * (x**2 + z**2)
            I33 += m * (x**2 + y**2)
            I12 += -m * x * y
            I13 += -m * x * z
            I23 += -m * y * z

        I = np.array([[I11, I12, I13],
                      [I12, I22, I23],
                      [I13, I23, I33]])

        evals, evecs = np.linalg.eigh(I)
        if vectors:
            return evals, evecs.transpose()
        else:
            return evals

    def get_angular_momentum(self):
        """Get total angular momentum with respect to the center of mass."""
        com = self.get_center_of_mass()
        positions = self.get_positions()
        positions -= com  # translate center of mass to origin
        return np.cross(positions, self.get_momenta()).sum(0)

    def rotate(self, v, a=None, center=(0, 0, 0), rotate_cell=False):
        """Rotate atoms.

        Rotate the angle *a* around the vector *v*.  If *a* is not
        given, the length of *v* is used as the angle.  If *a* is a
        vector, then *v* is rotated into *a*.  The point at *center*
        is fixed.  Use *center='COM'* to fix the center of mass.
        Vectors can also be strings: 'x', '-x', 'y', ... .

        Examples:

        Rotate 90 degrees around the z-axis, so that the x-axis is
        rotated into the y-axis:

        >>> a = pi / 2
        >>> atoms.rotate('z', a)
        >>> atoms.rotate((0, 0, 1), a)
        >>> atoms.rotate('-z', -a)
        >>> atoms.rotate((0, 0, a))
        >>> atoms.rotate('x', 'y')
        """

        norm = np.linalg.norm
        v = string2vector(v)
        if a is None:
            a = norm(v)
        if isinstance(a, (float, int)):
            v /= norm(v)
            c = cos(a)
            s = sin(a)
        else:
            v2 = string2vector(a)
            v /= norm(v)
            v2 /= norm(v2)
            c = np.dot(v, v2)
            v = np.cross(v, v2)
            s = norm(v)
            # In case *v* and *a* are parallel, np.cross(v, v2) vanish
            # and can't be used as a rotation axis. However, in this
            # case any rotation axis perpendicular to v2 will do.
            eps = 1e-7
            if s < eps:
                v = np.cross((0, 0, 1), v2)
            if norm(v) < eps:
                v = np.cross((1, 0, 0), v2)
            assert norm(v) >= eps
            if s > 0:
                v /= s
        
        if isinstance(center, str) and center.lower() == 'com':
            center = self.get_center_of_mass()

        p = self.arrays['positions'] - center
        self.arrays['positions'][:] = (c * p - 
                                       np.cross(p, s * v) + 
                                       np.outer(np.dot(p, v), (1.0 - c) * v) +
                                       center)
        if rotate_cell:
            rotcell = self.get_cell()
            rotcell[:] = (c * rotcell - 
                          np.cross(rotcell, s * v) + 
                          np.outer(np.dot(rotcell, v), (1.0 - c) * v))
            self.set_cell(rotcell)
                
    def rotate_euler(self, center=(0, 0, 0), phi=0.0, theta=0.0, psi=0.0):
        """Rotate atoms via Euler angles.
        
        See e.g http://mathworld.wolfram.com/EulerAngles.html for explanation.
        
        Parameters:
        
        center :
            The point to rotate about. A sequence of length 3 with the
            coordinates, or 'COM' to select the center of mass.
        phi :
            The 1st rotation angle around the z axis.
        theta :
            Rotation around the x axis.
        psi :
            2nd rotation around the z axis.
        
        """
        if isinstance(center, str) and center.lower() == 'com':
            center = self.get_center_of_mass()
        else:
            center = np.array(center)
        # First move the molecule to the origin In contrast to MATLAB,
        # numpy broadcasts the smaller array to the larger row-wise,
        # so there is no need to play with the Kronecker product.
        rcoords = self.positions - center
        # First Euler rotation about z in matrix form
        D = np.array(((cos(phi), sin(phi), 0.),
                      (-sin(phi), cos(phi), 0.),
                      (0., 0., 1.)))
        # Second Euler rotation about x:
        C = np.array(((1., 0., 0.),
                      (0., cos(theta), sin(theta)),
                      (0., -sin(theta), cos(theta))))
        # Third Euler rotation, 2nd rotation about z:
        B = np.array(((cos(psi), sin(psi), 0.),
                      (-sin(psi), cos(psi), 0.),
                      (0., 0., 1.)))
        # Total Euler rotation
        A = np.dot(B, np.dot(C, D))
        # Do the rotation
        rcoords = np.dot(A, np.transpose(rcoords))
        # Move back to the rotation point
        self.positions = np.transpose(rcoords) + center

    def get_dihedral(self, list):
        """Calculate dihedral angle.

        Calculate dihedral angle between the vectors list[0]->list[1]
        and list[2]->list[3], where list contains the atomic indexes
        in question.
        """

        # vector 0->1, 1->2, 2->3 and their normalized cross products:
        a = self.positions[list[1]] - self.positions[list[0]]
        b = self.positions[list[2]] - self.positions[list[1]]
        c = self.positions[list[3]] - self.positions[list[2]]
        bxa = np.cross(b, a)
        bxa /= np.linalg.norm(bxa)
        cxb = np.cross(c, b)
        cxb /= np.linalg.norm(cxb)
        angle = np.vdot(bxa, cxb)
        # check for numerical trouble due to finite precision:
        if angle < -1:
            angle = -1
        if angle > 1:
            angle = 1
        angle = np.arccos(angle)
        if np.vdot(bxa, c) > 0:
            angle = 2 * np.pi - angle
        return angle

    def _masked_rotate(self, center, axis, diff, mask):
        # do rotation of subgroup by copying it to temporary atoms object
        # and then rotating that
        #
        # recursive object definition might not be the most elegant thing,
        # more generally useful might be a rotation function with a mask?
        group = self.__class__()
        for i in range(len(self)):
            if mask[i]:
                group += self[i]
        group.translate(-center)
        group.rotate(axis, diff)
        group.translate(center)
        # set positions in original atoms object
        j = 0
        for i in range(len(self)):
            if mask[i]:
                self.positions[i] = group[j].get_position()
                j += 1

    def set_dihedral(self, list, angle, mask=None):
        """
        set the dihedral angle between vectors list[0]->list[1] and 
        list[2]->list[3] by changing the atom indexed by list[3]
        if mask is not None, all the atoms described in mask 
        (read: the entire subgroup) are moved
        
        example: the following defines a very crude 
        ethane-like molecule and twists one half of it by 30 degrees.

        >>> atoms = Atoms('HHCCHH', [[-1, 1, 0], [-1, -1, 0], [0, 0, 0],
                                     [1, 0, 0], [2, 1, 0], [2, -1, 0]])
        >>> atoms.set_dihedral([1,2,3,4],7*pi/6,mask=[0,0,0,1,1,1])
        """
        # if not provided, set mask to the last atom in the
        # dihedral description
        if mask is None:
            mask = np.zeros(len(self))
            mask[list[3]] = 1
        # compute necessary in dihedral change, from current value
        current = self.get_dihedral(list)
        diff = angle - current
        axis = self.positions[list[2]] - self.positions[list[1]]
        center = self.positions[list[2]]
        self._masked_rotate(center, axis, diff, mask)
        
    def rotate_dihedral(self, list, angle, mask=None):
        """Rotate dihedral angle.

        Complementing the two routines above: rotate a group by a
        predefined dihedral angle, starting from its current
        configuration
        """
        start = self.get_dihedral(list)
        self.set_dihedral(list, angle + start, mask)
    
    def get_angle(self, list):
        """Get angle formed by three atoms.
        
        calculate angle between the vectors list[0]->list[1] and
        list[1]->list[2], where list contains the atomic indexes in
        question."""
        # normalized vector 1->0, 1->2:
        v10 = self.positions[list[0]] - self.positions[list[1]]
        v12 = self.positions[list[2]] - self.positions[list[1]]
        v10 /= np.linalg.norm(v10)
        v12 /= np.linalg.norm(v12)
        angle = np.vdot(v10, v12)
        angle = np.arccos(angle)
        return angle

    def set_angle(self, list, angle, mask=None):
        """Set angle formed by three atoms.
        
        Sets the angle between vectors list[1]->list[0] and 
        list[1]->list[2].

        Same usage as in set_dihedral."""
        # If not provided, set mask to the last atom in the angle description
        if mask is None:
            mask = np.zeros(len(self))
            mask[list[2]] = 1
        # Compute necessary in angle change, from current value
        current = self.get_angle(list)
        diff = current - angle
        # Do rotation of subgroup by copying it to temporary atoms object and
        # then rotating that
        v10 = self.positions[list[0]] - self.positions[list[1]]
        v12 = self.positions[list[2]] - self.positions[list[1]]
        v10 /= np.linalg.norm(v10)
        v12 /= np.linalg.norm(v12)
        axis = np.cross(v10, v12)
        center = self.positions[list[1]]
        self._masked_rotate(center, axis, diff, mask)
    
    def rattle(self, stdev=0.001, seed=None):
        """Randomly displace atoms.

        This method adds random displacements to the atomic positions,
        taking a possible constraint into account.  The random numbers are
        drawn from a normal distribution of standard deviation stdev.

        For a parallel calculation, it is important to use the same
        seed on all processors!  """
        
        rs = np.random.RandomState(seed)
        positions = self.arrays['positions']
        self.set_positions(positions +
                           rs.normal(scale=stdev, size=positions.shape))
        
    def get_distance(self, a0, a1, mic=False):
        """Return distance between two atoms.

        Use mic=True to use the Minimum Image Convention.
        """

        R = self.arrays['positions']
        D = R[a1] - R[a0]
        if mic:
            Dr = np.linalg.solve(self._cell.T, D)
            D = np.dot(Dr - np.round(Dr) * self._pbc, self._cell)
        return np.linalg.norm(D)

    def set_distance(self, a0, a1, distance, fix=0.5):
        """Set the distance between two atoms.

        Set the distance between atoms *a0* and *a1* to *distance*.
        By default, the center of the two atoms will be fixed.  Use
        *fix=0* to fix the first atom, *fix=1* to fix the second
        atom and *fix=0.5* (default) to fix the center of the bond."""

        R = self.arrays['positions']
        D = R[a1] - R[a0]
        x = 1.0 - distance / np.linalg.norm(D)
        R[a0] += (x * fix) * D
        R[a1] -= (x * (1.0 - fix)) * D

    def get_scaled_positions(self):
        """Get positions relative to unit cell.

        Atoms outside the unit cell will be wrapped into the cell in
        those directions with periodic boundary conditions so that the
        scaled coordinates are between zero and one."""

        scaled = np.linalg.solve(self._cell.T, self.arrays['positions'].T).T
        for i in range(3):
            if self._pbc[i]:
                # Yes, we need to do it twice.
                # See the scaled_positions.py test
                scaled[:, i] %= 1.0
                scaled[:, i] %= 1.0
        return scaled

    def set_scaled_positions(self, scaled):
        """Set positions relative to unit cell."""
        self.arrays['positions'][:] = np.dot(scaled, self._cell)

    def get_temperature(self):
        """Get the temperature. in Kelvin"""
        ekin = self.get_kinetic_energy() / len(self)
        return ekin / (1.5 * units.kB)

    def get_isotropic_pressure(self, stress):
        """Get the current calculated pressure, assume isotropic medium.
            in Bar
        """
        if type(stress) == type(1.0) or type(stress) == type(1):
            return -stress * 1e-5 / units.Pascal
        elif stress.shape == (3, 3):
            return (-(stress[0, 0] + stress[1, 1] + stress[2, 2]) / 3.0) * \
                    1e-5 / units.Pascal
        elif stress.shape == (6,):
            return (-(stress[0] + stress[1] + stress[2]) / 3.0) * \
                   1e-5 / units.Pascal
        else:
            raise ValueError('The external stress has the wrong shape.')

    def __eq__(self, other):
        """Check for identity of two atoms objects.

        Identity means: same positions, atomic numbers, unit cell and
        periodic boundary conditions."""
        try:
            a = self.arrays
            b = other.arrays
            return (len(self) == len(other) and
                    (a['positions'] == b['positions']).all() and
                    (a['numbers'] == b['numbers']).all() and
                    (self._cell == other.cell).all() and
                    (self._pbc == other.pbc).all())
        except AttributeError:
            return NotImplemented

    def __ne__(self, other):
        eq = self.__eq__(other)
        if eq is NotImplemented:
            return eq
        else:
            return not eq

    __hash__ = None

    def get_volume(self):
        """Get volume of unit cell."""
        return abs(np.linalg.det(self._cell))
    
    def _get_positions(self):
        """Return reference to positions-array for in-place manipulations."""
        return self.arrays['positions']

    def _set_positions(self, pos):
        """Set positions directly, bypassing constraints."""
        self.arrays['positions'][:] = pos

    positions = property(_get_positions, _set_positions,
                         doc='Attribute for direct ' +
                         'manipulation of the positions.')

    def _get_atomic_numbers(self):
        """Return reference to atomic numbers for in-place 
        manipulations."""
        return self.arrays['numbers']

    numbers = property(_get_atomic_numbers, set_atomic_numbers,
                       doc='Attribute for direct ' +
                       'manipulation of the atomic numbers.')

    def _get_cell(self):
        """Return reference to unit cell for in-place manipulations."""
        return self._cell
    
    cell = property(_get_cell, set_cell, doc='Attribute for direct ' +
                       'manipulation of the unit cell.')

    def _get_pbc(self):
        """Return reference to pbc-flags for in-place manipulations."""
        return self._pbc
    
    pbc = property(_get_pbc, set_pbc,
                   doc='Attribute for direct manipulation ' +
                   'of the periodic boundary condition flags.')

    def get_name(self):
        """Return a name extracted from the elements."""
        elements = {}
        for a in self:
            try:
                elements[a.symbol] += 1
            except:
                elements[a.symbol] = 1
        name = ''
        for element in elements:
            name += element
            if elements[element] > 1:
                name += str(elements[element])
        return name

    def write(self, filename, format=None):

        if format == None:
            format = self.format
        if format == 'vasp':
            write_vasp(filename, self)
        elif format == 'xyz':
            write_xyz(filename, self)
        elif format == 'con':
            write_con(filename, self)
        else:
            raise Exception, "Unknown file format: %s" % format
        

def string2symbols(s):
    """Convert string to list of chemical symbols."""
    n = len(s)

    if n == 0:
        return []
    
    c = s[0]
    
    if c.isdigit():
        i = 1
        while i < n and s[i].isdigit():
            i += 1
        return int(s[:i]) * string2symbols(s[i:])

    if c == '(':
        p = 0
        for i, c in enumerate(s):
            if c == '(':
                p += 1
            elif c == ')':
                p -= 1
                if p == 0:
                    break
        j = i + 1
        while j < n and s[j].isdigit():
            j += 1
        if j > i + 1:
            m = int(s[i + 1:j])
        else:
            m = 1
        return m * string2symbols(s[1:i]) + string2symbols(s[j:])

    if c.isupper():
        i = 1
        if 1 < n and s[1].islower():
            i += 1
        j = i
        while j < n and s[j].isdigit():
            j += 1
        if j > i:
            m = int(s[i:j])
        else:
            m = 1
        return m * [s[:i]] + string2symbols(s[j:])
    else:
        raise ValueError


def symbols2numbers(symbols):
    if isinstance(symbols, str):
        symbols = string2symbols(symbols)
    numbers = []
    for s in symbols:
        if isinstance(s, str):
            numbers.append(atomic_numbers[s])
        else:
            numbers.append(s)
    return numbers


def string2vector(v):
    if isinstance(v, str):
        if v[0] == '-':
            return -string2vector(v[1:])
        w = np.zeros(3)
        w['xyz'.index(v)] = 1.0
        return w
    return np.array(v, float)


def default(data, dflt):
    """Helper function for setting default values."""
    if data is None:
        return None
    elif isinstance(data, (list, tuple)):
        newdata = []
        allnone = True
        for x in data:
            if x is None:
                newdata.append(dflt)
            else:
                newdata.append(x)
                allnone = False
        if allnone:
            return None
        return newdata
    else:
        return data

def slice2enlist(s):
    """Convert a slice object into a list of (new, old) tuples."""
    if isinstance(s, (list, tuple)):
        return enumerate(s)
    if s.step == None:
        step = 1
    else:
        step = s.step
    if s.start == None:
        start = 0
    else:
        start = s.start
    return enumerate(range(start, s.stop, step))


class FixConstraint:
    """Base class for classes that fix one or more atoms in some way."""

    def index_shuffle(self, ind):
        """Change the indices.

        When the ordering of the atoms in the Atoms object changes,
        this method can be called to shuffle the indices of the
        constraints.

        ind -- List or tuple of indices.

        """
        raise NotImplementedError

    def repeat(self, m, n):
        """ basic method to multiply by m, needs to know the length
        of the underlying atoms object for the assignment of
        multiplied constraints to work.
        """
        raise NotImplementedError


class FixConstraintSingle(FixConstraint):
    """Base class for classes that fix a single atom."""

    def index_shuffle(self, ind):
        """The atom index must be stored as self.a."""
        newa = -1   # Signal error
        for new, old in slice2enlist(ind):
            if old == self.a:
                newa = new
                break
        if newa == -1:
            raise IndexError('Constraint not part of slice')
        self.a = newa


class FixAtoms(FixConstraint):
    """Constraint object for fixing some chosen atoms."""
    def __init__(self, indices=None, mask=None):
        """Constrain chosen atoms.

        Parameters
        ----------
        indices : list of int
           Indices for those atoms that should be constrained.
        mask : list of bool
           One boolean per atom indicating if the atom should be
           constrained or not.

        Examples
        --------
        Fix all Copper atoms:

        >>> c = FixAtoms(mask=[s == 'Cu' for s in atoms.get_chemical_symbols()])
        >>> atoms.set_constraint(c)

        Fix all atoms with z-coordinate less than 1.0 Angstrom:

        >>> c = FixAtoms(mask=atoms.positions[:, 2] < 1.0)
        >>> atoms.set_constraint(c)
        """

        if indices is None and mask is None:
            raise ValueError('Use "indices" or "mask".')
        if indices is not None and mask is not None:
            raise ValueError('Use only one of "indices" and "mask".')

        if mask is not None:
            self.index = np.asarray(mask, bool)
        else:
            # Check for duplicates
            srt = np.sort(indices)
            for i in range(len(indices) - 1):
                if srt[i] == srt[i+1]:
                    raise ValueError(
                        'FixAtoms: The indices array contained duplicates. '
                        'Perhaps you wanted to specify a mask instead, but '
                        'forgot the mask= keyword.')
            self.index = np.asarray(indices, int)

        if self.index.ndim != 1:
            raise ValueError('Wrong argument to FixAtoms class!')

    def adjust_positions(self, old, new):
        new[self.index] = old[self.index]

    def adjust_forces(self, positions, forces):
        forces[self.index] = 0.0

    def index_shuffle(self, ind):
        # See docstring of superclass
        if self.index.dtype == bool:
            self.index = self.index[ind]
        else:
            index = []
            for new, old in slice2enlist(ind):
                if old in self.index:
                    index.append(new)
            if len(index) == 0:
                raise IndexError('All indices in FixAtoms not part of slice')
            self.index = np.asarray(index, int)

    def copy(self):
        if self.index.dtype == bool:
            return FixAtoms(mask=self.index.copy())
        else:
            return FixAtoms(indices=self.index.copy())

    def __repr__(self):
        if self.index.dtype == bool:
            return 'FixAtoms(mask=%s)' % ints2string(self.index.astype(int))
        return 'FixAtoms(indices=%s)' % ints2string(self.index)

    def repeat(self, m, n):
        i0 = 0
        l = len(self.index)
        natoms = 0
        if isinstance(m, int):
            m = (m, m, m)
        index_new = []
        for m2 in range(m[2]):
            for m1 in range(m[1]):
                for m0 in range(m[0]):
                    i1 = i0 + n
                    if self.index.dtype == bool:
                        index_new.extend(self.index)
                    else:
                        index_new += [i+natoms for i in self.index]
                    i0 = i1
                    natoms += n
        if self.index.dtype == bool:
            self.index = np.asarray(index_new, bool)
        else:
            self.index = np.asarray(index_new, int)
        return self

    def delete_atom(self, ind):
        """ Removes atom number ind from the index array, if present.
        Required for removing atoms with existing FixAtoms constraints.
        """
        if self.index.dtype == bool:
            self.index = np.delete(self.index, ind)
        else:
            if ind in self.index:
                i = list(self.index).index(ind)
                self.index = np.delete(self.index, i)
            for i in range(len(self.index)):
                if self.index[i] >= ind:
                    self.index[i] -= 1

def ints2string(x, threshold=10):
    """Convert ndarray of ints to string."""
    if len(x) <= threshold:
        return str(x.tolist())
    return str(x[:threshold].tolist())[:-1] + ', ...]'

class FixBondLengths(FixConstraint):
    def __init__(self, pairs, iterations=10):
        self.constraints = [FixBondLength(a1, a2)
                            for a1, a2 in pairs]
        self.iterations = iterations

    def adjust_positions(self, old, new):
        for i in range(self.iterations):
            for constraint in self.constraints:
                constraint.adjust_positions(old, new)

    def adjust_forces(self, positions, forces):
        for i in range(self.iterations):
            for constraint in self.constraints:
                constraint.adjust_forces(positions, forces)

    def copy(self):
        return FixBondLengths([constraint.indices
                               for constraint in self.constraints])

class FixBondLength(FixConstraint):
    """Constraint object for fixing a bond length."""
    def __init__(self, a1, a2):
        """Fix distance between atoms with indices a1 and a2."""
        self.indices = [a1, a2]

    def adjust_positions(self, old, new):
        p1, p2 = old[self.indices]
        d = p2 - p1
        p = sqrt(np.dot(d, d))
        q1, q2 = new[self.indices]
        d = q2 - q1
        q = sqrt(np.dot(d, d))
        d *= 0.5 * (p - q) / q
        new[self.indices] = (q1 - d, q2 + d)

    def adjust_forces(self, positions, forces):
        d = np.subtract.reduce(positions[self.indices])
        d2 = np.dot(d, d)
        d *= 0.5 * np.dot(np.subtract.reduce(forces[self.indices]), d) / d2
        forces[self.indices] += (-d, d)

    def index_shuffle(self, ind):
        'Shuffle the indices of the two atoms in this constraint'
        newa = [-1, -1] # Signal error
        for new, old in slice2enlist(ind):
            for i, a in enumerate(self.indices):
                if old == a:
                    newa[i] = new
        if newa[0] == -1 or newa[1] == -1:
            raise IndexError('Constraint not part of slice')
        self.indices = newa

    def copy(self):
        return FixBondLength(*self.indices)

    def __repr__(self):
        return 'FixBondLength(%d, %d)' % tuple(self.indices)

class FixedMode(FixConstraint):
    """Constrain atoms to move along directions orthogonal to
    a given mode only."""

    def __init__(self, indices, mode):
        if indices is None:
            raise ValueError('Use "indices".')
        if indices is not None:
            self.index = np.asarray(indices, int)
        self.mode = (np.asarray(mode) / np.sqrt((mode **2).sum())).reshape(-1)

    def adjust_positions(self, oldpositions, newpositions):
        newpositions = newpositions.ravel()
        oldpositions = oldpositions.ravel()
        step = newpositions - oldpositions
        newpositions -= self.mode * np.dot(step, self.mode)
        newpositions = newpositions.reshape(-1, 3)
        oldpositions = oldpositions.reshape(-1, 3)

    def adjust_forces(self, positions, forces):
        forces = forces.ravel()
        forces -= self.mode * np.dot(forces, self.mode)
        forces = forces.reshape(-1, 3)

    def copy(self):
        return FixedMode(self.index.copy(), self.mode)

    def __repr__(self):
        return 'FixedMode(%s, %s)' % (ints2string(self.index),
                                      self.mode.tolist())

class FixedPlane(FixConstraintSingle):
    """Constrain an atom *a* to move in a given plane only.

    The plane is defined by its normal: *direction*."""

    def __init__(self, a, direction):
        self.a = a
        self.dir = np.asarray(direction) / sqrt(np.dot(direction, direction))

    def adjust_positions(self, oldpositions, newpositions):
        step = newpositions[self.a] - oldpositions[self.a]
        newpositions[self.a] -= self.dir * np.dot(step, self.dir)

    def adjust_forces(self, positions, forces):
        forces[self.a] -= self.dir * np.dot(forces[self.a], self.dir)

    def copy(self):
        return FixedPlane(self.a, self.dir)

    def __repr__(self):
        return 'FixedPlane(%d, %s)' % (self.a, self.dir.tolist())


class FixedLine(FixConstraintSingle):
    """Constrain an atom *a* to move on a given line only.

    The line is defined by its *direction*."""

    def __init__(self, a, direction):
        self.a = a
        self.dir = np.asarray(direction) / sqrt(np.dot(direction, direction))

    def adjust_positions(self, oldpositions, newpositions):
        step = newpositions[self.a] - oldpositions[self.a]
        x = np.dot(step, self.dir)
        newpositions[self.a] = oldpositions[self.a] + x * self.dir

    def adjust_forces(self, positions, forces):
        forces[self.a] = self.dir * np.dot(forces[self.a], self.dir)

    def copy(self):
        return FixedLine(self.a, self.dir)

    def __repr__(self):
        return 'FixedLine(%d, %s)' % (self.a, self.dir.tolist())

class FixCartesian(FixConstraintSingle):
    "Fix an atom in the directions of the cartesian coordinates."
    def __init__(self, a, mask=(1, 1, 1)):
        self.a = a
        self.mask = -(np.array(mask) - 1)

    def adjust_positions(self, old, new):
        step = new[self.a] - old[self.a]
        step *= self.mask
        new[self.a] = old[self.a] + step

    def adjust_forces(self, positions, forces):
        forces[self.a] *= self.mask

    def copy(self):
        return FixCartesian(self.a, 1 - self.mask)

    def __repr__(self):
        return 'FixCartesian(indice=%s mask=%s)' % (self.a, self.mask)

class fix_cartesian(FixCartesian):
    "Backwards compatibility for FixCartesian."
    def __init__(self, a, mask=(1, 1, 1)):
        import warnings
        super(fix_cartesian, self).__init__(a, mask)
        warnings.warn('fix_cartesian is deprecated. Please use FixCartesian'
                      ' instead.', DeprecationWarning, stacklevel=2)

class FixScaled(FixConstraintSingle):
    "Fix an atom in the directions of the unit vectors."
    def __init__(self, cell, a, mask=(1, 1, 1)):
        self.cell = cell
        self.a = a
        self.mask = np.array(mask)

    def adjust_positions(self, old, new):
        scaled_old = np.linalg.solve(self.cell.T, old.T).T
        scaled_new = np.linalg.solve(self.cell.T, new.T).T
        for n in range(3):
            if self.mask[n]:
                scaled_new[self.a, n] = scaled_old[self.a, n]
        new[self.a] = np.dot(scaled_new, self.cell)[self.a]

    def adjust_forces(self, positions, forces):
        scaled_forces = np.linalg.solve(self.cell.T, forces.T).T
        scaled_forces[self.a] *= -(self.mask - 1)
        forces[self.a] = np.dot(scaled_forces, self.cell)[self.a]

    def copy(self):
        return FixScaled(self.cell, self.a, self.mask)

    def __repr__(self):
        return 'FixScaled(%s, %d, %s)' % (repr(self.cell),
                                          self.a,
                                          repr(self.mask))

class fix_scaled(FixScaled):
    "Backwards compatibility for FixScaled."
    def __init__(self, cell, a, mask=(1, 1, 1)):
        import warnings
        super(fix_scaled, self).__init__(cell, a, mask)
        warnings.warn('fix_scaled is deprecated. Please use FixScaled '
                      'instead.', DeprecationWarning, stacklevel=2)

class SinglePointCalculator:
    """Special calculator for a single configuration.

    Used to remember the energy, force and stress for a given
    configuration.  If the positions, atomic numbers, unit cell, or
    boundary conditions are changed, then asking for
    energy/forces/stress will raise an exception."""
    
    def __init__(self, energy, forces, stress, magmoms, atoms):
        """Save energy, forces and stresses for the current configuration."""
        self.energy = energy
        if forces is not None:
            forces = np.array(forces, float)
        self.forces = forces
        if stress is not None:
            stress = np.array(stress, float)
        self.stress = stress
        if magmoms is not None:
            magmoms = np.array(magmoms, float)
        self.magmoms = magmoms
        self.atoms = atoms.copy()

    def calculation_required(self, atoms, quantities):
        ok = self.atoms == atoms
        return ('forces' in quantities and (self.forces is None or not ok) or
                'energy' in quantities and (self.energy is None or not ok) or
                'stress' in quantities and (self.stress is None or not ok) or
                'magmoms' in quantities and (self.magmoms is None or not ok))

    def update(self, atoms):
        if self.atoms != atoms:
            raise RuntimeError('Energy, forces and stress no longer correct.')

    def get_potential_energy(self, atoms=None):
        if atoms is not None:
            self.update(atoms)
        if self.energy is None:
            raise RuntimeError('No energy.')
        return self.energy

    def get_forces(self, atoms):
        self.update(atoms)
        if self.forces is None:
            raise RuntimeError('No forces.')
        return self.forces

    def get_stress(self, atoms):
        self.update(atoms)
        if self.stress is None:
            raise NotImplementedError
        return self.stress

    def get_spin_polarized(self):
        return self.magmoms is not None and self.magmoms.any()

    def get_magnetic_moments(self, atoms=None):
        if atoms is not None:
            self.update(atoms)
        if self.magmoms is not None:
            return self.magmoms
        else:
            return np.zeros(len(self.positions))

class SinglePointKPoint:
    def __init__(self, kpt, spin):
        self.k = kpt
        self.s = spin
        self.eps_n = []
        self.f_n = []

class SinglePointDFTCalculator(SinglePointCalculator):
    def __init__(self, energy, forces, stress, magmoms, atoms,
                 eFermi=None):
        SinglePointCalculator.__init__(self, energy, forces, stress, 
                                       magmoms, atoms)
        if eFermi is not None:
            self.eFermi = eFermi
        self.kpts = None

    def get_fermi_level(self):
        """Return the Fermi-level(s)."""
        return self.eFermi

    def get_bz_k_points(self):
        """Return the k-points."""
        if self.kpts is not None:
            # we assume that only the gamma point is defined
            return np.zeros((1, 3))
        return None

    def get_number_of_spins(self):
        """Return the number of spins in the calculation.

        Spin-paired calculations: 1, spin-polarized calculation: 2."""
        if self.kpts is not None:
            # we assume that only the gamma point is defined
            return len(self.kpts)
        return None

    def get_spin_polarized(self):
        """Is it a spin-polarized calculation?"""
        nos = self.get_number_of_spins()
        if nos is not None:
            return nos == 2
        return None
    
    def get_ibz_k_points(self):
        """Return k-points in the irreducible part of the Brillouin zone."""
        return self.get_bz_k_points()

    def get_occupation_numbers(self, kpt=0, spin=0):
        """Return occupation number array."""
        # we assume that only the gamma point is defined
        assert(kpt == 0)
        if self.kpts is not None:
            for kpt in self.kpts:
                if kpt.s == spin:
                    return kpt.f_n
        return None

    def get_eigenvalues(self, kpt=0, spin=0):
        """Return eigenvalue array."""
        # we assume that only the gamma point is defined
        assert(kpt == 0)
        if self.kpts is not None:
            for kpt in self.kpts:
                if kpt.s == spin:
                    return kpt.eps_n
        return None
