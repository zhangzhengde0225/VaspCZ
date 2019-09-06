import glob
import os
import numpy
import sys
import re
import subprocess

class Kdb():
    def check_svn_version(self):
        #does not work
        svn_info = subprocess.check_output("svn info")
        print svn_info

    def check_version(self):
        if sys.version_info[:2] != (2,7):
            print "python 2.7.X required."
            return False
        return True


    def atomAtomPbcVector(self, atoms, a, b):
        if not hasattr(atoms, 'ibox'):
            atoms.ibox = numpy.linalg.inv(atoms.get_cell())
        if not hasattr(atoms, 'pbcVectors'):
            atoms.pbcVectors = {}
        if (a, b) not in atoms.pbcVectors or (b, a) not in atoms.pbcVectors:
            atoms.pbcVectors[(a, b)] = self.pbc(atoms.positions[b] - atoms.positions[a], atoms.get_cell(), atoms.ibox)
            atoms.pbcVectors[(b, a)] = -atoms.pbcVectors[(a, b)]
        return atoms.pbcVectors[(a, b)]
                

    def atomAtomPbcDistance(self, atoms, a, b):
        if not hasattr(atoms, 'pbcDistances'):
            atoms.pbcDistances = {}
        if (a, b) not in atoms.pbcDistances or (b, a) not in atoms.pbcDistances:
            atoms.pbcDistances[(a, b)] = numpy.linalg.norm(self.atomAtomPbcVector(atoms, a, b))
            atoms.pbcDistances[(b, a)] = atoms.pbcDistances[(a, b)]
        return atoms.pbcDistances[(a, b)]
            

    def atomAtomDistance(self, atoms, a, b):
        if not hasattr(atoms, 'distances'):
            atoms.distances = {}
        if (a, b) not in atoms.distances or (b, a) not in atoms.distances:
            atoms.distances[(a, b)] = numpy.linalg.norm(atoms.positions[a] - atoms.positions[b])
            atoms.distances[(b, a)] = atoms.distances[(a, b)]
        return atoms.distances[(a, b)]


    def getNameList(self, atoms):
        """
        Returns a sorted list of element names.
        """
        nl = []
        for name in atoms.get_chemical_symbols():
            if name not in nl:
                nl.append(name)
        return sorted(nl)


    def nameCount(self, atoms):
        counts = {}
        for name in atoms.get_chemical_symbols():
            if not name in counts:
                counts[name] = 0
            counts[name] += 1
        return counts
            

    def pbc(self, r, box, ibox = None):
        """
        Applies periodic boundary conditions.
        Parameters:
            r:      the vector the boundary conditions are applied to
            box:    the box that defines the boundary conditions
            ibox:   the inverse of the box. This will be calcluated if not provided.
        """
        #if ibox == None:    
        #if not hasattr(ibox, 'shape'):
        if type(ibox) != numpy.ndarray and type(ibox) != list and type(ibox) != tuple: #MJW fix
            ibox = numpy.linalg.inv(box)
        vdir = numpy.dot(r, ibox)
        vdir = (vdir % 1.0 + 1.5) % 1.0 - 0.5
        return numpy.dot(vdir, box)


    def per_atom_norm(self, v, box, ibox = None):
        '''
        Returns a length N numpy array containing per atom distance
            v:      an Nx3 numpy array
            box:    box matrix that defines the boundary conditions
            ibox:   the inverse of the box. will be calculated if not provided
        '''
        diff = self.pbc(v, box, ibox)
        return numpy.array([numpy.linalg.norm(d) for d in diff])


    def load_mode(self, modefilein):
        ''' 
        Reads a mode.dat file into an N by 3 numpy array
            modefilein: filename
        '''
        f = open(modefilein, 'r')
        lines = f.readlines()
        f.close()
        mode = []
        for line in lines:
            l = line.strip().split()
            for j in range(3):
                mode.append(float(l[j]))
        mode = numpy.array(mode)
        mode.resize(len(mode)/3, 3)
        return mode


    def save_mode(self, modefileout, displace_vector):
        '''
        Saves an Nx3 numpy array into a mode.dat file. 
            modefileout:     filename
            displace_vector: the mode (Nx3 numpy array)
        '''
        f = open(modefileout, 'w')
        for i in range(len(displace_vector)):
            f.write("%.3f %.3f %.3f\n" % (displace_vector[i][0], 
                displace_vector[i][1], displace_vector[i][2]))

    def list_element_combinations(self, kdbdir):
        combinations = [os.path.basename(i) for i in glob.glob(os.path.join(kdbdir, "*"))]
        return combinations

    def combo_split(self, combo):
        elements = []
        for i in range(len(combo)):
            if combo[i] == combo[i].lower():
                elements[-1] += combo[i]
            else:
                elements.append(combo[i])
        return elements

    def is_symbol_subset(self, a, b):
        for symbol in a:
            if symbol not in b:
                return False
        return True

    def query_has_all(self, kdbdir, symbols):
        result = []
        combinations = [os.path.basename(i) for i in glob.glob(os.path.join(kdbdir, "*"))]
        for combo in combinations:
            elements = self.combo_split(combo)
            if not is_symbol_subset(symbols, elements):
                continue
            for N in glob.glob(os.path.join(kdbdir, combo, '*')):
                result.append(N)
        return result

