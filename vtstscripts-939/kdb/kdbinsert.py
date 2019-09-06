#!/usr/bin/env python

import numpy

from aselite import elements
from aselite import FixAtoms
from kdb import Kdb


class KdbInsert(Kdb):
    def __init__(self):
        pass

    def coordination_numbers(self, p, nf):
        nl = []
        for a in range(len(p)):
            nl.append([])
            for b in range(len(p)):
                if b != a:
                    dist = numpy.linalg.norm(p.get_positions()[a] - p.get_positions()[b])        
                    if dist < (elements[p.get_chemical_symbols()[a]]["radius"] + 
                               elements[p.get_chemical_symbols()[b]]["radius"]) * (1.0 + nf):
                        nl[a].append(b)
        return [len(l) for l in nl]

    def getMappings(self, a, b, nf, dc, mappings = None):
        """ A recursive depth-first search for a complete set of mappings from atoms
            in configuration a to atoms in configuration b. Do not use the mappings
            argument, this is only used internally for recursion.
            
            Returns None if no mapping was found, or a dictionary mapping atom 
            indices a to atom indices b.
            
            Note: If a and b are mirror images, this function will still return a 
            mapping from a to b, even though it may not be possible to align them 
            through translation and rotation. """
        # If this is the top-level user call, create and loop through top-level
        # mappings.
        if mappings == None:
            # Find the least common coordination number in b.
            bCoordinations = self.coordination_numbers(b, nf)
            bCoordinationsCounts = {}
            for coordination in bCoordinations:
                if coordination in bCoordinationsCounts:
                    bCoordinationsCounts[coordination] += 1
                else:
                    bCoordinationsCounts[coordination] = 1
            bLeastCommonCoordination = bCoordinationsCounts.keys()[0]
            for coordination in bCoordinationsCounts.keys():
                if bCoordinationsCounts[coordination] < bCoordinationsCounts[bLeastCommonCoordination]:
                    bLeastCommonCoordination = coordination
            # Find one atom in a with the least common coordination number in b. 
            # If it does not exist, return None.
            aCoordinations = self.coordination_numbers(a, nf)
            try:
                aAtom = aCoordinations.index(bLeastCommonCoordination)
            except ValueError:
                return None
            # Create a mapping from the atom chosen from a to each of the atoms with
            # the least common coordination number in b, and recurse.
            for i in range(len(bCoordinations)):
                if bCoordinations[i] == bLeastCommonCoordination:
                    # Make sure the element types are the same.
                    if a.get_chemical_symbols()[aAtom] != b.get_chemical_symbols()[i]:
                        continue
                    mappings = self.getMappings(a, b, nf, dc, {aAtom:i})
                    # If the result is not none, then we found a successful mapping.
                    if mappings is not None:
                        return mappings
            # There were no mappings.        
            return None
        
        # This is a recursed invocation of this function.
        else:
            # Find an atom from a that has not yet been mapped.
            unmappedA = 0
            while unmappedA < len(a):
                if unmappedA not in mappings.keys():
                    break
                unmappedA += 1
            # Calculate the distances from unmappedA to all mapped a atoms.
            distances = {}
            for i in mappings.keys():
                distances[i] = self.atomAtomDistance(a, unmappedA, i)
            
            # Loop over each unmapped b atom. Compare the distances between it and 
            # the mapped b atoms to the corresponding distances between unmappedA 
            # and the mapped atoms. If everything is similar, create a new mapping
            # and recurse.
            for bAtom in range(len(b)):
                if bAtom not in mappings.values():
                    for aAtom in distances:
                        # Break if type check fails.
                        if b.get_chemical_symbols()[bAtom] != a.get_chemical_symbols()[unmappedA]:
                            break
                        # Break if distance check fails  
                        bDist = self.atomAtomDistance(b, bAtom, mappings[aAtom])
                        if abs(distances[aAtom] - bDist) > dc:
                            break
                    else:
                        # All distances were good, so create a new mapping.
                        newMappings = mappings.copy()
                        newMappings[unmappedA] = bAtom
                        # If this is now a complete mapping from a to b, return it.
                        if len(newMappings) == len(a):
                            return newMappings
                        # Otherwise, recurse.
                        newMappings = self.getMappings(a, b, nf, dc, newMappings)
                        # Pass any successful mapping up the recursion chain. 
                        if newMappings is not None:
                            return newMappings     
            # There were no mappings.   
            return None 


    def stripUnselectedAtoms(self, atoms, selected):
        """ Removes any atoms from atoms that are not in selected and returns a new
        structure and a mapping from atoms in the old structure to atoms in the new 
        structure. """
        src = atoms.copy()
        dest = atoms.copy()
        while len(dest) > 0:
            dest.pop()
        mapping = {}
        index = 0
        constraints = []
        for i in selected:
            mapping[i] = index
            index += 1
            if i in src.constraints[0].index: 
                constraints.append(index)
            dest.append(src[i])
        dest.set_constraint(FixAtoms(constraints))
        return dest, mapping
        
        
    def getProcessMobileAtoms(self, r, s, p, mac):
        """ Returns a list of atom indices that move more than mac 
        between reactant and saddle, saddle and product, or 
        reactant and product. If no atoms move more than mac, returns
        the atom that moves the most. """
        mobileAtoms = []
        reactant2saddle = self.per_atom_norm(s.positions - r.positions, s.get_cell())
        product2saddle = self.per_atom_norm(s.positions - p.positions, s.get_cell())
        reactant2product = self.per_atom_norm(p.positions - r.positions, s.get_cell())
        for i in range(len(s)):
            if max(reactant2saddle[i], product2saddle[i], reactant2product[i]) > mac:
                mobileAtoms.append(i)
        if len(mobileAtoms) == 0:
            mobileAtoms.append(list(reactant2product).index(max(reactant2product)))
        return mobileAtoms


    def getProcessNeighbors(self, mobileAtoms, r, s, p, nf):
        """ Given a list mobile atoms, a reactant, saddle, and product, 
        returns a list of neighboring atoms according to the nf (NEIGHBOR_FUDGE)
        paramter."""
        neighborAtoms = []
        for atom in mobileAtoms:
            r1 = elements[s.get_chemical_symbols()[atom]]["radius"]
            for i in range(len(s)):
                if i in mobileAtoms or i in neighborAtoms:
                    continue
                r2 = elements[s.get_chemical_symbols()[i]]["radius"]
                maxDist = (r1 + r2) * (1.0 + nf)
                if self.atomAtomPbcDistance(r, atom, i) < maxDist:
                    neighborAtoms.append(i)
                elif self.atomAtomPbcDistance(s, atom, i) < maxDist:
                    neighborAtoms.append(i)
                elif self.atomAtomPbcDistance(p, atom, i) < maxDist:
                    neighborAtoms.append(i)
        return neighborAtoms

    #function will be overridden in remote/local classes
    def insert_into_db(self, **args):
        print "function not yet overloaded"



    def insert(self, reactant, saddle, product, mode=None, nf=0.2, dc=0.3, mac=0.7, kdbname='kdb.db'):

        # Keep a copy of the original data
        original_reactant = reactant.copy()
        original_saddle = saddle.copy()
        original_product = product.copy()

        if mode is not None:
            original_mode = mode.copy()
        else:
    	    original_mode = None
        mobileAtoms = self.getProcessMobileAtoms(reactant, saddle, product, mac)
            
        selectedAtoms = mobileAtoms + self.getProcessNeighbors(mobileAtoms, reactant, product, saddle, nf)
        
        # Quit if not enough selected atoms.
        if len(selectedAtoms) < 2:
            print "kdbinsert abort: Too few atoms in process, or neighbor_fudge too small."
            return 0
                
        # Remove unselected atoms.
        reactant, mapping = self.stripUnselectedAtoms(reactant, selectedAtoms)
        saddle,   mapping = self.stripUnselectedAtoms(saddle,   selectedAtoms)
        product,  mapping = self.stripUnselectedAtoms(product,  selectedAtoms)

        # Update the mode.
        if mode is not None:
            newMode = numpy.zeros((len(selectedAtoms), 3))
            for m in mapping:
                newMode[mapping[m]] = mode[m]
            mode = newMode

        # Remove PBC's.
        temp = reactant.copy()
        undone = range(len(temp))
        working = [undone.pop()]        
        while len(undone) > 0:
            if len(working) == 0:
                print "kdbinsert abort: Dissociated reactant, or neighbor_fudge too small."
                return 0
            a = working.pop()
            for i in undone[:]:
                v = self.pbc(temp.positions[i] - temp.positions[a], temp.get_cell())
                d = numpy.linalg.norm(v)
                if d < (elements[temp.get_chemical_symbols()[a]]["radius"] + 
                        elements[temp.get_chemical_symbols()[i]]["radius"]) * (1.0 + nf):
                    temp[i].position = temp[a].position + v
                    working.append(i)
                    undone.remove(i)
        v1s = self.pbc(saddle.positions - reactant.positions, reactant.get_cell())
        v12 = self.pbc(product.positions - reactant.positions, reactant.get_cell())
        reactant = temp
        saddle.positions = reactant.positions + v1s
        product.positions = reactant.positions + v12
        
        # Find saddle center of coordinates.
        coc = numpy.zeros((1,3))
        for i in range(len(saddle)):
            coc += saddle[i].position
        coc = coc / len(saddle)
        
        # Shift all structures so that the saddle center of coordinates is at 
        # [0, 0, 0].
        reactant.positions = reactant.positions - coc    
        saddle.positions = saddle.positions - coc    
        product.positions = product.positions - coc    

        # Give all structures a huge box.
        # TODO: all references to boxes should be removed after PBCs are removed.
        reactant.cell = numpy.identity(3) * 1024
        saddle.cell = numpy.identity(3) * 1024
        product.cell = numpy.identity(3) * 1024

        # get mobile_list
        mob_list = []
        for atom in mobileAtoms:
            mob_list.append(mapping[atom])

        arg_dict = {'or': original_reactant, 'os': original_saddle, 'op': original_product, 'om': original_mode,
                    'r': reactant, 's': saddle, 'p': product, 'm': mode, 'ma': mob_list,
                    'kdbname': kdbname, 'nf': nf, 'dc': dc, 'mac': mac}
        # function is overloaded in either local_insert.py or remote_insert.py
        return self.insert_into_db(**arg_dict)
