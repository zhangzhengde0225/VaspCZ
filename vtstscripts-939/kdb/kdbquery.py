#!/usr/bin/env python


import os
import sys
import numpy 
import glob
import shutil
import math
import copy

from optparse import OptionParser

from kdb import Kdb
from config import *
from aselite import elements
from aselite import write_vasp

class KdbQuery(Kdb):
    def __init__(self):
        self.return_dict = {}

    def isDistance(self, pbcvector, target, box, dc):
        for x in [-1, 0, 1]:
            for y in [-1, 0, 1]:
                for z in [-1, 0, 1]:
                    temp = pbcvector.copy()
                    temp += x * box[0]
                    temp += y * box[1]
                    temp += z * box[2]
                    if abs(numpy.linalg.norm(temp) - target) < dc:
                        return True
        return False
        

    def centroid(self, a, which=None):
        if which == None:
            which = range(len(a))
        c = numpy.array([0.0, 0.0, 0.0])
        for i in which:
            c += a.positions[i]
        c /= len(which)
        return c


    def clump(self, c, atoms, nf):
        # Remove PBC's.
        temp = c.copy()
        undone = atoms[:]
        working = [undone.pop()]
        while len(undone) > 0:
            if len(working) == 0:
                print "Dissociated reactant, or neighbor_fudge too small."
                return
            a = working.pop()
            for i in undone[:]:
                v = self.pbc(temp.positions[i] - temp.positions[a], temp.cell)
                d = numpy.linalg.norm(v)
                if d < (elements[temp.get_chemical_symbols()[a]]["radius"] + elements[temp.get_chemical_symbols()[i]]["radius"]) * (1.0 + nf):
                    temp.positions[i] = temp.positions[a] + v
                    working.append(i)
                    undone.remove(i)
        return temp

    def query_db(self, **args):
        print "function not yet overloaded"

    #note this function gets overloaded when interacting with remote DB
    def output_query(self, outputdir, numMatches, suggestion, sugproduct, modeTemp=None):
        #create direcotry if none exists.
        if not os.path.isdir(outputdir):
            os.mkdir(outputdir)
        #write the files to the output
        write_vasp(outputdir + "/SADDLE_%d" % numMatches, suggestion)
        write_vasp(outputdir + "/PRODUCT_%d" % numMatches, sugproduct)
        if modeTemp is not None:
            self.save_mode(outputdir + "/MODE_%d" % numMatches, modeTemp)
        os.system("touch %s/.done_%d" % (outputdir, numMatches))

    def query(self, reactant, outputdir = "./kdbmatches", nf=0.2, dc=0.3, nodupes = False, kdbname = 'kdb.db'):
        # XXX: I think the best way forward to allow parallel processes
        # here is to make the query function return atoms objects instead
        # of writing them to file there.

        # Get the ibox to speed up pbcs.
        ibox = numpy.linalg.inv(reactant.cell)
        
        # Remove directory if kdbmatches is already there.
        if os.path.isdir(outputdir):
            shutil.rmtree(outputdir)
        # A list of unique saddles, used for duplicate removal.
        uniques = []
        
        # Get a list of kdb entries that match the query configuration elementally.
        entries, name = self.query_db(kdbname = kdbname, reactant = reactant)

        if len(entries) == 0:
            print "No entries for those elements."
            return

        # For each nonfrozen atom in reactant, create a list of neighboring element
        # types and the count of each type.
        # TODO: this can be made N^2/2 trivially.
        # TODO: this can use SAP for ortho boxes.
        reactantNeighbors = {}
        for i in range(len(reactant)):
            if i in reactant.constraints[0].index:
                continue
            r1 = elements[reactant.get_chemical_symbols()[i]]["radius"]
            reactantNeighbors[i] = {}
            for j in range(len(reactant)):
                if j == i:
                    continue
                r2 = elements[reactant.get_chemical_symbols()[j]]["radius"]
                d = numpy.linalg.norm(self.pbc(reactant.positions[i] - reactant.positions[j], reactant.cell, ibox))
                if d > (r1 + r2) * (1 + nf):
                    continue
                if reactant.get_chemical_symbols()[j] not in reactantNeighbors[i]:
                    reactantNeighbors[i][reactant.get_chemical_symbols()[j]] = 0
                reactantNeighbors[i][reactant.get_chemical_symbols()[j]] += 1
        
        # Create a list of element types and counts for the entire reactant. 
        reactantNameCount = self.nameCount(reactant)
        numMatches = 0
        
        ###########################################################################
        # (Main) Loop over each kdb entry.
        ###########################################################################
        for entry in entries:
        
            entryMatches = 0
        
            mirrored = "not mirrored"
            if entry["mirror"]:
                mirrored = "mirrored"
            # print "checking", name, "with id:", entry['id'], mirrored
            print "KDB checking entry:", entry['id'], "(",mirrored,")"

            # Load the minimum.
            kdbmin = copy.deepcopy(entry['minimum'])        

            # Make sure the reactant has at least as many atoms of each type as the
            # kdb configuration.
            passedNameCount = True
            kdbNameCount = self.nameCount(kdbmin)
            for name in kdbNameCount:
                if name not in reactantNameCount:
                    passedNameCount = False
                    break
                if kdbNameCount[name] > reactantNameCount[name]:
                    passedNameCount = False
                    break
            if not passedNameCount:
                print "%10d  name count fail" % entryMatches
                continue

            # Load the mobile atoms list.
            kdbmobile = copy.deepcopy(entry['mobile'])     

            # Mirror the minimum if the mirror flag is set for this entry.
            if entry["mirror"]:
                for i in range(len(kdbmin)):
                    kdbmin.positions[i] += 2.0 * (kdbmin.positions[0] - kdbmin.positions[i])
            
            # For each mobile atom in kdbmin, create a list of neighboring element
            # types and the count of each type.
            kdbNeighbors = {}
            for i in kdbmobile:
                r1 = elements[kdbmin.get_chemical_symbols()[i]]["radius"]
                kdbNeighbors[i] = {}
                for j in range(len(kdbmin)):
                    if j == i:
                        continue
                    r2 = elements[kdbmin.get_chemical_symbols()[j]]["radius"]
                    d = numpy.linalg.norm(kdbmin.positions[i] - kdbmin.positions[j])
                    if d > (r1 + r2) * (1 + nf):
                        continue
                    if kdbmin.get_chemical_symbols()[j] not in kdbNeighbors[i]:
                        kdbNeighbors[i][kdbmin.get_chemical_symbols()[j]] = 0
                    kdbNeighbors[i][kdbmin.get_chemical_symbols()[j]] += 1

            kdbUnmapped = range(len(kdbmin)) # Keep track of the kdb atoms that have been mapped.

            # Create the initial mappings.
            mappings = None
            db_a = kdbmobile[0] # This will be the selected mobile atom.
            for m in kdbmobile:
                mMappings = []
                for freeAtom in reactantNeighbors.keys():
                    for elementType in reactantNeighbors[freeAtom]:
                        if elementType not in kdbNeighbors[m]:
                            break
                        if kdbNeighbors[m][elementType] != reactantNeighbors[freeAtom][elementType]:
                            break
                    else:
                        mMappings.append({m:freeAtom})
                if mappings == None:
                    mappings = mMappings
                if len(mMappings) < len(mappings):
                    mappings = mMappings
                    db_a = m
            
            kdbUnmapped.remove(db_a)

            while len(kdbUnmapped) > 0 and len(mappings) > 0:
                # Create a list of new mappings that will replace mappings at the
                # end of this iteration.
                newMappings = []
                # Select an unmapped atom from kdbmin.
                kdbAtom = kdbUnmapped.pop()
                # Get the distance between kdbAtom and every other atom in the kdb
                # configuration.
                kdbDistances = {}
                for i in range(len(kdbmin)):
                    kdbDistances[i] = numpy.linalg.norm(kdbmin.positions[kdbAtom] - kdbmin.positions[i])
                # Loop over each mapping and try to place kdbAtom.
                for mapping in mappings:
                    # Loop over each atom in the reactant.
                    for reactantAtom in range(len(reactant)):
                        # Make sure it has not already been mapped.
                        if reactantAtom in mapping.values():
                            continue
                        # Loop over the atoms in mapping and see if the distance
                        # between reactantAtom and mapping.values() atoms is the same
                        # within dc (DISTANCE_CUTOFF) of the distance between kdbAtom
                        # and mapping.keys() atoms.
                        for DA in mapping.keys():
                            RA = mapping[DA]
                            pbcVector = self.atomAtomPbcVector(reactant, RA, reactantAtom)
                            if PBC_MAPPING_CHECK:
                                if not self.isDistance(pbcVector, kdbDistances[DA], reactant.cell, dc):
                                    break
                            else:
                                if abs(kdbDistances[DA] - self.atomAtomPbcDistance(reactant, RA, reactantAtom)) > dc:
                                    break
                        else:
                            newMapping = mapping.copy()
                            newMapping[kdbAtom] = reactantAtom
                            newMappings.append(newMapping)
                mappings = newMappings

            # Load the mode.
            mode = copy.deepcopy(entry['mode'])

            # Loop over each mapping and try to find a rotation that aligns the
            # kdb configuration with the query configuration.
            for mapping in mappings:
                
                #print "in mappings loop"
                #print "mapping values: ",mapping.values()

                reactantrot = self.clump(reactant, mapping.values(), nf)

                # If no neighbors found, go to next mapping
                if reactantrot is None:
                    continue
            
                #print "mapping values, after clump: ",mapping.values()
                #print "reactantrot: ",reactantrot

                # Make a copy of kdbmin for rotation and put it in the box.
                kdbrot = kdbmin.copy()
                kdbrot.cell = reactant.cell.copy()
                
                # Rotation Matrix calculation start
                tb = kdbrot.copy()
                tb.positions -= self.centroid(tb)
                ta = tb.copy()
                offset = self.centroid(reactantrot, mapping.values())
                i = 0
                for m in mapping:
                    ta.positions[i] = tb.positions[m] + self.pbc((reactantrot.positions[mapping[m]] - offset) - tb.positions[m], reactantrot.cell)
                    i += 1
                ta.positions -= self.centroid(ta)
                m = numpy.dot(tb.positions.transpose(), ta.positions)
                sxx = m[0][0]
                sxy = m[0][1]
                sxz = m[0][2]
                syx = m[1][0]
                syy = m[1][1]
                syz = m[1][2]
                szx = m[2][0]
                szy = m[2][1]
                szz = m[2][2]
                n = numpy.zeros((4,4))
                n[0][1] = syz - szy
                n[0][2] = szx - sxz
                n[0][3] = sxy - syx
                n[1][2] = sxy + syx
                n[1][3] = szx + sxz
                n[2][3] = syz + szy
                n += n.transpose()
                n[0][0] =  sxx + syy + szz
                n[1][1] =  sxx - syy - szz
                n[2][2] = -sxx + syy - szz
                n[3][3] = -sxx - syy + szz
                w, v = numpy.linalg.eig(n)
                maxw = 0
                maxv = 0
                for i in range(len(w)):
                    if w[i] > maxw:
                        maxw = w[i]
                        maxv = v[:,i]
                Rmat = numpy.zeros((3,3))
                aa = maxv[0]**2
                bb = maxv[1]**2
                cc = maxv[2]**2
                dd = maxv[3]**2
                ab = maxv[0]*maxv[1]
                ac = maxv[0]*maxv[2]
                ad = maxv[0]*maxv[3]
                bc = maxv[1]*maxv[2]
                bd = maxv[1]*maxv[3]
                cd = maxv[2]*maxv[3]
                Rmat[0][0] = aa + bb - cc - dd
                Rmat[0][1] = 2*(bc-ad) 
                Rmat[0][2] = 2*(bd+ac) 
                Rmat[1][0] = 2*(bc+ad) 
                Rmat[1][1] = aa - bb + cc - dd
                Rmat[1][2] = 2*(cd-ab) 
                Rmat[2][0] = 2*(bd-ac) 
                Rmat[2][1] = 2*(cd+ab) 
                Rmat[2][2] = aa - bb - cc + dd
                Rmat = Rmat.transpose()
                # Rotation Matrix calculation end

                translation1 = self.centroid(kdbrot)
                kdbrot.positions -= translation1
                kdbrot.positions = numpy.dot(kdbrot.positions, Rmat)
                
                translation2 = self.centroid(reactantrot, mapping.values())
                
                kdbrot.positions += translation2

                # Calculate a score for this mapping.
                score = max([numpy.linalg.norm(self.pbc(kdbrot.positions[m] - reactantrot.positions[mapping[m]], reactantrot.cell)) for m in mapping])
                
                if score > dc:
                    continue
                
                # Load the saddle from the database.
                kdbSaddle = copy.deepcopy(entry['saddle'])
                
                # Mirror the saddle if the mirror flag is set for this entry.
                if entry["mirror"]:
                    for i in range(len(kdbSaddle)):
                        kdbSaddle.positions[i] += 2.0 * (kdbmin.positions[0] - kdbSaddle.positions[i])

                # Load the product from the database.
                kdbProduct = copy.deepcopy(entry['product'])           

                # Mirror the product if the mirror flag is set for this entry.
                if entry["mirror"]:
                    for i in range(len(kdbProduct)):
                        kdbProduct.positions[i] += 2.0 * (kdbmin.positions[0] - kdbProduct.positions[i])

                # Map the mode.
                if mode is not None:
                    modeTemp = reactantrot.positions * 0.0
                    for m in mapping:
                        modeTemp[mapping[m]] = mode[m]
                    try:
                        modeTemp /= numpy.linalg.norm(modeTemp)
                    except FloatingPointError:
                        mode = None

                # Perform the saddle transformation.
                kdbSaddle.positions -= translation1
                kdbSaddle.positions = numpy.dot(kdbSaddle.positions, Rmat)
                kdbSaddle.positions += translation2
                
                # Perform the mode transformation.
                if mode is not None:
                    modeTemp = numpy.dot(modeTemp, Rmat)

                # Perform the product transformation.
                kdbProduct.positions -= translation1
                kdbProduct.positions = numpy.dot(kdbProduct.positions, Rmat)
                kdbProduct.positions += translation2
                
                # Create the suggestion.
                suggestion = reactant.copy()
                sugproduct = reactant.copy()
                for m in mapping:
                    if mapping[m] not in suggestion.constraints[0].index:
                        suggestion.positions[mapping[m]] = kdbSaddle.positions[m]
                    if mapping[m] not in sugproduct.constraints[0].index:
                        sugproduct.positions[mapping[m]] = kdbProduct.positions[m]
                
                # Check for duplicates.
                if nodupes:
                    isdupe = False
                    for unique in uniques:
                        pan = self.per_atom_norm(unique.positions - suggestion.positions, suggestion.cell, ibox)
                        if max(pan) <= dc:
                            isdupe = True
                            break
                    if isdupe:
                        continue
                    uniques.append(suggestion.copy())
                
                # Rebox.
                if REBOX_SUGGESTIONS:
                    suggestion.positions = self.pbc(suggestion.positions, suggestion.cell)
                    sugproduct.positions = self.pbc(sugproduct.positions, sugproduct.cell)
                    
                # Write suggestion.
                if mode is not None:
                    self.output_query(outputdir, numMatches, suggestion, sugproduct, modeTemp)
                else:
                    self.output_query(outputdir, numMatches, suggestion, sugproduct)
                
                entryMatches += 1
                numMatches += 1

            #print "%10d" % entryMatches
            print "KDB matches: ", entryMatches
    
