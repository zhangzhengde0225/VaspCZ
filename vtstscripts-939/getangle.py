#!/usr/bin/env python

import time
import sys
import numpy as np
from string import Template
import fnmatch
import os
import re
import datetime

## READ IN VASP CON FILE ##
filein = sys.argv[1]

with open(filein, 'r') as f:
    first_line = f.readline()
    second_line = f.readline()
    third_line = f.readline()
    fourth_line = f.readline()
    fifth_line = f.readline()

foo = first_line[0]
elementSymbol = first_line.split()
a = second_line.split()
foo_v1 = third_line.split()
foo_v2 = fourth_line.split()
foo_v3 = fifth_line.split()

#Assign lattice vectors to float arrays
v1 = np.array([float(foo_v1[0]), float(foo_v1[1]), float(foo_v1[2])])
v2 = np.array([float(foo_v2[0]), float(foo_v2[1]), float(foo_v2[2])])
v3 = np.array([float(foo_v3[0]), float(foo_v3[1]), float(foo_v3[2])])

#Calcualte cross product of v1, v2 and save as surface area 
vol12 = np.cross(v2, v1)
vol23 = np.cross(v2, v3)
vol13 = np.cross(v3, v1)

SA12 = np.linalg.norm(vol12)
SA23 = np.linalg.norm(vol23)
SA13 = np.linalg.norm(vol13)

#volume of cell
vol = np.dot(vol12, v3)*float(a[0])*float(a[0])*float(a[0])

cosang12 = np.dot(v1, v2)/np.linalg.norm(v1)/np.linalg.norm(v2)
cosang13 = np.dot(v1, v3)/np.linalg.norm(v1)/np.linalg.norm(v3)
cosang23 = np.dot(v2, v3)/np.linalg.norm(v2)/np.linalg.norm(v3)
rad12 = np.arccos(cosang12)
rad13 = np.arccos(cosang13)
rad23 = np.arccos(cosang23)
ang12 = np.degrees(rad12)
ang13 = np.degrees(rad13)
ang23 = np.degrees(rad23)

print "Angle between v1 and v2: %f" % ang12
print "Angle between v1 and v3: %f" % ang13
print "Angle between v2 and v3: %f" % ang23
print "Surface Area of side v1, v2: %f" % SA12
print "Surface Area of side v1, v3: %f" % SA13
print "Surface Area of side v2, v3: %f" % SA23
print "Volume: %f" % vol
