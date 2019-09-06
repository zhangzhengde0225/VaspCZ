#!/usr/bin/env python

from os import system
from os.path import dirname, abspath, join


vtst_path = dirname(abspath(__file__))

system("grep ut insout.dat | grep -v itr | cut -c 5-100 > o.u.t.t.e.m.p")
system("gnuplot %s" % join(vtst_path, 'insplot.gnu'))
system("rm o.u.t.t.e.m.p")
