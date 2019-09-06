#!/usr/bin/env python

import sys
import glob

ignore = ['kdbaddpr.pl', 'kdbaddprnew.pl', 'kdbaddvpr.pl', 'kdbquery.pl', 'kdbquerynew.pl', 'akmcreset.pl', 'akmcupdate.pl', 'akmc.pl', 'akmccleanjobs.pl', 'pos2con.pl']

if len(sys.argv) > 1:
    for fname in sys.argv[1:]:
        if fname in ignore:
            ignore.remove(fname)
        else:
            ignore.append(fname)
    f = open(__file__, 'r')
    lines = f.readlines()
    f.close()
    f = open(__file__, 'w')
    for line in lines:
        if line.startswith('ignore ='):
            line = 'ignore = ' + repr(ignore) + "\n"
        f.write(line)
    f.close()
    print 'The coverage.py ignore list has been updated:', ignore
    sys.exit()

ig = [p.replace('.pl', '') for p in ignore]
pl = sorted([p.replace('.pl', '') for p in glob.glob("*.pl")])
py = [p.replace('.py', '') for p in glob.glob("*.py")]

print

for p in pl:
    if p in ig:
        print "\033[02;32m% 16s.pl (not going to implement)\033[00m" % p
    elif p in py:
        print "\033[01;32m% 16s.pl %s.py\033[00m" % (p, p)
    else:
        print "\033[01;31m% 16s.pl\033[00m" % p

print    

