
import sys
import os
import shutil

cwd = os.getcwd()
os.chdir(os.path.dirname(__file__))
from kdb import kdb
from kdb import remote_client
# this is not ideal but to pickle the ASE.Atoms objects we need the class files in the same directory
# that the program is being called from.
shutil.copyfile('kdb/aselite.py', 'aselite.py')
import aselite
os.remove('aselite.py')
os.remove('aselite.pyc')
os.chdir(cwd)

def run(args):
    if not kdb.Kdb().check_version():
        sys.exit()

    if len(args) < 1:
        print "first parameter sohuld be either: insert or query"
        sys.exit()
    if args[0] == 'insert':
        if len(args) < 4:
            print "parameters for insert should include reactant, saddle, and product files."
            sys.exit()

        try:
            reactant = aselite.read_any(args[1])
            saddle   = aselite.read_any(args[2])
            product  = aselite.read_any(args[3])
        except IOError:
            print "One or more files could not be read."
            sys.exit()
        try:
            mode = kdb.Kdb().load_mode(args[3])
        except:
            mode = None
        remote_client.server_insert(reactant, saddle, product, mode)

    elif args[0] == 'query':
        if len(args) < 2:
            print "parameters for query should include a reactant file."
            sys.exit()
        try:
            reactant = aselite.read_any(args[1])
        except IOError:
            print "could not read reactant file."
            sys.exit()
        remote_client.server_query(reactant)
    else:
    	print "first parameter sohuld be either: insert or query"


if __name__ == "__main__":
	args = sys.argv[1:]
	run(args)