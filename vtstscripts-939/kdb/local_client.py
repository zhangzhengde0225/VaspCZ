# Main.py is the user entry point into the program.
import sys

from kdb import Kdb
import local_insert
import local_query
import local_db
from aselite import read_any
from config import *

def run(args):
    if not Kdb().check_version():
        sys.exit()

    if len(args) < 1:
        print "first parameter should be either: insert, query"
        sys.exit()
    if args[0] == 'insert':
        if len(args) < 4:
            print "parameters for insert should include reactant, saddle, and product files."
            sys.exit()
        #read files
        try:
            reactant = read_any(args[1])
            saddle   = read_any(args[2])
            product  = read_any(args[3])
        except IOError:
            print "One or more files could not be read."
            sys.exit()
        try:
            mode = Kdb().load_mode(args[3])
        except:
            mode = None
        #grab params
        db = local_db.LocalDB(KDB_NAME)
    	params = db.get_params()    
    	#insert
        local_insert.LocalInsert().insert(reactant, saddle, product, mode=mode, dc=params['dc'], nf=params['nf'], mac=params['mac'], kdbname=KDB_NAME)
    
    elif args[0] == 'query':
        if len(args) < 2:
            print "parameters for query should include a reactant file."
            sys.exit()
        #read file
        try:
        	reactant = read_any(args[1])
        except IOError:
        	print "reactant file could not be read."
        	sys.exit()
        #grab params
        db = local_db.LocalDB(KDB_NAME)
        params = db.get_params()
        #query
        local_query.LocalQuery().query(reactant, "./kdbmatches", dc=params['dc'], nf=params['nf'], kdbname=KDB_NAME)

if __name__ == "__main__":
    args = sys.argv[1:]
    run(args)