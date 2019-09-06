from kdbquery import KdbQuery
from optparse import OptionParser
import sys
from aselite import read_any
from ase.io import write
from config import *
from local_db import LocalDB

class LocalQuery(KdbQuery):
    # this function will override the default query_db function
    def query_db(self, **args):
        # initialize db class
        db = LocalDB(args['kdbname'])
        # get name of reactant and all entries with the same name
        name = db.get_name(args['reactant'].get_chemical_symbols())
        entries = db.get_process(name)
        return entries, name

if __name__ == "__main__":

    # Parse command line options.
    parser = OptionParser(usage = "%prog [options] reactant.con")
    parser.add_option("-a", "--alias", dest = "alias", 
                      help = "the name of the SQL kdb",
                      default = KDB_NAME)
    parser.add_option("--nodupes", dest = "nodupes", action="store_true",
                      help = "detect and remove duplicate suggestions (can be expensive)")
    options, args = parser.parse_args()

    # Make sure we get the reactant file name.
    if len(args) < 1:
        parser.print_help()
        sys.exit()
        
    # Load the reactant con file.
    reactant = read_any(args[0])

    # get database params
    db = LocalDB(options.alias)
    params = db.get_params()

    #create isntance of LocalQuery and run query()
    query_sub_class = LocalQuery()
    query_sub_class.query(reactant, "./kdbmatches", dc=params['dc'], nf=params['nf'], nodupes=options.nodupes, kdbname=options.alias)
    