from kdbquery import KdbQuery
from optparse import OptionParser
import sys
from aselite import read_any
from remote_db import RemoteDB
from server_config import *

class RemoteQuery(KdbQuery):
    def __init__(self):
        self.return_dict = {}
    # overloads KdbQuery.query_db()
    def query_db(self, **args):
        db = RemoteDB()
        name = db.get_name(args['reactant'].get_chemical_symbols())
        entries = db.get_process(name)
        return entries, name

    def output_query(self, outpurdir, numMatches, suggestion, sugproduct, modeTemp=None):
        self.return_dict[numMatches] = [suggestion, sugproduct]
        if modeTemp is not None:
            self.return_dict[numMatches].append(modeTemp)

if __name__ == "__main__":

    # Parse command line options.
    parser = OptionParser(usage = "%prog [options] reactant.con")
    parser.add_option("--nodupes", dest = "nodupes", action="store_true",
                      help = "detect and remove duplicate suggestions (can be expensive)")
    options, args = parser.parse_args()

    # Make sure we get the reactant file name.
    if len(args) < 1:
        parser.print_help()
        sys.exit()
        
    # Load the reactant con file.
    reactant = read_any(args[0])
    db = RemoteDB()
    params = db.get_params()
    query_sub_class = RemoteQuery()
    query_sub_class.query(reactant, "./kdbmatches_remote", dc=params['dc'], nf=params['nf'], nodupes=options.nodupes, kdbname=db_name)