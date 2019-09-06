from kdbinsert import KdbInsert
from optparse import OptionParser
import sys
from aselite import read_any
from remote_db import RemoteDB
from server_config import *

class RemoteInsert(KdbInsert):
    #email and password are set prior to running insert if email/pw combo is present
    def __init__(self):
        self.email = None
        self.password = None
    # overloads KdbInsert.insert_into_db
    def insert_into_db(self, **args):
        #create database instance
        db = RemoteDB()
        # test if process is already in database
        name = db.get_name(args['s'].get_chemical_symbols())
        saddle_list = db.get_saddles(name)
        for db_saddle in saddle_list:
            if len(args['s']) != len(db_saddle[0]):
                continue
            if self.getMappings(args['s'], db_saddle[0], args['nf'], args['dc']) is not None:
                #print "SQL duplicate of", name, "with id:", db_saddle[1]
                return "SQL duplicate of " +  name +  " with id: " +  str(db_saddle[1])
        # add process to db
        try:
            db.add_process(args['or'], args['os'], args['op'], args['om'], 
                           args['r'], args['s'], args['p'], args['m'], args['ma'], self.email, self.password)
        except TypeError:
            print "Account info in user_config.py is not valid. Try running remote_config.py to set up account"
            return
        # Indicate that the process was inserted successfully.
        #print "good"
        return "good"

if __name__ == "__main__":
    insert_sub_class = RemoteInsert()
    # Parse command line options.
    parser = OptionParser(usage = "%prog [options] reactant saddle product mode")
    parser.add_option("-o", "--mode", dest = "mode", 
                      help = "optional mode file",
                      default = None)
    options, args = parser.parse_args()

    # Make sure we get the reactant, saddle, product, and mode files.
    if len(args) < 3:
        parser.print_help()
        sys.exit()

    # Load the reactant, saddle, product, and mode files.
    reactant = read_any(args[0])
    saddle   = read_any(args[1])
    product  = read_any(args[2])
    mode = None
    if options.mode is not None:
        mode = insert_sub_class.load_mode(options.mode)

    # load previous params
    db = RemoteDB()
    params = db.get_params()
    nf = params['nf']
    dc = params['dc']
    mac = params['mac']

    insert_sub_class.insert(reactant, saddle, product, mode, nf=nf, dc=dc, mac=mac)
