from kdbinsert import KdbInsert
from optparse import OptionParser
import sys
from aselite import read_any
from config import *
from local_db import LocalDB


class LocalInsert(KdbInsert):
    def __init__(self):
        pass

    # This function will overload the default insert_into_db function
    def insert_into_db(self, **args):
        # create instance of database
        db = LocalDB(args['kdbname'], args['nf'], args['dc'], args['mac'])
        # test if process is already in database
        name = db.get_name(args['s'].get_chemical_symbols())
        saddle_list = db.get_saddles(name)
        for db_saddle in saddle_list:
            if len(args['s']) != len(db_saddle[0]):
                continue
            if self.getMappings(args['s'], db_saddle[0], args['nf'], args['dc']) is not None:
                print "SQL duplicate of", name, "with id:", db_saddle[1]
                return "SQL duplicate of " + name + " with id: " + str(db_saddle[1])
        # add process to db
        db.add_process(args['or'], args['os'], args['op'], args['om'],
                       args['r'], args['s'], args['p'], args['m'], args['ma'])
        # Indicate that the process was inserted successfully.
        #print "good"
        print "KDB insert success"
        return "good"


if __name__ == "__main__":
    insert_sub_class = LocalInsert()

    # Parse command line options.
    parser = OptionParser(usage="%prog [options] reactant saddle product mode")
    parser.add_option("-o", "--mode", dest="mode",
                      help="optional mode file",
                      default=None)
    parser.add_option("-n", "--nf", dest="nf", action="store", type="float",
                      help="neighbor fudge parameter",
                      default=None)
    parser.add_option("-c", "--dc", dest="dc", action="store", type="float",
                      help="distance cutoff parameter",
                      default=None)
    parser.add_option("-m", "--mac", dest="mac", action="store", type="float",
                      help="mobile atom cutoff parameter",
                      default=None)
    options, args = parser.parse_args()

    # Make sure we get the reactant, saddle, product, and mode files.
    if len(args) < 3:
        parser.print_help()
        sys.exit()

    # Load the reactant, saddle, product, and mode files.
    reactant = read_any(args[0])
    saddle = read_any(args[1])
    product = read_any(args[2])
    mode = None
    if options.mode is not None:
        mode = insert_sub_class.load_mode(options.mode)

    # load previous params
    db = LocalDB(KDB_NAME)
    params = db.get_params()
    if options.nf is None:
        options.nf = params['nf']
    if options.dc is None:
        options.dc = params['dc']
    if options.mac is None:
        options.mac = params['mac']

    # run the insert standard insert function.
    insert_sub_class.insert(reactant, saddle, product, mode=mode, nf=options.nf, dc=options.dc, mac=options.mac,
                            kdbname=KDB_NAME)
