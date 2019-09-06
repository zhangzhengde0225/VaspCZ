#!/usr/bin/env python

import local_db
from kdbinsert import KdbInsert

class LocalUpdate(KdbInsert):
	# overrides default insert_into_db function
    def insert_into_db(self, **args):
        # create instance of database
        db = local_db.LocalDB(args['kdbname'], args['nf'], args['dc'], args['mac'])
        # add process to db
        db.add_process(args['or'], args['os'], args['op'], args['om'], 
                         args['r'], args['s'], args['p'], args['m'], args['ma'])
        # Indicate that the process was inserted successfully.
        print "good update"
