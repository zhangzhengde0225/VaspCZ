from remote_db import RemoteDB
from kdbinsert import KdbInsert
from optparse import OptionParser
import pymysql
from server_config import *

class RemoteUpdate(KdbInsert):
    def __init__(self):
        self.current_o_pro_pk = None

    # overloaded insert function from KdbInsert class
    def insert_into_db(self, **args):
        # create instance of database
        db = RemoteDB()
        # add process to db
        db.add_updated_process(args['r'], args['s'], args['p'], args['m'], args['ma'], self.current_o_pro_pk)
        # Indicate that the process was inserted successfully.
        print "good update"

    # changes the nf, dc, and mac values stored in the database
    def change_params(self, nf, dc, mac):
        conn = self.connect_db(db=backup_db_name)
        conn.execute('''DELETE FROM Param''')
        values = [('nf', nf), ('dc', dc), ('mac', mac)]
        for value in values:
            conn.execute('''INSERT INTO Param Values ('%s','%f')''' % value)
        conn.execute('''COMMIT''')
        conn.close()
        print "Parameters have been changed."

    # remove all processes that are not backups
    def remove_all_process(self):
        conn = self.connect_db()
        conn.execute('''DELETE FROM Atom''')
        conn.execute('''DELETE FROM Mobile''')
        conn.execute('''DELETE FROM Process''')
        conn.execute('''DELETE FROM Atoms''')
        conn.execute('''COMMIT''')
        conn.close()

    # query the backup database for all processes
    # and reinsert processes into queryable database
    def populate_kdb(self):
        conn = self.connect_db(db=backup_db_name)
        conn.execute('''SELECT pro_id, reactant_id, saddle_id, product_id FROM Process''')
        process_list = conn.fetchall()
        conn.close()
        db = RemoteDB()
        params = db.get_params()
        if len(process_list) == 0:
            print "No items in database to update."
            return
        print "Updating", db_name, "database."
        for process in process_list:
            reactant = db.get_atoms(process[1], db=backup_db_name)
            saddle   = db.get_atoms(process[2], db=backup_db_name)
            product  = db.get_atoms(process[3], db=backup_db_name)
            mode     = db.get_mode(process[1], db=backup_db_name)
            self.current_o_pro_pk = process[0]
            self.insert(reactant, saddle, product, mode=mode, nf=params['nf'], dc=params['dc'], mac=params['mac'])


    def connect_db(self, db=db_name):
        return pymysql.connect(host=host, port=port, user=user, passwd=password, db=db).cursor()

if __name__ == "__main__":
    db = RemoteDB()
    params = db.get_params()
    # Parse command line options.
    parser = OptionParser(usage = "%prog [options] reactant saddle product mode")
    parser.add_option("-n", "--nf", dest = "nf", action="store", type="float", 
                      help = "neighbor fudge parameter",
                      default = params['nf'])
    parser.add_option("-c", "--dc", dest = "dc", action="store", type="float", 
                      help = "distance cutoff parameter",
                      default = params['dc'])
    parser.add_option("-m", "--mac", dest = "mac", action="store", type="float", 
                      help = "mobile atom cutoff parameter",
                      default = params['mac'])
    options, args = parser.parse_args()

    update_class = RemoteUpdate()
    update_class.change_params(options.nf, options.dc, options.mac)
    update_class.remove_all_process()
    update_class.populate_kdb()
    print "Update complete."