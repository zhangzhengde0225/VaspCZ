from aselite import Atoms
from aselite import FixAtoms
from aselite import write_vasp
import numpy as np
import pymysql
from server_config import *
import re
import copy
import sys

class RemoteDB():
    def __init__(self):
        self.db_name = db_name
        self.backup_db_name = backup_db_name
        self.user_db_name = user_db_name
        #dictionaries hold highest primary key values to reduce number of queries while inserting
        self.pk_dict = {'Atoms':   self.get_max('Atoms',   'atoms_id',  self.db_name),
                        'Atom':    self.get_max('Atom',    'atom_id',   self.db_name),
                        'Process': self.get_max('Process', 'pro_id',    self.db_name),
                        'Mobile':  self.get_max('Mobile',  'mob_id',    self.db_name)}

        self.b_pk_dict = {'Atoms':   self.get_max('Atoms',   'atoms_id',  self.backup_db_name),
                          'Atom':    self.get_max('Atom',    'atom_id',   self.backup_db_name),
                          'Process': self.get_max('Process', 'pro_id',    self.backup_db_name)}

        self.u_pk_dict = {'User': self.get_max('User', 'user_id', self.user_db_name)}

    #############################################################
    # Insertions ################################################
    #############################################################

    # adds a process, each parameter is an ase.Atom object 
    # EXCEPT mode/o_mode, mode is a numpy array or None 
    # and mobile_set is a list of numbers                                   
    def add_process(self, o_reactant, o_saddle, o_product, o_mode, 
                    reactant, saddle, product, mode, mobile_set, email, passwd):

        #check email/passwd combo is in database
        user_id = self.get_user_id(email, passwd)
        if user_id == None:
            print "email/password incorrect."
            return -1

        # get name and primary keys
        symbols = o_reactant.get_chemical_symbols()
        name = self.get_name(symbols)
        o_pro_pk = self.get_id('Process', db='backup')
        pro_pk = self.get_id('Process')

        # connect to backup DB
        conn = self.connect_db(db=self.backup_db_name)
        conn.execute('''BEGIN''')
        # insert atoms objects into backup DB
        or_pk = self.add_atoms(o_reactant, o_mode, conn)
        os_pk = self.add_atoms(o_saddle, o_mode, conn)
        op_pk = self.add_atoms(o_product, o_mode, conn)
        # insert process to backup DB
        conn.execute('''INSERT INTO Process 
                        VALUES ('%d','%s','%d','%d','%d', '%d')''' % 
                     (o_pro_pk, name, or_pk, os_pk, op_pk, user_id))
        conn.execute('COMMIT')
        conn.close()

        # connect to DB
        conn = self.connect_db(db=self.db_name)
        conn.execute('''BEGIN''')
        # insert atoms objects into DB
        r_pk  = self.add_atoms(reactant, mode, conn)
        s_pk  = self.add_atoms(saddle, mode, conn)
        p_pk  = self.add_atoms(product, mode, conn)
        # insert process to DB
        conn.execute('''INSERT INTO Process 
                        VALUES ('%d','%s','%d','%d','%d','%d')''' % 
                     (pro_pk, name, o_pro_pk, r_pk, s_pk, p_pk))
        # insert mobile atoms list into DB
        self.add_mobile(mobile_set, pro_pk, conn)
        conn.execute('COMMIT')
        conn.close()
        return pro_pk

    # adds a process to the kdb table only. Doesn't create a backup. This function
    # should only be used within kdbupdate.
    def add_updated_process(self, reactant, saddle, product, mode, mobile_set, o_pro_pk):
        symbols = reactant.get_chemical_symbols()
        name = self.get_name(symbols)
        pro_pk = self.get_id('Process')
        conn = self.connect_db()
        conn.execute('''BEGIN''')
        r_pk  = self.add_atoms(reactant, mode, conn)
        s_pk  = self.add_atoms(saddle, mode, conn)
        p_pk  = self.add_atoms(product, mode, conn)
        conn.execute('''INSERT INTO Process 
                        VALUES ('%d','%s','%d','%d','%d','%d')''' % 
                     (pro_pk, name, o_pro_pk, r_pk, s_pk, p_pk))
        self.add_mobile(mobile_set, pro_pk, conn)
        conn.execute('COMMIT')
        conn.close()
        return pro_pk


     
    # adds a specific ase.Atom object (SHOULD NOT BE USED)
    # this function is only to be used inside add_process()    
    def add_atoms(self, atoms, mode, conn):
        # if mode = None, create array of zeros.
        if mode == None:
            mode = np.zeros((len(atoms),3))
        # get important data from ase.Atoms instance
        symbols     = atoms.get_chemical_symbols()
        positions   = atoms.get_positions()
        cell        = atoms.get_cell()
        constraints = atoms._get_constraints()
        # get unique primary key for the correct database
        conn.execute('SELECT DATABASE()')
        db_name = conn.fetchone()[0]
        db = ''
        if db_name == self.backup_db_name:
            db = 'backup'
        atoms_pk = self.get_id('Atoms', db=db)
        # add data to Atoms table
        conn.execute('''INSERT INTO Atoms VALUES ('%d','%f','%f','%f','%f','%f','%f','%f','%f','%f')''' %
                     (atoms_pk, 
                      cell[0][0],cell[0][1],cell[0][2],
                      cell[1][0],cell[1][1],cell[1][2],
                      cell[2][0],cell[2][1],cell[2][2]))
        # for each atom in atoms add the atom to the database
        for i in range(len(atoms)):
            atom_pk  = self.get_id('Atom', db=db)
            index    = i
            sub_mode = mode[i]
            symbol   = symbols[i]
            position = positions[i]
            if i in constraints[0].index:
                fixed = 1
            else:
                fixed = 0
            conn.execute('''INSERT INTO Atom VALUES ('%d','%d','%d','%s','%f','%f','%f','%f','%f','%f','%d')''' %
                         (atom_pk, atoms_pk, index, symbol, 
                          position[0], position[1], position[2], 
                          sub_mode[0], sub_mode[1], sub_mode[2], fixed))
        return atoms_pk

    # mobile_set is a list of atom indicies, pro_pk is a unique identifier
    # of the process that contains the mobile atoms
    def add_mobile(self, mobile_set, pro_pk, conn):
        for mobile in mobile_set:
            mob_pk = self.get_id('Mobile')
            conn.execute('''INSERT INTO Mobile VALUES ('%d','%d','%d')''' % 
                         (mob_pk, mobile, pro_pk))

    # add a user to the kdb_user database.
    # not all emails are hashed using SHA.
    def add_user(self, f_name, l_name, email, passwd):
        if not self.is_email(email):
            print "please enter a valid email."
            return "please enter a valid email."
        user_pk = self.get_id('User', db='user')
        conn = self.connect_db(db=self.user_db_name)
        try:
            conn.execute('''INSERT INTO User VALUES ('%d', '%s', '%s', '%s', SHA('%s'))''' %
                        (user_pk, f_name, l_name, email, passwd))
            conn.execute('''COMMIT''')
            conn.close()
            print "account added"
            return "account added"
        except pymysql.err.IntegrityError:
            print "This email address is already in our database."
            return "This email address is already in our database."
            conn.close()



    #############################################################
    # Queries ###################################################
    #############################################################

    # query db for ase.Atoms object based off atoms_id        
    def get_atoms(self, atoms_id, db=db_name):
        conn = self.connect_db(db=db)
        # query the db for the specific atoms entry
        conn.execute('''SELECT * FROM Atoms WHERE atoms_id = '%d' ''' % atoms_id)
        atoms = conn.fetchone()
        cell = []
        # populate the 3D cell block
        count = 0
        temp = []
        for i in range(len(atoms)):
            if i !=0:
                temp.append(atoms[i])
                if count == 2:
                    cell.append(temp)
                    temp = []
                    count = 0
                else:
                    count += 1
        # query the db for all atom entries with the same atoms_id as given
        conn.execute('''SELECT * FROM Atom WHERE atoms_id = '%d' '''% atoms_id)
        atoms = conn.fetchall()
        conn.close()
        symbols = []
        positions = []
        fixed = []
        mode = []
        # loop over each atom from query and pull out data needed to build ase.Atoms instance
        for i in range(len(atoms)):
            symbols.append(str(atoms[i][3]))
            positions.append([atoms[i][4],atoms[i][5],atoms[i][6]])
            mode.append([atoms[i][7],atoms[i][8],atoms[i][9]])
            fixed.append(atoms[i][10])
        # create FixAtoms mask for ase.Atoms' constraint
        constraint = FixAtoms(mask=fixed)
        # create ase.Atoms instance
        atoms = Atoms(symbols = symbols, positions = positions, cell = cell, constraint = constraint)
        #self.write_atoms('test.con', atoms)
        return atoms
    
    #generate a mode list from a given atoms_id
    def get_mode(self, atoms_id, db=db_name):
        conn = self.connect_db(db=db)
        conn.execute(''' SELECT mode0, mode1, mode2 FROM Atom WHERE atoms_id = '%d' ORDER BY num''' % atoms_id)
        mode_list = conn.fetchall()
        conn.close()
        mode = np.array(mode_list)
        return mode
            
    # name is the chemical name for the stucture, IE: Al or CuO
    def get_process(self, name):
        # query db for all processes with same name as given in params
        conn = self.connect_db()
        return_list = []
        conn.execute('''SELECT * FROM Process WHERE name = '%s' ''' % name)
        process_list = conn.fetchall()
        # iterate through all process matches and extract important data
        for process in process_list:
            # dictionary to hold important data
            pro_dict = {'minimum': None, 
                        'saddle' : None, 
                        'product': None, 
                        'mobile' : [], 
                        'mode'   : None, 
                        'mirror' : False,
                        'id'     : None}
            # query database for mobile atoms list of the current process
            conn.execute('''SELECT * FROM Mobile WHERE pro_id = '%d' ''' % process[0])
            mobile_list = conn.fetchall()
            # add the mobile atoms to dictionary
            for mobile in mobile_list:
                pro_dict['mobile'].append(mobile[1])
            # convert atoms_ids to ase.Atoms() objects
            pro_dict['minimum'] = self.get_atoms(process[3])
            pro_dict['saddle']  = self.get_atoms(process[4])
            pro_dict['product'] = self.get_atoms(process[5])
            pro_dict['mode']    = self.get_mode(process[3])
            pro_dict['id']      = process[0]
            return_list.append(pro_dict)
            # swap min and product and re add dictionary
            swap_dict = copy.deepcopy(pro_dict)
            swap_dict['minimum'], swap_dict['product'] = swap_dict['product'], swap_dict['minimum']
            return_list.append(swap_dict)
            #more lists to append for mirror changes
            mirror_dict = copy.deepcopy(pro_dict)
            mirror_dict['mirror'] = True
            return_list.append(mirror_dict)
            mirror_dict2 = copy.deepcopy(swap_dict)
            mirror_dict2['mirror'] = True
            return_list.append(mirror_dict2)
        conn.close()
        return return_list

    # function grabs all params in the param table and returns them
    def get_params(self):
        conn = self.connect_db(db=self.backup_db_name)
        conn.execute('''SELECT * FROM Param''')
        param_list = conn.fetchall()
        conn.close()
        return_dict = {}
        for tup in param_list:
            return_dict[tup[0]] = tup[1]
        return return_dict

    # function used to check for duplicate entries
    def get_saddles(self, name):
        # query db for all processes
        conn = self.connect_db()
        conn.execute('''SELECT saddle_id FROM Process WHERE name = '%s' ''' % name)
        sid_list = conn.fetchall()
        saddle_list = []
        for sid in sid_list:
            saddle = self.get_atoms(sid[0])
            saddle_list.append([saddle, sid[0]])
        return saddle_list

    # function grabs the user_id given an email/passwd combo
    def get_user_id(self, email, passwd):
        if not self.is_user(email, passwd):
            return None
        conn = self.connect_db(db=self.user_db_name)
        conn.execute('''SELECT user_id FROM User WHERE email = '%s' && password = SHA('%s')''' % (email, passwd))
        return conn.fetchone()[0]

    def is_user(self, email, passwd):
        conn = self.connect_db(db=self.user_db_name)
        conn.execute('''SELECT user_id FROM User WHERE email = '%s' && password = SHA('%s')''' % (email, passwd))
        if conn.fetchone() is None:
            return False
        else:
            return True


    #############################################################
    # Helper Functions ##########################################
    #############################################################
            
    # helper function to remove a process from the database
    #outdated do not use
    def remove_process(self, pro_pk):
        conn = self.connect_db()
        conn.execute('''SELECT * FROM Process
                                  WHERE pro_id = '%d' ''' % pro_pk)
        process = conn.fetchone()
        conn.execute('''DELETE FROM Mobile WHERE pro_id = '%d' ''' % pro_pk)
        conn.execute('''DELETE FROM Atom 
                        WHERE atoms_id in ('%d','%d','%d','%d','%d','%d')''' % (process[3],process[4],process[5],process[6],process[7],process[8]))
        conn.execute('''DELETE FROM Atoms
                        WHERE atoms_id in ('%d','%d','%d','%d','%d','%d')''' % (process[3],process[4],process[5],process[6],process[7],process[8]))
        conn.execute('''DELETE FROM Process WHERE pro_id = '%d' ''' % pro_pk)
        conn.execute('COMMIT')
        conn.close()

    def remove_process_user(self, email):
        #get user_id
        conn = self.connect_db(db=self.user_db_name)
        conn.execute('''SELECT user_id FROM User WHERE email = '%s' ''' % email)
        results = conn.fetchone()
        conn.close()
        if results == None:
            print "%s is not in the database." % email
            sys.exit()
        user_id = results[0]

        #get all backup processes from the user_id
        conn = self.connect_db(db=self.backup_db_name)
        conn.execute('''SELECT pro_id FROM Process WHERE user_id = '%d' ''' % user_id)
        results = conn.fetchall()
        conn.close()
        if results == ():
            print "%s has not added any processes." % email
            sys.exit()
        pro_ids = []
        for x in results:
            pro_ids.append(x[0])

        #loop through proceses and delete them
        for o_pro_id in pro_ids:
            #remove from querable table
            conn = self.connect_db()
            conn.execute('''SELECT pro_id from Process WHERE original_pro_id = '%d' ''' % o_pro_id)
            result = conn.fetchone()
            if result is not None:
                pro_id = result[0]
                conn.execute('''SELECT * FROM Process WHERE pro_id = '%d' ''' % pro_id)
                process = conn.fetchone()
                conn.execute('''DELETE FROM Mobile WHERE pro_id = '%d' ''' % pro_id)
                conn.execute('''DELETE FROM  Atom WHERE atoms_id in ('%d', '%d', '%d')''' % (process[2],process[3],process[4]))
                conn.execute('''DELETE FROM Process WHERE pro_id = '%d' ''' % pro_id)
                conn.execute('''DELETE FROM  Atoms WHERE atoms_id in ('%d', '%d', '%d')''' % (process[2],process[3],process[4]))
                conn.execute('COMMIT')
                conn.close()
            #remove from nonquerable table
        for o_pro_id in pro_ids:
            conn = self.connect_db(db=self.backup_db_name)
            conn.execute('''SELECT * FROM Process WHERE pro_id = '%d' ''' % o_pro_id)
            process = conn.fetchone()
            conn.execute('''DELETE FROM  Atom WHERE atoms_id in ('%d', '%d', '%d')''' % (process[2],process[3],process[4]))
            conn.execute('''DELETE FROM Process WHERE pro_id = '%d' ''' % o_pro_id)
            conn.execute('''DELETE FROM  Atoms WHERE atoms_id in ('%d', '%d', '%d')''' % (process[2],process[3],process[4]))
            conn.execute('COMMIT')
            conn.close()




    # writes an atoms object to a file as a .con file
    def write_atoms(self, filename, atoms):
        write_vasp(filename, atoms)

    # helper function to get name from list of symbols
    # IE: symbols = ['Cu', 'Cu', 'O'] -> name = "CuO"
    def get_name(self, symbols):
        name = ""
        for symbol in symbols:
            if symbol not in name:
                name = name + symbol
        return name

    # helper function to connect to the database
    def connect_db(self, db=db_name):
        return pymysql.connect(host=host, port=port, user=user, passwd=password, db=db).cursor()

    # helper function to get largest id for a specific table
    def get_max(self, table_name, id_name, db):
        conn = self.connect_db(db=db)
        # query the given table for the largest primary key
        conn.execute('''SELECT max({}) FROM {}'''.format(id_name, table_name))
        max_id = conn.fetchone()[0]
        conn.close()
        # if there are no entries in the database max_id = None 
        if not max_id:
            max_id = 0
        return max_id

    # helper function to get unique primary key number for a specific table
    def get_id(self, table_name, db=''):
        if db == 'backup':
            self.b_pk_dict[table_name] += 1
            return self.b_pk_dict[table_name]
        elif db == 'user':
            self.u_pk_dict[table_name] += 1
            return self.u_pk_dict[table_name]
        else:
            self.pk_dict[table_name] += 1
            return self.pk_dict[table_name]

    # helper function to valid email add_process
    def is_email(self, email):
        if not re.match(r'(.+)@(.+)\.(.{2,4})', email):
            return False
        return True