import sqlite3
from aselite import Atoms
from aselite import FixAtoms
from aselite import write_vasp
import local_update
import numpy
import sys
import copy

class LocalDB():
    def __init__(self, db_name, nf=None, dc=None, mac=None):
        self.db_name = db_name
        self.nf = nf
        self.dc = dc
        self.mac = mac
        # check if db_name is an existing database
        if not self.check_tables():
            # if not initialize the database
            self.create_tables()
        # create a dictionary to use for assaigning unique identifiers
        self.pk_dict = {'Atoms':   self.get_max('Atoms', 'atoms_id'),
                        'Atom':    self.get_max('Atom', 'atom_id'),
                        'Process': self.get_max('Process', 'pro_id'),
                        'Mobile':  self.get_max('Mobile', 'mob_id')
                       }
        # check if parameters are different
        params = self.get_params()
        if self.nf == None:
            self.nf = params['nf']
        if self.dc == None:
            self.dc = params['dc']
        if self.mac == None:
            self.mac = params['mac']
        if self.check_params():
            print "Parameters are different from what is currently in the database."
            answer = raw_input("Do you want to update database entries to reflect the new parameters? (yes, no)")
            if 'y' in answer.lower():
                print "Updating parameters DO NOT exit until update it complete or database could be corrupted."
                self.update_params()
                print "Update complete. Please retry previous command to insert new process."
                sys.exit()
            else:
                print "Structures in the database will be left untouched."

    #database initialization
    def create_tables(self):
        conn = sqlite3.connect(self.db_name)
        # atoms_id is a unique identifier
        # atoms_name is the structure name IE: 'Al' or 'CuO'
        # atoms_cellXX are bounds for the 3D cell IE: 00 = top left 02 = top right 
        conn.execute('''CREATE TABLE Atoms(atoms_id INT PRIMARY KEY UNIQUE NOT NULL, 
                                           atoms_cell00 REAL NOT NULL,
                                           atoms_cell01 REAL NOT NULL,
                                           atoms_cell02 REAL NOT NULL,
                                           atoms_cell10 REAL NOT NULL,
                                           atoms_cell11 REAL NOT NULL,
                                           atoms_cell12 REAL NOT NULL,
                                           atoms_cell20 REAL NOT NULL,
                                           atoms_cell21 REAL NOT NULL,
                                           atoms_cell22 REAL NOT NULL)''')
        # atom_id is a unique identifer
        # atoms_id is a reference value to identify what collection of atoms it belongs in
        # num is the ase.atom.index value
        # symbol is the atom's chemical symbol IE: 'Al' 'Cu'
        # x,y,z_coord is the x,y,z position of the atom
        # fixed is a boolean contraint if the atom is fixed or not
        conn.execute('''CREATE TABLE Atom(atom_id  INT PRIMARY KEY UNIQUE NOT NULL, 
                                          atoms_id INT NOT NULL,
                                          num      INT NOT NULL,  
                                          symbol  TEXT NOT NULL, 
                                          x_coord REAL NOT NULL, 
                                          y_coord REAL NOT NULL, 
                                          z_coord REAL NOT NULL, 
                                          mode0   REAL NOT NULL,
                                          mode1   REAL NOT NULL,
                                          mode2   REAL NOT NULL,
                                          fixed    INT NOT NULL, 
                                          FOREIGN KEY(atoms_id) REFERENCES Atoms(atoms_id))''')
        # pro_id is a unique identifer
        # all the rest of the values are references to atoms_ids
        conn.execute('''CREATE TABLE Process(pro_id INT PRIMARY KEY UNIQUE NOT NULL,
                                             name                TEXT NOT NULL,
                                             is_used              INT NOT NULL,
                                             original_reactant_id INT NOT NULL,
                                             original_saddle_id   INT NOT NULL,
                                             original_product_id  INT NOT NULL,
                                             reactant_id          INT,
                                             saddle_id            INT,
                                             product_id           INT,
                                             FOREIGN KEY (original_reactant_id) REFERENCES Atoms(atoms_id),
                                             FOREIGN KEY (original_saddle_id)   REFERENCES Atoms(atoms_id),
                                             FOREIGN KEY (original_product_id)  REFERENCES Atoms(atoms_id),
                                             FOREIGN KEY(reactant_id)           REFERENCES Atoms(atoms_id),
                                             FOREIGN KEY(saddle_id)             REFERENCES Atoms(atoms_id),
                                             FOREIGN KEY(product_id)            REFERENCES Atoms(atoms_id))''')
        # mob_id is a unqiue identifer
        # num is the atom number 
        # pro_id is a reference to the process the atom belongs to
        conn.execute('''CREATE TABLE Mobile(mob_id INT PRIMARY KEY UNIQUE NOT NULL,
                                            num    INT NOT NULL,
                                            pro_id INT NOT NULL,
                                            FOREIGN KEY (pro_id) REFERENCES Process(pro_id))''')

        # config_option is the name of the configuration, IE: nf, dc, or mac
        # config_value is the value of the corresponding corfigutation option
        conn.execute('''CREATE TABLE Param(config_option TEXT UNIQUE NOT NULL,
                                            config_value   INT NOT NULL)''')
        
        # add default parameters for Params table
        default_values = [('nf', 0.2), ('dc', 0.3), ('mac', .7)]
        conn.executemany('''INSERT INTO Param VALUES (?,?)''', default_values)

        # commit the data and close the database
        conn.commit()
        conn.close()

    #############################################################
    # Insertions ################################################
    #############################################################

    # adds a process, each parameter is an ase.Atom object 
    # EXCEPT mode/o_mode, mode is a numpy array or None 
    # and mobile_set is a list of numbers                                   
    def add_process(self, o_reactant, o_saddle, o_product, o_mode, 
                    reactant, saddle, product, mode, mobile_set):
        conn = self.connect_db()
        conn.execute('''BEGIN''')
        or_pk = self.add_atoms(o_reactant, o_mode, conn)
        os_pk = self.add_atoms(o_saddle, o_mode, conn)
        op_pk = self.add_atoms(o_product, o_mode, conn)
        if reactant:
            r_pk  = self.add_atoms(reactant, mode, conn)
            s_pk  = self.add_atoms(saddle, mode, conn)
            p_pk  = self.add_atoms(product, mode, conn)
            is_used = 1
        else:
            r_pk = None
            s_pk = None
            p_pk = None
            is_used = 0    
        symbols = o_reactant.get_chemical_symbols()
        name = self.get_name(symbols)
        pro_pk = self.get_id('Process')
        conn.execute('''INSERT INTO Process VALUES (?,?,?,?,?,?,?,?,?)''', 
                     (pro_pk, name, is_used, or_pk, os_pk, op_pk, r_pk, s_pk, p_pk)
                    )
        self.add_mobile(mobile_set, pro_pk, conn)
        conn.commit()
        conn.close()
        return pro_pk
     
    # adds a specific ase.Atom object (SHOULD NOT BE USED)
    # this function is only to be used inside add_process()    
    def add_atoms(self, atoms, mode, conn):
        # if mode = None, create array of zeros.
        #if mode == None:
        #if not hasattr(mode, 'shape'):
        if type(mode) != numpy.ndarray and type(mode) != list and type(mode) != tuple: # MJW fix
            mode = numpy.zeros((len(atoms),3))
        # get important data from ase.Atoms instance
        symbols     = atoms.get_chemical_symbols()
        positions   = atoms.get_positions()
        cell        = atoms.get_cell()
        constraints = atoms._get_constraints()
        # get unique primary key
        atoms_pk    = self.get_id('Atoms')
        # add data to Atoms table
        conn.execute('''INSERT INTO Atoms VALUES (?,?,?,?,?,?,?,?,?,?)''',
                     (atoms_pk, 
                      cell[0][0],cell[0][1],cell[0][2],
                      cell[1][0],cell[1][1],cell[1][2],
                      cell[2][0],cell[2][1],cell[2][2]))
        for i in range(len(atoms)):
            atom_pk  = self.get_id('Atom')
            index    = i
            sub_mode = mode[i]
            symbol   = symbols[i]
            position = positions[i]
            if i in constraints[0].index:
                fixed = 1
            else:
                fixed = 0
            conn.execute('''INSERT INTO Atom VALUES (?,?,?,?,?,?,?,?,?,?,?)''',
                         (atom_pk, atoms_pk, index,  symbol, 
                          position[0], position[1], position[2], 
                          sub_mode[0], sub_mode[1], sub_mode[2], fixed))
        return atoms_pk

    # mobile_set is a list of atom indicies, pro_pk is a unique identifier
    # of the process that contains the mobile atoms
    def add_mobile(self, mobile_set, pro_pk, conn):
        for mobile in mobile_set:
            mob_pk = self.get_id('Mobile')
            conn.execute('''INSERT INTO Mobile VALUES (?,?,?)''', 
                         (mob_pk, mobile, pro_pk))
            
    #############################################################
    # Queries ###################################################
    #############################################################

    # query db for ase.Atoms object based off atoms_id        
    def get_atoms(self, atoms_id):
        conn = self.connect_db()
        # query the db for the specific atoms entry
        atoms = conn.execute('SELECT * FROM Atoms WHERE atoms_id = ?', (str(atoms_id),) ).fetchone()
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
        atoms = conn.execute('SELECT * FROM Atom WHERE atoms_id = ? ORDER BY num', (str(atoms_id),) ).fetchall()
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
        cell = [[1,0,0], [0,1,0], [0,0,1]]
        #atoms = Atoms(symbols = symbols, positions = positions, cell = cell, constraint = constraint)
        atoms = Atoms(symbols = symbols, positions = positions, cell = cell)
        #self.write_atoms('test.con', atoms)
        return atoms
    
    #generate a mode list from a given atoms_id
    def get_mode(self, atoms_id):
        conn = self.connect_db()
        mode_list = conn.execute(''' SELECT mode0, mode1, mode2 FROM Atom WHERE atoms_id = ? ORDER BY num''', (atoms_id,)).fetchall()
        conn.close()
        mode = numpy.array(mode_list)
        return mode
            
    # name is the chemical name for the stucture, IE: Al or CuO
    def get_process(self, name):
        # query db for all processes with same name as given in params
        conn = self.connect_db()
        return_list = []
        process_list = conn.execute('''SELECT * FROM Process WHERE name = ?''', (name,)).fetchall()
        # iterate through all process matches and extract important data
        for process in process_list:
            #check if the process is useable
            if not process[2]:
                continue
            # dictionary to hold important data
            pro_dict = {'minimum': None, 
                        'saddle' : None, 
                        'product': None, 
                        'mobile' : [], 
                        'mode'   : None, 
                        'mirror' : False,
                        'id'     : None}
            # query database for mobile atoms list of the current process
            mobile_list =  conn.execute('''SELECT * FROM Mobile WHERE pro_id = ?''', (process[0],)).fetchall()
            # add the mobile atoms to dictionary
            for mobile in mobile_list:
                pro_dict['mobile'].append(mobile[1])
            # convert atoms_ids to ase.Atoms() objects
            pro_dict['minimum'] = self.get_atoms(process[6])
            pro_dict['saddle']  = self.get_atoms(process[7])
            pro_dict['product'] = self.get_atoms(process[8])
            pro_dict['mode']    = self.get_mode(process[6])
            pro_dict['id']      = process[0]
            return_list.append(pro_dict)
            #pro_dict with mirror set to true
            mirror_dict = copy.deepcopy(pro_dict)
            mirror_dict['mirror'] = True
            return_list.append(mirror_dict)
            # swap min and product and re add dictionary
            swap_dict = copy.deepcopy(pro_dict)
            swap_dict['minimum'], swap_dict['product'] = swap_dict['product'], swap_dict['minimum']
            return_list.append(swap_dict)
            # swap_dict with mirror set to true
            mirror_dict2 = copy.deepcopy(swap_dict)
            mirror_dict2['mirror'] = True
            return_list.append(mirror_dict2)
        conn.close()
        return return_list 

    # This function is used to check for duplicate entries
    def get_saddles(self, name):
        # query db for all processes
        conn = self.connect_db()
        sid_list = conn.execute('''SELECT saddle_id FROM PROCESS WHERE name = ?''', (name,) ).fetchall() 
        conn.close()
        saddle_list = []
        for sid in sid_list:
            saddle = self.get_atoms(sid[0])
            saddle_list.append([saddle, sid[0]])
        return saddle_list

    # Retrieves the MAC, DC, and NF parameters.
    def get_params(self):
        conn = self.connect_db()
        params = conn.execute('''SELECT * FROM Param''').fetchall()
        conn.close()
        param_dict = {}
        for param in params:
            param_dict[str(param[0])] = param[1]
        return param_dict


    #############################################################
    # Helper Functions ##########################################
    #############################################################

    # updates the params to define mobile atom's neighbors and
    # what is considered a mobile atom
    def update_params(self):
        # update Params table
        conn = self.connect_db()
        conn.execute('''DELETE FROM Param''')
        values = [('nf', self.nf), ('dc', self.dc), ('mac', self.mac)]
        conn.executemany('''INSERT INTO Param VALUES (?,?)''', values)
        conn.commit()
        # query for all processes
        process_list = conn.execute('''SELECT * FROM Process''').fetchall()
        # loop over all processes from the query
        for process in process_list:
            o_reactant    = self.get_atoms(process[3])
            o_saddle      = self.get_atoms(process[4])
            o_product     = self.get_atoms(process[5])
            o_mode        = self.get_mode(process[3])
            kdbupdate = local_update.LocalUpdate()
            num = kdbupdate.insert(o_reactant, o_saddle, o_product, mode=o_mode, nf=self.nf, dc=self.dc, mac=self.mac, kdbname=self.db_name)
            if num == 0:
                print "process number", process[0], "is now not useable for queries."
                print "However it will not be deleted from the database."
                # remove reactant, saddle, product, mobile but leave the original
                # reactant, saddle, and product
                conn.execute('''UPDATE Process 
                                SET is_used = 0, reactant_id = NULL, saddle_id = NULL, product_id = NULL
                                WHERE pro_id = ?''', (process[0],) )
                conn.execute('''DELETE FROM Atoms 
                                WHERE atoms_id in (?,?,?)''', (process[6],process[7],process[8]))
                conn.execute('''DELETE FROM Atom
                                WHERE atoms_id in (?,?,?)''', (process[6],process[7],process[8]))
                conn.execute('''DELETE FROM Mobile
                                WHERE pro_id = ?''', (process[0],) )
                conn.commit()
            else:
                self.remove_process(process[0])
                conn.commit()
        conn.close()
            
    # helper function to remove a process from the database
    def remove_process(self, pro_pk):
        conn = self.connect_db()
        process = conn.execute('''SELECT * FROM Process
                                  WHERE pro_id = ?''', (pro_pk,)).fetchone()
        conn.execute('''DELETE FROM Mobile WHERE pro_id = ?''', (pro_pk,))
        conn.execute('''DELETE FROM Atom 
                        WHERE atoms_id in (?,?,?,?,?,?)''', (process[3],process[4],process[5],process[6],process[7],process[8]))
        conn.execute('''DELETE FROM Atoms
                        WHERE atoms_id in (?,?,?,?,?,?)''', (process[3],process[4],process[5],process[6],process[7],process[8]))
        conn.execute('''DELETE FROM Process WHERE pro_id = ?''', (pro_pk,))
        conn.commit()
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
    def connect_db(self):
        # check if tables exist, if false create tables
        return sqlite3.connect(self.db_name)
        
    # helper function to decide if database needs to create tables or not
    def check_tables(self):
        conn = sqlite3.connect(self.db_name)
        # query master table (created by sqlite3) for other tables
        tables = conn.execute(''' SELECT name FROM sqlite_master WHERE type='table' ''').fetchall()
        conn.close()
        # tables is a list and will return false if empty
        if tables:
            return True
        return False

    # helper function to get largest id for a specific table
    def get_max(self, table_name, id_name):
        conn = self.connect_db()
        # query the given table for the largest primary key
        max_id = conn.execute('''SELECT max({}) FROM {}'''.format(id_name, table_name)).fetchone()[0]
        conn.close()
        # if there are no entries in the database max_id = None 
        if not max_id:
            max_id = 0
        return max_id

    # helper function to get unique primary key number for a specific table
    def get_id(self, table_name):
        self.pk_dict[table_name] += 1
        return self.pk_dict[table_name]

    # helper function to check if given params are different from params in database.
    def check_params(self):
        conn = self.connect_db()
        params = conn.execute('''SELECT * FROM Param''').fetchall()
        for param in params:
            if 'nf'  == param[0].lower() and self.nf  != param[1]:
                conn.close()
                return True
            if 'dc'  == param[0].lower() and self.dc  != param[1]:
                conn.close()
                return True
            if 'mac' == param[0].lower() and self.mac != param[1]:
                conn.close()
                return True
            conn.close()    
            return False

