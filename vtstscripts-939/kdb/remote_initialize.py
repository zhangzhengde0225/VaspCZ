import pymysql
import sys
from server_config import *

# Try to create the databases without any tables.
def create_databases():
    if user == '' or password == '' or port == None or host == '':
        print "Cannot initialize Database until server_config.py is properly setup."
        return
    print "creating databases"
    conn = pymysql.connect(host=host, port=port, user=user, passwd=password).cursor()
    try:
        conn.execute('CREATE DATABASE %s' % db_name)
        conn.execute('CREATE DATABASE %s' % backup_db_name)
        conn.execute('CREATE DATABASE %s' % user_db_name)
    except pymysql.err.ProgrammingError:
        print "databases were previoulsy created. Exiting now."
        conn.close()
        sys.exit()
    conn.close()

# creates user table
def create_user_table():
    conn = connect_db(db=user_db_name)
    # user_id is a unique identifier
    # email is unique but not an identifier
    # password is MD5 hashed
    conn.execute('''CREATE TABLE User(user_id INT PRIMARY KEY UNIQUE NOT NULL,
                                      first_name VARCHAR(20) NOT NULL,
                                      last_name VARCHAR(20) NOT NULL,
                                      email VARCHAR(255) NOT NULL UNIQUE,
                                      password VARCHAR(100) NOT NULL)''')
    conn.execute('''UPDATE User SET password = MD5(password)''')
    conn.execute('''COMMIT''')
    conn.close()

# create tables for database that holds queryable Atoms
def create_tables():
    print "creating tables"
    conn = connect_db(db=db_name)
    # atoms_id is a unique identifier
    # atoms_name is the structure name IE: 'Al' or 'CuO'
    # atoms_cellXX are bounds for the 3D cell IE: 00 = top left, 02 = top right 
    conn.execute('''CREATE TABLE Atoms(atoms_id      INT PRIMARY KEY UNIQUE NOT NULL, 
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
    conn.execute('''CREATE TABLE Process(pro_id          INT PRIMARY KEY UNIQUE NOT NULL,
                                         name           TEXT NOT NULL,
                                         original_pro_id INT NOT NULL,
                                         reactant_id     INT NOT NULL,
                                         saddle_id       INT NOT NULL,
                                         product_id      INT NOT NULL,
                                         FOREIGN KEY(reactant_id)     REFERENCES Atoms(atoms_id),
                                         FOREIGN KEY(saddle_id)       REFERENCES Atoms(atoms_id),
                                         FOREIGN KEY(product_id)      REFERENCES Atoms(atoms_id))''')
    # this foreign key refers to a Process table in the backup database
    conn.execute('''ALTER TABLE Process
                    ADD FOREIGN KEY fk_name(original_pro_id) REFERENCES %s.Process(pro_id)
                    ON DELETE CASCADE''' % backup_db_name)

    # mob_id is a unqiue identifer
    # num is the atom number 
    # pro_id is a reference to the process the atom belongs to
    conn.execute('''CREATE TABLE Mobile(mob_id INT PRIMARY KEY UNIQUE NOT NULL,
                                        num    INT NOT NULL,
                                        pro_id INT NOT NULL,
                                        FOREIGN KEY (pro_id) REFERENCES Process(pro_id))''')

    # commit and close the connection
    conn.execute('COMMIT')
    conn.close()

# creates tables to hold backups. NOT queryable objects.
def create_backup_tables():
    print "creating  backup tables"
    conn = connect_db(db=backup_db_name)
    
    # atoms_id is a unique identifier
    # atoms_name is the structure name IE: 'Al' or 'CuO'
    # atoms_cellXX are bounds for the 3D cell IE: 00 = top left, 02 = top right 
    conn.execute('''CREATE TABLE Atoms(atoms_id      INT PRIMARY KEY UNIQUE NOT NULL, 
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
    conn.execute('''CREATE TABLE Process(pro_id          INT PRIMARY KEY UNIQUE NOT NULL,
                                         name           TEXT NOT NULL,
                                         reactant_id     INT NOT NULL,
                                         saddle_id       INT NOT NULL,
                                         product_id      INT NOT NULL,
                                         user_id         INT NOT NULL,
                                         FOREIGN KEY(reactant_id) REFERENCES Atoms(atoms_id),
                                         FOREIGN KEY(saddle_id)   REFERENCES Atoms(atoms_id),
                                         FOREIGN KEY(product_id)  REFERENCES Atoms(atoms_id))''')

    conn.execute('''ALTER TABLE Process
                    ADD FOREIGN KEY fk_name(user_id) REFERENCES %s.User(user_id)
                    ON DELETE CASCADE''' % user_db_name)

    # config_option is the name of the configuration, IE: nf, dc, or mac
    # config_value is the value of the corresponding corfigutation option
    conn.execute('''CREATE TABLE Param(config_option varchar(25) UNIQUE NOT NULL,
                                       config_value  REAL NOT NULL)''')
    
    # add default parameters for Params table
    default_values = ('nf', .2, 'dc', .3, 'mac', .7)
    conn.execute('''INSERT INTO Param (config_option, config_value)
                    VALUES ('%s','%f'), ('%s','%f'), ('%s','%f')''' % default_values)

    # commit and close the connection
    conn.execute('COMMIT')
    conn.close()

# function used to completely purge the database.
# generally should only be used for testing.
def purge_db():
    print "purging database"
    conn = pymysql.connect(host=host, port=port, user=user, passwd=password).cursor()
    conn.execute('DROP DATABASE %s' % db_name)
    conn.execute('DROP DATABASE %s' % backup_db_name)
    conn.execute('DROP DATABASE %s' % user_db_name)
    conn.close()

# standard connect to database function.
# this should maybe just be put in kdb.py (but it would have to be overloaded for remote/local configuations)
def connect_db(db=db_name):
    return pymysql.connect(host=host, port=port, user=user, passwd=password, db=db).cursor()


if __name__ == "__main__":
    args = sys.argv
    if len(args) > 1:
        if 'userpurge' in args:
            import remote_db
            db = remote_db.RemoteDB()
            if len(args)>2:
                db.remove_process_user(args[2])
                print "removed"
            else:
                print "You need to specify email."
        elif 'purge' in args:
            purge_db()
            create_databases()
            create_user_table()
            create_backup_tables()
            create_tables()
    else:   
        create_databases()
        create_user_table()
        create_backup_tables()
        create_tables()