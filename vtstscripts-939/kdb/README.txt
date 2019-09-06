------------------------------------------------------------------------------
About:
------------------------------------------------------------------------------
KDB is split into two sections: local and remote. Both sections use the same 
kdbinsert and kdbquery classes to do all calculations on atoms objects. However
when kdbinsert inserts or when kdbquery queries, specific local and remote 
functions are called to ensure the data is inserted/retrieved from the correct
database implementation.

The local database is implemented with SQLite3 and primarly to be 
used with eon to assist with AKMC calculations. The remote database is 
implemented with MySQL and can be used more broadly. The MySQL implementation 
will allow users to retrieve saddles that other user's submissions. 

------------------------------------------------------------------------------
How To:
------------------------------------------------------------------------------
-install/setup clientside:
	-You can grab KDB with svn using this command:
	 svn co http://theory.cm.utexas.edu/svn/kdb
	-If you are using KDB through eon or vtsttools kdb comes with these packages
-Stand alone:
	-local:
		-insert:
			Ensure reactant, saddle, and product files are located in your directory
			python local_client.py insert reactant saddle product
			note reactant, saddle, and product are the specific file names
		-query:
			Much like insert:
			python local_client.py query reactant
	-remote:
		-insert:
			python remote_client.py insert reactant saddle product
			note you will be asked to create an account on your first insert
			just follow the prompts.
		-query:
			python remote_client.py query reactant
			note you do not need an account to query
-eon:
	The only useable version of kdb in eon is the local database. To use kdb
	to assist your AKMC calculations edit your config.ini file with the
	following text:
	[KDB]
	use_kdb = true

-install server:
	If you plan on hosting your own server or are in charge of maintaining the
	DB on theory this is the right place to read. 
	-install MySQL and remember the root username and password for mysql
	-install screen (or any similiar program that will keep processes running
	 after exiting terminal)
	-on the server machine open server_config.py in an editor and set the
	 root username/password and port you are running mysql on (3306 is default)
	-create the databases (python remote_initialize.py)
	-run bottle to accept HTTP requests (python routes.py)
	-you should see info about bottle and now the server is running
-server maintenance:
	-note server maintenance can not be done on the client.
	-purge the database. This will completely remove all kdb related data
	 including kdb accounts. (python remote_initialize.py purge)
	-remove all processes by username with this command:
	 python remove_initialize.py userpurge email 
	 where email is the users email address

------------------------------------------------------------------------------
Development Notes:
------------------------------------------------------------------------------ 
KDB is written so that changes to the math behind KDB inserts and queries only
has to be changed once to be reflected across both the local and remote
versions. The class structure is as follows: 

-kdb is the base class that holds functions for kdbinsert and kdbquery. 
	-kdbinsert is a child of kdb and does all the atoms manipulation before
	 the object is inserted into a database.
	 	-remote_insert is a child of kdbinsert and overloads the
	 	 insert_into_db function so that it inserts into the MySQL database.
	 	-local_insert is a child of kdbinsert and overloads the 
	 	 insert_into_db function so that it inserts into the SQLite database.
	-kdbquery is a child of kdb and does all the atoms manipulation before
	 the object is queried from a database.
	 	-remote_query is a child of kdbquery and overloads the query_db
	 	 function so that it queries the MySQL database.
	 	-local_query is a child of kdbquery and overloads the query_db
	 	 function so that is queries the SQLite database.

Insert functions should only be called from either LocalInsert or RemoteInsert
classes. If insert is called from the KdbInsert class it will return 
an error message because insert_into_db has not been overloaded yet.

Similary query functions should only be called from either LocalQuery or 
RemoteQuery classes to ensure query_db has been overloaded.

All MySQL server side interactions are handled with the pymysql classes contained in
the remote folder. pymysql is an open-source project that uses 100% python
to make MySQL interactions.

aselite is a stripped down version of the ASE classes. aselite contains 
file I/O for poscar/con files, the Atoms class, and various other ASE classes
used in the KDB.

All remote interactions are done through HTTP requests. The server is running
bottle and listens for HTTP requests on port 8080. 

Entry points to the software are local_client.py and remote_client.py.

------------------------------------------------------------------------------
Test Cases:
------------------------------------------------------------------------------
In the kdb directory is a test directory. Here test.py hosts a few test cases.
python test.py local
will test the local database's ability to insert, check for duplicates, query,
check math on insert, and math on query.

python test.py remote
will test insert and query on the remote database. Note you do not need to 
check the math again because it is the same for both local and remote.
If you want to test remote you first need to create an account using the
remote_client.py file then copy the .kdb file that is created into the test
folder 

python test.py local remote
will test both databases.

------------------------------------------------------------------------------
Future Improvements:
------------------------------------------------------------------------------
- Using SQLAlchemy could replace PyMySql and reduce duplicated code between the
  local_db and remote_db files.
- Instead of running bottle on the server end to interact with MySQL we should
  handle all HTTP requests sever side with Apache.
- reactant comparison in kdbquery could be done more efficiently with another
  implementation. Possibly convert the atoms objects to graphs and run a
  sub-graph isomorphism algorithm to check for similiarites (VF2). However
  this approach has been known to run extremely slow at times.
  
------------------------------------------------------------------------------
Future Improvements from Rye:
------------------------------------------------------------------------------
- should try to use ASE's PBC functions
- getting reactantNeighbors should use memoized distance functions in kdb.py 