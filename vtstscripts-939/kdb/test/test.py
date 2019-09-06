# NOTE: to test MySQL/anything in remote, you must first 
#       install MySQL-server (atleast until I setup the DB on Theory server)
#       then run remote_initialize.py

import os
import sys
import commands
import filecmp

#get args
args = sys.argv

if len(args) <= 1:
    print "test.py requires atleast one arguement. Please 'local' and or 'remote' as parameters."

if "local" in args:
    print "testing local"
    #remove any old queries and or sqlite3 databases
    try:
        os.system("rm -rf kdbmatches")
        os.system('rm kdb.db')
    except:
        pass
    # sqlite3 tests
    # test insert function with sqlite3 database
    print "Testing insert."
    out1 = commands.getoutput("python ../local_client.py insert test_vars/reactant.con test_vars/saddle.con test_vars/product.con --mode test_vars/mode.dat")
    print out1
    print ""
    print ""
    # test check for duplicates with sqlite3 database
    print "Testing duplicate check."
    out2 = commands.getoutput("python ../local_client.py insert test_vars/reactant.con test_vars/saddle.con test_vars/product.con --mode test_vars/mode.dat")
    print out2
    print ""
    print ""
    # test query function with sqlite3 database
    print "Testing query."
    out3 = commands.getoutput("python ../local_client.py query test_vars/reactant.con")
    print out3
    print ""
    print ""

    print "Testing insert math"
    if filecmp.cmp('kdb.db', 'test_vars/kdb.db'):
        print "math good"
    else:
        print "insert math different from expected"

    print "Testing query math"
    common_files = ['PRODUCT_0','PRODUCT_1','PRODUCT_2','PRODUCT_3','PRODUCT_4','PRODUCT_5','PRODUCT_6','PRODUCT_7',
                    'SADDLE_0','SADDLE_1','SADDLE_2','SADDLE_3','SADDLE_4','SADDLE_5','SADDLE_6','SADDLE_7']
    result = filecmp.cmpfiles('kdbmatches', 'test_vars/kdbmatches', common_files)
    if not result[1] and not result[2]:
        print "math good"
    elif not result[1] and result[2]:
        print "some errors occured:"
        print result[2]
    else:
        print "the list of files that are different."
        print result[1]


if "remote" in args:
    print "testing insert"
    out1 = commands.getoutput("python ../remote_client.py insert test_vars/reactant.con test_vars/saddle.con test_vars/product.con --mode test_vars/mode.dat")
    print out1
    print ""
    print ""

    print "testing query"
    out2 = commands.getoutput("python ../remote_client.py query test_vars/reactant.con")
    print out2
    print ""
    print ""