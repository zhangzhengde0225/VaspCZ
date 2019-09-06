import urllib
import httplib
import pickle
from aselite import read_any, write_vasp
import sys
import os
import shutil
import getpass
import remote_config
from kdb import Kdb

host = 'theory.cm.utexas.edu'
port = 8080

def server_create_account():
    #grab info from user
    info = remote_config.RemoteConfig().config()
    #create/populate dictionary
    params = {}
    params['first']    = info[0]
    params['last']     = info[1]
    params['email']    = info[2]
    params['password'] = info[3]
    #format for http post
    params  = urllib.urlencode(params)
    headers = {'Content-type': 'application/x-www-form-urlencoded', 'Accept': 'text/plain'}
    conn    = httplib.HTTPConnection(host=host, port=port)
    #send http POST request
    conn.request('POST', '/account_create', params, headers)
    #grab results
    response = conn.getresponse()
    #print response.status, response.reason
    data = response.read()
    #if account created store email/password for later use
    if data == "account added":
        print data
        set_account_info([info[2],info[3]])
        return 1
    else:
        print data
        answer = raw_input('Would you like to try again? [y/n] ')
        if 'y' in answer.lower()    :
            out = server_create_account()
        else:
            return 0
    return out


def get_account_info():
    with open('.kdb', 'rb') as infile:
        output = pickle.load(infile)
    return output

def set_account_info(info):
    with open('.kdb', 'w') as outfile:
            pickle.dump([info[0], info[1]], outfile)

def server_insert(reactant, saddle, product, mode):
    try:
        email_pass = get_account_info()
    except IOError:
        answer = raw_input('no account information found. Do you have an account? [y/n] ')
        if 'y' in answer.lower():
            email_pass = []
            email_pass.append(raw_input('email: '))
            email_pass.append(getpass.getpass('password: '))
            set_account_info(email_pass)
            print "Attempting to insert process."
        else:
            print "No problem, lets create one."
            if server_create_account():
                email_pass = get_account_info()
            else:
                return 0

    params = {}
    params['email']    = email_pass[0]
    params['password'] = email_pass[1]
    params['reactant'] = pickle.dumps(reactant)
    params['saddle']   = pickle.dumps(saddle)
    params['product']  = pickle.dumps(product)
    if mode is not None:
        params['mode'] = pickle.dumps(mode)
    params  = urllib.urlencode(params)
    headers = {'Content-type': 'application/x-www-form-urlencoded', 'Accept': 'text/plain'}
    conn    = httplib.HTTPConnection(host=host, port=port)
    conn.request('POST', '/insert', params, headers)
    response = conn.getresponse()
    #print response.status, response.reason
    data = response.read()
    print data
    if data == "invalid account info":
        os.remove('.kdb')

def server_query(reactant):
    params = {}
    params['reactant'] = pickle.dumps(reactant)
    params  = urllib.urlencode(params)
    headers = {'Content-type': 'application/x-www-form-urlencoded', 'Accept': 'text/plain'}
    conn = httplib.HTTPConnection(host=host, port=port)
    conn.request('POST', '/query', params, headers)
    response = conn.getresponse()
    #print response.status, response.reason
    data = response.read()
    suggestion_dict = pickle.loads(data)
    try:
        os.mkdir('kdbmatches')
    except:
        shutil.rmtree('kdbmatches')
        os.mkdir('kdbmatches')
    for key in suggestion_dict:
        write_vasp('kdbmatches' + "/SADDLE_%d" % key, suggestion_dict[key][0])
        write_vasp('kdbmatches' + "/PRODUCT_%d" % key, suggestion_dict[key][1])
        try:
            Kdb().save_mode('kdbmatches' + "/MODE_%d" % key, suggestion_dict[key][2])
        except:
            pass
    print "done, output now in kdbmatches/"

def run(args):
    if not Kdb().check_version():
        sys.exit()
    if len(args) < 1:
        print "first parameter sohuld be either: inert or query"
        sys.exit()
    if args[0] == 'insert':
        if len(args) < 4:
            print "parameters for insert should include reactant, saddle, and product files."
            sys.exit()
        try:
            reactant = read_any(args[1])
            saddle   = read_any(args[2])
            product  = read_any(args[3])
        except IOError:
            print "One or more files could not be read."
            sys.exit()
        try:
            mode = Kdb().load_mode(args[3])
        except:
            mode = None

        server_insert(reactant, saddle, product, mode)
    elif args[0] == 'query':
        if len(args) < 2:
            print "parameters for query should include a reactant file."
            sys.exit()
        try:
            reactant = read_any(args[1])
        except IOError:
            print "could not read reactant file."
            sys.exit()
        server_query(reactant)


if __name__ == "__main__":
    args = sys.argv[1:]
    run(args)
