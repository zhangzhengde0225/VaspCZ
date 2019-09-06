from bottle import route, post, run, request
import remote_insert, remote_query, remote_db
import pickle

host = '0.0.0.0'
port = 8080

@post('/account_create')
def create_account():
	email    = request.forms.get('email')
	password = request.forms.get('password')
	first    = request.forms.get('first')
	last     = request.forms.get('last')

	db = remote_db.RemoteDB()
	output = db.add_user(first, last, email, password)
	return output

@post('/insert')
def insert():
	db = remote_db.RemoteDB()
	if not db.is_user(request.forms.get('email'),request.forms.get('password')):
		return "invalid account info"
	reactant = pickle.loads(request.forms.get('reactant'))
	saddle = pickle.loads(request.forms.get('saddle'))
	product = pickle.loads(request.forms.get('product'))
	try:
		mode = pickle.loads(request.forms.get('mode'))
	except:
		mode = None
	insert_class = remote_insert.RemoteInsert()
	insert_class.email = request.forms.get('email')
	insert_class.password = request.forms.get('password')
	output = insert_class.insert(reactant, saddle, product, mode)
	return str(output)

@post('/query')
def query():
	reactant = pickle.loads(request.forms.get('reactant'))
	query_class = remote_query.RemoteQuery()
	query_class.query(reactant)
	output = pickle.dumps(query_class.return_dict)
	return output

def start():
	run(host=host, port=port)

if __name__ == "__main__":
	start()


