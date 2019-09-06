import getpass

# remote config is used to create new accounts on the database
class RemoteConfig():
	def config(self):
		f_name = raw_input('Please enter your first name: ')
		l_name = raw_input('Please enter your last name: ')
		email = raw_input('Please enter your email address: ')
		passwd = getpass.getpass('Please enter your password (text will not show up): ')
		return [f_name, l_name, email, passwd]

	def is_yes(self, string):
		if 'y' in string.lower():
			return True
		else:
			return False