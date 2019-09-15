import os


def run():
	maintain_files = 'ini,fin'.split(',')
	files = os.listdir('./')
	del_files = []
	for file in files:
		if file in maintain_files:
			pass
		else:
			del_files.append(file)
	ipt = input(f'Those file/dir will be deleted: {del_files}\nconfirm ([y]es/no): ')
	if ipt in ['y', 'YES', 'yes', 'Y']:
		for file in del_files:
			os.system(f'rm -rf {file}')
	else:
		print(f'Did not del any file.')


if __name__ == '__main__':
	run()