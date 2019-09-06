import os, sys
import install

# 获取安装目录
home = os.path.expanduser('~')
ins_path = install.install_path if install.install_path != '~' else home
install_path = ins_path if '~/' not in ins_path else os.path.join(home, ins_path.split('~/')[-1])
# 获取lib安装目录
lib_path = None
for path in sys.path:
	if os.path.basename(path) == 'site-packages':
		lib_path = path
		break
if lib_path is None:
	raise NameError('Did not found python lib path when uninstall VaspCZ lib')

print(f'{"install path":<25}{install_path}')
print(f'{"lib path":<25}{lib_path}')


def uninstall():
	print(f'{"":-<20}{"VaspCZ Uninstalling...:^20"}{"":-<}')
	# 卸载VaspCZ
	current_path = os.getcwd()
	os.chdir(install_path)
	if os.path.isdir('VaspCZ'):
		os.system(f'rm -rf VaspCZ')
	if os.path.isdir('vtst'):
		os.system(f'rm -rf vtst')
	os.chdir(lib_path)
	if os.path.isdir('VaspCZ'):
		os.system(f'rm -rf VaspCZ')
	print(f'{"":-<20}{"VaspCZ Uninstalled:^20"}{"":-<}')


if __name__ == '__main__':
	uninstall()
