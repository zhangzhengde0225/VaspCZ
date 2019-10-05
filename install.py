import os, sys
import argparse
import subprocess

# Default Configuration
install_path = '~/bin'  # VaspCZ的安装目录，默认在用户/bin下
sh_path = '~'  # 程序中需要提交Vasp计算程序时，qsub Vasp.sh命令中Vasp.sh的目录，默认是用户根目录
pseudopotential_path = '~'
shortcut = 'vcz'  # 程序快捷键
vtst = True  # 是否安装VTST tools工具

__version__ = '1.0.2'


def write_path_to_bashrc(path):
	with open(os.path.expanduser('~')+'/.bashrc', 'r') as f:
		data = f.readlines()
	flag = 0
	for line in data:
		if path in line:
			flag = 1
			print(f'path:{path} already exists in .bashrc')
			break
	if flag == 0:
		print(f'write:{path} to .bashrc')
		data.append('\n')
		data.append(f'#{path.split("/")[-1]}\n')
		data.append(f'export PATH={path}:$PATH\n')
	with open(os.path.expanduser('~') + '/.bashrc', 'w') as f:
		f.writelines(data)


def install(prefix, shortcut):
	print('installing VaspCZ software...')
	if not os.path.isdir(prefix):
		os.makedirs(prefix)
	file_path = os.getcwd()  # 获取安装文件的路径

	os.chdir(prefix)  # 进到要安装的目录
	print(f'install path: {prefix}')
	if os.path.isdir('VaspCZ'):
		os.system('rm -rf VaspCZ')
	os.mkdir('VaspCZ')
	os.chdir('VaspCZ')

	# 源文件
	if os.path.isdir('sourcecode'):
		os.system(f'rm -rf sourcecode')
	os.mkdir('sourcecode')
	os.system('cp -rf '+file_path+'/sourcecode/* ./sourcecode/')  # 拷贝源文件
	# 例子

	os.system(f'cp -rf {file_path}/examples.zip .')
	os.system(f'unzip examples.zip')
	# if os.path.isdir('examples'):
	# 	os.system(f'rm -rf examples')
	# os.mkdir(f'examples')
	# os.system(f'cp -rf {file_path}/examples/* ./examples/')

	# 修改Main.py的表头
	with open(f'sourcecode/VaspCZ{__version__}.py', 'r') as f:
		data = f.readlines()
	if '#!' in data[0]:
		data[0] = f'#!{sys.executable}\n'  # 修改表头为当前系统的python位置
	else:
		data.insert(0, f'#!{sys.executable}\n')  # 如果原来的文件没有，就添加表头为当前系统的python位置
	with open(f'sourcecode/VaspCZ{__version__}.py', 'w') as f:
		f.writelines(data)

	# 做软链接
	for shortc in shortcut.split(','):
		os.system(f'ln -sf sourcecode/VaspCZ{__version__}.py '+shortc)

	# 路径到bashrc里面
	path = os.path.join(prefix, 'VaspCZ')
	write_path_to_bashrc(path)

	# 在安装目录创建文件build-in_data, 目前存储Vasp.sh的路径
	global Vaspsh_path
	global Vasp_Pseudopotential_path
	with open('sourcecode/build-in_data.txt', 'w') as f:
		f.writelines(f'VaspCZ_path={prefix}/VaspCZ\n')
		f.writelines(f'Vaspsh_path={Vaspsh_path}\n')
		f.writelines(f'Vasp_Pseudopotential_path={Vasp_Pseudopotential_path}')
	os.chdir(file_path)
	print('VaspCZ software installed successfully.\n')
	print(f'安装成功，使用完全功能请手动完成以下配置：\n1. 将适合当前计算平台的PBS提交任务脚本拷贝到目录"{Vaspsh_path}"下并命名为Vasp.sh')
	print(f'2. 将你的贋势文件拷贝到目录"{Vasp_Pseudopotential_path}"下并命名：PseudoPotential/[贋势类型]')


def side_vtst(prefix):
	print(f'installing vtst tools...')
	file_path = os.getcwd()  # 安装文件的目录
	os.chdir(prefix)
	print(f'install path: {prefix}')
	if not os.path.isdir('vtst'):
		os.mkdir('vtst')
	path = os.path.join(file_path, 'vtstscripts-939')
	os.system(f'cp -rf {path}/* ./vtst')  # 拷贝vtst
	# 写路径到bashrc
	path = os.path.join(prefix, 'vtst')
	write_path_to_bashrc(path)
	os.chdir(file_path)
	print(f'vtst tools installed successfully.\n')


def install_lib(prefix):
	# 获取当前python的依赖库的路径
	print(f'installing VaspCZ python lib...')
	lib_path = None
	for path in sys.path:
		if os.path.basename(path) == 'site-packages':
			lib_path = path
			break
	if lib_path is None:
		raise NameError('Did not found python lib path when install VaspCZ lib')
	file_path = os.getcwd()  # 安装文件的目录
	print(f'lib path: {lib_path}')
	os.chdir(lib_path)
	if os.path.isdir('VaspCZ'):
		os.system('rm -rf VaspCZ')
	os.mkdir('VaspCZ')
	os.chdir('VaspCZ')
	os.system(f"cp -rf {file_path}/sourcecode/__init__.py .")
	os.system(f"cp -rf {file_path}/sourcecode/zzdlib.py .")
	if os.path.isdir('sourcecode'):
		os.system('rm -rf sourcecode')
	os.mkdir('sourcecode')
	os.system(f"cp -rf {prefix}/VaspCZ/sourcecode/* ./sourcecode/")
	os.chdir(file_path)
	print(f'VaspCZ python lib installed successfully.\n')



if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('-p', '--prefix', default=install_path, help='installation path')  # 安装路径
	parser.add_argument('-s', '--shortcut', default=shortcut, help='Shortcat of the VaspCZ')  # 安装后的快捷键
	parser.add_argument('-v', '--vtst', default=vtst, help='install vtst or not, bool, default=True')
	parser.add_argument('-sh_path', '--Vaspsh_path', default=sh_path, help='the sample vasp.sh path')
	parser.add_argument('-psp_path', '--Vasp_Pseudopotential_path', default=pseudopotential_path)
	args = parser.parse_args()
	home = os.path.expanduser('~')
	Vaspsh_path = os.path.expanduser('~') if args.Vaspsh_path == '~' else args.Vaspsh_path
	Vaspsh_path = Vaspsh_path if '~/' not in Vaspsh_path else os.path.join(home, Vaspsh_path.split('~/')[-1])
	psp_path = os.path.expanduser('~') if args.Vasp_Pseudopotential_path == '~' else args.Vasp_Pseudopotential_path
	Vasp_Pseudopotential_path = psp_path if '~/' not in psp_path else os.path.join(home, psp_path.split('~/')[-1])
	prefix = os.path.expanduser('~') if args.prefix == '~' else args.prefix
	prefix = prefix if '~/' not in prefix else os.path.join(home, prefix.split('~/')[-1])
	print(f'{"":-<20}{"install info":^20}{"":-<20}')
	print(f'{"install path:":<25}{prefix}\n{"shortcut:":<25}{args.shortcut}\n{"install_vtst:":<25}{args.vtst}')
	print(f'{"Vasp.sh path:":<25}{Vaspsh_path}\n{"Pseudopotential_path:":<25}{Vasp_Pseudopotential_path}')
	# 安装程序
	print(f'{"":-<20}{"installing":^20}{"":-<20}')
	install(prefix, args.shortcut)
	# 安装VaspCZ库
	install_lib(prefix)
	# 安装vtsttool
	if args.vtst:
		side_vtst(prefix)
	os.chdir(os.path.expanduser('~'))
	# os.system('source .bashrc')
	subprocess.call(f'source {os.path.expanduser("~")}/.bashrc', shell=True)
	print(f'{"":-<20}{"install finished":^20}{"":-<20}')

