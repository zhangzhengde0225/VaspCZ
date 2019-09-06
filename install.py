import os, sys
import argparse
import subprocess


__version__ = '1.0.1'

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
	if not os.path.isdir(prefix):
		os.makedirs(prefix)
	file_path = os.getcwd()  # 获取安装文件的路径
	os.chdir(prefix)  # 进到要安装的目录
	if not os.path.isdir('VaspCZ'):
		os.mkdir('VaspCZ')
	os.chdir('VaspCZ')

	if os.path.isdir('sourcecode'):
		os.system(f'rm -rf sourcecode')
	os.mkdir('sourcecode')
	os.system('cp -rf '+file_path+'/sourcecode/* ./sourcecode/')  # 拷贝源文件

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
	with open('sourcecode/build-in_data.txt', 'w') as f:
		f.writelines(f'Vaspsh_path={Vaspsh_path}\n')

	os.chdir(file_path)


def side_vtst(prefix):
	file_path = os.getcwd()  # 安装文件的目录
	os.chdir(prefix)
	if not os.path.isdir('vtst'):
		os.mkdir('vtst')
	path = os.path.join(file_path, 'vtstscripts-939')
	os.system(f'cp -rf {path}/* ./vtst')  # 拷贝vtst
	# 写路径到bashrc
	path = os.path.join(prefix, 'vtst')
	write_path_to_bashrc(path)
	os.chdir(file_path)


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('-p', '--prefix', default=os.path.join(os.path.expanduser('~'), 'bin'), help='installation path')  # 安装路径
	parser.add_argument('-s', '--shortcut', default='v, V', help='Shortcat of the VaspCZ')  # 安装后的快捷键
	parser.add_argument('-v', '--vtst', default=True, help='install vtst or not, bool, default=True')
	parser.add_argument('-sh_path', '--Vaspsh_path', default='~', help='the sample vasp.sh path')
	args = parser.parse_args()
	Vaspsh_path = args.Vaspsh_path
	print(f'prefx:{args.prefix} shortcut:{args.shortcut} install_vtst:{args.vtst}')
	install(args.prefix, args.shortcut)
	# 安装vtsttool
	if args.vtst:
		side_vtst(args.prefix)
	os.chdir(os.path.expanduser('~'))
	# os.system('source .bashrc')
	subprocess.call(f'source {os.path.expanduser("~")}/.bashrc', shell=True)

