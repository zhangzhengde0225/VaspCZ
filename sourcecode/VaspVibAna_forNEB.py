"""
NEB的振动分析
程序运行在NEB目录，需要NEB完成
读取ini/CONTCAR fin/CONTCAR  鞍点(一般02)/CONTCAR
对比ini/CONTCAR fin/CONTCAR 的原子位置，确定是第几个是扩散元素。
对ini/CONTCAR saddle/CONTCAR 扩散元素T T T，非扩散元素F F F
修改INCAR进行振动分析
"""

import os, sys
import subprocess
import zzdlib
import argparse
import numpy as np

python = sys.executable
current_py_folder = os.path.dirname(os.path.abspath(__file__))
VaspCZ_path = [os.path.dirname(current_py_folder) if 'sourcecode' in current_py_folder else current_py_folder][0] + '/sourcecode'


def modify_POSCAR(data, index):
	"""
	根据输入的数据和索引修改POSCAR，添加Selective Dynamics, 索引所在的位置设置为T T T, 其他位置设置为 F F F
	:return:
	"""
	POSCAR_data = []
	direct, direct_index = zzdlib.File.getLine(data, keywords='Direct')
	decoded_data = zzdlib.Vasp.decode_POSCAR(data)
	number_of_atom = decoded_data[2]
	for i in range(len(data)):
		if i < direct_index:
			POSCAR_data.append(data[i])  # 前面的部分
		elif i == direct_index:  # Direct部分，要加一个Selective Dynamics
			POSCAR_data.append('Selective Dynamics\n')
			POSCAR_data.append(data[i])
		elif direct_index < i <= direct_index+np.sum(number_of_atom):  # 原子位置部分
			if i != direct_index+index+1:
				POSCAR_data.append(data[i].strip('\n') + ' F F F\n')
			else:
				POSCAR_data.append(data[i].strip('\n') + ' T T T\n')
		else:  # 最后的部分
			POSCAR_data.append(data[i])
	return POSCAR_data


def get_saddle_image():
	if not os.path.isfile('neb.dat'):
		os.system('nebbarrier.pl')
	with open('neb.dat', 'r') as f:
		data = f.readlines()
	# 解码成n行，5列的numpy数组
	data = [line.split() for line in data]
	data = np.array(data).astype(float)
	# print(data, data.shape)
	index = np.argmax(data, axis=0)[2]
	saddle_image = data[index][0]
	saddle_image = f'0{int(saddle_image)}'
	return saddle_image


def run(include_fin=False):
	# 检查NEB是否完成
	out = subprocess.check_output([python, f'{VaspCZ_path}/NEBCheck1.1.py'], shell=False, stderr=subprocess.STDOUT).decode('utf-8')
	isDone = bool([0 if out.find(f'Path:{"./":<40}  NEB计算完成!')==-1 else 1][0])
	# isDone = ['Done' if ("Path:./" in line and "NEB计算完成" in line) for line in out.split('\n') else None]
	# print(out, type(out), isDone, type(isDone))
	if isDone==False:
		print(f'当前目录下NEB计算并未完成，退出程序')
		exit()
	else:
		print(f'检查：当前目录{os.getcwd()}  NEB计算完成\n开始振动分析')

	# 读取文件ini/Opt和fin/Opt下的POSCAR文件，最终获取扩散原子的索引
	with open('./ini/Opt/POSCAR', 'r') as f:
		iniPOS = f.readlines()
	with open('./fin/Opt/POSCAR', 'r') as f:
		finPOS = f.readlines()
	# print(iniCONT)
	ini_result = zzdlib.Vasp.decode_POSCAR(iniPOS)  # 结果是：vector, elements, number_of_atom, position
	fin_result = zzdlib.Vasp.decode_POSCAR(finPOS)
	ini_position = ini_result[3]
	fin_position = fin_result[3]
	# print(ini_position.shape)
	distance = np.sqrt(np.sum(np.square(ini_position-fin_position), axis=1))
	# print(distance)
	index = distance.astype(bool).tolist().index(True)  # 获取了扩散原子的索引号，这个索引号的原子设置为T T T其他为F F F
	# print(index)

	# 把ini的CONTCAR读取为POSCAR分析振动
	with open('ini/CONTCAR', 'r') as f:
		data = f.readlines()
	POSCAR_for_iV = modify_POSCAR(data, index)
	# fin/CONTCAR 读取变为POSCAR分析振动
	with open('fin/CONTCAR', 'r') as f:
		data = f.readlines()
	POSCAR_for_fV = modify_POSCAR(data, index)
	# 鞍点CONTCAR 读取分析振动
	saddle_image = get_saddle_image()
	with open(f'{saddle_image}/CONTCAR', 'r') as f:
		data = f.readlines()
	POSCAR_for_sV = modify_POSCAR(data, index)

	# 开始振动分析
	if not os.path.isdir('vib_analysis'):
		os.mkdir('./vib_analysis')
	os.chdir('./vib_analysis')

	# 初态振动分析
	vib_analysis(f'ini_state', POSCAR_data=POSCAR_for_iV)
	vib_analysis(f'sad_state', POSCAR_data=POSCAR_for_sV)
	if include_fin:
		vib_analysis(f'fin_state', POSCAR_data=POSCAR_for_fV)

	os.chdir('..')


def vib_analysis(folder, POSCAR_data):
	if not os.path.isdir(folder):
		os.mkdir(folder)
	os.chdir(folder)
	os.system(f'cp ../../ini/Opt/INCAR .')
	os.system(f'cp ../../ini/Opt/POTCAR .')
	os.system(f'cp ../../ini/Opt/KPOINTS .')
	Vaspsh_path = zzdlib.File.Vaspsh_path()
	os.system(f'cp {Vaspsh_path}/Vasp.sh .')
	with open('POSCAR', 'w') as f:
		f.writelines(POSCAR_data)

	# print(os.getcwd())
	# 修改INCAR
	data_INCAR = zzdlib.File.openFile('INCAR', 'r')
	data_INCAR = zzdlib.File.substituteData(data_INCAR, keywords='SYSTEM', newline='SYSTEM=Vib\n')  # 修改
	data_INCAR = zzdlib.File.substituteData(data_INCAR, keywords='NSW', newline='NSW=1\n')  # 修改
	data_INCAR = zzdlib.File.substituteData(data_INCAR, keywords='POTIM', newline='POTIM=0.03\n')  # 修改
	data_INCAR = zzdlib.File.substituteData(data_INCAR, keywords='IBRION', newline='IBRION=5\n')  # 修改
	data_INCAR.append('NFREE=2\n')  # 添加
	data_INCAR.append('ISYM=0\n')  # 添加
	data_INCAR.append('PREC=Accurate\n')  # 添加
	data_INCAR = zzdlib.File.substituteData(data_INCAR, keywords='NPAR', newline='\n')  # 删除
	data_INCAR = zzdlib.File.substituteData(data_INCAR, keywords='NCORE', newline='\n')  # 删除
	zzdlib.File.openFile('INCAR', 'w', data=data_INCAR)

	# POTCAR 不修改
	# KPOINS 不修改

	# Vasp.sh
	data_VaspSh = zzdlib.File.openFile('../../ini/Opt/Vasp.sh', 'r')
	jobname_old = zzdlib.File.getLine(data_VaspSh, keywords='#PBS -N')[0].split()[2]
	jobname = jobname_old[:-2]+folder[0]+'V'  # 在原来的基础上去掉最后两个字母，加folder第一个字符，i f s ， 加上V代表振动分析
	global nodes
	global ppn
	zzdlib.Vasp.modify_Vasp_sh(jobname=jobname, nodes=nodes, ppn=ppn)
	zzdlib.Vasp.check_and_qsub()

	os.chdir('..')



if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='manual to this script')
	parser.add_argument('-nd', '--nodes', type=str, default='1')
	parser.add_argument('-np', '--ppn', type=str, default='8')
	parser.add_argument('-fin', '--include_fin', type=str, default='False')
	args = parser.parse_args()
	nodes = args.nodes
	ppn = args.ppn
	include_fin = bool(args.include_fin)
	run(include_fin=include_fin)