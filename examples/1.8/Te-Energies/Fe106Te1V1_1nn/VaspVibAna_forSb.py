"""
振动分析.
仅仅支持fcc 3x3x3 超胞共108Fe，掺杂溶质solute和空位vacancy时两种振动频率的计算。
运行在Fe106Te1V1 或 Fe107Te1 文件夹下
"""

import os, sys
import subprocess
import zzdlib
import argparse
import numpy as np

python = sys.executable
current_py_folder = os.path.dirname(os.path.abspath(__file__))
VaspCZ_path = [os.path.dirname(current_py_folder) if 'sourcecode' in current_py_folder else current_py_folder][0] + '/sourcecode'

def modify_POSCAR(data, indexes):
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
			tmp_i = i - direct_index - 1
			if tmp_i not in indexes:
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
	"""
	默认溶质位置为0.33 0.50 0.50  空位位置为：0.50 0.50 0.66
	需要进行振动分析的位置为：
		sol 0.33 0.50 0.50
		1nn 0.33 0.33 0.66
		2nn 0.66 0.50 0.50
		3nn 0.66 0.33 0.66
		4nn 0.66 0.50 0.83
	"""
	# 读取文件ini/Opt和fin/Opt下的POSCAR文件，最终获取扩散原子的索引
	with open('./POSCAR', 'r') as f:
		data = f.readlines()
	pos_result = zzdlib.Vasp.decode_POSCAR(data)  # 结果是：vector, elements, number_of_atom, position
	pos_position = pos_result[3]

	# 获取所有索引，TTT 其他FFF
	valuable_position = [
		('0.3333', '0.5000', '0.5000'), ('0.3333', '0.3333', '0.6666'), ('0.6666', '0.5000', '0.5000'),
		('0.6666', '0.3333', '0.6666'), ('0.6666', '0.5000', '0.8333')]
	indexes = []
	for i in range(len(pos_position)):
		x,y,z = pos_position[i].tolist()
		x = f'{x:<.4f}'
		y = f'{y:<.4f}'
		z = f'{z:<.4f}'
		xyz = (x, y, z)
		if xyz in valuable_position:
			# print(xyz, i)
			indexes.append(i)
	print(indexes)

	# 把ini的CONTCAR读取为POSCAR分析振动
	with open('./CONTCAR', 'r') as f:
		data = f.readlines()
	POSCAR_for_vib = modify_POSCAR(data, indexes)

	# 获取任务名
	jobname = f'fc{os.path.basename(os.getcwd())[5::]}_vib'
	# 开始振动分析
	vib_analysis('vib_analysis', POSCAR_data = POSCAR_for_vib, jobname=jobname)

def vib_analysis(folder, POSCAR_data, jobname):
	if not os.path.isdir(folder):
		os.mkdir(folder)
	os.chdir(folder)
	os.system(f'cp ../INCAR .')
	os.system(f'cp ../POTCAR .')
	os.system(f'cp ../KPOINTS .')
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

	global nodes
	global ppn
	zzdlib.Vasp.modify_Vasp_sh(jobname=jobname, nodes=nodes, ppn=ppn)
	zzdlib.Vasp.check_and_qsub()

	os.chdir('..')



if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='manual to this script')
	parser.add_argument('-nd', '--nodes', type=str, default='1')
	parser.add_argument('-np', '--ppn', type=str, default='8')
	args = parser.parse_args()
	nodes = args.nodes
	ppn = args.ppn
	run()
