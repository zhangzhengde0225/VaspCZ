"""
振动分析.
仅仅支持fcc 3x3x3 超胞共108Fe，掺杂溶质solute和空位vacancy时两种振动频率的计算。
运行在fcc_withele 文件夹下，提交2个振动分析：Fe106Te1V1 或 Fe107Te1

2019.8.7 更新
原来只考虑1NN原子振动的变化，且VASP结果默认不区分各个原子的振动频率
现在考虑1NN和2NN原子振动的变化。
对于Fe107Te1, Te位置为0.33 0.50 0.50， 1NN 12个 代表位置为0.50 0.50 0.66, 2NN 6个 代表位置为 0.66 0.50 0.50， 共3个原子的振动。
对于Fe106Te1V1, 考虑空位的1NN 2NN和溶质的1NN 2NN的所有原子的并集，共有18+10=28个原子的振动。
"""

import os, sys
import subprocess
import zzdlib
import argparse
import numpy as np

python = sys.executable
current_py_folder = os.path.dirname(os.path.abspath(__file__))
VaspCZ_path = [os.path.dirname(current_py_folder) if 'sourcecode' in current_py_folder else current_py_folder][0] + '/sourcecode'

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


def run(nodes, ppn, ele):
	print(f'开始进行fcc_with{ele}振动分析')
	os.chdir(f'{ele}-Energies')


	cases = ['Fe107ele1', 'Fe106ele1V1']
	pos_dict = {
		'Fe107': ['505066', '665050', '335050'],
		'Fe106': [
			'005050', '166650', '165033', '165066', '163350', '336633', '333333', '335016', '338350',
			'331650', '333366', '335050', '335083', '336666', '505000', '505033', '501666', '503350', '503383',
			'506650', '506683', '508366', '663366', '666666', '665050', '665083', '835066']
	}

	for i in range(len(cases)):
		case = cases[i].replace('ele', ele)
		if not os.path.isdir(case):
			print(f'case {case} 不存在，退出程序')
			exit()
		os.chdir(case)
		if os.path.isdir('vib_analysis'):
			print(f'{ele} {case} vib_analysis 文件夹已存在，退出程序')
		else:
			os.mkdir('vib_analysis')
		os.chdir('vib_analysis')

		# 读取文件ini/Opt和fin/Opt下的POSCAR文件，最终获取扩散原子的索引
		POSCAR_data = zzdlib.File.openFile('../POSCAR', 'r')
		pos_result = zzdlib.Vasp.decode_POSCAR(POSCAR_data)  # 结果是：vector, elements, number_of_atom, position
		pos_position = pos_result[3]

		pos_list = pos_dict[case[:5]]
		for j in range(len(pos_list)):
			pos_tmp = pos_list[j]  # pos现在一个字符串
			x_tmp = f'0.{pos_tmp[:2]}{pos_tmp[1]}{pos_tmp[1]}'
			y_tmp = f'0.{pos_tmp[2:4]}{pos_tmp[3]}{pos_tmp[3]}'
			z_tmp = f'0.{pos_tmp[4::]}{pos_tmp[5]}{pos_tmp[5]}'
			pos = (x_tmp, y_tmp, z_tmp)
			indexes = []
			for k in range(len(pos_position)):
				x, y, z = pos_position[k].tolist()
				x = f'{x:<.4f}'
				y = f'{y:<.4f}'
				z = f'{z:<.4f}'
				xyz = (x, y, z)
				if xyz == pos:
					# print(xyz, i)
					indexes.append(k)

			# 每个索引对应要计算振动的一个原子，每一次不同的Poslist就是计算不同的原子的振动
			dir_name = f'{pos[0][2:4]}{pos[1][2:4]}{pos[2][2:4]}'  # 如 505066
			os.mkdir(dir_name)
			os.chdir(dir_name)

			# INCAR 需要修改
			os.system(f'cp ../../INCAR .')
			zzdlib.Vasp.modify_INCAR_for_vibration_analysis()

			# POSCAR 从上级从是CONTCAR拷贝，需要修改
			os.system(f'cp ../../CONTCAR POSCAR')

			data_POS = zzdlib.File.openFile('POSCAR', 'r')
			POSCAR_data = zzdlib.Vasp.modify_POSCAR_Selective_Dynamics(data=data_POS, indexes=indexes)
			zzdlib.File.openFile('POSCAR', 'w', data=POSCAR_data)

			# KPOINTS 无序修改
			os.system(f'cp ../../KPOINTS .')
			# POTCAR 无序修改
			os.system(f'cp ../../POTCAR .')

			# Vasp.sh 需要修改
			Vaspsh_path = zzdlib.File.Vaspsh_path()
			os.system(f'cp {Vaspsh_path}/Vasp.sh .')
			jobname = f'fc{case[3:5]}{ele}{dir_name}V'
			zzdlib.Vasp.modify_Vasp_sh(jobname, nodes=nodes, ppn=ppn)

			zzdlib.Vasp.check_and_qsub(need_input=True)
			os.chdir('..')

		os.chdir('../..')
	os.chdir('..')


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='manual to this script')
	parser.add_argument('-nd', '--nodes', type=str, default='1')
	parser.add_argument('-np', '--ppn', type=str, default='8')
	args = parser.parse_args()
	nodes = args.nodes
	ppn = args.ppn
	ele = os.path.basename(os.getcwd()).split('with')[-1]
	print(f'输入参数： nodes:{nodes} ppn:{ppn} ele:{ele}')
	run(nodes, ppn, ele)
