"""
之前所有算的bcc_withele，需要补全一些计算，运行在bcc_withele文件夹下：
补全内容：
	1.ele-V下的1nn 2nn Fe_12 Fe_13 Fe_14下的振动分析
	2.ele-V下的w5任务的提交。
	3.ele-Energies下的 53Fe1Al的振动分析，52Fe1Al1V的振动分析
"""

import os
import argparse
import sys
sys.path.append(os.path.join(os.environ['HOME'], 'bin/VaspCZ/sourcecode'))
import zzdlib
import subprocess

pys_path = os.path.join(os.environ['HOME'], 'bin/VaspCZ/sourcecode')


def ele_V_vib(ele, nodes, ppn):
	print(f'执行ele-V 1nn 2nn Fe_12 Fe_13 Fe_14 振动分析')
	os.chdir(f'{ele}-V')
	cases = ['1nn', '2nn', 'Fe_12', 'Fe_13', 'Fe_14']
	include_fin = [False, False, True, True, True]
	for i in range(len(cases)):
		case = cases[i]
		os.chdir(case)
		code = [
			sys.executable, f'{pys_path}/VaspVibAna_forNEB.py', f'-nd={nodes}', f'-np={ppn}',
			f'--include_fin={include_fin[i]}']
		print(code)
		subprocess.call(code, shell=False)
		os.chdir('..')
	os.chdir('..')


def ele_V_w5(ele, nodes, ppn, w5_path):
	print(f'执行ele-V w5 ini fin 的结构优化')
	os.chdir(f'{ele}-V')
	os.mkdir('w5')
	os.chdir('w5')

	w5_path = os.path.join(os.environ['HOME'], w5_path)
	os.system(f'cp -rf {w5_path}/* .')

	infi_list = ['ini', 'fin']
	for i in range(len(infi_list)):
		infi = infi_list[i]
		os.chdir(f'{infi}/Opt')
		# INCAR 无需修改
		# KPOINT 无序修改
		# POSCAR 修改元素
		zzdlib.Vasp.modify_POSCAR_ele(oldele='Te', new_ele=ele)
		# POTCAR 匹配
		zzdlib.Vasp.gennerate_POTCAR()  # 从当前目录下的POSCAR 读取元素并生成PBE贋势
		# Vasp.sh 匹配平台，修改名称
		Vaspsh_path = zzdlib.File.Vaspsh_path()
		os.system(f'cp {Vaspsh_path}/Vasp.sh .')
		jobname = f'bc{ele}w5{infi[0]}O'
		zzdlib.Vasp.modify_Vasp_sh(jobname, nodes=nodes, ppn=ppn)
		# 检查和提交
		zzdlib.Vasp.check_and_qsub()
		os.chdir('../..')
	os.chdir('../..')

def ele_Energies_vib(ele, nodes, ppn):
	print(f'执行ele-Energies Fe531{ele} Fe52{ele}1V1的振动分析')
	os.chdir(f'{ele}-Energies/{ele}-V')

	# cases = [f'53Fe1{ele}', f'52Fe1{ele}1V']
	cases = [f'53Fe1{ele}']
	for i in range(len(cases)):
		case = cases[i]
		os.chdir(case)
		if os.path.isdir('vib_analysis'):
			print(f'{ele} {case} vib_analysis已存在，退出')
			exit()
		else:
			os.mkdir('vib_analysis')
		os.chdir('vib_analysis')

		# 获取所有索引，TTT 其他FFF
		solute = ('0.3333', '0.3333', '0.3333')
		nn1 = ('0.5000', '0.5000', '0.5000')
		nn2 = ('0.3333', '0.3333', '0.6666')
		nn3 = ('0.3333', '0.6666', '0.6666')
		nn5 = ('0.6666', '0.6666', '0.6666')
		valuable_position = [
			('0.3333', '0.3333', '0.6666'), ('0.3333', '0.6666', '0.6666'), ('0.6666', '0.6666', '0.6666'),
			('0.5000', '0.5000', '0.5000'), ('0.3333', '0.3333', '0.3333')]
		xnn_list = ['solute', 'nn1', 'nn2', 'nn3', 'nn5']
		if case == f'52Fe1{ele}1V':
			xnn_list.remove('nn1')
		for j in range(len(xnn_list)):
			xnn = eval(xnn_list[j])

			# 创建
			os.mkdir(xnn_list[j])
			os.chdir(xnn_list[j])
			# 拷贝文件
			# INCAR 需要修改
			os.system(f'cp ../../INCAR .')
			zzdlib.Vasp.modify_INCAR_for_vibration_analysis()

			# POSCAR 从上级从是CONTCAR拷贝，需要修改
			os.system(f'cp ../../CONTCAR POSCAR')

			with open('../../POSCAR', 'r') as f:  # 上一级的POSCAR
				data = f.readlines()
			pos_result = zzdlib.Vasp.decode_POSCAR(data)  # 结果是：vector, elements, number_of_atom, position
			pos_position = pos_result[3]

			indexes = []
			for k in range(len(pos_position)):
				x, y, z = pos_position[k].tolist()
				x = f'{x:<.4f}'
				y = f'{y:<.4f}'
				z = f'{z:<.4f}'
				xyz = (x, y, z)
				if xyz == xnn:
					# print(xyz, i)
					indexes.append(k)

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
			jobname = f'bc{case[4::]}{xnn_list[j][-1]}_vib'
			zzdlib.Vasp.modify_Vasp_sh(jobname, nodes=nodes, ppn=ppn)

			zzdlib.Vasp.check_and_qsub()
			os.chdir('..')

		os.chdir('..')
		os.chdir('..')
	os.chdir('../..')

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('-nd', '--nodes', default='1', type=str)
	parser.add_argument('-np', '--ppn', default='8', type=str)
	parser.add_argument('-w5p', '--w5_path', default='zhangzd/bcc_54atom/w5_prepare', type=str)
	args = parser.parse_args()
	nodes = args.nodes
	ppn = args.ppn
	w5_path = args.w5_path
	ele = os.path.basename(os.getcwd()).split('with')[-1]
	ipt = input(f'功能：1.{ele}-Energies_vib  2.{ele}-V_vib   3.{ele}-V_w5_Opt\n请输入需要的功能(默认all): ')
	if ipt in ['1', 'all', 'A', 'a', '']:
		ele_Energies_vib(ele, nodes, ppn)
	if ipt in ['2', 'all', 'A', 'a', '']:
		ele_V_vib(ele, nodes, ppn)
	if ipt in ['3', 'all', 'A', 'a', '']:
		ele_V_w5(ele, nodes, ppn, w5_path)


