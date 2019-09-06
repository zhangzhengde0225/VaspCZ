"""
运行在fcc_withele目录下，进入ele-Energies/，获取Fe106ele1V1文件夹下的输入文件，把POSCAR改为对应的2nn-6nn，然后进行计算
"""

import os
import subprocess
import sys
sys.path.append(os.path.join(os.environ['HOME'], 'bin/VaspCZ/sourcecode'))
import zzdlib
import numpy as np
import argparse

import VaspCZ

pys_path = os.path.join(os.environ['HOME'], 'bin/VaspCZ/sourcecode')

def run(nodes, ppn, ele):
	os.chdir(f'{ele}-Energies')
	caselist = '1nn,2nn,3nn,4nn,5nn,6nn'.split(',')
	vaclist = '505066,665050,666633,668350,835066,668316'.split(',')
	for i in range(len(caselist)):
		case = caselist[i]
		vac = vaclist[i]
		dirname = f'Fe106{ele}1V1_{case}'
		if os.path.isdir(dirname):
			print(f'{dirname} 已存在，退出')
			exit()
		os.mkdir(dirname)
		os.chdir(dirname)
		if case == '1nn':
			os.system(f'cp ../Fe106{ele}1V1/* .')
		else:
			os.system(f'cp ../Fe106{ele}1V1/INCAR .')
			os.system(f'cp ../Fe106{ele}1V1/POTCAR .')
			os.system(f'cp ../Fe106{ele}1V1/KPOINTS .')
			# 获取和修改POSCAR
			os.system(f'cp ../Fe106{ele}1V1/POSCAR .')
			pos_data = zzdlib.File.openFile('./POSCAR', 'r')
			new_pos_data = []
			for j in range(len(pos_data)):
				line = pos_data[j]
				if j < 8:  # 前8行是其他信息，到direct
					new_pos_data.append(line)
				else:
					try:
						x, y, z = line.strip('\n').split()
					except Exception as e:
						print(j, line)
						x, y, z = '0', '0', '0'
					vacx = f'0.{vac[:2]}'
					vacy = f'0.{vac[2:4]}'
					vacz = f'0.{vac[4::]}'
					if x[:4] == vacx and y[:4] == vacy and z[:4] == vacz:
						thatline = f'0.500000000 0.500000000 0.666666666\n'
						new_pos_data.append(thatline)
					else:
						new_pos_data.append(line)
			zzdlib.File.openFile('./POSCAR', 'w', data=new_pos_data)

			# Vasp.sh
			jobname = f'fc06{ele}{case}Opt'
			Vaspsh_path = zzdlib.File.Vaspsh_path()
			os.system(f'cp {Vaspsh_path}/Vasp.sh .')
			zzdlib.Vasp.modify_Vasp_sh(jobname, nodes=nodes, ppn=ppn)
			zzdlib.Vasp.check_and_qsub(need_input=True)
		os.chdir('..')
	os.chdir('..')


def check_and_reqsub():
	elelist = 'Al,Cr,Cu,Mo,Nb,Ru,Sb,Tc'.split(',')
	for ele in elelist:
		first_dir = f'fcc_with{ele}/{ele}-Energies'
		os.chdir(first_dir)
		xnnlist = '1nn,2nn,3nn,4nn,5nn,6nn'.split(',')
		for xnn in xnnlist:
			second_dir = f'Fe106{ele}1V1_{xnn}'
			os.chdir(second_dir)
			data_log = zzdlib.File.openFile('log', 'r')
			res = zzdlib.File.getLine(data=data_log, keywords='node8')
			if res[0] != 'Not Match':  # 匹配到Node8了
				zzdlib.Vasp.check_and_qsub(need_input=True)
			else:
				print(f'fcc_with{ele} {xnn} 未检测到node8, 跳过')
			os.chdir('..')
		os.chdir('../..')


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='manual to this script')
	parser.add_argument('-nd', '--nodes', type=str, default='1')
	parser.add_argument('-np', '--ppn', type=str, default='12')
	parser.add_argument('-func', '--function', type=str, default='cal')
	args = parser.parse_args()
	nodes = args.nodes
	ppn = args.ppn
	func = args.function
	ele = os.path.basename(os.getcwd()).split('with')[-1]
	print(f'输入参数： nodes:{nodes} ppn:{ppn} ele:{ele} func:{func}')
	if func == 'cal':
		run(nodes, ppn, ele)
	else:
		check_and_reqsub()


