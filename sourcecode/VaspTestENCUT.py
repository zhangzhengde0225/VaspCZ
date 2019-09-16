#!/home/zhangzhengde/bin/bin/python3
# coding=utf-8

import os
import argparse
import VaspCZ.zzdlib as zzd


def modify_vasp_sh(jobname, nodes, ppn):
	with open('./Vasp.sh', 'r') as f:
		data = f.readlines()
	new_data = []
	for line in data:
		if ' #PBS -N' in line:
			new_data.append(f' #PBS -N {jobname}\n')
		elif ' #PBS -l nodes' in line:
			new_data.append(f' #PBS -l nodes={nodes}:ppn={ppn}\n')
		else:
			new_data.append(line)
	with open('./Vasp.sh', 'w') as f:
		f.writelines(new_data)


def run(jobname, nodes, ppn, encut):
	input_files = 'INCAR,POSCAR,POTCAR,KPOINTS'.split(',')
	for i in input_files:
		if i not in os.listdir():
			raise NameError(f'ENCUT Test: input file "{i}" missing in current dir.')
	if os.path.isdir(encut):  # 有目录什么也不做
		print(f'ENCUT:{encut} already exists, do nothing.')
		pass
	else:
		os.system('mkdir '+encut)  # 创建目录
		for file in input_files:
			if os.path.isfile(file):
				os.system(f'cp {file} {encut}')# 拷贝输入文件
		os.chdir(encut)  # 进入创建的目录
		vasp_sh_path = zzd.File.Vaspsh_path()
		os.system(f'cp {vasp_sh_path}/Vasp.sh .')
		# 需修改INCAR
		data_INCAR = zzd.File.openFile('INCAR', 'r')
		data_new = zzd.File.substituteData(data_INCAR, keywords='ENCUT', newline=f'ENCUT={encut}\n')
		zzd.File.openFile('INCAR', 'w', data=data_new)
		# 无需修改POTCAR
		# 无需修改POSCAR
		# 无需修改KPOINTS
		# 修改Vasp.sh，指定任务和任务名，修改，提交任务
		modify_vasp_sh(f'{jobname}_{encut}', nodes, ppn)
		# 测试代码，打印
		#os.system('cat KPOINTS')
		#os.system('cat Vasp.sh')
		zzd.Vasp.check_and_qsub()
		os.chdir('..')


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('-jb', '--jobname_prefix', default='k_test', type=str)
	parser.add_argument('-nd', '--nodes', default='1', type=str)
	parser.add_argument('-np', '--ppn', default='8', type=str)
	parser.add_argument('-EN', '--ENCUTs', default='200,250,300,350,400,450,500,550,600,650,700', type=str)
	args = parser.parse_args()
	jobname = args.jobname_prefix
	nodes = args.nodes
	ENCUTs = args.ENCUTs.split(',')
	ppn = args.ppn
	print(f'running k_point test \n parameter: \njobname_prefix:{jobname} nodes:{nodes} ppn:{ppn} \nENCUTs:{ENCUTs}')

	inp = input('confirm run ([y]es/no):  ')
	if inp in ['', 'y', 'yes', 'Y', 'Yes', 'YES']:
		for encut in ENCUTs:
			run(jobname, nodes, ppn, encut)
	else:
		print('Did not run.')