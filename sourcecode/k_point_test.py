#!/home/zhangzhengde/bin/bin/python3
#coding=utf-8

import os
import argparse

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


def run(jobname, nodes, ppn, K):
	if os.path.isdir(K):  # 有目录什么也不做
		print(f'k_mesh:{K} already exists, do nothing.')
		pass
	else:
		os.system('mkdir '+K)  # 创建目录
		for file in os.listdir():
			if os.path.isfile(file):
				os.system(f'cp {file} {K}')# 拷贝输入文件
		os.chdir(K)  # 进入创建的目录
		# 无需修改INCAT
		# 无需修改POTCAR
		# 无需修改POSCAR
		# 修改KPOINTS
		with open('./KPOINTS', 'r') as f:
			data = f.readlines()
		data[3] = f'{K[0]} {K[1]} {K[2]}\n'
		with open('./KPOINTS', 'w') as f:
			f.writelines(data)
		# 修改Vasp.sh，指定任务和任务名，修改，提交任务
		modify_vasp_sh(f'{jobname}_{K}', nodes, ppn)
		# 测试代码，打印
		#os.system('cat KPOINTS')
		#os.system('cat Vasp.sh')
		os.system('qsub Vasp.sh')  # 提交任务
		os.chdir('..')


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('-jb', '--jobname_prefix', default='k_test', type=str)
	parser.add_argument('-nd', '--nodes', default='1', type=str)
	parser.add_argument('-np', '--ppn', default='8', type=str)
	parser.add_argument('-k', '--k_mesh', default='111,333,555,777,999', type=str)
	args = parser.parse_args()
	jobname = args.jobname_prefix
	nodes = args.nodes
	k_mesh = args.k_mesh.split(',')
	ppn = args.ppn
	print(f'running k_point test \n parameter: \njobname_prefix:{jobname} nodes:{nodes} ppn:{ppn} \nk_mesh:{k_mesh}')

	inp = input('confirm run (yes/no) (default:yes): ')
	if inp in ['', 'y', 'yes', 'Y', 'Yes', 'YES']:
		for K in k_mesh:
			run(jobname, nodes, ppn, K)
	else:
		print('did not run.')