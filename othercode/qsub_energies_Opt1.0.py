"""
在ele-Energies文件夹下运行
1.0
	进入两个文件夹，提交结构优化。[Fe106Cr1v1, Fe107Cr1]
	Cr是old-ele, 替换成新的ele
"""

import os
import argparse
import sys
sys.path.append(os.path.join(os.environ['HOME'], 'bin/VaspCZ/sourcecode'))
import zzdlib


def run(nodes, ppn, ele, oele):
	dirs = [f'Fe106{oele}1V1', f'Fe107{oele}1']
	for dir in dirs:
		new_dir = dir.replace(oele, ele)
		os.system(f'mv {dir} {new_dir}')
		os.chdir(new_dir)

		# INCAR无需修改
		# POSCAR需要修改
		zzdlib.Vasp.modify_POSCAR_ele(oele, ele)

		# POTCAR 需要修改
		zzdlib.Vasp.gennerate_POTCAR()

		# KPOINT无需修改
		# Vasp.sh 修改nodes 和ppn
		jobname = f'fc{new_dir[5::]}_Opt'
		zzdlib.Vasp.modify_Vasp_sh(jobname, nodes, ppn)
		zzdlib.Vasp.check_and_qsub()
		os.chdir('..')


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('-nd', '--nodes', default='1', type=str)
	parser.add_argument('-np', '--ppn', default='8', type=str)
	parser.add_argument('-ele', '--element', default='Cr', type=str)
	parser.add_argument('-oele', '--old_element', default='Te', type=str)
	args = parser.parse_args()
	nodes = args.nodes
	ppn = args.ppn
	ele = args.element
	oele = args.old_element
	run(nodes, ppn, ele, oele)