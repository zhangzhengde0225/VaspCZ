import os
import subprocess
import sys
sys.path.append(os.path.join(os.environ['HOME'], 'bin/VaspCZ/sourcecode'))
import zzdlib
import numpy as np

import VaspCZ

pys_path = os.path.join(os.environ['HOME'], 'bin/VaspCZ/sourcecode')


def run_Hf_and_Hb(ele):
	# 4个文件夹下的结构优化的WARNING 检查和能量检查。
	E_pure = zzdlib.Vasp.check_WARNING_and_Energy(path=f'{ele}-Energies/Fe108')
	E_ele1 = zzdlib.Vasp.check_WARNING_and_Energy(path=f'{ele}-Energies/Fe107{ele}1')
	E_V1 = zzdlib.Vasp.check_WARNING_and_Energy(path=f'{ele}-Energies/Fe107V1')

	xnn_list = '1nn,2nn,3nn,4nn,5nn,6nn'.split(',')
	for xnn in xnn_list:
		E_ele1V1 = zzdlib.Vasp.check_WARNING_and_Energy(path=f'{ele}-Energies/Fe106{ele}1V1_{xnn}')

		Hf = E_V1 - 107/108*E_pure
		Hb = E_ele1V1 + E_pure - E_ele1 - E_V1
		print(f'xnn: {xnn} Hf: {Hf:.4f}  Hb: {Hb:.4f}')


if __name__ == '__main__':
	os.system(f'source ~/.bashrc')
	ele = os.path.basename(os.getcwd()).split('with')[-1]
	print(f'ele: {ele}')
	run_Hf_and_Hb(ele)