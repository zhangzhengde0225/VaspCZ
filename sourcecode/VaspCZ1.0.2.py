#!/home/zhangzhengde/bin/bin/python3
#coding=utf-8
"""
这是应用程序Linux用户界面。
"""

import os, sys
import numpy as np
import time
import VaspGenerate
import VaspCZ.zzdlib as zzd
import utils
import subprocess

python = sys.executable
current_py_folder = os.path.dirname(os.path.abspath(__file__))
VaspCZ_path = [
	os.path.dirname(current_py_folder) if 'sourcecode' in current_py_folder else current_py_folder][0] + '/sourcecode'

__version__ = '1.0.2'

def run():
	# print(VaspCZ_path)
	i = 0
	while True:
		content = {0: 'Exit', 1: 'Opt and Sta module', 2: 'NEB module', 3: 'Test module'}
		ipt = input(utils.gui_string(
			title='VaspCZ interface',
			content=content,
			footnote=f'by: Zhengde Zhang (drivener@163.com) version: {__version__}'))
		try:
			ipt = content[int(ipt)]
		except Exception as e:
			print(f'VaspCZ 功能{e}选择错误，请正确输入.')
			pass
		# print(ipt)
		if ipt == 'Opt and Sta module':
			print(f'OS module selected')
			content_os = utils.zip_content([
				'Back', 'Generate inputs (example)', 'Generate INCAR for Sta', 'Generate POTCAR',
				'Generate KPOINTS', 'Generate Vasp.sh',
				'Vasp Keep Inputs', 'Vasp Pre-check and Qsub', 'Check Results'])
			while True:
				ipt1 = input(utils.gui_string(
					'Optimization and Static calculation', content=content_os))
				try:
					ipt1 = content_os[int(ipt1)]
				except Exception as e:
					print(f'OS module 功能{e}选择错误，请正确输入.')
				if ipt1 == 'Generate inputs (example)':
					VaspGenerate.generate_inputs(examp='fcc_Fe_3x3x3')
					exit()
				elif ipt1 == 'Generate INCAR for Sta':
					VaspGenerate.generate_INCAR_for_Sta()
					exit()
				elif ipt1 == 'Generate POTCAR':
					utils.deal_with_gen_pot()
					exit()
				elif ipt1 == 'Generate KPOINTS':
					utils.deal_with_gen_kpoints()
					exit()
				elif ipt1 == 'Generate Vasp.sh':
					utils.deal_with_gen_vasp_sh()
					exit()
				elif ipt1 == 'Vasp Keep Inputs':
					utils.deal_with_vasp_keep_inputs()
					exit()
				elif ipt1 == 'Vasp Pre-check and Qsub':
					zzd.Vasp.check_and_qsub(need_input=True)
					exit()
				elif ipt1 == 'Check Results':
					utils.deal_with_check_results()
					exit()
				elif ipt1 == 'Back':
					break
		elif ipt == 'NEB module':
			print(f'NEB module selected')
			content_neb = utils.zip_content([
				'Back', 'NEB Opt-Sta', 'NEB Sta-NEB', 'NEB Vibration Analysis', 'NEB Keep INFI/Opt Inputs', 'NEB Keep Inputs',
				'NEB Check RMS', 'NEB Check Dist', 'NEB Check Results', 'NEB Check Vibration Results'
			])
			while True:
				ipt2 = input(utils.gui_string(
					title='NEB calculation', content=content_neb
				))
				try:
					ipt2 = content_neb[int(ipt2)]
				except Exception as e:
					print(f'NEB module 功能{e}选择错误，请正确输入.')
				# print(ipt2)
				if ipt2 == 'Back':
					break
				elif ipt2 == 'NEB Opt-Sta':
					utils.deal_with_neb_opt_sta()
				elif ipt2 == 'NEB Sta-NEB':
					utils.deal_with_neb_sta_neb()
					exit()
				elif ipt2 == 'NEB Vibration Analysis':
					utils.deal_with_neb_vibration_analysis()
					exit()
				elif ipt2 == 'NEB Keep INFI/Opt Inputs':
					subprocess.call(f'{python} {VaspCZ_path}/VaspNEBKeepINFI_OptInputs.py', shell=True)
					exit()
				elif ipt2 == 'NEB Keep Inputs':
					subprocess.call(f'{python} {VaspCZ_path}/VaspNEBKeepInputs.py', shell=True)
					exit()
				elif ipt2 == 'NEB Check RMS':
					subprocess.call(f'{python} {VaspCZ_path}/VaspNEBCheckRMS.py', shell=True)
					exit()
				elif ipt2 == 'NEB Check Dist':
					ipt2_nebcd = input('Check ([POS]/CONT):  ')
					pos_or_cont = 'POSCAR' if (ipt2_nebcd == 'POS' or ipt2_nebcd == '') else 'CONTCAR'
					subprocess.call(f'{python} {VaspCZ_path}/VaspNEBCheckDist.py --POSorCONT={pos_or_cont}', shell=True)
					exit()
				elif ipt2 == 'NEB Check Results':
					utils.deal_with_neb_check_results()
				elif ipt2 == 'NEB Check Vibration Results':
					code = [sys.executable, f'{VaspCZ_path}/VaspVibAna_Result.py', '--isprint=True']
					subprocess.call(code, shell=False)
					exit()
				else:
					pass
		elif ipt == 'Test module':
			print(f'Tese module selected')
			content_test = utils.zip_content([
				'Back',
				'ENCUT Test',
				'KPOINTS Mesh Test'
			])
			while True:
				ipt3 = input(utils.gui_string(title='Vasp Test Module', content=content_test))
				try:
					ipt3 = content_test[int(ipt3)]
				except Exception as e:
					print(f'Test module 功能{e}选择错误，请正确输入.')
				if ipt3 == 'Back':
					break
				elif ipt3 == 'ENCUT Test':
					utils.deal_with_test_encut()
					exit()
				elif ipt3 == 'KPOINTS Mesh Test':
					utils.deal_with_test_kpoints()
					exit()
				else:

					pass
		elif ipt == 'Exit':
			break
		else:
			pass



if __name__ == '__main__':
	run()