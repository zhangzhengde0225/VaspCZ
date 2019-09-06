import os
import subprocess
import sys
sys.path.append(os.path.join(os.environ['HOME'], 'bin/VaspCZ/sourcecode'))
import zzdlib
import numpy as np

import VaspCZ

pys_path = os.path.join(os.environ['HOME'], 'bin/VaspCZ/sourcecode')


def decode_vib_frequency(data):
	# data 是shell代码出来的结果，是一个列表
	frequncies = []
	for i in range(len(data)):
		line = data[i]
		if 'THz' in line:
			frequncy = line.split('=')[-1].split()[0]
			frequncy = float(f'{float(frequncy):.4f}')  # 保留4位有效数字
			frequncies.append(frequncy)
	frequncies_np = np.array(frequncies)
	return frequncies_np


def grep_fre(path):
	code = f'grep THz {path}/OUTCAR'
	result = zzdlib.getshellResult(code)
	fre = decode_vib_frequency(result)
	return fre


def run_Sb(ele):
	vac_1nn = grep_fre(f'{ele}-Energies/Fe107V1/vib_analysis/1nn')
	vac_2nn = grep_fre(f'{ele}-Energies/Fe107V1/vib_analysis/2nn')
	print(f'vac_1nn: {vac_1nn}  vac_2nn: {vac_2nn}')
	item1 = np.power(np.prod(vac_1nn), 12) * np.power(np.prod(vac_2nn), 6)

	pure = grep_fre(f'{ele}-Energies/Fe108/vib_analysis')
	print(f'pure: {pure}')
	item2 = np.power(np.prod(pure), 10)  # 几次方

	sol = grep_fre(f'{ele}-Energies/Fe107{ele}1/vib_analysis/335050')  # 溶质位置335050
	sol_1nn = grep_fre(f'{ele}-Energies/Fe107{ele}1/vib_analysis/505066')
	sol_2nn = grep_fre(f'{ele}-Energies/Fe107{ele}1/vib_analysis/665050')
	print(f'sol: {sol}  sol_1nn: {sol_1nn}  sol_2nn: {sol_2nn}')
	item3 = np.prod(sol) * np.power(np.prod(sol_1nn), 12) * np.power(np.prod(sol_2nn), 6)

	vs = grep_fre(f'{ele}-Energies/Fe106{ele}1V1/vib_analysis/*')
	print(vs.shape)
	item4 = np.prod(vs)

	inner = item1 * item3 / item2 / item4
	Sb = np.log(inner)

	print(f'item1: {item1}\nitem2: {item2} \nitem3: {item3} \nitem4: {item4}\ninner: {inner:.4f} Sb: {Sb:.4f}')


def run_Hf_and_Hb(ele):
	# 4个文件夹下的结构优化的WARNING 检查和能量检查。
	E_ele1V1 = zzdlib.Vasp.check_WARNING_and_Energy(path=f'{ele}-Energies/Fe106{ele}1V1')
	E_pure = zzdlib.Vasp.check_WARNING_and_Energy(path=f'{ele}-Energies/Fe108')
	E_ele1 = zzdlib.Vasp.check_WARNING_and_Energy(path=f'{ele}-Energies/Fe107{ele}1')
	E_V1 = zzdlib.Vasp.check_WARNING_and_Energy(path=f'{ele}-Energies/Fe107V1')

	Hf = E_V1 - 107/108*E_pure
	Hb = E_ele1V1 + E_pure - E_ele1 - E_V1
	print(f'Hf: {Hf:.4f}  Hb: {Hb:.4f}')

def run_Sf(ele):
	vac_1nn = grep_fre(f'{ele}-Energies/Fe107V1/vib_analysis/1nn')
	vac_2nn = grep_fre(f'{ele}-Energies/Fe107V1/vib_analysis/2nn')
	# print(f'vac_1nn: {vac_1nn}  vac_2nn: {vac_2nn}')
	item1 = np.power(np.prod(vac_1nn), 12) * np.power(np.prod(vac_2nn), 6)

	pure = grep_fre(f'{ele}-Energies/Fe108/vib_analysis')
	# print(f'pure: {pure}')
	item2 = np.power(np.prod(pure), 18)  # 几次方

	inner = item2/item1
	Sf = np.log(inner)
	print(f'inner: {inner:.4f} Sf: {Sf:.4f}')


def get_w(Hm, v, T=1000):
	"""
	w = v*exp(-Hm/(kB*T))
	:return:
	"""
	# 假定T=1000
	kB = 1.380649 * 10 ** -23
	Hm_unit_j = Hm * 1.602 * 10 ** -19
	w = v * 10**12 * np.exp(-Hm_unit_j/(kB*T))
	return w

def get_f2(Hm_and_v, T=1000):
	"""
	f2 = (2*w1 + 7*w3*F)/(2*w1 + 2*w2+ 7*w3*F)
	:param Hm_and_v:
	:return:
	"""
	data = Hm_and_v
	w1 = get_w(Hm_and_v[0,2], Hm_and_v[1,2], T=T)
	w3_2nn = get_w(data[0,6], data[1,6], T=T)
	w3_3nn = get_w(data[0,8], data[1,8], T=T)
	w3_4nn = get_w(data[0,10], data[1,10], T=T)
	w3 = (2*w3_2nn + 4*w3_3nn + w3_4nn)/7
	w2 = get_w(data[0,4], data[1,4], T=T)

	# 获取zeta
	w0 = get_w(data[0,0], data[1,0], T=T)
	w4_2nn = get_w(data[0,7], data[1,7], T=T)
	w4_3nn = get_w(data[0,9], data[1,9], T=T)
	w4_4nn = get_w(data[0,11], data[1,11], T=T)
	w4 = (2*w4_2nn + 4*w4_3nn + w4_4nn)/7
	# print(w4, w0)
	zeta = w4/w0
	fenzi = (10*np.power(zeta,4) + 180.5*np.power(zeta,3) + 927*np.power(zeta,2) + 1341*zeta)
	fenmu = 7*(2*np.power(zeta,4) + 40.2*np.power(zeta,3) + 254*np.power(zeta,2) + 597*zeta + 436)
	F = 1 - fenzi/fenmu
	print(f'zeta: {zeta}')
	print(f'F: {F}')
	f2 = (2* w1 + 7 * w3 * F) / (2 * w1 + 2 * w2 + 7 * w3 * F)
	return f2

def get_Hm_and_v_from_raw():
	cases = 'w0_3nn,w1,w2,w3_2nn,w3_3nn,w3_4nn'.split(',')
	Hm_list = []
	v_list = []
	for i in range(len(cases)):
		case = cases[i]
		os.chdir(case)

		# 获取Hm
		# print(pys_path, os.getcwd())
		res = zzdlib.getshellResult(f'python3 {pys_path}/NEBCheck.py --func=2')
		# print(res)
		line, index = zzdlib.File.getLine(res, 'IMAGE')
		IMAGE = []
		Barrier = []
		# print(res)
		for j in range(index + 1, len(res)):
			if res[j] == '\n':
				break
			img = int(res[j].split()[0])
			barrier = float(res[j].split()[-1])
			IMAGE.append(img)
			Barrier.append(barrier)
		Hm_foreward = np.max(Barrier)
		Hm_backward = np.max(Barrier) - Barrier[-1]
		Hm_list.append(Hm_foreward)
		Hm_list.append(Hm_backward)

		# 获取v
		code = f'python3 {pys_path}/VaspVibAna_Result.py --isprint=True'
		res = zzdlib.getshellResult(code)
		v_foreward = zzdlib.File.getLine(res, keywords='foreward')[0].split()[-1]
		v_backward = zzdlib.File.getLine(res, keywords='backward')[0].split()[-1]
		if v_backward == 'Match':
			v_backward = '0'
		v_list.append(float(v_foreward))
		v_list.append(float(v_backward))
		os.chdir('..')

	Hm_and_v = np.concatenate((np.array(Hm_list).reshape(1, -1), np.array(v_list).reshape(1, -1)), axis=0)
	return Hm_and_v


def run_f2(ele, T=1000):
	# 先获取Hm和v的数据
	os.chdir(f'{ele}-V')

	if os.path.isfile('.Hm_and_v'):
		with open('.Hm_and_v') as f:
			data = f.read()
			Hm_and_v = np.array(eval(data))

	else:
		Hm_and_v = get_Hm_and_v_from_raw()
		with open('.Hm_and_v', 'w') as f:
			f.write(str(Hm_and_v.tolist()))

	os.chdir('..')

	# 算f0
	Hm_self = [1.3943 for i in range(12)]
	v_self = [6.2805 for i in range(12)]

	Hm_and_v_self = np.concatenate((np.array(Hm_self).reshape(1, -1), np.array(v_self).reshape(1, -1)), axis=0)
	# print(Hm_and_v_self)
	f0 = get_f2(Hm_and_v_self, T=T)
	print(f'f0: {f0}')

	# 计算
	f2 = get_f2(Hm_and_v, T=T)
	print(f'f2: {f2}')
	return f2


if __name__ == '__main__':
	os.system(f'source ~/.bashrc')
	ele = os.path.basename(os.getcwd()).split('with')[-1]
	print(f'ele: {ele}')

	run_Sb(ele)
	run_Sf(ele)
	run_Hf_and_Hb(ele)
	run_f2(ele, T=1000)



