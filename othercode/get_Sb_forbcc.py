import os
import subprocess
import sys
sys.path.append(os.path.join(os.environ['HOME'], 'bin/VaspCZ/sourcecode'))
import zzdlib
import numpy as np


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


def run(ele, correct=False):
	code1 = f'grep THz {ele}-Energies/{ele}-V/52Fe1{ele}1V/vib_analysis/*/OUTCAR'
	result1 = zzdlib.getshellResult(code1)
	fre1 = decode_vib_frequency(result1)

	code2 = f'grep THz {ele}-Energies/{ele}-V/53Fe1{ele}/vib_analysis/*/OUTCAR'
	result2 = zzdlib.getshellResult(code2)
	fre2 = decode_vib_frequency(result2)

	print(fre1.shape, fre2.shape)

	# 主要是因为bcc_withCr 中 52Fe1Cr1V的atom46原子振动计算没有完成，fre1缺少一个原子的振动数据，用其他的原子的振动的平均值补充上去
	if correct:
		number = fre1.shape[0]
		row = int(number/3)
		print(f'总{number} 行{row}')
		fre1 = fre1.reshape(row, 3)
		fre1_mean = np.mean(fre1, axis=0)
		fre1_mean = fre1_mean.reshape(1, 3)
		fre1 = np.concatenate((fre1, fre1_mean), axis=0)
		fre1 = fre1.reshape(number+3)

	vac_1nn = np.array([6.6160, 6.6107, 5.5804])
	vac_2nn = np.array([7.0517, 7.0461, 6.6759])
	item1 = np.power(np.prod(vac_1nn), 8) * np.power(np.prod(vac_2nn), 6)

	pure = np.array([7.0748, 7.0727, 7.0708])
	item2 = np.power(np.prod(pure), 15)

	item4 = np.prod(fre1)
	item3 = np.prod(fre2)

	inner = item1 * item3 / item2 / item4
	Sb = np.log(inner)

	print(f'item1: {item1}\nitem2: {item2} \nitem3: {item3} \nitem4: {item4}\n inner: {inner} Sb: {Sb}')


def run2(ele):
	atoms = [53, 16, 12, 25, 14, 17, 23, 26, 39, 40, 37, 42, 31, 48, 27, 36, 30, 28, 13, 11, 5]
	fre1 = np.array([])
	for i in range(len(atoms)):
		atom = str(atoms[i])
		code = f'grep THz {ele}-Energies/{ele}-V/52Fe1{ele}1V/vib_analysis/atom{atom}/OUTCAR'
		result = zzdlib.getshellResult(code)
		fre = decode_vib_frequency(result)
		fre1 = np.concatenate((fre1, fre), axis=0)
	print(fre1.shape)

	vac_1nn = np.array([6.6160, 6.6107, 5.5804])
	vac_2nn = np.array([7.0517, 7.0461, 6.6759])
	item1 = np.power(np.prod(vac_1nn), 8) * np.power(np.prod(vac_2nn), 6)

	pure = np.array([7.0748, 7.0727, 7.0708])
	item2 = np.power(np.prod(pure), 8)

	sol = np.array([6.4545, 6.4515, 6.4482])
	sol_1nn = np.array([7.0020, 6.9931, 6.8190])
	sol_2nn = np.array([7.1250, 7.0793, 7.0716])
	item3 = np.prod(sol) * np.power(np.prod(sol_1nn), 8) * np.power(np.prod(sol_2nn), 6)

	item4 = np.prod(fre1)

	inner = item1 * item3 / item2 / item4
	Sb = np.log(inner)

	print(f'item1: {item1}\nitem2: {item2} \nitem3: {item3} \nitem4: {item4}\n inner: {inner} Sb: {Sb}')



if __name__ == '__main__':
	ele = 'Cr'
	ipt = input(f'选择功能(默认1)：1.所有原子  2.部分原子')
	if ipt == '1' or ipt == '':
		run(ele, correct=True)
	else:
		run2(ele)


