"""
振动分析结束后，获取振动分析结果，计算出当前态下的振动振动频率。
公式为：v = (求积vi_equilirium)/(求积vi_saddle)
平衡态一般取初态ini的结果，如果有反向过程，如bcc九频模型中w4是w3的反过程，取w3过渡态的fin作为平衡态，即可就出w4的振动频率
"""

import zzdlib
import os
import argparse
import numpy as np


def run_walk(isprint=True):
	run_path = os.getcwd()
	for dirpath, dirnames, filenames in os.walk('./'):
		if 'vib_analysis' in dirnames:
			os.chdir(f'{dirpath}/vib_analysis')
			run(isprint=isprint)
			os.chdir(run_path)

def run(isprint=True):
	# 先检查振动分析都没有完成
	dirnames = os.listdir()
	Donelist = []
	for state in ['ini', 'sad', 'fin']:
		if f'{state}_state' in dirnames:
			data =zzdlib.File.openFile(f'{state}_state/log', 'r')
			Total_list = zzdlib.File.getAllline(data, keywords='Total')
			if int(Total_list[0].split('/')[-1].strip('\n')) == len(Total_list):
				isDone=True
			else:
				isDone=False
		else:
			isDone=False
		Donelist.append(isDone)
	tmp_path = f'{os.getcwd()}'
	print(f'{tmp_path:<40}', Donelist)

	# 计算振动频率
	if Donelist[0] and Donelist[1]:  #
		# print(f'开始计算正向过程振动频率')
		cal('ini_state', 'sad_state', isprint=isprint, title='foreward')
	if Donelist[1] and Donelist[2]:
		# print(f'开始计算反向过程振动频率')
		cal('fin_state', 'sad_state', isprint=isprint, title='backward')

def cal(path1, path2, isprint=True, title=''):
	data1 = zzdlib.File.openFile(f'{path1}/OUTCAR', 'r')
	data2 = zzdlib.File.openFile(f'{path2}/OUTCAR', 'r')
	frequncy1 = zzdlib.File.getAllline(data1, keywords='THz')
	frequncy2 = zzdlib.File.getAllline(data2, keywords='THz')
	if isprint:
		print(f'{os.getcwd()}')
		for line in frequncy1:
			print(line.strip('\n'))
		for line in frequncy2:
			print(line.strip('\n'))
	fre1_value = [fre.split('=')[1].split('THz')[0].strip() for fre in frequncy1]
	fre2_value = [fre.split('=')[1].split('THz')[0].strip() for fre in frequncy2]
	fre1 = np.array(fre1_value).astype(float)
	fre2 = np.array(fre2_value).astype(float)
	fre = np.prod(fre1)/np.prod(fre2[:-1])  # 这一这里去掉了最后一个振动频率，但是不能确定最后一个振动频率就要去掉的振动频率
	print(f'{os.getcwd():<40} {title:8} freqency: {fre:.4f}')



if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='manual to this script')
	parser.add_argument('-p', '--isprint', type=bool, default=True)
	args = parser.parse_args()
	isprint = args.isprint
	print(f'{"":-<30}{"开始检查振动分析结果":^30}{"":<30}')
	run_walk(isprint=isprint)
	print(f'{"":-<30}{"检查结束":^30}{"":<30}')
