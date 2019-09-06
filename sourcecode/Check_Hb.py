"""
针对fcc获得Hb
遍历当前文件夹下的文件下的文件夹，fcc_withele存在的时候，就计算一个Hb
"""

import os

def cal_Hb():
	Hb = item1 + item2 - item3 -item4
	return Hb

def run():

	for path, dirnames, filenames in os.walk():
		for dirname in dirnames:
			if dirname[:-2] == 'fcc_with':
				print(path, dirname)

if __name__ == '__main__':
	run()
