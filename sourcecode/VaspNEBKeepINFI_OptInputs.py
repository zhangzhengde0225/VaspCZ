#!/home/zhangzhengde/bin/bin/python3
#coding=utf-8


"""
更新，包含了提示要删除文件。
"""

import os
import sys

alldir = os.listdir()

# 只保留ini fin
del_files = []
for dir in alldir:
	if dir == 'ini' or dir == 'fin':
		continue
	else:
		del_files.append(dir)

del_ini_files = os.listdir('ini')
del_ini_files.remove('Opt')
del_fin_files = os.listdir('fin')
del_fin_files.remove('Opt')

del_ini_Opt_files = os.listdir('ini/Opt')
del_fin_Opt_files = os.listdir('fin/Opt')
for file in ['INCAR', 'POSCAR', 'POTCAR', 'KPOINTS', 'Vasp.sh']:
	del_ini_Opt_files.remove(file)
	del_fin_Opt_files.remove(file)

ipt = input(f'即将删除文件和文件夹\n当前目录下: {del_files}\nini/下: {del_ini_files}\nini/Opt/下: {del_ini_Opt_files}\nfin/下: {del_fin_files}\nfin/Opt/下: {del_fin_Opt_files}\n是否确定([y]es/no): ')

if ipt in ['y', 'yes', 'Y', 'YES', '']:
	for file in del_files:
		os.system(f'rm -rf {file}')
	for infi in ['ini', 'fin']:
		os.chdir(infi)
		need_del_files = eval(f'del_{infi}_files')
		for file in need_del_files:
			os.system(f'rm -rf {file}')
		os.chdir('Opt')
		need_del_files = eval(f'del_{infi}_Opt_files')
		for file in need_del_files:
			os.system(f'rm -rf {file}')
		os.chdir('../..')
else:
	print(f'未删除任何东西')

