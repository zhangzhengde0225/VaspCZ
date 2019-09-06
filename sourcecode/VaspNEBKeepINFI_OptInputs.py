#!/home/zhangzhengde/bin/bin/python3
#coding=utf-8


import os
import sys



alldir = os.listdir()


#只保留ini fin
for dir in alldir:
	if dir == 'ini' or dir == 'fin':
		continue
	else:
		os.system('rm -rf '+dir)

for infi in ['ini','fin']:
	os.chdir(infi)
	os.system('rm *')
	os.chdir('Opt')

	Flist = 'INCAR,POSCAR,POTCAR,KPOINTS,Vasp.sh'.split(',')
	dir = 'KeepInputsDir'
	os.system('mkdir '+dir)
	for i in range(len(Flist)):
		os.system('cp ./'+Flist[i]+' ./'+dir+'/')
	os.system('rm *')
	os.system('cp ./'+dir+'/* ./')
	os.system('rm -rf '+dir)
	os.chdir('../..')

