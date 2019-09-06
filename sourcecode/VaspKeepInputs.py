#!/home/zhangzhengde/bin/bin/python3
#coding=utf-8


import os
import sys


Files = 'INCAR,POSCAR,POTCAR,KPOINTS,Vasp.sh'
Flist = Files.split(',')
for i in range(1,len(sys.argv)):
	Flist.append(sys.argv[i])
print(Flist)


dir = 'KeepInputsDir'
os.system('mkdir '+dir)
for i in range(len(Flist)):
	os.system('cp ./'+Flist[i]+' ./'+dir+'/')
os.system('rm *')
os.system('cp ./'+dir+'/* ./')
os.system('rm -rf '+dir)


