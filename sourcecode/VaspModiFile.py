#!/home/zhangzhengde/bin/bin/python3
#coding=utf-8

import sys
sys.path.append('/home/zhangzhengde/bin/pythonlib')
import zzdlib
import os


data_INCAR = zzdlib.File.openFile('./INCAR','r')
data_INCAR_old = data_INCAR.copy()
data_INCAR = zzdlib.File.substituteData(data_INCAR,'ISIF','ISIF=2')
data_INCAR = zzdlib.File.substituteData(data_INCAR,'NELM','NELM=400','default')
if zzdlib.File.getLine(data_INCAR,'NELMDL') == 'Not Match':
	data_INCAR.append('NELMDL=5\n')
	data_INCAR.append('LMAXMIX=4\n')
	data_INCAR.append('AMIX=0.2\n')
	data_INCAR.append('BMIX=0.0001\n')
	data_INCAR.append('AMIX_MAG=0.8\n')
	data_INCAR.append('BMIX_MAG=0.0001\n')

zzdlib.File.openFile('./INCAR','w',data=data_INCAR)

for i in range(len(data_INCAR)):
	try:
		print('{}    {}'.format(data_INCAR_old[i].strip('\n'),data_INCAR[i].strip('\n')))
	except Exception as e:
		print('{}    {}'.format(' ',data_INCAR[i].strip('\n')))

