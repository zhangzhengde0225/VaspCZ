#!/home/zhangzhengde/bin/bin/python3
#coding=utf-8

import sys
sys.path.append('/home/zhangzhengde/bin/pythonlib/')
import zzdlib
import os

data_log = zzdlib.File.openFile('./log','r')
print('{:^5}  {}'.format(' ',data_log[12].strip('\n')))
for i in range(len(data_log)):
	elec = data_log[i].split()[0].strip(':')
	if elec == 'RMM' or elec == 'DAV':  # 该行开头是RMM DAV而下一行不是：
		try:
			elec2 = data_log[i+1].split()[0].strip(':')
			if (elec2 == 'RMM' or elec2 == 'DAV') == False:
				ions = data_log[i+1].split()[0] +data_log[i+1].split()[1].strip('=')
				print('{:^5}  {}'.format(ions,data_log[i].strip('\n')))
		except:
			print('{:^5}  {}'.format('run',data_log[i].strip('\n')))




