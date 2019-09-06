#!/home/zhangzhengde/bin/bin/python3
#coding=utf-8

import sys
sys.path.append('/home/zhangzhengde/bin/pythonlib')
import zzdlib as zzd
import os



data_Sh = zzd.File.openFile('./Vasp.sh','r')
#传入第1个参数是节点数，第2个参数是核数，第3个参数的任务名，第4个参数是是否加急
defnc = zzd.File.getLine(data_Sh,'#PBS -l nodes')[0].split()[-1]
defnodes = defnc.split(':')[0].split('=')[-1]
defppn = defnc.split(':')[1].split('=')[-1]
defEMER = zzd.File.getLine(data_Sh,'#PBS -q EMERGENCY')
try:
	para1 = eval(sys.argv[3])
	if type(para1) == int:
		print('错误的任务名，退出')
		exit()
except Exception as e:
	para1 = '1'

for i in range(1,len(sys.argv)):
	if i ==3:
		if sys.argv[i] !='def':
			data_Sh = zzd.File.substituteData(data_Sh,'#PBS -N ',' #PBS -N '+sys.argv[3])
	if i ==2: #2 3 一起了
		if sys.argv[1] != defnodes and sys.argv[2] == defppn:
			data_Sh = zzd.File.substituteData(data_Sh,'#PBS -l nodes',' #PBS -l nodes='+sys.argv[1]+':ppn='+defppn)
		if sys.argv[2] != defppn and sys.argv[1] == defnodes:
			data_Sh = zzd.File.substituteData(data_Sh,'#PBS -l nodes',' #PBS -l nodes='+defnodes+':ppn='+sys.argv[2])
		if sys.argv[1] != defnodes and sys.argv[2] != defppn:
			data_Sh = zzd.File.substituteData(data_Sh,'#PBS -l nodes',' #PBS -l nodes='+sys.argv[1]+':ppn='+sys.argv[2])
	if i ==4:
		if defEMER == 'Not Match':#没找到说明默认是不加急的
			if sys.argv[i] == 'yes' or sys.argv[i] == 'y':
				data_Sh.insert(4,' #PBS -q EMERGENCY\n')
		else: #默认加急
			if sys.argv[i] == 'no' or sys.argv[i] =='n':
				data_Sh.remove(' #PBS -q EMERGENCY\n')
zzd.File.openFile('./Vasp.sh','w',data=data_Sh)

	

if zzd.Vasp.checkInputs():
	usrsel=input('Parepare done! Would you like to submit the job? (default=yes)：')
	if usrsel  ==  'no' or usrsel =='n':
		print('未提交任务')
	else:
		os.system('qsub Vasp.sh')
else:
	print('前检查有问题，请人工检查')



