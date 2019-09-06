#!/home/zhangzhengde/bin/bin/python3

import sys
sys.path.append('/home/zhangzhengde/bin/pythonlib')
import os
import zzdlib
import argparse


inl = ['ini','fin']

parser = argparse.ArgumentParser(description = 'manual to this script')
parser.add_argument('--nc',type=str,default='0,0')
parser.add_argument('--EMER',type=str,default='default')
args = parser.parse_args()
print('脚本运行，参数:--nc={}  --EMER={}'.format(args.nc,args.EMER))


for i in inl:
	os.system('cp ./'+i+'/Opt/INCAR ./'+i)
	os.system('cp ./'+i+'/Opt/CONTCAR ./'+i+'/POSCAR')
	os.system('cp ./'+i+'/Opt/POTCAR ./'+i)
	os.system('cp ./'+i+'/Opt/KPOINTS ./'+i)
	os.system('cp ./'+i+'/Opt/Vasp.sh ./'+i)	
	os.chdir(i)
	with open('./INCAR','r') as f:
		data_INCAR = f.readlines()
		for nl in range(len(data_INCAR)):
			if 'SYSTEM' in data_INCAR[nl]:
				data_INCAR[nl] = 'SYSTEM=Static\n' #修改表头
			if 'NSW' in data_INCAR[nl]:
				data_INCAR[nl] = 'NSW=1\n' #修改NSW
			if 'IBRION' in data_INCAR[nl]:
				data_INCAR[nl] = 'IBRION=-1\n' #修改IBRION
			if 'EDIFFG' in data_INCAR[nl]:
				data_INCAR[nl] = '#'+data_INCAR[nl] #去掉EDIFFG
	with open('./INCAR','w') as f:
		f.writelines(data_INCAR)
		f.close()
		print("修改"+i+"静态INCAR完成")
	data_Sh = zzdlib.File.openFile('./Vasp.sh','r')
	#修改任务名
	oldname = zzdlib.File.getLine(data_Sh,'#PBS -N')[0].strip('\n').split()[-1]
	jobname = oldname[:-1]+'S'
	data_Sh = zzdlib.File.substituteData(data_Sh,'#PBS -N',' #PBS -N '+jobname)
	#修改nodes
	ndAndnc = zzdlib.File.getLine(data_Sh,'#PBS -l nodes')[0].strip('\n').split()[-1]
	nd = ndAndnc.split(':')[0].split('=')[-1]
	nc = ndAndnc.split(':')[1].split('=')[-1]
	if args.nc.split(',')[0] != '0':
		nd = args.nc.split(',')[0]
	if args.nc.split(',')[1] != '0':
		nc = args.nc.split(',')[1]
	data_Sh = zzdlib.File.substituteData(data_Sh,'#PBS -l nodes',' #PBS -l nodes='+nd+':ppn='+nc)
	#修改是否紧急
	EMER = zzdlib.File.getLine(data_Sh,'#PBS -q EMERGENCY')[0]#获取默认加急状态
	if EMER == 'Not Match':
		if args.EMER == 'yes' or args.EMER == 'y':
			data_Sh.insert(4,' #PBS -q EMERGENCY\n')
	else:
		if args.EMER == 'no' or args.EMER == 'n':
			data_Sh.remove(' #PBS -q EMERGENCY\n')
	zzdlib.File.openFile('./Vasp.sh','w',data_Sh)

	if zzdlib.Vasp.checkInputs():
		usrsel=input('前检查无问题，是否要提交任务(默认yes)：')
		if usrsel  ==  'no' or usrsel =='n':
			print('未提交任务')
		else:
			os.system('qsub Vasp.sh')
	else:
		print('前检查有问题，请人工检查')
	os.chdir('..')

	



