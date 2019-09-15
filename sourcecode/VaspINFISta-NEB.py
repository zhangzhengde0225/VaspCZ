#!/home/zhangzhengde/bin/bin/python3

import sys
sys.path.append('/home/zhangzhengde/bin/pythonlib')
import os
import VaspCZ.zzdlib as zzdlib
import argparse


parser = argparse.ArgumentParser(description='manual to this script')
parser.add_argument('--nc', type=str, default='0,0')
parser.add_argument('--EMER', type=str, default='default')
parser.add_argument('--jobname', type=str, default='vaspNEB')
args = parser.parse_args()
print(f'脚本运行，参数:--nc={args.nc}  --EMER={args.EMER}  --jobname={args.jobname}')


dist = zzdlib.getshellResult('dist.pl ini/CONTCAR fin/CONTCAR')
dist = eval(dist[-1])
print('ini和fin中CONTCAR的dist为：{}'.format(dist))
if dist >= 9:
	print('dist过大，请检查')
	exit()
else:  # 向下取整数，如果是偶数则加一，如果是奇数直接用。0-1.9输入1，2-3.9属于3，4-5.9属于5
	image = int(dist/0.8)
	if image % 2 == 0:  # 是偶数
		image = image + 1
	os.system('nebmake.pl ini/CONTCAR fin/CONTCAR '+str(image))  # 插入中间的IMAGE
	os.system('cp ini/OUTCAR 00/')
	os.system('cp fin/OUTCAR 0'+str(image+1)+'/')
	# 拷贝输入文件
	os.system('cp ini/Opt/INCAR .')
	os.system('cp ini/Opt/KPOINTS .')
	os.system('cp ini/Opt/POTCAR .')
	# os.system('cp ini/Opt/Vasp.sh .')
	# 修改INCAR
	data_INCAR = zzdlib.File.openFile(path='./INCAR',mode='r')
	if data_INCAR[0] == 'SYSTEM=Opt\n':
		data_INCAR[0] = 'SYSTEM=NEB\n'
		data_INCAR = zzdlib.File.substituteData(data=data_INCAR,keywords='IBRION', newline='IBRION=1')
		data_INCAR.append('NFREE=2\n')
		data_INCAR.append('#neb\n')
		data_INCAR.append('IMAGES='+str(image)+'\n')
		data_INCAR.append('SPRING=-5\n')
		data_INCAR.append('LCLIMB=.TRUE.\n')
		data_INCAR.append('ICHAIN=0\n')
	zzdlib.File.openFile(path='./INCAR' ,mode='w', data=data_INCAR)
	# 修改Vasp.sh
	vaspsh_path = zzdlib.File.Vaspsh_path()
	os.system(f'cp {vaspsh_path}/Vasp.sh .')
	data_Sh = zzdlib.File.openFile('./Vasp.sh','r')
	# oldname = zzdlib.File.getLine(data_Sh,'#PBS -N')[0].strip('\n').split()[-1] #获取旧任务名
	# jobname = oldname[:-2] +'NEB'
	jobname = args.jobname
	data_Sh = zzdlib.File.substituteData(data_Sh, '#PBS -N', ' #PBS -N '+jobname)
	if args.nc.split(',')[0] == '0':
		nodes = str(image)
	else:
		nodes = args.nc.split(',')[0]
	if args.nc.split(',')[1] == '0':
		ppn = '8'
	else:
		ppn = args.nc.split(',')[1]
	data_Sh = zzdlib.File.substituteData(data_Sh, '#PBS -l nodes', ' #PBS -l nodes='+nodes+':ppn='+ppn)
	EMER = zzdlib.File.getLine(data_Sh,'#PBS -q EMERGENCY')[0]  # 获取默认加急状态
	if EMER == 'Not Match':
		isEMER = False
		if args.EMER == 'yes' or args.EMER == 'y':
			data_Sh.insert(4, ' #PBS -q EMERGENCY\n')
	else:
		isEMER = True
		if args.EMER == 'no' or args.EMER == 'n':
			data_Sh.remove(' #PBS -q EMERGENCY\n')
	zzdlib.File.openFile('./Vasp.sh', 'w', data_Sh)
	# 提交任务
	# zzdlib.Vasp.checkInputs()
	zzdlib.Vasp.check_and_qsub()
	# if zzdlib.Vasp.checkInputs():
	# 	usrsel=input('前检查无问题，是否要提交任务(默认yes)：')
	# 	if usrsel == 'no' or usrsel == 'n':
	# 		print('未提交任务')
	# 	else:
	# 		os.system('qsub Vasp.sh')
	# else:
	# 	print('前检查有问题，请人工检查')
	
