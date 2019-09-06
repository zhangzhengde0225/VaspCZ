#!/home/zhuguifeng/BIN/bin/python3
# -*- coding: utf-8-*-

import sys
sys.path.append('/home/zhangzhengde/bin/pythonlib')
import zzdlib as zzd
import os
import subprocess
import time
import argparse

def isNEB(dirpath):
	with open(dirpath+'/log') as file_log:
		flag = 0
		data_log = file_log.readlines()
		for i in range(6):
			if 'each image running on' in data_log[i]:
				flag = 1
	if flag  == 1:
		return True
	elif flag == 0:
		return False
def openFile(filepath):
	with open(filepath) as file:
		data = file.readlines()
	return data

def isTaskrun(dirpath):#提供路径，判断该路径下Vasp.sh中的任务名是否在mjb中正在计算，是返回True和mjb中该行信息，否返回False和空字符串
	data_Vaspsh  = openFile(dirpath+'/Vasp.sh')
	for nl in range(len(data_Vaspsh)): #从Vasp.sh中获取任务名
		if '#PBS -N' in data_Vaspsh[nl]:
			jobname = data_Vaspsh[nl].split('-N')[1].strip('\n').strip(' ')
	usrname = os.getcwd().split('home/')[1].split('/')[0]  
	ob_mjb = subprocess.Popen(['qstat','-x','-u',usrname],stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE) 
	data_mjb = ob_mjb.stdout.readlines()
	for nl in list(reversed(range(len(data_mjb)))):
		data_mjb[nl] = bytes.decode(data_mjb[nl]) 
		if jobname in data_mjb[nl]: #如果任务名存在
			line = data_mjb[nl].strip('\n')
			if line.split()[9] == 'R': #且任务正在运行
				isrun = True
				line = line
				break
			else:
				isrun = False
				line = 'job exist do not running'
		else:
			isrun = False
			line = 'job not exist'
	return isrun,line


	#1.检查是否有错误
def CheckWARNING(directory,reQsub):	
	for dirpath, dirnames, filenames in os.walk(directory):
		warnflag = 0
		doneflag = 0
		forceflag = 0
		step =0
		if 'log' in filenames:
			try:
				if isNEB(dirpath) == True:
					data_INCAR = openFile(dirpath+'/INCAR')
					for i in range(len(data_INCAR)):
						if 'NSW' in data_INCAR[i]:
							NSW = data_INCAR[i].split('=')[1].strip('\n').strip(' ')
						if 'IMAGES=' in data_INCAR[i]:
							IMAGES = data_INCAR[i].split('=')[1].strip('\n').strip(' ')
							break
					os.system('cp '+dirpath+'/log '+dirpath+'/01/stdout')#拷贝neb目录的log到01文件夹中
					for i in range(1,eval(IMAGES)+1):
						data_log = openFile(dirpath+'/0'+str(i)+'/stdout')
						for j in range(len(data_log)):
							if 'reached required accuracy' in data_log[j]:#判断是否完成计算
								#print('Path:{}  计算完成。'.format(dirpath))
								doneflag = doneflag+1
							if 'WARNING' in data_log[j]:
								print('Path:{}  stdout WARNING:'.format(dirpath+'/0'+str(i)+'/'))
								print(data_log[j])
								warnflag = 1
								with open(dirpath+'/0'+str(i)+'/OUTCAR') as file_OUTCAR:#log中存在WARNING再检测OUTCAR中的
									data_OUTCAR = file_OUTCAR.readlines()
									for nl in range(len(data_OUTCAR)):
										if 'WARNING' in data_OUTCAR[nl]:
											print('Path:{}  OUTCAR WARNING:'.format(dirpath+'/0'+str(i)+'/')) 
											print(data_OUTCAR[nl])
											warnflag=1
											break
									break
							if 'Ctrl-C caught... cleaning up processes' in data_log[j]:
								forceflag =1
								break
					data_rootdirlog = openFile(dirpath+'/log') #获取Log中的数据
					if len(data_rootdirlog) == 0:
						loglen = 'log without data'
					else:
						if 'running on' in data_rootdirlog[0]:
							loglen = 'log with right data'
						for nl in list(reversed(range(len(data_rootdirlog)))):
							if 'F=' in data_rootdirlog[nl]:
								step = data_rootdirlog[nl].split('=')[0].strip('F').strip(' ')
								break
					if doneflag <eval(IMAGES):
						isrun,line = isTaskrun(dirpath)
						if isrun == True:
							if forceflag == 1:
								print('Path:{:<40}  NEB强制结束计算    {}F'.format(dirpath,step))
							elif step == NSW:
								print('Path:{:<40}  NEB计算达到设定步数    {}F'.format(dirpath,step))
							else:
								print('Path:{:<40}  NEB正在计算...    {}F'.format(dirpath,step))
						else:
							print('Path:{:<40}  NEB计算结束       {}F'.format(dirpath,step))
					elif doneflag == eval(IMAGES):
						print('Path:{:<40}  NEB计算完成!'.format(dirpath))
			except Exception as e:
				print('Path:{:<40}  错误类型:{}'.format(dirpath,e))
				#判断错误，是否重新提交任务。
				data_rootdirlog = openFile(dirpath+'/log')
				if len(data_rootdirlog) == 0: #log中无数据
					loglen = 'log without data'
					isrun,line = isTaskrun(dirpath)
					if isrun == True:
						runtime = line.split()[10]
						runpath = os.getcwd()
						if int(line.split()[10].split(':')[0]) >=1 or int(line.split()[10].split(':')[1]) >= 2: #且任务运行时间大于1分钟
							#提示运行，删除并且重新提交任务
							print('发现路径{}/log无值且程序运行时间为{}'.format(runpath+dirpath.strip('.'),runtime))
							choose = reQsub
							if choose == 'yes' or choose=='y':
								jobID = line.split()[0].split('.')[0]
								os.chdir(runpath+dirpath.strip('.'))
								os.system('qdel '+jobID)
								print('已经删除任务{}.msvr1'.format(jobID))	
								time.sleep(1)
								os.system('qsub Vasp.sh')
								print('已经重新提交{}NEB计算'.format(runpath+dirpath.strip('.')))
								os.chdir(runpath)
							elif choose == '' or choose=='no' or choose=='n':
								continue
							else:
								print('输入错误')
								continue
						else:
							print('发现路径{}/log无值且程序运行时间为{},程序刚刚运行'.format(runpath+dirpath.strip('.'),runtime))
				else:
					#print(data_rootdirlog)
					for mml in range(len(data_rootdirlog)):
						if data_rootdirlog[mml] == 'Ctrl-C caught... cleaning up processes\n':
							print('Path:{:<40}   手动退出'.format(dirpath))
							break


	if warnflag == 0:
		print('log和OUTCAR中无警告')
		

def CheckEnergy(directory):
	#2.打印能量等信息
	runpath = os.getcwd()
	for dirpath, dirnames, filenames in os.walk(directory):
		if 'log' in filenames:
			try:
				if isNEB(dirpath) == True:
					os.chdir(dirpath)
					obj = subprocess.Popen([f'{os.environ["HOME"]}/bin/vtst/nebef.pl'],stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
					nebef = str(obj.stdout.read()).split("'")[1].split('\\n')
					print('Path:{}'.format(runpath+dirpath))
					print(' IMAGE          RMS          Energy            Barrier')
					for i in range(len(nebef)):
						print(nebef[i].strip('\n'))
					os.chdir(runpath)
			except Exception as e:
				print('Path:{} 错误类型:{}'.format(runpath+dirpath,e))

def getBarrier(IS='-1'):
	#IS代表输入的Ionstep，根据IS获取该步时候的扩散势垒
	Path = []
	Barrier = []
	Step = []
	Image = []
	runpath = os.getcwd()
	for dirpath, dirnames, filenames in os.walk('./'):
		if 'INCAR' in filenames and 'POTCAR' in filenames and 'KPOINTS' in filenames:
			data_INCAR=zzd.File.openFile(dirpath+'/INCAR','r')
			#print('dir',data_INCAR[0].split('=')[-1])
			if 'NEB' in data_INCAR[0].split('=')[-1]:
				#print('NEB')
				data_Sh = zzd.File.openFile(dirpath+'/Vasp.sh','r')
				jobname = zzd.File.getLine(data_Sh,'#PBS -N')[0].split()[-1][:10]  #获取任务名，且只取前10个字符
				stat = zzd.Vasp.checkJobstatus(jobname)
				if stat == 'R' or stat == 'Q':
					print('路径{} NEB计算正在{}'.format(dirpath,stat))
				elif 'log' in os.listdir(dirpath):#未提交或者已经算完：
					data_log = zzd.File.openFile(dirpath+'/log','r')
					isRA = zzd.File.getLine(data_log,'reached required accuracy')[0]
					images = zzd.File.getLine(data_INCAR,'IMAGES')[0].split('=')[-1]
					Fstep = zzd.File.getLine(data_INCAR,'NSW')[0].split('=')[-1]
					ionstep = zzd.File.getAllline(data_log,'F=')[-1].split()[0]
					if 'reached required accuracy' in isRA or Fstep == ionstep:
						#print('RA')
						barriers = []
						os.chdir(dirpath)
						if IS == '-1':#默认用nebef.pl获取势垒
							nebef = zzd.getshellResult('nebef.pl')
							for inn in nebef:
								barriers.append(float(inn.split()[3]))
						else:#获取输出的步数IS，计算该步数下的势垒。
							try:	
								data_OUTini = zzd.File.openFile('./00/OUTCAR','r')
								Eini = zzd.File.getAllline(data_OUTini,'energy without entropy')[-1].split()[-1]
								data_OUTfin = zzd.File.openFile('./0'+str(int(images)+1)+'/OUTCAR','r')
								Efin = zzd.File.getAllline(data_OUTini,'energy without entropy')[-1].split()[-1]
								barriers.append(float(Eini)-float(Eini))
								for ii in range(1,int(images)+1):
									data = zzd.File.openFile('./0'+str(ii)+'/stdout','r')
									Ener = zzd.File.getLine(data,IS+' F=')[0].split()[4]
									barriers.append(float(Ener)-float(Eini))
							except Exception as e:
								print('获取第{}步势垒出现问题，提示为：{}'.format(IS,e))
						Barrier.append(max(barriers))
						Step.append(ionstep)
						Path.append(dirpath[:18])
						Image.append(images)
						os.chdir(runpath)
					else:
						print('{} 未达到收敛标准或者计算步数未满，退出'.format(dirpath))
				else:
					print('出现了奇怪的问题')
	print('{:<20}{:<10}{:<5}{:<3}'.format('path','barrier','step','images'))
	for i in range(len(Barrier)):
		print('{:<20}{:<10.7f}{:<5}{:<3}'.format(Path[i],Barrier[i],Step[i],Image[i]))


if __name__ == '__main__':
	print('脚本运行，功能列表： 1.检查计算是否完成(默认)  2.获取NEB计算结果  3.以上全部  --func=1 --reQsub=no')
	parser = argparse.ArgumentParser(description ='Manual')
	parser.add_argument('--func',type=str,default='1')
	parser.add_argument('--reQsub',type=str,default='no')
	args = parser.parse_args()

	choose = args.func
	print('NEBCheckBegins...   --func={}'.format(choose))
	if choose == '1' or choose =='3':
		CheckWARNING('./',args.reQsub)
	if choose == '2' or choose =='3':
		CheckEnergy('./')
	if choose == '4':
		IS = input('检查NEB第几步的势垒(默认最后一步)： ')
		if IS == '':
			IS='-1'
		getBarrier(IS)
	print('NEBCheck done!')

					
			
