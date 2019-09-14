#!/home/zhangzhengde/bin/bin/python3
#coding=utf-8

import os, sys
import zzdlib as zzd
import subprocess
import string
import shutil

def CheckWARNING(usript):
	#usrinput为1时，只打印是否完成，为2时候，只打印能量信息，为3时，打印全部
	Path = []
	Energy =[]
	Step =[]
	mag =[]
	dist =[]
	RMS =[]
	warnflag = 0
	jstat = ['Running...','Done','SC Running...','SC done','Stoped']
	for dirpath,dirnames,filenames in os.walk('./'):
		if 'INCAR' in filenames and 'POSCAR' in filenames and 'Vasp.sh' in filenames:
			data_Sh = zzd.File.openFile(dirpath+'/Vasp.sh','r')
			data_INCAR = zzd.File.openFile(dirpath+'/INCAR','r')
			SYSTEM = zzd.File.getLine(data_INCAR,'SYSTEM')[0].split('=')[-1]
			jobname = zzd.File.getLine(data_Sh,'#PBS -N')[0].strip('\n').split()[-1]
			jobstatus = zzd.Vasp.checkJobstatus(jobname)
			#print('Path:{}  {}   {}'.format(dirpath,jobstatus,jobname))
			if jobstatus == 'Q':
				print('Path:{:<60}任务正在排队...   '.format(dirpath))
			elif jobstatus == 'R' or 'log' in filenames: #正在计算和算完了用相同的检测方法
				if usript == '1' or usript =='3':
					data_log = zzd.File.openFile(dirpath+'/log','r')
					reachRA = zzd.File.getLine(data_log,'reached required accuracy')[0]
					try:
						RMM = zzd.File.getAllline(data_log,'RMM:')[-1].split()[1]
						RMM = RMM+' RMM'
					except:
						RMM = '0 RMM'
					try:
						DAV = zzd.File.getAllline(data_log,'DAV:')[-1].split()[1]
						DAV = DAV+' DAV'
					except:
						DAV = '0 DAV'
					termination = zzd.File.getLine(data_log,'Ctrl-C caught... cleaning up processes')[0]
					try:
						ionstep = zzd.File.getAllline(data_log,'F=')[-1].split()[0]
					except:
						ionstep = '0'	
					#print('Path:{}  {}   {}'.format(dirpath,reachRA,v1F))
					if 'reached required accuracy' in reachRA:
						print('Path:{:<60}{:<15}{:>3}F{:>9}'.format(dirpath,jstat[1],ionstep,RMM))
					elif SYSTEM == 'Self' or SYSTEM =='Static' or SYSTEM == 'self' or SYSTEM =='static':
						vRMM = zzd.File.getLine(data_log,'RMM:')[0]
						vDAV = zzd.File.getLine(data_log,'DAV:')[0]
						if ('DAV:' in vDAV or 'RMM:' in vRMM) and ionstep =='0':
							print('Path:{:<60}{:<15}{:>3}F{:>9}{:>7}{:>5}'.format(dirpath,jstat[2],ionstep,RMM,DAV,'*'))
						elif ionstep == '1': #没有完成，但有1F了
							print('Path:{:<60}{:<15}{:>3}F{:>9}'.format(dirpath,jstat[3],ionstep,RMM))
					elif 'cleaning up processes' in termination:
						print('Path:{:<60}{:<15}{:>3}F{:>9}{:>7}'.format(dirpath,jstat[4],ionstep,RMM,DAV))
					else:
						print('Path:{:<60}{:<15}{:>3}F{:>9}{:>7}{:>5}'.format(dirpath,jstat[0],ionstep,RMM,DAV,'*'))
				
					#检查警告
					try:
						data_OUTCAR = zzd.File.openFile(dirpath+'/OUTCAR','r')
						WARNING_log = zzd.File.getAllline(data_log,'WARNING')
						WARNING_OUT = zzd.File.getAllline(data_OUTCAR,'WARNING')
						ERROR_log = zzd.File.getAllline(data_log,'ERROR')
						if WARNING_log !=[] or WARNING_OUT !=[]:
							print('Path:{:<60}    出现警告'.format(dirpath,ionstep))
							zzd.File.printData(WARNING_log)
							zzd.File.printData(WARNING_OUT)
							warnflag = 1
						if ERROR_log !=[]:
							zzd.File.printData(ERROR_log)
							warnflag =1
					except Exception as e:
						print('path:{} {}'.format(dirpath,e))
				
				if usript == '2' or usript == '3':
					#打印能量信息
					data_log = zzd.File.openFile(dirpath+'/log','r')
					try:
						log_lastF = zzd.File.getAllline(data_log,'F=')[-1]
					except:
						log_lastF = '0 F= 0 E0= 0  d E 0  mag=    0' #在存在log，Log中连1步都没算完的时候
					if len(dirpath) >15:
						path = '...'+dirpath[-15:]
					else:
						path = dirpath
					Path.append(path)
					Energy.append(log_lastF.split()[4])
					Step.append(log_lastF.split()[0])
					try:
						mag.append(log_lastF.split()[9])
					except:
						mag.append('0')
					dist.append(zzd.getshellResult('dist.pl '+dirpath+'/POSCAR '+dirpath+'/CONTCAR')[-1].strip('\n'))
					try:
						RMS.append(zzd.getshellResult('grep RMS '+dirpath+'/OUTCAR')[-1].split()[4])
					except:#如果出现没算完，grep返回一个空的列表的时候
						RMS.append('0')

	if (usript == '1' or usript == '3') and warnflag == 0:
		print('无警告或错误')
	if (usript == '2' or usript == '3'):
		print('{:<18}{:<11}{:<5}{:^10}{:^10}{:^10}'.format('路径','能量','步数','mag','dist','RMS'))
		for i in range(len(Path)):
			print('{:<20}{:<13.4f}{:<7}{:<12.4f}{:<10.4f}{:<10.4f}'.format(Path[i],eval(Energy[i]),Step[i],eval(mag[i]),eval(dist[i]),eval(RMS[i])))
	
def VaspNEBCheckDist(POSorCONT):
	if '00' in os.listdir() and 'INCAR' in os.listdir():
		data_log = zzd.File.openFile('./INCAR')
		image = zzd.File.getLine(data_log,'IMAGES')[0].split('=')[-1].strip('\n')
		if int(image) <=9:
			#os.system('cp ini/CONTCAR 00/CONTCAR')
			#os.system('cp fin/CONTCAR 0'+str(int(image)+1)+'/CONTCAR')
			for i in range(0,int(image)+1):
				if i == 0:
					dist = zzd.getshellResult('dist.pl ./0'+str(i)+'/POSCAR ./0'+str(i+1)+'/'+POSorCONT)
				elif i == int(image):
					dist = zzd.getshellResult('dist.pl ./0'+str(i)+'/'+POSorCONT+' ./0'+str(i+1)+'/POSCAR')
				else:
					dist = zzd.getshellResult('dist.pl ./0'+str(i)+'/'+POSorCONT+' ./0'+str(i+1)+'/'+POSorCONT)
				print('{}  0{}-0{}  {}'.format(POSorCONT,i,i+1,dist[0].strip('\n')))
		else:
			print('image too large')
	else:
		print('当前不在NEB目录，退出程序')
		exit()
def VaspNEBCheckRMS():
	print('CheckNEBRMS is running...')
	if '00' in os.listdir() and 'INCAR' in os.listdir():
		data_INCAR = zzd.File.openFile('./INCAR')
		data_log = zzd.File.openFile('./log')
		image = zzd.File.getLine(data_INCAR,'IMAGES')[0].split('=')[-1].strip('\n')
		#print(zzd.File.getAllline(data_log,'F='))
		ionstep = zzd.File.getAllline(data_log,'F=')[-1].split()[0]	
		if int(image) <=9:
			#print('aa')
			#data00 = zzd.getshellResult('grep RMS 00/OUTCAR')
			data01 = zzd.getshellResult('grep RMS 01/OUTCAR')
			data02 = zzd.getshellResult('grep RMS 02/OUTCAR')
			data03 = zzd.getshellResult('grep RMS 03/OUTCAR')
			#data04 = zzd.getshellResult('grep RMS 04/OUTCAR')
			#print(data00)
			print('{:^3}{:^10}{:^10}{:^10}{:^10}'.format('Step','01-RMS','02-RMS','03-RMS','01+02+03'))
			for i in range(int(ionstep)):
				print('{:>3}{:>10}{:>10}{:>10}{:>10.6f}'.format(i+1,data01[i].split()[4],data02[i].split()[4],data03[i].split()[4],float(data01[i].split()[4])+float(data02[i].split()[4])+float(data03[i].split()[4])))
	else:
		print('不在NEB目录，退出程序')
		exit()



if __name__ =='__main__':
	VaspCZ_path = os.path.dirname(os.path.abspath(__file__))+ '/sourcecode'
	while True:
		ipt = input('''
|============================================================|
|                        VASP TOOLS                          |
|------------------------------------------------------------|
|        (1)Check                                            |
|        (2)NEBCheck                                         |
|        (3)VaspOpt-Sta                                      |
|        (4)VaspINFISta-NEB                                  |
|        (5)VaspKeepInputs                                   |
|        (6)VaspModiFile                                     |
|        (7)VaspGetNELM                                      |
|        (8)VaspQsub                                         |
|        (9)VaspMiniTools                                    |
|        (0)Quit                                             |
|------------------------------------------------------------|
|               by: Zhengde Zhang (zhangzhengde@sinap.ac.cn) |
|============================================================|
Input(default=1):  ''')
		if ipt == '1' or ipt == '':
			while True:
				ipt1 = input('''
|============================================================|
|                       Check Parameter                      |
|------------------------------------------------------------|
|        (1)Only Job Status                                  |
|        (2)Current Results                                  |
|        (3)All                                              |
|        (0)Back                                             |
|============================================================|
Input(default=3):  ''') 
				if ipt1 == '':
					ipt1 = '3'
				if ipt1 != '0':
					print('\nRun Vasp Check...')
					CheckWARNING(ipt1)
					exit()
				else:
					break
				
		elif ipt == '2':
			while True:
				ipt2 = input('''
|============================================================|
|                    Check NEB Parameter                     |
|------------------------------------------------------------|
|        (1)Only Job Status                                  |
|        (2)NEB Results                                      |
|        (3)Above                                            |
|        (4)NEB Barrier                                      |
|        (0)Back                                             |
|============================================================|
Input(default=4):  ''')
				if ipt2 == '1':
					subprocess.call(f'{sys.executable} {VaspCZ_path}/NEBCheck1.1.py --func=1',shell=True)
				elif ipt2 == '2':
					subprocess.call(f'{sys.executable} {VaspCZ_path}/NEBCheck.py --func=2',shell=True)
				elif ipt2 == '3':
					subprocess.call(f'{sys.executable} {VaspCZ_path}/NEBCheck.py --func=3',shell=True)
				elif ipt2 == '4' or ipt2 == '':
					subprocess.call(f'{sys.executable} {VaspCZ_path}/NEBCheck.py --func=4',shell=True)
				elif ipt2 == '0':
					break
				else:
					continue
				exit()
		elif ipt == '3':
			ipt3plus = input('''
|============================================================|
|                   Vasp  Opt-Sta                            |
|------------------------------------------------------------|
|        (1) Opt-Sta (Current Folder)                        |
|        (2) INI FIN Opt-Sta                                 |
|        (0) back                                            | 
|============================================================|
Input:(defalut=1) ''')
			while True:
				if ipt3plus == '0':
					break
				elif ipt3plus == '1' or ipt3plus == '':
					subprocess.call(f'{sys.executable} {VaspCZ_path}/VaspOpt-Sta.py')
					exit()
				elif ipt3plus == '2':
					ipt3 = input('''
|============================================================|	
|                   Vasp INI FIN Opt-Sta                     |
|------------------------------------------------------------|
|  This tool is used for Running self-consistent calculations|
|when structure optimization are finished in ini/Opt and fin/|
|Opt respectively.                                           |
|------------------------------------------------------------|
|  Parameter $1 $2 $3 represent nodes,ppn and isEMERGENCY re-|
|spectively.                                                 |
|  i.e input '2 12 yes' means nodes=2,ppn=12 and run on EMER-|
|GENCY node. Default setting is from ini/Opt/Vasp.sh and fin-|
|/Opt/Vasp.sh                                                |
|        (1) 0 0 default
|        (0) back
|============================================================|
Input:  ''')
					while True:
						if ipt3 == '0':
							break
						else:	
							try:
								nc = ipt3.split()[0]+','+ipt3.split()[1]
							except:
								nc = '0,0'
							try:
								EMER = ipt3.split()[2]
							except:
								EMER='default'
							subprocess.call(f'{sys.executable} {VaspCZ_path}/VaspINFIOpt-Sta.py  --nc='+nc+' --EMER='+EMER,shell=True)
							exit()

		elif ipt == '4':
			ipt4 = input('''
|============================================================|
|                   Vasp INI FIN Sta-NEB                     |
|------------------------------------------------------------|
|  This tool is used for Running NEB calculations when self- |
|onsistent are finished in ini/ and fin/ respectively.       |
|------------------------------------------------------------|
|  Parameter $1 $2 $3 represent nodes,ppn and isEMERGENCY re-|
|spectively.                                                 |
|  i.e input '2 12 yes' means nodes=2,ppn=12 and run on EMER-|
|GENCY node. Default setting is from ini/Opt/Vasp.sh and fin-|
|/Opt/Vasp.sh                                                |
|        (1) 0 0 default                                     |
|        (0) back                                            |
|============================================================|
Input:  ''')
			while True:
				if ipt4 == '0':
					break
				else:
					try:
						nc = ipt4.split()[0]+','+ipt4.split()[1]
					except:
						nc = '0,0'
					try:
						EMER = ipt4.split()[2]
					except:
						EMER ='default'
					subprocess.call(f'{sys.executable} {VaspCZ_path}/VaspINFISta-NEB.py  --nc='+nc+' --EMER='+EMER,shell=True)
					exit()
		elif ipt == '5':
			subprocess.call(f'{sys.executable} {VaspCZ_path}/VaspKeepInputs.py',shell=True)
			break
		elif ipt == '6':
			subprocess.call(f'{sys.executable} {VaspCZ_path}/VaspModiFile.py',shell=True)
			break
		elif ipt == '7':
			subprocess.call(f'{sys.executable} {VaspCZ_path}/VaspGetNELM.py',shell=True)
			break
		elif ipt == '8':
			subprocess.call(f'{sys.executable} {VaspCZ_path}/VaspQsub.py',shell=True)
			break
		elif ipt == '9':
			while True:
				ipt9 = input('''
|============================================================|
|                     VASP Mini Tools                        |
|------------------------------------------------------------|
|        (1)VaspNEBCheckdist                                 |
|        (2)VaspNEBCheckRMS                                  |
|        (3)VaspNEBKeepINFI_Optinputs                        |
|        (4)NEBKeepInputs                                    |
|        (5)k-point mesh test                                |
|        (0)back                                             |
|============================================================|
Input(default=2):  ''')
				if ipt9 == '0':
					break
				elif ipt9 == '1':
					ipt91 = input('Check POS or CONT(default=POS):  ')
					if ipt91 == 'CONT':
						VaspNEBCheckDist('CONTCAR')
					else:
						VaspNEBCheckDist('POSCAR')	
				elif ipt9 == '' or ipt9 == '2':
					VaspNEBCheckRMS()
				elif ipt9 == '3':
					subprocess.call(f'{sys.executable} {VaspCZ_path}/VaspNEBKeepINFI_OptInputs.py',shell=True)
				elif ipt9 == '4':
					subprocess.call(f'{sys.executable} {VaspCZ_path}/NEBKeepInputs.py',shell=True)
				elif ipt9 == '5':
					while True:
						ipt95 = input('''
|============================================================|
|                  VASP k-point mesh test                    |
|------------------------------------------------------------|
|  This is the Vasp k-point mesh test tool. Please prepare t-|
|he input files (INCAR, POSCAR, POTCAR, KPOINTS, Vasp.sh) of |
|VASP in current directory before run the tool.              |
|  Please input the jobname_prefix, nodes, ppn and the k_mesh|
|you want to test. (separate with one white space)           |
|  For example: --jobname_prefix=k_test --nodes=1 --ppn=8    |
| --k_mesh=111,333,555,777,999                               | 
|         -----------------------------------------          |    	
|             parameter              default                 |
|         -----------------------------------------          |
|            jobname_prefix           k_test                 |
|    	         nodes                   1                   |
|                 ppn                    8                   |
|                k_mesh         111,333,555,777,999          |
|         -----------------------------------------          |
|  Go back while input 0                                     |
|============================================================|
Input:  ''')
						if ipt95 == '0':
							break
						else:
							try:
								inputs = ipt95.split()
								dp = {
									'jobname_prefix':'k_test', 'nodes':'1', 'ppn':'8',
									'k_mesh': '111,333,555,777,999'}
								for item in inputs:
									para, value = item.strip('--').split('=')
									print(f'输入的参数为{para}:{value}')
									dp[para] = value
								# print(dp)
								code = f'{sys.executable} {VaspCZ_path}/k_point_test.py --jobname_prefix={dp["jobname_prefix"]} --nodes={dp["nodes"]} --ppn={dp["ppn"]} --k_mesh={dp["k_mesh"]}'
								print(code)
								subprocess.call(code, shell=True)
								exit()
							except Exception as e:
								print(f'Vasp k-mesh test error: {e}')
								continue
				else:
					continue
		elif ipt == '0':
			exit()
		else:
			print('Error Input')
			continue






