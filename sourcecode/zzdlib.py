#!/home/zhangzhengde/bin/bin/python
#coding=utf-8

import os
import subprocess
import time
import numpy as np

user = 'zhangzhengde'
current_py_folder = os.path.dirname(os.path.abspath(__file__))
VaspCZ_path = [os.path.dirname(current_py_folder) if 'sourcecode' in current_py_folder else current_py_folder][0] + '/sourcecode'

def setUser(inputusername):
	global user
	user = inputusername

def getshellResult(code): #获取shell代码的输出，每一行加一个\n结尾，作为字符串组成列表。最后一行如果为空行则删除
	result = subprocess.check_output(code,stderr=subprocess.STDOUT,shell=True)
	result = result.decode('utf-8').split('\n')
	for i in range(len(result)):
		result[i] = result[i]+'\n'	
	if result[-1] == '\n':
		result = result[:-1]
	return result

def listSum(numlist):#输入是以数字或数字字符串组成的List,返回和，类型是Int或者float
	sum = 0
	flag = 0 
	for i in range(len(numlist)):
		try:
			sum  = sum + eval(numlist[i])
			flag = flag +1
		except:
			sum = sum + numlist[i] 
			flag = flag+1
	if flag != len(numlist):
		print('zzdlib函数listSum错误')
	return sum

class File():
	def openFile(path,mode='r',data=None):
		if mode == 'r':
			with open(path,mode) as f:
				data = f.readlines()
				return data
		elif mode == 'w':
			with open(path,mode) as f:
				f.writelines(data)
	def substituteData(data,keywords,newline,mode='default'):  #给出关键词和新行，默认情形搜索出现第一次出现关键词的行并替换，mode等于别的可以设置替换全部出现关键字的行
		for i in range(len(data)):
			if keywords in data[i]:
				if newline[-1] == '\n':
					data[i] = newline
				else:
					data[i] = newline + '\n'
				if mode == 'default':
					break
		return data
	def getLine(data,keywords):#给出关键词，返回有关键词的第一行并返回，返回为字符串和所在的行号
		for i in range(len(data)):
			if keywords in data[i]:
				return data[i].strip('\n'),i
		return 'Not Match',0


	def getAllline(data,keywords):  #给出关键词，返回所有带有关键词的行，返回为列表
		result = []
		for i in range(len(data)):
			if keywords in data[i]:
				result.append(data[i])
		return result
	def getNullline(data):
		data2 = data.copy()
		result = []
		for i in range(len(data2)):
			data2[i] = data2[i].strip(' ')
			if data2[i] == '\n':
				result.append(i)
		return result
	
	def printData(data):
		for i in range(len(data)):
			print(data[i].strip('\n'))

	def Vaspsh_path():
		with open(f'{VaspCZ_path}/build-in_data.txt', 'r') as f:
			data =f.readlines()
		Vaspsh_path = data[1].split('=')[1].strip('\n')
		return Vaspsh_path


class Vasp():
	def checkInputs():#在提交任务之前，检查Vasp的各项输入：
		#检查INCAR
		data_INCAR = File.openFile('./INCAR','r')
		SYSTEM = File.getLine(data_INCAR,'SYSTEM')[0].split('=')[-1]
		ENCUT = File.getLine(data_INCAR,'ENCUT')[0].split('=')[-1]
		ISIF = File.getLine(data_INCAR,'ISIF')[0].split('=')[-1]
		IBRION = File.getLine(data_INCAR,'IBRION')[0].split('=')[-1]
		if IBRION == '-1':
			IBRION = 'No update'
		elif IBRION == '0':
			IBRION = 'Molecular dynamics'
		elif IBRION == '1':
			IBRION = 'qussi-Newton'
		elif IBRION == '2':
			IBRION = 'conjugate-gradient'
		ISPIN = File.getLine(data_INCAR,'ISPIN')[0].split('=')[-1]
		if ISPIN =='1' or ISPIN == 'Not Match':
			ISPIN = 'non spin'
			spinNum = 9999
		elif ISPIN == '2':
			ISPIN = File.getLine(data_INCAR,'MAGMOM')[0].split('=')[-1]
			kind = ISPIN.split()
			spinNum = 0
			for kl in range(len(kind)):#kind举例：['23*5','23*-5','1']
				spinNum = spinNum + eval(kind[kl].split('*')[0])
				
		EDIFF = File.getLine(data_INCAR,'EDIFF')[0].split('=')[-1]
		EDIFFG = File.getLine(data_INCAR,'EDIFFG')[0].split('=')[-1]
		
		#检查POSCAR
		if SYSTEM == 'NEB':
			data_POSCAR = File.openFile('./00/POSCAR','r')
		else:
			data_POSCAR = File.openFile('./POSCAR','r')
		elelist = data_POSCAR[5].split()
		numlist = data_POSCAR[6].split()
		number = listSum(numlist)
		POS_out = ''
		for i in range(len(elelist)):
			POS_out = POS_out + elelist[i]+numlist[i]+' '
		
		#POTCAR
		data_POT = File.openFile('./POTCAR','r')
		PAW = File.getAllline(data_POT,'PAW')
		TITEL = File.getAllline(PAW,'TITEL')
		POTelelist = []
		POT_out = ''
		for tl in range(len(TITEL)):
			POTelelist.append(TITEL[tl].split()[3])
			POT_out = POT_out + TITEL[tl].split()[3] +' '
		
		#POINTS
		data_KP = File.openFile('./KPOINTS','r')
		method = data_KP[2].strip('\n')
		grid = data_KP[3].strip('\n')
	
		#Vasp.sh
		data_Sh  = File.openFile('./Vasp.sh','r')
		jobname = File.getLine(data_Sh,'#PBS -N')[0].strip('\n').split()[-1]
		jobnodes = File.getLine(data_Sh,'#PBS -l')[0].strip('\n').split()[-1]
		EMER = File.getLine(data_Sh,'#PBS -q')[0].strip('\n').split()[-1]
		if EMER == 'Match':
			EMER = 'No'
		elif EMER == 'EMERGENCY':
			EMER = 'Yes'

		print('Vasp前检查:\n路径：{}\n计算任务:{}  截断能:{}  ISIF:{}  离子更新:{}  磁性:{}  电子收敛:{}  离子收敛:{}'.format(os.getcwd(),SYSTEM,ENCUT,ISIF,IBRION,ISPIN,EDIFF,EDIFFG))
		print('POSCAR原子:{}种共计{}个 {}  POTCAR原子:{}  KPOINTS方法:{} 网格:{}  任务名:{}  节点与核:{}  加急:{}'.format(len(elelist),number,POS_out,POT_out,method,grid,jobname,jobnodes,EMER))
		
		
		if POTelelist == elelist or spinNum == number or spinNum == 9999:
			return True
		else:
			return False

	def check_and_qsub(need_input=True):
		'''
		need_input: 是否需要输入提交任务。默认是需要，设置为False不需要，检查完成后无错误直接提交
		:return:
		'''
		if Vasp.checkInputs():
			if need_input:
				ipt = input(f'是否提交任务(默认yes)：')
				if ipt in ['yes', 'y', '', 'Yes', 'YES', 'Y']:
					print('已提交任务')
					os.system('qsub Vasp.sh')
				else:
					print('未提交任务')
			else:
				print('已提交任务')
				os.system('qsub Vasp.sh')
		else:
			print(f'Vasp输入文件有误，请检查, path: {os.getcwd()}')
			exit()

	def checkJobstatus(name,includeE=True):
		'''
		name可以是jobname或者jobid或者qstat -x -u后该行的任意信息。
		只会查询那些状态为R Q E 的任务，F的不管。
		'''
		mjb = []
		mjb = getshellResult('qstat -x -u '+user)
		for i in range(len(mjb)):
			if name in mjb[i]:
				status = mjb[i].split()[9]
				runtime = mjb[i].split()[10]
				if status == 'R' or status == 'Q':
					return status,runtime
				if includeE == True:
					if status == 'E':
						return status,runtime
		return 'Not Found', '0'

	def keepInputs(addfile=[],workdir='./'):
		Flist = 'INCAR,POSCAR,POTCAR,KPOINTS,Vasp.sh'.split(',')
		for File in addfile:
			Flist.append(File)
		print(Flist)
		
		dir = 'KeepInputsDir'
		os.system('mkdir '+dir)
		for File in Flist:
			os.system('cp '+File+' ./'+dir)
		os.system('rm *')
		os.system('cp ./'+dir+'/* ./')
		os.system('rm -rf '+dir)
			
	def checkNEBperiod():
		'''
		遍历当前路径下的所有文件夹，如果发现有neb计算，判断ini 和fin分别的计算周期，并返回
		return: 列表，每个NEB作为一个元素，每个元素也是一个列表，列表下元素分别为NEB的路径、NEB阶段和状态，初态阶段和状态(在计算、未提交、算完)，
				[[NEB1path,NEBperiod,iniperiod,,finperiod],...]    3/3 2/3 1/3 Done NotDone 
		'''
		result = []
		currentpath = os.getcwd()
		#print(currentpath)
		for dirpath,dirnames,filenames in os.walk('./'):
			if ('ini' in dirnames) and ('fin' in dirnames):
				NEBpath = currentpath+dirpath.strip('.')  #返回的路径是完整路径
				#3/3
				if 'log' in filenames:
					#print(NEBpath)
					if Vasp.checkisDone(dirpath):
						NEBperiod = '3/3 Done'
					else:
						NEBperiod = '3/3 NotDone'
					iniperiod = '3/3 Done'
					finperiod = '3/3 Done'
				else: #可能是2/3 或者1/3
					#开始判断下一级
					NEBperiod = 'Not Match'
					if 'log' in os.listdir(dirpath+'/ini'):
						#print('in ini')
						if Vasp.checkisDone(dirpath+'/ini',isSelf=True):
							iniperiod = '2/3 Done'
						else:
							iniperiod = '2/3 NotDone'
					else:
						#print('ini without log')
						if 'log' in os.listdir(dirpath+'/ini/Opt'):
							if Vasp.checkisDone(dirpath+'/ini/Opt'):
								iniperiod = '1/3 Done'
							else:
								iniperiod = '1/3 NotDone'
						else:
							iniperiod ='0/3 NotInit'
					
					if 'log' in os.listdir(dirpath+'/fin'):
						if Vasp.checkisDone(dirpath+'/fin',isSelf=True):
							finperiod = '2/3 Done'
						else:
							finperiod = '2/3 NotDone'
					else:
						if 'log' in os.listdir(dirpath+'/fin/Opt'):
							if Vasp.checkisDone(dirpath+'/fin/Opt'):
								finperiod = '1/3 Done'
							else:
								finperiod = '1/3 NotDone'
						else:
							finperiod = '0/3 NotInit'
				#print(NEBpath)
				#print(NEBperiod)
				#print(iniperiod)
				#print(finperiod)
				result.append([NEBpath,NEBperiod,iniperiod,finperiod])

		return result
	
	def checkisDone(path,isSelf=False):
		'''
		传入路径和是否为自洽计算，通过Log判断该计算是否完成，返回True或False
		path:传入路径
		'''
		data_log = File.openFile(path+'/log')
		#print(data_log)
		#print(len(data_log))
		#有时候会出现log中无数据，报错
		if len(data_log) == 0:
			print('当前路径{}的log中无数据，请检查'.format(os.getcwd()+path.strip('.')))

		if isSelf==True:
			F = File.getLine(data_log,'F=')[0]
			#print('F',F)
			if '1 F=' in F:
				return True
			else:
				return False
		else:
			for i in range(len(data_log)-10,len(data_log)):
				if 'reached required accuracy' in data_log[i]:
					return True
			return False

	def decode_POSCAR(POSCAR):
		"""
		解码POSCAR，返回一个基矢、原子种类、原子数目、每个原子的位置（取前4位）
		:param POSCAR:
		:return:
		"""
		scale = float(POSCAR[1])  # 缩放系数
		a = np.array([float(tmp) for tmp in POSCAR[2].split()]) * scale  # 基矢
		b = np.array([float(tmp) for tmp in POSCAR[3].split()]) * scale  # 基矢
		c = np.array([float(tmp) for tmp in POSCAR[4].split()]) * scale  # 基矢
		vector = np.concatenate([a, b, c]).reshape(3, 3)
		elements = POSCAR[5].split()
		number_of_atom = np.array([int(tmp) for tmp in POSCAR[6].split()])
		direct, index = File.getLine(POSCAR, keywords='Direct')  # 获取Direct所在的索引，它后面的n行就是各个原子的位置。
		number_of_atom_sum = int(number_of_atom.sum())
		position = []
		for i in range(number_of_atom_sum):
			x, y, z = POSCAR[index+1+i].split()[:3]
			x, y, z = x[:6], y[:6], z[:6]
			position.append(np.array([float(x), float(y), float(z)]))
		position = np.array(position)

		# print(vector, vector.shape)
		# print(elements, number_of_atom, number_of_atom_sum)
		# print(direct, index)
		# print(position, position.shape)

		return vector, elements, number_of_atom, position

	def modify_Vasp_sh(jobname, nodes, ppn):
		with open('./Vasp.sh', 'r') as f:
			data = f.readlines()
		new_data = []
		for line in data:
			if ' #PBS -N' in line:
				new_data.append(f' #PBS -N {jobname}\n')
			elif ' #PBS -l nodes' in line:
				new_data.append(f' #PBS -l nodes={nodes}:ppn={ppn}\n')
			else:
				new_data.append(line)
		with open('./Vasp.sh', 'w') as f:
			f.writelines(new_data)

	def gennerate_POTCAR(elements=None, pseudotype='PBE'):
		if elements is None:
			with open('POSCAR', 'r') as f:
				data = f.readlines()
			res = Vasp.decode_POSCAR(data)
			elements = res[1]
		path1 = os.path.join(os.environ['HOME'], 'PseudoPotential')
		path2 = pseudotype
		os.system('rm POTCAR')
		for i in range(len(elements)):
			path3 = elements[i]
			path = os.path.join(path1, path2, path3)
			if os.path.isfile(f'{path}/POTCAR.Z'):
				if i == 0:
					os.system(f'zcat {path}/POTCAR.Z >POTCAR')
				else:
					os.system(f'zcat {path}/POTCAR.Z >>POTCAR')
			elif os.path.isfile(f'{path}/POTCAR'):
				code = f'cat {path}/POTCAR >POTCAR' if i==0 else f'cat {path}/POTCAR >>POTCAR'
				os.system(code)
			else:
				print(f'gennerate POTCAR error, element {elements[i]} not found')
				exit()

	def modify_POSCAR_ele(oldele, new_ele):
		"""
		读取并修改当前POSCAR的元素
		:param new_ele:
		:return:
		"""
		with open('POSCAR', 'r') as f:
			data = f.readlines()
		new_data = []
		for line in data:
			if oldele in line:
				new_data.append(line.replace(oldele, new_ele))
			else:
				new_data.append(line)
		with open('POSCAR', 'w') as f:
			f.writelines(new_data)

	def modify_POSCAR_Selective_Dynamics(data, indexes):
		"""
		根据输入的数据和索引修改POSCAR，添加Selective Dynamics, 索引所在的位置设置为T T T, 其他位置设置为 F F F
		注意：indexes以POSCAR中一个原子所在位置为初始0
		:return:
		"""
		POSCAR_data = []
		direct, direct_index = File.getLine(data, keywords='Direct')
		decoded_data = Vasp.decode_POSCAR(data)
		number_of_atom = decoded_data[2]
		for i in range(len(data)):
			if i < direct_index:
				POSCAR_data.append(data[i])  # 前面的部分
			elif i == direct_index:  # Direct部分，要加一个Selective Dynamics
				POSCAR_data.append('Selective Dynamics\n')
				POSCAR_data.append(data[i])
			elif direct_index < i <= direct_index + np.sum(number_of_atom):  # 原子位置部分
				tmp_i = i - direct_index - 1
				if tmp_i not in indexes:
					POSCAR_data.append(data[i].strip('\n') + ' F F F\n')
				else:
					POSCAR_data.append(data[i].strip('\n') + ' T T T\n')
			else:  # 最后的部分
				POSCAR_data.append(data[i])
		return POSCAR_data

	def modify_INCAR_for_vibration_analysis():
		"""
		修改当前目录的INCAR为振动分析的INCAR并保存
		注意：在Opt基础上的INCAR进行修改
		:return:
		"""
		with open('INCAR', 'r') as f:
			data_INCAR = f.readlines()
		data_INCAR = File.substituteData(data_INCAR, keywords='SYSTEM', newline='SYSTEM=Vib\n')  # 修改
		data_INCAR = File.substituteData(data_INCAR, keywords='NSW', newline='NSW=1\n')  # 修改
		data_INCAR = File.substituteData(data_INCAR, keywords='POTIM', newline='POTIM=0.03\n')  # 修改
		data_INCAR = File.substituteData(data_INCAR, keywords='IBRION', newline='IBRION=5\n')  # 修改
		data_INCAR = File.substituteData(data_INCAR, keywords='NFREE', newline='\n')  # 删除
		data_INCAR.append('NFREE=2\n')  # 添加
		data_INCAR = File.substituteData(data_INCAR, keywords='ISYM', newline='\n')  # 删除
		data_INCAR.append('ISYM=0\n')  # 添加
		data_INCAR = File.substituteData(data_INCAR, keywords='PREC', newline='\n')  # 删除
		data_INCAR.append('PREC=Accurate\n')  # 添加
		data_INCAR = File.substituteData(data_INCAR, keywords='NPAR', newline='\n')  # 删除
		data_INCAR = File.substituteData(data_INCAR, keywords='NCORE', newline='\n')  # 删除
		with open('INCAR', 'w') as f:
			f.writelines(data_INCAR)

	def check_WARNING_and_Energy(path='.'):
		"""
		检查路径下的结构优化是否完成，有无WARNING，返回能量
		:return:
		"""
		if not Vasp.checkisDone(path):
			print(f'当前路径{os.getcwd()} Vasp计算没有完成，退出程序')
			exit()

		OUTCAR_file = os.path.join(path, 'OUTCAR')
		OUT_data = File.openFile(OUTCAR_file, 'r')
		flag = 0
		WARNING_list = []
		for i in range(len(OUT_data)):
			line = OUT_data[i]
			if 'WARNING' in line:
				flag = 1
				WARNING_list.append(line)

		log_file = os.path.join(path, 'log')
		log_data = File.openFile(log_file, 'r')
		for i in range(len(log_data)):
			line = log_data[i]
			if 'WARNING' in line:
				flag = 1
				WARNING_list.append(line)

		if flag == 1:
			print(f'路径{os.getcwd()} 有警告！！！ 警告内容：')
			for WARN in WARNING_list:
				print(WARN.strip('\n'))

		F_list = File.getAllline(log_data, keywords='F=')
		energy = F_list[-1].split('F=')[-1].split()[0]
		energy = float(energy)
		return energy




	








