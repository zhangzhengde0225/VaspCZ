import os
import VaspCZ.zzdlib as zzd
import argparse

def CheckWARNING(usript):
	# usrinput为1时，只打印是否完成，为2时候，只打印能量信息，为3时，打印全部
	Path = []
	Energy = []
	Step = []
	mag = []
	dist = []
	RMS = []
	warnflag = 0
	jstat = ['Running...', 'Done', 'SC Running...', 'SC done', 'Stoped']
	for dirpath, dirnames, filenames in os.walk('./'):
		if 'INCAR' in filenames and 'POSCAR' in filenames and 'Vasp.sh' in filenames:
			data_Sh = zzd.File.openFile(dirpath + '/Vasp.sh', 'r')
			data_INCAR = zzd.File.openFile(dirpath + '/INCAR', 'r')
			SYSTEM = zzd.File.getLine(data_INCAR, 'SYSTEM')[0].split('=')[-1]
			jobname = zzd.File.getLine(data_Sh, '#PBS -N')[0].strip('\n').split()[-1]
			# jobstatus = zzd.Vasp.checkJobstatus(jobname)
			# # print('Path:{}  {}   {}'.format(dirpath,jobstatus,jobname))
			# if jobstatus == 'Q':
			# 	print('Path:{:<60}任务正在排队...   '.format(dirpath))
			# elif jobstatus == 'R' or 'log' in filenames:  # 正在计算和算完了用相同的检测方法
			if 'log' in filenames:
				if usript == '1' or usript == '3':
					data_log = zzd.File.openFile(dirpath + '/log', 'r')
					reachRA = zzd.File.getLine(data_log, 'reached required accuracy')[0]
					try:
						RMM = zzd.File.getAllline(data_log, 'RMM:')[-1].split()[1]
						RMM = RMM + ' RMM'
					except:
						RMM = '0 RMM'
					try:
						DAV = zzd.File.getAllline(data_log, 'DAV:')[-1].split()[1]
						DAV = DAV + ' DAV'
					except:
						DAV = '0 DAV'
					termination = zzd.File.getLine(data_log, 'Ctrl-C caught... cleaning up processes')[0]
					try:
						ionstep = zzd.File.getAllline(data_log, 'F=')[-1].split()[0]
					except:
						ionstep = '0'
					# print('Path:{}  {}   {}'.format(dirpath,reachRA,v1F))
					if 'reached required accuracy' in reachRA:
						print('Path:{:<60}{:<15}{:>3}F{:>9}'.format(dirpath, jstat[1], ionstep, RMM))
					elif SYSTEM == 'Self' or SYSTEM == 'Static' or SYSTEM == 'self' or SYSTEM == 'static':
						vRMM = zzd.File.getLine(data_log, 'RMM:')[0]
						vDAV = zzd.File.getLine(data_log, 'DAV:')[0]
						if ('DAV:' in vDAV or 'RMM:' in vRMM) and ionstep == '0':
							print('Path:{:<60}{:<15}{:>3}F{:>9}{:>7}{:>5}'.format(dirpath, jstat[2], ionstep, RMM, DAV,
																				  '*'))
						elif ionstep == '1':  # 没有完成，但有1F了
							print('Path:{:<60}{:<15}{:>3}F{:>9}'.format(dirpath, jstat[3], ionstep, RMM))
					elif 'cleaning up processes' in termination:
						print('Path:{:<60}{:<15}{:>3}F{:>9}{:>7}'.format(dirpath, jstat[4], ionstep, RMM, DAV))
					else:
						print(
							'Path:{:<60}{:<15}{:>3}F{:>9}{:>7}{:>5}'.format(dirpath, jstat[0], ionstep, RMM, DAV, '*'))

					# 检查警告
					try:
						data_OUTCAR = zzd.File.openFile(dirpath + '/OUTCAR', 'r')
						WARNING_log = zzd.File.getAllline(data_log, 'WARNING')
						WARNING_OUT = zzd.File.getAllline(data_OUTCAR, 'WARNING')
						ERROR_log = zzd.File.getAllline(data_log, 'ERROR')
						if WARNING_log != [] or WARNING_OUT != []:
							print('Path:{:<60}    出现警告'.format(dirpath, ionstep))
							zzd.File.printData(WARNING_log)
							zzd.File.printData(WARNING_OUT)
							warnflag = 1
						if ERROR_log != []:
							zzd.File.printData(ERROR_log)
							warnflag = 1
					except Exception as e:
						print('path:{} {}'.format(dirpath, e))

				if usript == '2' or usript == '3':
					# 打印能量信息
					data_log = zzd.File.openFile(dirpath + '/log', 'r')
					try:
						log_lastF = zzd.File.getAllline(data_log, 'F=')[-1]
					except:
						log_lastF = '0 F= 0 E0= 0  d E 0  mag=    0'  # 在存在log，Log中连1步都没算完的时候
					if len(dirpath) > 15:
						path = '...' + dirpath[-15:]
					else:
						path = dirpath
					Path.append(path)
					Energy.append(log_lastF.split()[4])
					Step.append(log_lastF.split()[0])
					try:
						mag.append(log_lastF.split()[9])
					except:
						mag.append('0')
					dist.append(
						zzd.getshellResult('dist.pl ' + dirpath + '/POSCAR ' + dirpath + '/CONTCAR')[-1].strip('\n'))
					try:
						RMS.append(zzd.getshellResult('grep RMS ' + dirpath + '/OUTCAR')[-1].split()[4])
					except:  # 如果出现没算完，grep返回一个空的列表的时候
						RMS.append('0')

	if (usript == '1' or usript == '3') and warnflag == 0:
		print('无警告或错误')
	if (usript == '2' or usript == '3'):
		print('{:<18}{:<11}{:<5}{:^10}{:^10}{:^10}'.format('路径', '能量', '步数', 'mag', 'dist', 'RMS'))
		for i in range(len(Path)):
			print('{:<20}{:<13.4f}{:<7}{:<12.4f}{:<10.4f}{:<10.4f}'.format(Path[i], eval(Energy[i]), Step[i],
																		   eval(mag[i]), eval(dist[i]), eval(RMS[i])))


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('-func', '--function', default='3', type=str)
	args = parser.parse_args()
	CheckWARNING(args.function)