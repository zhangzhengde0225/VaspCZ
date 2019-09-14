import os
import VaspCZ.zzdlib as zzd


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


if __name__ == '__main__':
	VaspNEBCheckRMS()
