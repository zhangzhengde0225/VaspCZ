import os
import argparse
import VaspCZ.zzdlib as zzd


def VaspNEBCheckDist(POSorCONT):
	if '00' in os.listdir() and 'INCAR' in os.listdir():
		data_log = zzd.File.openFile('./INCAR')
		image = zzd.File.getLine(data_log,'IMAGES')[0].split('=')[-1].strip('\n')
		if int(image) <=9:
			# os.system('cp ini/CONTCAR 00/CONTCAR')
			# os.system('cp fin/CONTCAR 0'+str(int(image)+1)+'/CONTCAR')
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


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='manual to this script')
	parser.add_argument('--POSorCONT', type=str, default='POSCAR')
	args = parser.parse_args()
	VaspNEBCheckDist(args.POSorCONT)
