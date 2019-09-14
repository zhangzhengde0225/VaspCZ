import os
import VaspCZ.zzdlib as zzd



def generate_inputs(examp='fcc_Fe_3x3x3'):
	"""
	从example下拷贝例子
	:return:
	"""
	path = zzd.File.VaspCZ_software_path()
	for file in 'INCAR,POSCAR,POTCAR,KPOINTS'.split(','):
		os.system(f'cp {path}/examples/{examp}/{file} .')
		print(f'生成文件：{file} (example)')


def generate_INCAR_for_Sta():
	print(f'修改当前结构优化INCAR为静态计算INCAR')
	try:
		data_INCAR = zzd.File.openFile('INCAR', 'r')
	except Exception as e:
		raise NameError(f'{e} 当前路径无INCAR文件')
	for nl in range(len(data_INCAR)):
		if 'SYSTEM' in data_INCAR[nl]:
			data_INCAR[nl] = 'SYSTEM=Static\n'  # 修改表头
		if 'NSW' in data_INCAR[nl]:
			data_INCAR[nl] = 'NSW=1\n'  # 修改NSW
		if 'IBRION' in data_INCAR[nl]:
			data_INCAR[nl] = 'IBRION=-1\n'  # 修改IBRION
		if 'EDIFFG' in data_INCAR[nl]:
			data_INCAR[nl] = '#' + data_INCAR[nl]  # 去掉EDIFFG
	zzd.File.openFile('INCAR', 'w', data=data_INCAR)
	print(f'修改为静态计算INCAR完成')

def run():
	a = zzd.File.VaspCZ_software_path()
	print(a)


if __name__ == '__main__':
	run()
