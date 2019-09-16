import os,sys
import VaspCZ.zzdlib as zzd
import subprocess


python = sys.executable
current_py_folder = os.path.dirname(os.path.abspath(__file__))
VaspCZ_path = [
	os.path.dirname(current_py_folder) if 'sourcecode' in current_py_folder else current_py_folder][0] + '/sourcecode'

def gui_string(title, content, footnote=None, iptnote='', mode='select', isprint=False):
	"""
	生成GUI界面，默认文本宽度为62字符
	:param title: 标题
	:param content: 内容
	:param footnote: 注脚，如果有
	:param iptnote: 最后一行Input后面的提示字符串，可输入默认值
	:param mode: 等于select时内容格式为数字加功能，选择，等于其他时，内容格式为说明字符串
	:return:
	"""
	gui = f'|{"":=<60}|\n'  # 第一行
	gui += f'|{title:^60}|\n'  # 标题行
	gui += f'|{"":-<60}|\n'  # 第三行
	if mode == 'select':  # 选择性内容
		for key in list(content.keys())[1::]:  # 内容行
			cont = content[key]
			gui += f'|{"":<8}({key}) {cont:<48}|\n'
	else:  # 说明性内容
		for key in list(content.keys())[1::]:  # 内容行
			cont = content[key]
			len_list = [56, 60, 60, 60, 60, 60, 60, 60]
			for j in range(len(len_list)):
				line = cont[0: len_list[j]]
				cont = cont[len_list[j]::]
				if isprint:
					print('line:', line)
					print('cont:', cont)
				if j == 0:
					line = f'{line:<56}' if len(line) < 56 else line
				else:
					line = f'{line:<60}' if len(line) < 60 else line
				if line[:56] == f'{"":<56}':
					break
				elif line[-1] != ' ':  # 处理最后一个字符不是空的情况
					if isprint:
						print('linexx', line, len(line), line[-1])
					if line[-2] == ' ':
						cont = line[-1] + cont
						line = line[:-1]
					else:
						cont = ' ' if cont == '' else cont
						if cont[0] == ' ':
							pass
						else:
							cont = line[-1] + cont
							line = line[:-1] + '-'
				if j == 0:

					gui += f'|{"":<4}{line:<56}|\n'
				else:

					gui += f'|{line:<60}|\n'
	key = list(content.keys())[0]
	gui += f'|{"":<8}({key}) {content[key]:<48}|\n'
	if footnote is not None:  # 如果有注脚行
		gui += f'|{"":-<60}|\n'  # 注脚行前一行
		gui += f'|{footnote:>60}|\n'

	gui += f'|{"":=<60}|\n'  # 最后一行
	gui += f'Input{iptnote}:  '  # 输入行
	return gui


def zip_content(content):
	return dict(zip(range(len(content)), content))


def deal_with_gen_pot():
	while True:
		elements = zzd.Vasp.decode_POSCAR(zzd.File.openFile('POSCAR', 'r'))[1]
		path = zzd.File.Vasp_pseudo_path()+'/PseudoPotential'
		content_os_gen_pot = zip_content([
			'Exit',
			f'The POTCAR will be generated according to the elements in POSCAR (current {elements}) from path: "{path}", default type is "PBE".',
			f'Change settings by input like: [ele1, ele2] PBE'])
		ipt11 = input(gui_string(
			title='Generate POTCAR', content=content_os_gen_pot, mode='string'))
		if ipt11 == '0':
			break
		elif ipt11 == '':
			elements = elements
			psudotype = 'PBE'
		else:
			try:
				elements = ipt11.split(']')[0].strip('[').split(',')
				psudotype = ipt11.split(']')[1].strip()
			except Exception as e:
				raise NameError(f'deal_with_gen_pot: 输入为{ipt11}. 输入格式错误。')
		psudotype = 'PBE' if psudotype == '' else psudotype
		print(f'生成元素：{elements}的贋势，类型：{psudotype}')
		zzd.Vasp.generate_POTCAR(elements=elements, pseudotype=psudotype)
		break


def deal_with_gen_kpoints():
	while True:
		content_os_gen_kp = zip_content([
			'Exit',
			f'The KPOINTS will be generated in vector: "5 5 5" with Monkhorst type.',
			f'Change settings by input like: 5 5 5 M'])
		ipt = input(gui_string(
			title='Generate KPOINTS', content=content_os_gen_kp, mode='string'))
		if ipt == '0':
			break
		elif ipt == '':
			vector = '5 5 5'
			kptype = 'Monkhorst'
		else:
			try:
				vector = f'{ipt.split()[0]} {ipt.split()[1]} {ipt.split()[2]}'
				kptype = ipt.split()[3]
			except Exception as e:
				raise NameError(f'deal_with_gen_pot: 输入为{ipt}. 输入格式错误。')
		kptype = 'Monkhorst' if kptype[0] == 'M' else kptype
		kptype = 'Gamma' if kptype[0] == 'G' else kptype
		print(f'生成网格："{vector}" 的贋势，类型：{kptype}')
		zzd.Vasp.generate_KPOINTS(vector=vector, kptype=kptype)
		break

def deal_with_gen_vasp_sh():
	path = zzd.File.Vaspsh_path()+'/Vasp.sh'
	content = zip_content([
		'Exit',
		f'The Vasp.sh file will be generated from template: {path}',
		f'{"Default nodes:":<20}{"1":>10}',
		f'{"Default ppn:":<20}{"12":>10}',
		f'{"Default job name:":<20}{"vaspjob":>10}',
		f'Change settings by input like: 1 12 jobname'
	])
	ipt = input(gui_string(title='Generate Vasp.sh', content=content, mode='string'))
	if ipt == '0':
		return None
	elif ipt == '':
		nodes, ppn, jobname = ('1', '12', 'vaspjob')
	else:
		try:
			nodes, ppn, jobname = ipt.split()
		except Exception as e:
			raise NameError(f'{e} deal_with_Vasp_sh error, 输入错误')
	print(f'生成Vasp.sh nodes: {nodes} ppn: {ppn} jobname: {jobname}')
	os.system(f'cp {path} .')
	data = zzd.File.openFile('Vasp.sh', 'r')
	data = zzd.File.substituteData(data=data, keywords='#PBS -N', newline=f' #PBS -N {jobname}\n')
	data = zzd.File.substituteData(data=data, keywords='#PBS -l nodes', newline=f' #PBS -l nodes={nodes}:ppn={ppn}\n')
	zzd.File.openFile('Vasp.sh', 'w', data=data)


def deal_with_vasp_keep_inputs():
	content = zip_content([
		'Exit',
		'The INCAR, POSCAR, POTCAR, KPOINTS and Vasp.sh will be kept in current directory while other files will be removed',
		'Add files need to keep input like: file1, file2'
	])
	ipt = input(gui_string(title='Vasp Keep Inputs', content=content, mode='string'))
	if ipt == '0':
		return None
	elif ipt == '':
		addfile = []
	else:
		addfile = ipt.strip().split(',')
	zzd.Vasp.keepInputs(addfile=addfile, workdir='./', need_confirm=True)

def deal_with_check_results():
	content = zip_content([
		'Exit',
		'Only Job Status',
		'Current Results',
		'All'
	])
	ipt = input(gui_string('Vasp Check Results', content=content, mode='select', iptnote='(default=3)'))
	if ipt == '0':
		return None
	elif ipt == '':
		ipt = '3'
	print(VaspCZ_path)
	subprocess.call(f'{python} {VaspCZ_path}/VaspCheckResults.py -func={ipt}', shell=True)

def deal_with_neb_opt_sta():
	content = zip_content([
		'Back',
		'Opt-Sta (Current directory)',
		'INI FIN Opt-Sta',
	])
	while True:
		ipt = input(gui_string('NEB Optimization-Static Calculation', content=content, iptnote='(default=1)'))
		if ipt == '0':
			break
		elif ipt == '1' or ipt == '':
			subprocess.call(f'{python} {VaspCZ_path}/VaspOpt-Sta.py', shell=True)
			exit()
		elif ipt == '2':
			data_Vaspsh = zzd.File.openFile('./ini/Opt/Vasp.sh', 'r')
			d_jobname = zzd.File.getLine(data_Vaspsh, '#PBS -N')[0].strip('\n').split()[-1]
			d_jobname = f'{d_jobname[:-1]}S'
			d_nodes = zzd.File.getLine(data_Vaspsh, '#PBS -l nodes')[0].strip('\n').split()[-1].split(':')[0].split('=')[-1]
			d_ppn = zzd.File.getLine(data_Vaspsh, '#PBS -l nodes')[0].strip('\n').split()[-1].split(':')[-1].split('=')[-1]
			content2 = zip_content([
				'Exit',
				'The static calculations in sub directories ini/ and fin/ will be performed when optimizations in ini/Opt and fin/Opt are finished.',
				'Default nodes and ppn from ini/Opt/Vasp.sh',
				f'Default nodes:               {d_nodes:>15}',
				f'Default ppn:                 {d_ppn:>15}',
				f'Default jobname:             {d_jobname:>15}',
				'Change settings by input like: nodes ppn jobname'
			])
			ipt2 = input(gui_string(
				title='Vasp NEB INI FIN Opt-Sta', content=content2, mode='string'))
			if ipt2 == '0':
				exit()
			elif ipt2 == '':
				nc = f'{d_nodes},{d_ppn}'
				jobname = d_jobname
			else:
				try:
					nc = f'{ipt2.split()[0]},{ipt2.split()[1]}'
					jobname = ipt2.split()[2]
				except Exception as e:
					raise NameError(f'{e} deal_with_neb_opt_sta error, 输入错误')
			subprocess.call(f'{python} {VaspCZ_path}/VaspINFIOpt-Sta.py --nc={nc} --jobname={jobname}', shell=True)
			exit()
		else:
			pass

def deal_with_neb_sta_neb():
	data_Vaspsh = zzd.File.openFile('./ini/Opt/Vasp.sh', 'r')
	d_jobname = zzd.File.getLine(data_Vaspsh, '#PBS -N')[0].strip('\n').split()[-1]
	d_jobname = f'{d_jobname[:-1]}N'
	d_nodes = zzd.File.getLine(data_Vaspsh, '#PBS -l nodes')[0].strip('\n').split()[-1].split(':')[0].split('=')[-1]
	d_ppn = zzd.File.getLine(data_Vaspsh, '#PBS -l nodes')[0].strip('\n').split()[-1].split(':')[-1].split('=')[-1]

	dist = zzd.getshellResult('dist.pl ini/CONTCAR fin/CONTCAR')
	dist = eval(dist[-1])
	print('ini和fin中CONTCAR的dist为：{}'.format(dist))
	if dist >= 9:
		print('dist过大，请检查')
		raise NameError(f'初态和末态距离太大，插入态数目大于9，不合理，请检查。')
	else:  # 向下取整数，如果是偶数则加一，如果是奇数直接用。0-1.9输入1，2-3.9属于3，4-5.9属于5
		image = int(dist / 0.8)
		if image % 2 == 0:  # 是偶数
			image = image + 1

	content = zip_content([
		'Exit',
		'The NEB calculation will be performed when static calculations in ini/ and fin/ are finished.',
		'The INCAR of NEB will be generated automatically fitted to NEB calculation based on ini/Opt/INCAR, The images is approximately equal to (distance between ini/CONTCAR and fin/CONTCAR)/0.8',
		'Default nodes and ppn from ini/Opt/Vasp.sh',
		f'Default nodes:                 {image:>15}',
		f'Defalut ppn:                   {d_ppn:>15}',
		f'Defalut jobname:               {d_jobname:>15}',
		'Change settings by input like: nodes ppn jobname'
	])
	ipt = input(gui_string(title='Vasp NEB Sta-NEB', content=content, mode='string'))
	if ipt == '0':
		return None
	elif ipt == '':
		nc = f'{image},{d_ppn}'
		jobname = d_jobname
	else:
		try:
			nc = f'{ipt.split()[0]},{ipt.split()[1]}'
			jobname = ipt.split()[2]
		except Exception as e:
			raise NameError(f'{e} deal_with_neb_sta_neb error, 输入错误')
	subprocess.call(f'{python} {VaspCZ_path}/VaspINFISta-NEB.py  --nc={nc} --jobname={jobname}', shell=True)


def deal_with_neb_vibration_analysis():
	content = zip_content([
		'Exit',
		'The vibration analysis will be preformed after NEB calculation is Done.',
		'The attempt frequency of migration atom in initial, saddle and finnal state will be calculated.',
		'Default nodes:          1',
		'Default ppn:            8',
		f'Defacult cal fin:  False (do not include finnal state)',
		'Change settings by input like: nodes ppn True(or False)'
	])
	ipt = input(gui_string(
		title='NEB Vibration Analysis',
		content=content, mode='string'
	))
	if ipt == '0':
		return None
	elif ipt == '':
		nodes, ppn, include = ('1', '8', 'False')
	else:
		try:
			nodes, ppn, include = ipt.split()[0], ipt.split()[1], ipt.split()[2]
		except Exception as e:
			raise NameError(f'{e} deal_with_neb_vibration_analysis error, 输入错误')
	subprocess.call(
		f'{python} {VaspCZ_path}/VaspVibAna_forNEB.py --nodes={nodes} --ppn={ppn} --include_fin={include}', shell=True)


def deal_with_neb_check_results():
	content = zip_content([
		'Back',
		'Only NEB Status',
		'NEB Results',
		'All',
		'NEB Barrier in typical step'
	])
	while True:
		ipt = input(gui_string(title='NEB Check Results', content=content, mode='select', iptnote="(1/2/[3]/4)"))
		if ipt == '0':
			break
		elif ipt == '':
			ipt = '3'
			subprocess.call(f'{python} {VaspCZ_path}/NEBCheck1.1.py --func={ipt}', shell=True)
			exit()
		elif ipt in ['1', '2', '3', '4']:
			subprocess.call(f'{python} {VaspCZ_path}/NEBCheck1.1.py --func={ipt}', shell=True)
			exit()
		else:
			print(f'输入：{ipt}，选择功能错误，请正确输入')
			pass


def deal_with_test_kpoints():
	content = zip_content([
		'Exit',
		'The KPOINT test will be performed when the input files in current dir (INCAR, POSCAR, POTCAR, KPOINTS)',
		'Default setting jobname prefix:      ktest_',
		'Default setting nodes:                    1',
		'Default setting ppn:                      8',
		'Default setting k_mesh: 111,333,555,777,999',
		'Change settings by input like: ktest_ 1 8 111,333,555,777,999'
	])
	ipt = input(gui_string(title='Vasp KPOINTS Test', content=content, mode='string'))
	if ipt == '0':
		return None
	elif ipt == '':
		prefix, nodes, ppn, k_mesh = ('ktest_', '1', '8', '111,333,555,777,999')
	else:
		try:
			prefix, nodes, ppn, k_mesh = ipt.split()[0], ipt.split()[1], ipt.split()[2], ipt.split()[3]
		except Exception as e:
			raise NameError(f'{e} deal_with_test_kpoints error, 输入错误')
	subprocess.call(f'{python} {VaspCZ_path}/VaspTestKPOINT.py --jobname_prefix={prefix} --nodes={nodes} --ppn={ppn} --k_mesh={k_mesh}', shell=True)


def deal_with_test_encut():
	content = zip_content([
		'Exit',
		'The ENCUT test will be performed when the five input in current dir (INCAR, POSCAR, POTCAR, KPOINTS)',
		'Default setting jobname prefix:      ENCUT_',
		'Default setting nodes:                    1',
		'Default setting ppn:                      8',
		'Default setting ENCUTs: 200,250,300,350,400,450,500,550,600,650,700',
		'Change settings by input like: ENCUT_ 1 8 200,250,300,350,400,450,500,550,600,650,700'
	])
	ipt = input(gui_string(title='Vasp ENCUT Test', content=content, mode='string'))
	default_params = ('k_test', '1', '8', '200,250,300,350,400,450,500,550,600,650,700')
	if ipt == '0':
		return None
	elif ipt == '':
		params = default_params
	else:
		params = [ipt.split(i) for i in range(4) if i < len(ipt.split())]  # 获取前几位输入参数
		params += [default_params[i] for i in range(len(ipt.split(), 4))]  # 获取默认参数
		try:
			a = eval(params[1])
			a = eval(params[2])
		except Exception as e:
			raise NameError(f'{e} deal_with_test_encut error, 输入错误')
	prefix, nodes, ppn, ENCUTs = params
	subprocess.call(
		f'{python} {VaspCZ_path}/VaspTestENCUT.py -jb={prefix} -nd={nodes} -np={ppn} -EN={ENCUTs}', shell=True)
