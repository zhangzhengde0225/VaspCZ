# written by Michael Waters

def interpolate_images(image_list, num_new_images, kind = 'linear', use_image_distance_in_spline = False):
	''' This function interpolates a list of ASE "Atoms" to a new list of images'''
	nimages = len(image_list)
	natoms = atoms.positions.shape[0]

	if nimages == 2:
		print ('Only 2 images, kind will be linear')
		kind = 'linear'
	elif nimages < 2:
		print('YOU NEED AT LEAST 2 IMAGES FOR INTERPOLATION!')

	from ase.geometry.geometry import find_mic
	from numpy import zeros, linspace, sqrt
	from scipy.interpolate import interp1d



	####################
	distance_seq = zeros(nimages)
	position_collection = zeros((natoms,3,nimages))
	image_index = 0 # for the first image, we don't need distances, just the orginal positions
	for atom_index in range(0,natoms):
		for dim_index in range(3):
			position_collection[atom_index, dim_index, image_index] = image_list[image_index].positions[atom_index, dim_index]

	for image_index in range(1,nimages):
		D = image_list[image_index].positions - image_list[image_index-1].positions
		D_min, D_min_len =  find_mic(D, cell =  image_list[image_index].get_cell() )
		# D_min is list of minimum image vectors
		distance =  sqrt((D_min**2).sum())
		distance_seq[image_index] = distance_seq[image_index-1] + distance
		for atom_index in range(0,natoms):
			for dim_index in range(3):
				position_collection[atom_index, dim_index, image_index] = \
					position_collection[atom_index, dim_index, image_index-1] + D_min[atom_index][dim_index]

	seq = linspace(0,1, nimages)
	# splines have a coordinate output and input that is image number scaled from 0 to 1.
	# we could try using the RMS distance/Frobenius distance/L2 norm along the path then scale it:
	if use_image_distance_in_spline:
		seq = distance_seq/distance_seq.max()

	# builds a spline for every atom's x,y,z coordinates
	spline_func_collection = []
	for atom_index in range(0,natoms):
		spline_func_collection.append([])
		for dim_index in range(3):
			func = interp1d(seq, position_collection[atom_index,dim_index] ,kind=kind)
			spline_func_collection[atom_index].append(func)


	################
	mag_collection = zeros((natoms,nimages))
	for image_index in range(0,nimages):
		mag_mom = image_list[image_index].get_initial_magnetic_moments()
		#print(mag_mom)
		for atom_index in range(0,natoms):
			mag_collection[atom_index, image_index ] = mag_mom[atom_index]

	#print(mag_collection)
	mag_spline_func_collection = []
	for atom_index in range(0,natoms):
		func = interp1d(seq, mag_collection[atom_index], kind=kind)
		mag_spline_func_collection.append(func)

	#############################
	from copy import deepcopy

	new_image_list = []
	new_seq = linspace(0,1,  num_new_images )
	new_mag_mom = zeros(natoms)
	for new_image_index in range(0, num_new_images):

		new_image = deepcopy(image_list[0])
		# I'm initializing new 'Atoms' objects with deepcopy, there must be a better
		# way which will also handle lattice vector changes
		pos = new_seq[new_image_index]

		for atom_index in range(0,natoms):
			for dim_index in range(3):
				new_image.positions[atom_index, dim_index] = spline_func_collection[atom_index][dim_index](pos)

		for atom_index in range(0,natoms):
			new_mag_mom[atom_index] = mag_spline_func_collection[atom_index](pos)
		new_image.set_initial_magnetic_moments(new_mag_mom)
		new_image_list.append(new_image)

	return new_image_list


def rms_distance(imageA,imageB):
	from numpy import sqrt
	from ase.geometry.geometry import find_mic
	D = imageB.positions-imageA.positions # 2d arrays
	D_min, D_min_len = find_mic( D, imageB.cell )
	distance =  sqrt((D_min**2).sum())
	return distance

def compute_image_rms_distances(image_list):
	from numpy import zeros
	distances = zeros(len(image_list)-1)
	for i in range(0, len(image_list)-1):
		distances[i] = rms_distance(image_list[i+1], image_list[i])
	return distances
############## These functions are meant to make working with VASP easier


def try_mkdir(direct):
	from os import mkdir
	from os.path import isdir
	if isdir(direct) == False:
		mkdir(direct)



def get_nimages(directory = ''):
	from os.path import isfile
	images = 1
	while isfile(directory +"%02d/CONTCAR"%images):
	    images+=1
	images-=2

	print(images,"Images Found")
	return images

def read_mag_cols(fname='OUTCAR'):

	fid = open(fname,'r')
	lines = fid.readlines()

	mag_line = -1

	for i in range(len(lines)):
		if " magnetization (x)" in lines[i]:
			mag_line = i

	#print (lines[mag_line+4:mag_line+4+n_atoms])
	mag_cols = [[],[],[],[],[]]
	line_index = mag_line+4
	while '---' not in lines[line_index]:

		#sline = line.split()
		sline = lines[line_index].split()
		mag_cols[0].append(int(sline[0]))
		for icol in range(1,5):
			mag_cols[icol].append(float(sline[icol]))

		line_index += 1
	#the last column, mag_cols[-1] has the total magnetic moment
	return mag_cols

###################### Test the function here

if __name__=='__main__':


	num_new_images = 5# this number matches the IMAGES tag in VASP
	use_image_distance_in_spline = True

	# handy function for getting the number of images in VASP format in this directory
	nimages= get_nimages()


	from ase import io
	from numpy import array
	image_list =[]
	for image in range(0,nimages+2):

		atoms = io.read('%02d/CONTCAR'%image)
		mag_cols = read_mag_cols('%02d/OUTCAR'%image)
		atoms.set_initial_magnetic_moments(mag_cols[-1])

		image_list.append(atoms)




	### now that the images are read, we can use the interpolating function
	# the +2 is because vasp doesn't count the first and last images in the IMAGES tag
	interpolated_image_list =  interpolate_images(image_list,  num_new_images+2,
					kind = 'cubic', use_image_distance_in_spline = use_image_distance_in_spline)

	if use_image_distance_in_spline:
		print('Compare image spacing before and after interpolation:')
		print(compute_image_rms_distances(image_list))
		print(compute_image_rms_distances(interpolated_image_list))



	# with the interpolated images, we can write them to a subdirectory
	sub_dir = 'interpolated_images/'
	try_mkdir(sub_dir)
	for new_image_index in range(0, num_new_images+2):

		imdir = sub_dir+'%02d/'%new_image_index
		try_mkdir(imdir)

		atoms = interpolated_image_list[new_image_index]

		fname = imdir+'POSCAR'
		io.write(fname, atoms, format='vasp')
		#fname = sub_dir+imdir+'CONTCAR'
		#io.write(fname, atoms, format='vasp')

		######## this part makes a MAGMOM line for our INCAR file
		magmom_name = imdir+'MAGMOM'
		mag_mom = atoms.get_initial_magnetic_moments()

		fid = open(magmom_name,'w')
		fid.write('MAGMOM =')

		for atom_index in range(mag_mom.shape[0]):
			fid.write(' %.2f'%( mag_mom[atom_index] ))
		fid.close()
