import setuptools

__version__ = '1.0.2'

with open('README.md', 'r') as f:
	long_description = f.read()

setuptools.setup(
	name='VaspCZ',
	version=__version__,
	description='Vasp Check by Zzd',
	long_description=long_description,
	author='Zhengde Zhang',
	author_email='drivener@163.com',
	url='https://github.com/zhangzhengde0225/VaspCZ',
	keywords='django markdown editor editormd',
	packages=setuptools.find_packages(),
	zip_safe=False,
	include_package_data=True,
	install_requires=[
		'numpy'
	],
	classifiers=(
		"Programming Language :: Python",
		"Development Status :: 4 - Beta",
		"Environment :: Console",
		"Framework :: Django",
		"License :: OSI Approved :: MIT License",
		"Natural Language :: Chinese (Simplified)",
		"Operating System :: POSIX :: Linux",
		"Programming Language :: Python :: 3.6"
	)
)