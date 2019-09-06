#!/public/software/apps/python/3.7.1/bin/python3
# coding=utf-8

import os
import argparse
import subprocess


def getshellResult(code):
	res = subprocess.check_output(code, stderr=subprocess.STDOUT, shell=True)
	res = res.decode('utf-8')
	res = res.split('\n')
	return res


def run(code, usr):
	res = getshellResult(code)
	used_core = 0
	for i in range(len(res)):
		line = res[i]
		if usr in line and 'R' in line:
			job_core = line.split()[6]
			used_core += int(job_core)
	print(f'已使用核：{used_core}')


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('-cd', '--code', default='qstat -a', help='code to show the job queue')
	parser.add_argument('-usr', '--user', default='yangyuqi')
	args = parser.parse_args()
	code = args.code
	usr = args.user
	run(code, usr)