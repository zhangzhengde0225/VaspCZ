import os
import subprocess


def dd():

	print('mdmdmd')

res = subprocess.check_output(dd, stderr=subprocess.STDOUT, shell=True)
print(res)