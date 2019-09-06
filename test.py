import sys
# print(sys.modules)

import VaspCZ.zzdlib as zzd

b = zzd.File.openFile('install.py', 'r')
print(b)