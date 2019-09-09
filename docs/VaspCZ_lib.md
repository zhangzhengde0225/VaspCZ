# VaspCZ python API 文档
## 一.安装和使用
### 1. 安装
请查看项目的[README.md](https://github.com/zhangzhengde0225/VaspCZ)
### 2. 使用
```
import VaspCZ.zzdlib as zzd
```
### 3. 简单示例
1. 获取linux shell代码的返回值：
```
res = zzd.getshellResult('pwd')
print(res)
```
输出：
```angular2html
['/home/zhangzhengde\n']
```
返回shell命令的控制台输出的结果，由每一行组成一个元素的列表。
2. 读取文件内容
```angular2html
data = zzd.File.open('')
```
## 二. API文档
### 1. shell模块
```
VaspCZ.zzdlib.getshellResult(code)
```

返回shell命令控制台输出的结果，由每一行组成一个元素的列表。

+ **code** : str
    
    shell命令，如：ls, qstat等  
    
+ **return** : list
    
    返回列表，其中每个元素是控制台输出的一行，保留换行符\n
    
### 2. File模块
```angular2html
class VaspCZ.zzdlib.File()
```
集成文件处理功能。

------
```angular2html
VaspCZ.zzdlib.File.openFile(path, [mode='r', data=None])
```
读取文件或保存文件。小括号()内为必须参数，中括号[]内为可选参数。

+ **path** : str
    
    文件所在路径。
    
+ **mode** : str, r or w
    
    可选，模式，r是读取，w是写入。

+ **data** : list

    写入模式时必须参数，由每一行字符串作为一个元素组成的列表，每一行以换行符\n结尾。
   
+ **return** : list or None

    返回值。读取模式返回该文件每一行组成的列表，写入模式无返回

---
```
VaspCZ.zzdlib.File.substitudeData(data, keywords, newline, [mode='default'])
```

传入文件数据，给出关键词和新行，默认情形搜索出现第一次出现关键词的行并替换，mode不等于default是替换全部出现关键字的行，返回替换后的数据。

+ **data** : list

    文件每一行组成的列表
    
+ **keyword** : str

    搜索文件内容的关键词
    
+ **newline** : str

    替换该行的内容，注意需要以\n结尾
    
+ **mode** : str

    可选，模式，默认为default, 只替换第一次出现关键词的行，修改为其他任意字符串，替换所有出现关键词的行
    
+ **return** : str

    返回值。返回替换后的数据
   
例子：读取当前文件夹下的INCAR文件并将截断能设置改为400
```
import VaspCZ.zzdlib as zzd

data = zzd.File.openFile('INCAR', 'r')
newdata = zzd.File.substitudeData(data, 'ENCUT', 'ENCUT=400\n')
zzd.File.openFile('INCAR', 'w', data=newdata)
```

---------------

```
VaspCZ.zzdlib.File.getLine(data,keywords)
```

给出关键词，招傲有关键词的第一行并返回，返回为字符串和所在的行索引。该功能用于获取特定想信息或者用于判断。

+ **data** : list

    文件每一行组成的列表
    
+  **keyword** : str

    搜索关键词
    
+ **return** : tuple

    返回值。元组，第一个元素是该行字符串，不带换行符\n，第二个元素是该行所在的索引。如为找到匹配的关键词，则返回为('Not Match', 0)
  
------

```angular2html
VaspCZ.zzdlib.getAllline(data, keywords)
```

给出关键词，返回所有带有关键词的所有行，返回为列表。该功能用于选择性获得文件特定行。

+ **data** : list

    文件每一行组成的列表
   
+ **keyword** : str

    搜索关键词
   
+ **return** : list

    返回值。所有有关键词的行组成的列表，注意每个行的换行符\n保留。    
    
---

```angular2html
VaspCZ.zzdlib.getNullline(data)
```

获取文件数据中是空位的索引。

+ **data** : list

    文件每一行组成的列表
   
+ **data** : list

    返回值。所有空行存在的的索引组成的列表。
   
---

```angular2html
VaspCZ.zzdlib.Vaspsh_path()
```

获取VaspCZ软件默认的PBS提交任务脚本Vasp.sh所在的文件路径。

- **return** : str

    返回值。字符串。Vasp.sh文件路径，默认路径为用户根目录，如/home/zhangzhengde，安装软件时候可在install.py中更改位置。
   
### 3. Vasp模块
Vasp计算类。
```angular2html
class Vasp()
```
类名。

-----
3.1
```angular2html
VaspCZ.zzdlib.Vasp.decode_POSCAR(POSCAR):
```
解码POSCAR，返回一个基矢、原子种类、原子数目、每个原子的位置（取前4位）

- **POSCAR** : list

    由File模块读取的POSCAR文件数据。

- **return** : tuple

    返回值。4个元素的元组。第一个元素为晶格基矢，3\3的numpy数组。第二个元素为原子种类，由元素名称符串组成的列表。第三个元素为原子数目，由整数组成的一维numpy数组。第四个元素为原子位置，n*3的numpy数组，n是原子总数，每行是一个原子，一到三列分别为该原子x, y, z坐标。

例子：读取和解码当前目录的POSCAR文件
```angular2html
import VaspCZ.zzdlib as zzd

data = zzd.File.openFile('./POSCAR', 'r')  # 读取POSCAR文件数据
res = zzd.Vasp.decode_POSCAR(data)  # 解码POSCAR
print(f'{res[0]}\n{res[1]}\n{res[2]}\n{res[3]}')  # 打印解码后的数据
```
原POSCAR文件为：
```angular2html
Fe
1.0
  3.5239999294         0.0000000000         0.0000000000
  0.0000000000         3.5239999294         0.0000000000
  0.0000000000         0.0000000000         3.5239999294
  Fe
  4
Direct
0.000000000         0.000000000         0.000000000
0.000000000         0.500000000         0.500000000
0.500000000         0.000000000         0.500000000
0.500000000         0.500000000         0.000000000
```

输出为：
```angular2html
res[0]:
array([[3.52399993 0.         0.        ]
 [0.         3.52399993 0.        ]
 [0.         0.         3.52399993]])
 
res[1]:
['Fe']

res[2]:
array([4])

res[3]:
array([[0.  0.  0. ]
 [0.  0.5 0.5]
 [0.5 0.  0.5]
 [0.5 0.5 0. ]])
```

-----
3.2
```angular2html
VaspCZ.zzdlib.Vasp.modify_POSCAR_ele(oldele, new_ele):
```

修改当前路径下POSCAR的原子种类，适合批量修改。

- **oldele** : str
    
    原来的元素名称。
    
- **new_ele** : str

    新的元素名称。
    
- **return** : None

----
3.3
```angular2html
VaspCZ.zzdlib.Vasp.gennerate_POTCAR([elements=None, pseudotype='PBE']):
```

在当前路径生成POTCAR文件，需要在安装中正确是指贋势文件目录，默认贋势文件目录为用户根目录。贋势目录名为：PseudoPotential。

- **elements** : list

    可选。默认为None，从当前POSCAR文件读取，如输入由元素名称组成的列表，则按列表生成
    
- **pseudotype** : str

    可选。默认为PBE，从PseudoPotential/PBE文件夹下读取元素贋势。
    
-----
3.4
```angular2html
modify_POSCAR_Selective_Dynamics(data, indexes)
```
根据输入的数据和索引修改POSCAR，添加Selective Dynamics, 索引所在的位置设置为T T T, 其他位置设置为 F F F

- **data** : list

    由File模块读取的POSCAR文件数据
    
- **indexes** : list

    由索引组成的列表。
    注意：indexes以POSCAR中一个原子所在位置为初始0

- **return** : list

    返回值。修改后的POSCAR文件数据。

---
3.5
```angular2html
modify_INCAR_for_vibration_analysis():
```
修改当前目录的INCAR为振动分析的INCAR并保存。
修改内容包括：SYSTEM=Vib, NSW=1, POTIM=0.3, IBRION=5, NFREE=2, ISYM=0, PREC=Accurate。
注意：在Opt的INCAR上进行修改。

- **return** : None

-----
3.6
```angular2html
VaspCZ.zzdlib.Vasp.checkInputs()
```
Vasp前检查。提交计算任务前，检查当前目录Vasp的各项输入文件，将计算信息打印到控制台，包含：计算路径、SYSTEM、截断能、ISIF、离子更新方法、是否有磁性、电子收敛标准、离子收敛标准、原子种类个数、POTCAR原子类型、KPOINTS方法、网格大小、任务名、节点数与核数、是否加急。

并检查INCAR设置磁性时原子数目与POSCAR原子数目是否一致，POSCAR和POTCAR的原子种类是否一致。


- **return** : bool

    返回值。检查通过返回True, 未通过返回False

---
3.7
```angular2html
VaspCZ.zzdlib.Vasp.check_and_qsub([need_input=True]):
```
检查前检查并提交任务。内部集成了上一个检查输入文件函数，使用中推荐该函数。

- **need_input** : bool

    可选，默认为True，检查通过后会询问用户是否确认提交任务，确认后提交，若为False，检查通过后自动提交任务。
   
- **return** : None
    
注意：该函数需要保证PBS提交任务的脚本命名为Vasp.sh

例子：进入到examples/fcc_Fe_primitive文件夹，准备好输入文件INCAR, POSCAR, POTCAR, KPOINT和Vasp.sh后。
```angular2html
import VaspCZ.zzdlib as zzd

zzd.Vasp.check_and_qsub()
```
控制台输出检查结果并询问：
```angular2html
Vasp前检查:
路径：/Users/tanmenglu/PycharmProjects/VaspCZ/examples/fcc_Fe_primitive
计算任务:Opt  截断能:400  ISIF:3  离子更新:conjugate-gradient  磁性:non spin  电子收敛:1E-6  离子收敛:-0.01
POSCAR原子:1种共计4个 Fe4   POTCAR原子:Fe   KPOINTS方法:Monkhorst 网格:4 4 4  任务名:fcc111  节点与核:nodes=1:ppn=24  加急:batch
是否提交任务(默认yes)：
```
按回车即可提交任务，输入n或no，回车，不提交任务。

--------
3.8
```angular2html
VaspCZ.zzdlib.Vasp.keepInputs([addfile=[], workdir='./']):
```
删除工作目录下的文件，仅保留输入文件。默认保留文件为：INCAR，POSCAR，POTCAR， KPOINTS和Vasp.sh

- **addfile** : list

    可选，默认为空，输入额外需要保留的文件名。
  
- **workdir** : str

    可选，默认为当前目录。
      
- **return** : None

---
3.9
```angular2html
VaspCZ.zzdlib.Vasp.checkNEBperiod():
```

遍历当前路径下的所有文件夹，如果发现有neb计算，判断ini和fin分别的计算周期，并返回
		
- **retrun** : list

    返回值。列表，元素为内层列表，每一个NEB计算路径的结果作为内层列表，每个内层列表的元素为：NEB的路径、NEB阶段和状态、ini计算状态和fin计算状态。

---
```angular2html

```


