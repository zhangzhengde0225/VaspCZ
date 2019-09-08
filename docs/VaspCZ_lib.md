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
---
```angular2html
VaspCZ.zzdlib.File.openFile(path, [mode='r', data=None])
```
读取文件或保存文件。

+ **path** : str
    
    文件所在路径。
    
+ **mode** : str, r or w
    
    模式，r是读取，w是写入。

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

    模式，默认为default, 只替换第一次出现关键词的行，修改为其他任意字符串，替换所有出现关键词的行
    
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

给出关键词，返回有关键词的第一行并返回，返回为字符串和所在的行号
    
    


