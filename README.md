### 求github star中，如果你觉得本项目不错，烦请点击项目右上角Star，感谢！~

### 反馈调试中，如使用中遇到问题，敬请上报到drivener@163.com，再次感谢！~

# 一、简介

VaspCZ(Vasp Check by Zzd)是作者在读博期间为提高科研效率而开发的Vasp辅助程序。该程序包含软件部分和API部分。

软件部分提供了Linux字符串用户界面，用于在超算平台中快捷提交任务和检查结果。包含三个模块：结构优化和静态计算(OS)模块、过渡态计算(NEB)模块和测试(Test)模块。

API部分为软件部分的底层，是自己写的一个python库。为有python基础的研究者提供了调用相关功能的接口，可以实现自定义计算和编写上层应用。库名：VaspCZ.zzdlib，包含三个模块：shell模块、File模块和Vasp模块。[API说明文档](https://github.com/zhangzhengde0225/VaspCZ/blob/master/docs/VaspCZ_lib.md)。

## 软件部分

超算平台linux(centos, ubuntu, redhat等)用户界面，无需python基础，会基本的linux命令即可。

+ OS模块(结构优化和静态计算模块)
    - 1.1 产生Vasp输入文件(示例)
    - 1.2 修改INCAR为静态计算INCAR
    - 1.3 产生贋势文件POTCAR
    - 1.4 产生网格文件KPOINTS
    - 1.5 产生提交任务脚本Vasp.sh
    - 1.6 仅保留Vasp输入文件
    - 1.7 一键前检查并提交任务
        * 自动检查有无输入错误，批量提交任务
    - 1.8 一键后检查并打印计算结果
        * 自动检查有无错误，自动检查所有路径下OS计算结果
    
+ NEB模块(过渡态计算模块)
    - 2.1 一键结构优化到静态计算
        * 自动检查有无输入错误，批量提交任务
    - 2.2 一键静态计算到过渡态计算
        * 自动检查有无输入错误，批量提交任务
    - 2.3 过渡态振动分析
    - 2.4 仅保留结构优化输入文件
    - 2.5 仅保留过渡态输入数据
    - 2.6 检查过渡态受力情况
    - 2.7 检查过渡态各态原子距离
    - 2.8 一键检查过渡态计算结果
        * 自动检查有无错误，自动检查所有路径下NEB计算结果
    - 2.9 一键检查过渡态振动分析结果
        * 自动检查有无错误，自动检查所有路径下振动分析结果
     
+ Test模块(截断能测试和K点测试模块)
    - 3.1 截断能测试
        * 一键截断能测试
    - 3.2 K点测试
        * 一键K点测试
        
## API部分

软件部分底层，python的API库，库名：VaspCZ.zzdlib，需有python基础，知道如何调用库。

+ shell模块
    - 1.1 快捷获取控制台输出
+ File模块
    - 2.1 一行代码读取和写入文件数据
    - 2.2 快捷修改文件内容
    - 2.3 快捷获取文件内容
    - 2.4 快捷获取文件内容(多行)
    - 2.5 获取文件空行(特定需求)
    - 2.6 获取路径
+ Vasp模块
    - 3.1 快捷解码POSCAR文件
    - 3.2 快捷修改POSCAR文件
    - 3.3 快捷生产POTCAR文件
    - 3.4 快速固定原子
    - 3.5 快捷修改振动分析INCAR文件
    - 3.6 检查输入文件，提交任务
    - 3.7 一键检查输入文件，提交任务(推荐)
    - 3.8 仅保留输入文件
    - 3.9 一键获取NEB计算周期

---------
# 二、安装
## 1. 安装和卸载
代码下载：
```angular2html
git clone https://github.com/zhangzhengde0225/VaspCZ.git
```
或者访问[github VaspCZ 网址](https://github.com/zhangzhengde0225/VaspCZ)下载，下载后解压：
```
unzip VaspCZ.zip
```
进入安装程序目录：
```angular2html
cd VaspCZ
```

默认VaspCZ安装配置为：

程序|安装路径|说明
:-----:|:-----:|:-----:
VaspCZ软件 | 用户根目录/bin | linux主程序
vtst | 用户根目录/bin | [VTST过渡态工具](http://theory.cm.utexas.edu/vtsttools/)
VaspCZ.zzdlib | 相应python3的site-packages目录| python API 接口
Vasp.sh | 用户根目录 | 超算平台PBS系统提交任务的脚本，需要自行准备拷贝到该路径下并命名为Vasp.sh
PseudoPotential | 用户根目录 | 生成POTCAR所需的贋势文件的文件夹，需要自行准备并拷贝到该路径下，命名方式为~/PseudoPotential/[贋势名]
vcz | 用户根目录/bin/VaspCZ | VaspCZ主程序快捷键

如需修改安装配置，请修改install.py第6-10行对应设置再进行安装。

安装：
```angular2html
python3 install.py
```
#### 注意：请将python更新至3.6及以上版本，原因：源代码中使用了f-strings功能，需要3.6及以上版本支持。
输入快捷键运行程序：
```angular2html
vcz
```
程序界面如图：

<img src="https://github.com/zhangzhengde0225/VaspCZ/blob/master/figs/VaspCZ_mainface.png" width="500" align=center>

如出现程序界面，说明安装成功。

如更新版本，安装前请先卸载：
```angular2html
python3 uninstall.py
```

## 2. 错误提示
(1) 权限不足

如果安装时提示：
```angular2html
PermissionError: [Errno 13] Permission denied: 'VaspCZ'
``` 
请使用管理员账号用以下命令安装。
```
sudo python3 install.py
```
如无管理员账号，请给当前用户安装独立的python后再安装VaspCZ

[源码安装用户独立python3教程](https://github.com/zhangzhengde0225/VaspCZ/blob/master/docs/python3_install_tutorial.md)

(2) 缺少python库

VaspCZ运行需要的库有：
```angular2html
numpy
```
如果提示：
```angular2html
ModuleNotFoundError: No module named 'numpy'
```
使用pip3安装相应库即可：
```angular2html
pip3 install numpy
```
-------
# 三、使用和示例
本章主要描述软件部分的使用方法和示例，Python API接口部分仅描述功能，接口详细信息参见[API说明文档](https://github.com/zhangzhengde0225/VaspCZ/blob/master/docs/VaspCZ_lib.md)。
## 1. VaspCZ 软件部分 (主程序)
软件部分提供了Linux字符串用户界面，用于在超算平台中快捷提交任务和检查结果。包含三个模块：结构优化和静态计算(OS)模块、过渡态计算(NEB)模块和截断能K点测试(Test)模块。

成功安装后输入快捷键即可进入用户界面：
```angular2html
vcz
```

<img src="https://github.com/zhangzhengde0225/VaspCZ/blob/master/figs/VaspCZ_mainface.png" width="500" align=center>

输入模块对应的选项即可进入相应模块。

### (1) Opt and Sta 模块

该模块提供了用于快捷进行结构优化(Optimization)计算和静态(Static)计算的功能。

用户界面：

<img src="https://github.com/zhangzhengde0225/VaspCZ/blob/master/figs/VaspCZ_OS_module.png" width="500" align="center">

包含功能：

功能标签|功能名称|
:-----:|:-----:
1.1|产生Vasp输入文件(示例)
1.2|修改INCAR为静态计算INCAR
1.3|产生贋势文件POTCAR
1.4|产生网格文件KPOINTS
1.5|产生提交任务脚本Vasp.sh
1.6|仅保留Vasp输入文件
1.7|前检查并提交任务
1.8|后检查并打印计算结果


#### OS模块功能示例：
进入到项目自带的examples：(请将"/home/zhangzd/bin"替换你的VaspCZ安装路径)
```angular2html
cd /home/zhangzd/bin/VaspCZ/examples/
```
------
OS模块下1.1-1.7功能示例：
进入1.1-1.7
```angular2html
cd 1.1-1.7
```
该文件夹为空文件夹。
输入：
```angular2html
vcz
1
```
而后输入1-7数字可以执行相应功能

#### 1.1 产生Vasp输入文件(示例)

会在该目录下产生Vasp的5个输入文件的示例：INCAR、POSCAR、POTCAR、KPOINTS和Vasp.sh

注意：生成Vasp.sh文件需要配置：Vasp.sh为PBS系统提交任务的脚本，因不同平台的脚本内容会有所不同，请将适合该平台的脚本正确拷贝到安装目录下，默认为：用户根目录，目录结构如下所示：
```angular2html
用户根目录(或配置的Vasp.sh路径)
|
|   Vasp.sh
|   ...(files)
```

#### 1.2 修改INCAR为静态计算的INCAR

在当前路径的结构优化INCAR上修改为静态计算的INCAR。

修改项目：

```angular2html
SYSTEM=Static
IBRION=-1
NSW=1
# EDIFFG=-0.01
```

#### 1.3 产生POTCAR

输入元素列表和贋势类型产生POTCAR。

默认产生适配当前目录下的POSCAR内的元素的POTCAR，默认贋势类型为PBE。

注意：将从安装VaspCZ时配置的贋势路径下读取数据，默认为用户根目录。使用该功能请将贋势文件夹命名为PseudoPotential并按如下目录安装。

```
用户根目录(或配置贋势安装路径)
|   
+---PseudoPotential
    |
    +---PBE
    |   |
    |   +---H
    |   +---He
    |   +---...(dirs)
    | 
    +---PW91
    +---LDA
    +---US_LDA_GGA
    +---...(dirs)
```

#### 1.4 产生KPOINTS

输入网格和方法产生KPOINTS文件。

默认网格为：5 5 5。

默认方法为：Monkhorst 。与Vasp官网一致，方法可只输入开头的字母如：M，可选方法有：M(Monkhorst)，A(Auto)

#### 1.5 产生Vasp.sh

输入任务所需节点数、核数和任务名产生提交任务脚本Vasp.sh。

默认：节点数：1 核数：12 任务名：jobname

注意：将从VaspCZ安装时候配置的Vasp.sh路径下读取数据，默认为用户根目录。使用该功能前请正确安装Vasp.sh。

#### 1.6 保留Vasp输入文件

删除其他所有文件和文件夹，仅保留Vasp的5个输入文件(INCAR、POSCAR、POTCAR、KPOINTS和Vasp.sh)，用于计算出现问题，重新算。

选择该功能后可输入文件名添加需要额外保留的文件。

#### 1.7 前检查并提交任务

准备好输入文件后，进行前检查，检查INCAR、POSCAR和POTCAR是否匹配，检查通过后将打印检查信息，并提示是否提交任务。

-----
OS模块下1.8功能示例：
退出1.1-1.7并进入1.8
```angular2html
cd ..
cd 1.8
```
该文件夹为计算好的Fe-Te体系不同情形下的结构优化结果。

输入运行1.8功能：
```angular2html
vcz
1
8
```

#### 1.8 检查结果

检查当前目录及所有子目录下的结构优化和静态计算的结果，如OUTCAR或者log中有错误(ERROR)或警告(WARNING)或提示所在位置。

输出如图所示：

<img src="https://github.com/zhangzhengde0225/VaspCZ/blob/master/figs/VaspCZ_1.8_output.png" width="700" align=center>

检查所有路径计算是否完成，输出当前路径、完成状态、离子步数和电子步数。

检查完后，输出当前路径、能量、离子步数、磁矩、POSCAR和CONTCAR原子之间的距离、原子最大受力。

------
### (2) NEB 模块
该模块提供了便捷的NEB方法计算过渡态的功能。

用户界面：

<img src="https://github.com/zhangzhengde0225/VaspCZ/blob/master/figs/VaspCZ_NEB_module.png" width="500" align=center>

包含功能：

功能标签|功能名称|
:-----:|:-----:
2.1|一键结构优化到静态计算
2.2|一键静态计算到过渡态计算
2.3|过渡态振动分析
2.4|仅保留结构优化输入文件
2.5|仅保留过渡态输入数据
2.6|检查过渡态受力情况
2.7|检查过渡态各态原子距离
2.8|检查过渡态计算结果
2.9|检查过渡态振动分析结果

过渡态计算的一般过程：先做结构优化，而后静态计算，最后过渡态计算，如需再振动分析。目录结果如下：
```angular2html
NEB计算目录
| ...(files)
+---ini
|   | ...(files)
|   +---Opt
|       | ...(files)
+---fin
    | ...(files)
    +---Opt
        | ...(files)
```
在准备进行过渡态计算的目录下，创建文件夹ini和fin分别代表初态和末态，在它们之下再分别创建Opt文件夹。

计算步骤如下：

1. ini/Opt/下进行初态的结构优化。

2. fin/Opt/下进行末态的结构优化。

3. ini/下在结构优化完成后进行静态计算以获得更准确的能量。

4. fin/下在结构优化完成后末态静态计算。

5. 当前路径下在两个静态计算完成后进行过渡态计算。

6. 如需，过渡态完成后当前路径下进行振动分析。


#### NEB模块功能示例：
进入到VaspCZ安装目录examples文件夹下：
```angular2html
cd /home/zhangzd/bin/VaspCZ/examples 
```

#### 2.1 一键结构优化到静态计算
如前过渡态的一般过程所示，结构优化完成后，自动进行初末态的静态计算。

进入2.1：
```angular2html
cd 2.1
```
该文件夹下包含一般性的过渡态计算结构，且ini/Opt和fin/Opt下计算已完成。(可用OS模块的1.8功能检查结果)

调用vcz，选择功能2.1：
```angular2html
vcz
2
1
```
此时：选择1为当前文件夹下的静态计算到结构优化，选择2为一键提交ini/和fin/文件下下的静态计算。

输入节点数、核数和文件名提交任务。默认为：

参数|默认值|
:---:|:---:
节点数|ini/Opt/Vasp.sh中读取
核数|ini/Opt/Vasp.sh中读取
任务名|ini/Opt/Vasp.sh中的最后一位改为S

#### 2.2 一键静态计算到过渡态计算
如前过渡态的一般过程所示，静态计算完成后，自动进行过渡态计算。

进入2.2文件夹:
```angular2html
cd 2.2
```
该文件夹下包含一般性的过渡态计算结构，且ini/Opt、fin/Opt、ini/和fin/下计算已完成。(可用OS模块的1.8功能检查结果)

输入vcz调用程序选择功能2.2即可实现自动提交过渡态计算任务。

输入节点数、核数和文件名提交任务。

默认参数为：

参数|默认值
:---:|:---:
节点数|~初末态结构原子距离和/0.8，取奇数
核数|ini/Opt/Vasp.sh中读取
任务名|ini/Opt/Vasp.sh中的最后一位改为N

#### 2.3 过渡态振动分析

过渡态完成后，计算迁移原子在初态、过渡态和末态中三个自由度上的尝试频率。

使用初态和过渡态的尝试频率可以计算该迁移过程的有效率。

计算方法为：该原子在初态时三个自由度上的尝试频率之积 比 该原子过渡态时的两个自由度上(共三个自由度，其中一个是虚频)的频率之积。

进入2.3文件夹：
```angular2html
cd 2.3
```

该文件夹下包含已经计算好的过渡态文件。(可用NEB模块的2.8功能检查结果)

调用vcz，并选择功能2.3

输入任务节点数、核数和是否包含末态振动提交任务。

默认参数为：

参数|默认值
:---:|:---:
节点数|1
核数|8
是否包含末态|False

提交任务后会创建vib_analysis文件夹，内再创建ini_state, sad_state和fin_state，计算不同结构中迁移原子的振动频率(尝试频率)。

#### 2.4 仅保留结构优化输入文件

删除当前目录下的所有文件和文件夹，仅保留ini/Opt/下和fin/Opt下的5个输入文件(INCAR, POSCAR, POTCAR, KPOINTS和Vasp.sh)。

该功能用于过渡态计算错误时回滚到结构优化重新计算。

#### 2.5 仅保留过渡态输入数据

删除当前目录下的文件和文件夹，仅保留ini/和fin/文件夹下所有内容。

该功能英语过渡态计算错误时回滚到过渡态重新计算。删除后调用NEB模块的2.2功能即可重新提交NEB任务。

#### 2.6 检查过渡态受力情况

NEB计算完成或正在计算中，检查每一离子步，每个IMAGE下的受力状况。

例如：进入examples/2.6-2.9文件夹，调用vcz 2.6 功能：
```angular2html
cd 2.6-2.9
vcz
2
6
```
输出如图所示：

<img src="https://github.com/zhangzhengde0225/VaspCZ/blob/master/figs/VaspCZ_2.6_output.png" width="500" align=center>

第一列为离子步，第二到四列为插入态IMAGE01、IMAGE02和IMAGE03在对应离子步下该结构中原子所受的最大力，第五列为前面的二到四列之和。如数据所示，第8步时所有插入态原子最大受力小于0.01 eV/Å，达到INCAR中的收敛要求。

该功能用于检查过渡态计算不收敛时较为合理的结构。例：假如INCAR中设置NSW=100，计算达100步未收敛，通常第100步并非合理的结构。借助此功能可找到最大受力和最小的步数，将该步的结构取出进行进一步分析和计算。


#### 2.7 检查过渡态各态原子距离

NEB计算完成或NEB计算生成插入态后，检查每个态之间原子的距离和。

例如：进入examples/2.6-2.9文件夹，调用vcz 2.7功能：
```angular2html
cd 2.6-2.9
vcz
2
7
```
选择需要检查的结构，默认为POS，代表POSCAR，可选为CONT，代表CONTCAR。

输出如图所示：

<img src="https://github.com/zhangzhengde0225/VaspCZ/blob/master/figs/VaspCZ_2.7_output.png" width="300" align=center>

第一列是POSCAR或CONTCAR，第二列是IMAGE，第三列是原子距离和。其值来自于vtst工具，如第一行的值为：
```angular2html
dist.pl 00/POSCAR 01/POSCAR
```
计算前检查POSCAR，用于确保插入过渡态准备，线性插入时各态距离和应相等。

计算后检查CONTCAR，用于查看过渡态中是否有某个态弛豫到不可预测的结构导致过渡态不收敛。

#### 2.8 检查过渡态计算结果

NEB计算完成后或计算中，检查当前目录及所有子目录下的NEB计算结果(忽略静态计算和结构优化)，如OUTCAR或者log中有错误(ERROR)或警告(WARNING)或提示所在位置，检查完成后输出结果。

例如：进入examples/2.6-2.9文件夹，调用vcz 2.8功能：
```angular2html
cd 2.6-2.9
vcz
2
8
```
输出如图所示：

<img src="https://github.com/zhangzhengde0225/VaspCZ/blob/master/figs/VaspCZ_2.8_output.png" width="700" align=center>

每一个有NEB计算的路径都会输出计算结果。第一列为不同的IMAGE，第二列为原子最大受力，第三列为该IMAGE总能，第四列为以IMAGE00作为参考原点是的能量差，最大能量差即为势垒，对应的IMAGE为鞍点。

如数据所示，该扩散过程(fcc Fe的自扩散)的扩散势垒为1.39 eV.

#### 2.9 检查过渡态振动分析结果

NEB振动分析结束后，检查当前目录及所有子目录下的原子振动频率(尝试频率)结果并计算有效频率。

例如：进入examples/2.6-2.9文件夹，调用vcz 2.9功能：
```angular2html
cd 2.6-2.9
vcz
2
9
```
输出如图所示：

<img src="https://github.com/zhangzhengde0225/VaspCZ/blob/master/figs/VaspCZ_2.9_output.png" width="700" align=center>

如数据所示，[True, True, False] 说明该子目录下包含初态和鞍点态振动分析，不包含末态。

第一个1f 2f 3f为迁移原子在初态结构中三个方向的振动，振动频率分别为6.59, 6.17和4.99 THz。第二个1f 2f 3f为迁移原子在鞍点态结构中三个方向的振动，振动频率分别为6.92, 4.67和5.81 THz，其中f/i表示第三个方向上为虚频。

该扩散过程原子的有效频率为：初态三个振动之积比鞍点态两个振动之积(排除虚频)，结果为: 6.28     THz。

本例是fcc Fe的自扩散，扩散前后结构等价，初态和末态相同，因此无需算末态振动。

通常体系中对称性不高，如有2个以上缺陷时，初态和末态是不等价的，此时反向方扩散的势垒就是以末态能量为原点时鞍点的能量，对应的有效频率为末态三个振动之积比鞍点态两个振动之积(排除虚频)。(在NEB模块2.3功能中输入参数包含末态时，2.9功能会自动计算反方向扩散的有效频率。)  

------
### (3) Test 模块
通常，一个体系在大规模进行计算和分析之前，需要进行截断能测试和K点测试确定合适的ENCUT设置和KPOINS设置。

该模块提供了快捷的Vasp截断能测试和K点测试功能。

用户界面：

<img src="https://github.com/zhangzhengde0225/VaspCZ/blob/master/figs/VaspCZ_Test_module.png" width="500" align=center>

包含功能：

功能标签|功能名称|
:-----:|:-----:
3.1|截断能测试
3.2|K点测试

#### Test模块功能示例：
进入到VaspCZ安装目录examples文件夹下：
```angular2html
cd /home/zhangzd/bin/VaspCZ/examples 
```

#### 3.1 截断能测试

做截断能测试的目的是选取一个合适的截断能，截断能决定了Vasp计算过程中被作为贋势处理的电子波函数的范围。截断能太小，计算得到的体系总能不可信，截断能太大，计算中迭代需要花费大量资源。

准备好输入文件(INCAR，POSCAR，POTCAR，KPOINTS和Vasp.sh)后，输入参数即可快捷提交截断能测试任务。

例如：进入examples/3.1文件夹，调用vcz 3.1功能：
```angular2html
cd 3.1
vcz
3
1
```
输入参数有：任务名前缀、节点数、核数和截断能列表。

默认参数为：

参数|默认值
:---:|:---:
任务名前缀|ENCUT_
节点数|1
核数|8
截断能列表|200,250,300,350,400,450,500,550,600,650,700

注意：截断能列表以英文逗号隔开。

提交任务后会以截断能为名创建文件夹，在每个文件夹内修改INCAR文件中的ENCUT为对应值，而后提交结构优化任务，任务名为任务名前缀+截断能。

计算完成后，可以使用OS模块的1.8功能检查各截断能时体系的总能，体系总能之差小于0.001 eV时，该截断能可选为合适的截断能。

#### 3.2 K点测试

做K点测试的目的是选取一个KPOINS设置，K点决定了Vasp计算过程中倒空间的网格分隔点数，体系越大，合适的K点网格一般越小。

准备好输入文件(INCAR，POSCAR，POTCAR，KPOINTS和Vasp.sh)后，输入参数即可快捷提交K点测试任务。

例如：进入examples/3.2文件夹，调用vcz 3.2功能：
```angular2html
cd 3.2
vcz
3
2
```
输入参数有：任务名前缀、节点数、核数和K点列表。

默认参数为：

参数|默认值
:---:|:---:
任务名前缀|ktest_
节点数|1
核数|8
K点列表|111,333,555,777,999

注意：K点列表以英文逗号隔开。

提交任务后会以K点为名创建文件夹，在每个文件夹内KPOINTS文件中的网格为K点，而后提交结构优化任务，任务名为任务名前缀+K点。

计算完成后，可以使用OS模块的1.8功能检查各K点时体系的总能。

-----
## 2. VaspCZ python API
python API部分为有python基础的研究者提供了本项目同通用功能的接口。通过库便捷调用相关功能，以实现自定义高通量计算。库名：VaspCZ.zzdlib，包含：shell模块，File模块和Vasp模块

### 安装和导入
安装软件时自动安装库，安装说明见本说明第三章。

导入：进入python3交互界面或在.py文件中导入库：
```angular2html
import VaspCZ.zzdlib as zzd
```
---
此处只列出各模块功能，详细接口说明见[API文档](https://github.com/zhangzhengde0225/VaspCZ/blob/master/docs/VaspCZ_lib.md)
### (1) shell模块
标签|代码|功能
:---:|:---:|:---:
1.1|VaspCZ.zzdlib.getshellResult(code)|返回shell命令控制台输出的结果，由每一行组成一个元素的列表。

### (2) File模块
标签|代码|功能
:---:|:---:|:---:
2.1|VaspCZ.zzdlib.File.openFile(path, [mode='r', data=None])|读取文件或保存文件
2.2|VaspCZ.zzdlib.File.substitudeData(data, keywords, newline, [mode='default'])|传入文件数据，给出关键词和新行，默认情形搜索出现第一次出现关键词的行并替换，mode不等于default是替换全部出现关键字的行，返回替换后的数据。
2.3|VaspCZ.zzdlib.File.getLine(data,keywords)|给出关键词，招傲有关键词的第一行并返回，返回为字符串和所在的行索引。该功能用于获取特定想信息或者用于判断。
2.4|VaspCZ.zzdlib.getAllline(data, keywords)|给出关键词，返回所有带有关键词的所有行，返回为列表。该功能用于选择性获得文件特定行。
2.5|VaspCZ.zzdlib.getNullline(data)|获取文件数据中是空位的索引。
2.6|VaspCZ.zzdlib.Vaspsh_path()|获取VaspCZ软件默认的PBS提交任务脚本Vasp.sh所在的文件路径。

### (3) Vasp模块
标签|代码|功能
:---|:---:|:---:
3.1|VaspCZ.zzdlib.Vasp.decode_POSCAR(POSCAR)|解码POSCAR，返回一个基矢、原子种类、原子数目、每个原子的位置（取前4位）
3.2|VaspCZ.zzdlib.Vasp.modify_POSCAR_ele(oldele, new_ele)|修改当前路径下POSCAR的原子种类，适合批量修改。
3.3|VaspCZ.zzdlib.Vasp.gennerate_POTCAR([elements=None, pseudotype='PBE'])|在当前路径生成POTCAR文件，需要在安装中正确是指贋势文件目录，默认贋势文件目录为用户根目录。贋势目录名为：PseudoPotential。
3.4|modify_POSCAR_Selective_Dynamics(data, indexes)|根据输入的数据和索引修改POSCAR，添加Selective Dynamics, 索引所在的位置设置为T T T, 其他位置设置为 F F F
3.5|modify_INCAR_for_vibration_analysis()|修改当前目录的INCAR为振动分析的INCAR并保存.
3.6|VaspCZ.zzdlib.Vasp.checkInputs()|Vasp前检查。提交计算任务前，检查当前目录Vasp的各项输入文件，将计算信息打印到控制台，包含：计算路径、SYSTEM、截断能、ISIF、离子更新方法、是否有磁性、电子收敛标准、离子收敛标准、原子种类个数、POTCAR原子类型、KPOINTS方法、网格大小、任务名、节点数与核数、是否加急。
3.7|VaspCZ.zzdlib.Vasp.check_and_qsub([need_input=True])|检查前检查并提交任务。内部集成了上一个检查输入文件函数，使用中推荐该函数。
3.8|VaspCZ.zzdlib.Vasp.keepInputs([addfile=[], workdir='./'])|删除工作目录下的文件，仅保留输入文件。默认保留文件为：INCAR，POSCAR，POTCAR， KPOINTS和Vasp.sh
3.9|VaspCZ.zzdlib.Vasp.checkNEBperiod()|遍历当前路径下的所有文件夹，如果发现有neb计算，判断ini和fin分别的计算周期，并返回

-----
# 四、其他说明
该项目已免费开源，[开源许可](https://github.com/zhangzhengde0225/VaspCZ/blob/master/LICENSE.md)。

欢迎开发和补充，如用于商业用途请注明出处。

如遇bug，敬请将说明、提示代码、截图等信息上报到drivener@163.com。

如对程序有疑问，请联系drivener@163.com。

作者水平有限，代码有诸多不足之处，还望斧正。


当前版本: 1.0.1

