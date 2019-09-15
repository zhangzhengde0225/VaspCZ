# 一、简介：

VaspCZ(Vasp Check by Zzd)是作者在读博期间为提高科研效率而开发的Vasp辅助程序。该程序包换软件部分和API部分。

软件部分提供了Linux字符串用户界面，用于在超算平台中快捷提交任务和检查结果。包含三个模块：结构优化和静态计算(OS)模块、过渡态计算(NEB)模块和测试(Test)模块。

API部分为有python基础的研究者开发，以python库的形式便捷调用相关功能，以实现自定义高通量计算。[API说明文档](https://github.com/zhangzhengde0225/VaspCZ/blob/master/docs/VaspCZ_lib.md)

--------
# 二、程序框架

xxxx
---------
# 三、安装
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
# 四、使用和示例
## 1. VaspCZ 主程序
提供了Linux字符串用户界面，用于在超算平台中快捷提交任务和检查结果。包含三个模块：结构优化和静态计算(OS)模块、过渡态计算(NEB)模块和截断能K点测试(Test)模块。

成功安装后输入快捷键即可进入用户界面：
```angular2html
vcz
```

<img src="https://github.com/zhangzhengde0225/VaspCZ/blob/master/figs/VaspCZ_mainface.png" width="500" align=center>

输入模块对应的选择即可进入相应模块。

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


#### 示例：
进入到项目自带的examples：(如安装中改变了，请将"/home/zhangzd/bin"替换为你配置的主程序安装路径)
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

会在该目录下产生Vasp的5个输入文件的例子：INCAR、POSCAR、POTCAR、KPOINTS和Vasp.sh

其中，Vasp.sh为PBS系统提交任务的脚本，因不同平台的脚本内容会有所不同，请将适合该平台的脚本正确拷贝到安装目录下，默认为：用户根目录，目录结构如下所示：
```angular2html
用户根目录(或配置的Vasp.sh路径)
|
|   Vasp.sh
|   ...
```

#### 1.2 修改INCAR为静态计算的INCAR

在当前路径的结构优化INCAR上修改为静态计算的INCAR。修改项目：

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
    |   +---...
    | 
    +---PW91
    +---LDA
    +---US_LDA_GGA
    +---...
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

检查所有路径计算是否完成，输出当前路径、离子步数和电子步数。

检查完后，输入当前路径、能量、离子步数、磁矩、POSCAR和CONTCAR原子之间的距离、原子最大受力。

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
2.7|检查过渡态原子移动情况
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
在准备进行过渡态计算的目录下，创建文件夹ini和fin分别代表初态和末态，在它们之下再分别创建Opt文件夹。计算步骤如下：

```markflow
st=>start:Start
e=>End
op=>operation: My Operation
cond=>condition: Yes or No?
inp=>Input

st->op->cond
cond(yes)->e
cond(no)->inp->op
```
```flow
st=>start: index
op=>operation: 申请
op2=>operation: 结果页
op3=>operation: 查询本地
i1=>inputoutput: bid入库
i2=>inputoutput: 填写个人信息
c1=>condition: 检查登录
c2=>condition: 登录
c3=>condition: 查询本地记录
c4=>condition: 检测状态
c5=>operation: 风控审核
e=>end

st->op->c1()
c1(no)->c2(yes)->op()
c1(yes)->c3(no)->i1(right)->i2(right)->c5()->op2->e
c1(yes)->c3(yes)->c4(no)->i2
c1(yes)->c3(yes)->c4(yes)->op3->op2
c3()->e
```

```flow
sequenceDiagram
    Alice ->> Bob: Hello Bob, how are you?
    Bob-->>John: How about you John?
    Bob--x Alice: I am good thanks!
    Bob-x John: I am good thanks!
    Note right of John: Bob thinks a long<br/>long time, so long<br/>that the text does<br/>not fit on a row.

    Bob-->Alice: Checking with John...
    Alice->John: Yes... John, how are you?
```
1. ini/Opt/下进行初态的结构优化

2. fin/Opt/下进行末态的结构优化

3. ini/下在结构优化完成后进行静态计算以获得更准确的能量

4. fin/下在结构优化完成后末态静态计算。

5. 当前路径下在两个静态计算完成后进行过渡态计算

6. 如需，过渡态完成后当前路径下进行振动分析。


#### 功能示例：
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

输入：
```angular2html
vcz
```
选择功能2.1

再选择1为当前文件夹下的静态计算到结构优化。

选择2为一键提交ini/和fin/文件下下的静态计算。

可输入节点数、核数和文件名提交任务。默认为：

参数|默认值|
:---:|:---:
节点数|ini/Opt/Vasp.sh中读取
核数|ini/Opt/Vasp.sh中读取
任务名|ini/Opt/Vasp.sh中的最后一位改为S

#### 2.2 一键静态计算到过渡态计算
如前过渡态的一般过程所示，静态计算完成后，自动进行过渡态计算。

进入2.2:
```angular2html
cd 2.2
```
该文件夹下包含一般性的过渡态计算结构，且ini/Opt、fin/Opt、ini/和fin/下计算已完成。(可用OS模块的1.8功能检查结果)

输入vcz调用程序选择功能2.2即可实现自动提交过渡态计算任务。

可输入节点数、核数和文件名提交任务，默认为：

参数|默认值
:---:|:---:
节点数|~初末态结构原子距离和/0.8，取奇数
核数|ini/Opt/Vasp.sh中读取
任务名|ini/Opt/Vasp.sh中的最后一位改为N




####
## 2. VaspCZ python API


当前版本: 1.0.1

