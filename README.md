# 一、简介：

VaspCZ(Vasp Check by Zzd)是作者在读博期间为提高科研效率而开发的Vasp辅助程序。该程序包换软件部分和API部分。

软件部分提供了Linux字符串用户界面，用于在超算平台中快捷提交任务和检查结果。包含三个模块：结构优化和静态计算(OS)模块、过渡态计算(NEB)模块和测试(Test)模块。

API部分为有python基础的研究者开发，以python库的形式便捷调用相关功能，以实现自定义高通量计算。[API说明文档](https://github.com/zhangzhengde0225/VaspCZ/blob/master/docs/VaspCZ_lib.md)

# 二、程序框架

xxxx

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





卸载：
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

[python官网](https://www.python.org)

[源码安装python教程](https://github.com/zhangzhengde0225/VaspCZ/blob/master/docs/python3_install_tutorial.md)

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

# 四、使用
## 1. VaspCZ linux 主程序
提供了Linux字符串用户界面，用于在超算平台中快捷提交任务和检查结果。包含三个模块：结构优化和静态计算(OS)模块、过渡态计算(NEB)模块和测试(Test)模块。

成功安装后输入快捷键即可进入用户界面：
```angular2html
vcz
```


### (1) Opt and Sta 模块
    

##【各个模块功能介绍】

## 2. VaspCZ python API


当前版本: 1.0.1

