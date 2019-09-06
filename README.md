简介：

VaspCZ(Vasp Check by Zzd)是作者在读博期间为提高科研效率而开发的Vasp辅助程序，该程序包含三个模块：结构优化和静态计算(OS)模块、过渡态计算(NEB)模块和测试(Test)模块，并提供了Linux用户界面。

【程序框架图】

【各个模块功能介绍】

当前版本: 1.0.1

安装和使用(终端)

    git clone https://github.com/zhangzhengde0225/VaspCZ.git

    cd VaspCZ
    
    python3 install.py
    
默认VaspCZ软件安装路径为：用户根目录/bin下，默认PBS系统提交任务的脚本Vasp.sh存在于用户根目录下，默认贋势文件夹Pseudopotential存在于用户根目录下，默认安装完成后程序快捷键为vcz，默认会安装VTST过渡态计算工具。如需修改默认设置，请修改install.py第6-10行对应设置再进行安装。

    
    
