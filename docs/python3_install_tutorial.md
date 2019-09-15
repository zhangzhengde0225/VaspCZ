# Python3.6.1 源码安装教程
用户独立python

## 1. 下载和安装
代码下载：
```angular2html
wget http://python.org/ftp/python/3.7.1/Python-3.7.1.tar.xz
```
或者官网下载最新版：[地址](https://www.python.org/downloads/source/)

解压并进入
```angular2html
tar xf Python-3.6.1.tar.xz
cd Python-3.6.1
```
配置安装路径：(--prefix后输入自己的安装目录)
```angular2html
./configure --prefix=/home/zhangzd/bin/python3
```
编译并安装
```angular2html
make
make install
```

## 2. 添加环境变量
安装后的python不在环境变量中不能直接调起，需要添加环境变量。

添加python3和pip3的软链接：(注：/home/zhangzd/bin/python3 应换为你的安装路径)
```angular2html
ln -s /home/zhangzd/bin/python3/bin/python3 /home/zhangzd/bin/python3/python3
ln -s /home/zhangzd/bin/python3/bin/pip3 /home/zhangzd/bin/python3/pip3
```
需要把软链路径写入到~/.bashrc文件中。
打开~/.bashrc文件：
```angular2html
vi ~/.bashrc
```
写入的内容为(写入一个路径同时包含python3和pip3)：(注：/home/zhangzd/bin/python3 应换为你的软链接的路径不带最后的文件名)
```angular2html
export PATH=/home/zhangzd/bin/python3/:$PATH
```
而后更新环境变量：
```angular2html
source ~/.bashrc
```
检查用户独立python3是否安装成功
```angular2html
which python3
which pip3
```
如果输出
```angular2html
~/bin/python3/python3
~/bin/python3/pip3
```
则表示成功安装。~代表你的用户根目录，你对该python3拥有完全权限。

用户独立的python3安装成功后，可进入VaspCZ安装目录继续安装。