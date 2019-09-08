
# VaspCZ python API 文档

## 1.安装和使用
### (1) 安装
请查看项目的[README.md](https://github.com/zhangzhengde0225/VaspCZ)
### (2) 使用
```
import VaspCZ.zzdlib as zzd
```
### (3) 简单示例
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

## 2. API文档

### 1. shell模块
```
VaspCZ.zzdlib.getshellResult(code)
```
返回shell命令控制台输出的结果，由每一行组成一个元素的列表。

+ **code** : str
    
    shell命令，如：ls, qstat等  
    
+ **return** : list
    
    返回列表，其中每个元素是控制台输出的一行，保留换行符\n
    
    



> code


 
**1\. 查询指定项目属性**
###### 接口功能
> 获取制定项目的分类信息

###### URL
> [http://www.api.com/index.php](www.api.com/index.php)

###### 支持格式
> JSON

###### HTTP请求方式
> GET

###### 请求参数
|参数|必选|类型|说明|
|:-----  |:-------|:-----|-----                               |
|name    |ture    |string|请求的项目名                          |
|type    |true    |int   |请求项目的类型。1：类型一；2：类型二 。|

###### 返回字段
> |返回字段|字段类型|说明                              |
|:-----   |:------|:-----------------------------   |
|status   |int    |返回结果状态。0：正常；1：错误。   |
|company  |string | 所属公司名                      |
|category |string |所属类型                         |

###### 接口示例
> 地址：[http://www.api.com/index.php?name="可口可乐"&type=1](http://www.api.com/index.php?name="可口可乐"&type=1)
``` javascript
{
    "statue": 0,
    "company": "可口可乐",
    "category": "饮料",
}