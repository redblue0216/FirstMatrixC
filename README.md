# FirstMatrixC矩阵计算库

![shields_version](/static/shields_version.svg)  ![shields_license](/static/shields_license.svg)  ![shields_author](/static/shields_author.svg)  ![shiedls_gcc](/static/shields_gcc.svg) 
![FirstMatrixCsymbol](/static/FirstMatrixCsymbol.JPG)

## 介绍
+ FirstMatrixC是一个基于C语言实现的矩阵计算库，主要功能包括矩阵基本运算、矩阵分解运算、矩阵变换运算和矩阵特殊运算，主要技术包括二级架构的模块化编程、动态内存管理、条件编译、防御性编程和新建矩阵数据结构。
+ FirstMatrixC是一个练习项目，主要用于提供矩阵计算相关知识的介绍说明使用，相关文档可见Doc文件夹，相关视频后期会上传到个人网站和B站，敬请期待。

## 安装
+ FirstMatrixC提供了编译好的静态链接库，如果需要自己编译可以使用src文件夹下的makefile自行编译。

## 设计
+ FirstMatrixC实现语言：C语言
+ FirstMatrixC采用模块化设计思想，设计为两层架构。底层设计为矩阵存储相关操作，顶层设计为具体的矩阵运算相关操作，顶层依赖于底层矩阵存储模块。
+ FirstMatrixC预定义了数据类型，方便在不同机器平台上调试程序。
+ FirstMatrix还采用了防御性编程，提高程序稳定性。
+ 以下是FirstMatrixC支持的矩阵计算功能（加粗为已实现功能，括号内包括功能接口函数）
	+ **矩阵存储(matrix_store)**
		+ **创建矩阵结构(CreateMatrix)**
		+ **设置矩阵数据(SetMatrixData)**
		+ **转换矩阵索引(TransformerIndex)**
		+ **获取矩阵元素(GetMatrixElement)**
		+ **设置矩阵归零(SetMatrixZero)**
	+ 矩阵基本运算
		+ **矩阵加法(MatrixAdd)**
		+ **矩阵标量乘法(MatrixScalarMultiply)**
		+ **矩阵乘法(MatrixMultiply)**
		+ **矩阵转置(MatrixTranspose)**
		+ **矩阵求逆(MatrixInverse)**
		+ **矩阵行列式(MatrixDeterminant)**
		+ 矩阵范数
	+ 矩阵分解运算
		+ **矩阵LU分解**
		+ **矩阵Cholesky分解**
		+ **矩阵QR分解**
		+ 特征值分解
		+ **奇异值分解**
	+ 矩阵变换运算
		+ 矩阵Gram-Schmidt变换
		+ 矩阵HouseHolder变换
		+ 矩阵Givens变换
	+ 矩阵特殊运算
		+ 矩阵Kronecker积
