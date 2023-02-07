//模块名称: basic_data_type
//设计：施华
//编码：施华
//简介：这是first_matrix_c矩阵库的基础数据类型预定义模块
//备注：无



/*
**功能介绍
**
**预定义基础数据类型
**
*********************************************************
**
**技术介绍
**
**typedef
*/



/*
**基本标准库
*/
#include <stdio.h>
#include <stdlib.h>



//使用条件编译是为了防止include的循环引用问题
#ifndef __BASIC_DATA_TYPE_H__
#define __BASIC_DATA_TYPE_H__
/*
**基础数据类型
*/
//错误编码
typedef unsigned int ERROR_ID;
//索引
typedef int INDEX;
//字符串
typedef char* STRING;
//无返回值
typedef void VOID;
//整数
typedef int INTEGER;
//实数
typedef float REAL;



//结束条件编译
#endif
/*******************************************************/
/*******************************************************/


