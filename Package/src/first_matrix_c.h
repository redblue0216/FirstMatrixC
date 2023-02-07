//模块名称: first_matrix_c
//设计：施华
//编码：施华
//简介：这是first_matrix_c矩阵库第一层接口头文件
//备注：无



/*
**功能介绍
**
**矩阵库接口整理
**（1）矩阵存储接口
**（2）矩阵基本运算接口
**（3）矩阵分解运算接口
**（4）矩阵交换运算接口
**（5）矩阵特别运算接口
**
*********************************************************
**
**技术介绍
**
**include
*/


//使用条件编译是为了防止include的循环引用问题
#ifndef FIRST_MATRIX_C_H__
#define FIRST_MATRIX_C_H__
/*
**基本标准库
*/
#include <stdio.h>
#include <stdlib.h>


//此处以下都注释掉是为了防止include的循环引用问题
/*
**基础矩阵数据类型
*/
#include "basic_data_type.h"


/*
**基础错误类型
*/
#include "basic_error_id.h"


/*
**基础矩阵数据结构
*/
#include "basic_data_structure.h"


/*
**矩阵存储模块
*/
#include "matrix_store.h"


/*
**矩阵基础运算模块
*/
#include "matrix_basic_operation.h"


/*
**矩阵分解运算模块
*/
#include "matrix_decomposition_operation.h"


/*
**矩阵变换运算模块
*/
#include "matrix_transformation_operation.h"


/*
**矩阵特征运算模块
*/
#include "matrix_special_operation.h"


//结束条件编译
#endif
/*******************************************************/
/*******************************************************/


