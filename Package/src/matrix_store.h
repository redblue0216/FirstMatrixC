//模块名称: matrix_store
//设计：施华
//编码：施华
//简介：这是first_matrix_c矩阵库的矩阵存储模块，提供存储相关操作函数声明
//备注：无



/*
**功能介绍
**
**提供矩阵存储相关操作，包括
**（1）创建矩阵结构
**（2）设置矩阵数据
**（3）转换矩阵索引(二维变一维)
**（4）获取矩阵元素
**（5）设置矩阵归零
**
*********************************************************
**
**技术介绍
**
**声明函数
*/



//使用条件编译是为了防止include的循环引用问题
#ifndef __MATRIX_STORE_H__
#define __MATRIX_STORE_H__
/*
**基本标准库
*/
#include <stdio.h>
#include <stdlib.h>


/*
**基础数据类型
*/
#include "basic_data_type.h"


/*
**基础错误类型
*/
#include "basic_error_id.h"


/*
**基础数据结构
*/
#include "basic_data_structure.h"



/*
**矩阵存储相关操作函数声明
*/
// 创建矩阵结构
Matrix CreateMatrix(INDEX row,INDEX column);
// 设置矩阵数据
Matrix SetMatrixData(Matrix matrix,const REAL *array);
// 转换矩阵索引
INDEX TransformerIndex(INDEX i,INDEX j,INDEX matrix_column);
// 获取矩阵元素
REAL GetMatrixElement(Matrix matrix,INDEX i,INDEX j);
// 设置矩阵元素的具体值
VOID SetMatrixElement(Matrix matrix,INDEX i,INDEX j,REAL element);
// 设置矩阵归零
VOID SetMatrixZero(Matrix matrix);
// 求取最大值宏
#define MAX(a, b) ((a) > (b) ? (a) : (b))


//结束条件编译
#endif
/*******************************************************/
/*******************************************************/


