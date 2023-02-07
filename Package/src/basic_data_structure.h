//模块名称: basic_data_structrue
//设计：施华
//编码：施华
//简介：这是first_matrix_c矩阵库的基础数据结构模块
//备注：无



/*
**功能介绍
**
**矩阵基础数据结构
**
*********************************************************
**
**技术介绍
**
**struct
*/



//使用条件编译是为了防止include的循环引用问题
#ifndef __BASIC_DATA_STRUCTURE_H__
#define __BASIC_DATA_STRUCTURE_H__
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
**矩阵数据结构Matrix
*/
typedef struct
{
    INTEGER row,column;
    REAL *data;  // 使用一维数组存储矩阵元素数据
}Matrix;


/*
**余子式数据结构Cofactor
*/
typedef struct
{
    INTEGER order;
    REAL *element;  // 使用一维数组存储余子式对应元素
    REAL *reverse_order_number;  // 使用一维数组存储逆序数
    REAL *data;  // 使用一维数组存储余子式元素数据
}Cofactor;


/*
**矩阵序列数据结构
*/
typedef struct 
{
    INTEGER matrix_number;
    INTEGER order;
    REAL *element;  // 使用一维数组存储余子式对应元素
    REAL *reverse_order_number;  // 使用一维数组存储逆序数
    Matrix *data; // 使用一维数组存储Matrix数据类型数据
}MatrixSeries;


/*
**矩阵序列的序列数据结构
*/
typedef struct
{
    INTEGER series_number;
    MatrixSeries *data; // 使用一维数组存储MatrixSeries数据类型数据
}SeriesSeries;


/*
**伴随矩阵数据结构
*/
typedef struct
{
    REAL *reverse_order_number;  // 使用一维数组存储逆序数
    REAL *data;  // 存储MatrixMatrix数据类型数据
    INTEGER order;
    REAL matrix_determinant;  // 使用实数存储矩阵行列式的值  
}AdjointMatrix;



//结束条件编译
#endif
/*******************************************************/
/*******************************************************/


