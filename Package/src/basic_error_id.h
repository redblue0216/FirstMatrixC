//模块名称: basic_error_id
//设计：施华
//编码：施华
//简介：这是first_matrix_c矩阵库的基础错误类型预定义模块
//备注：无



/*
**功能介绍
**
**预定义基础错误类型
**
*********************************************************
**
**技术介绍
**
**define
*/



//使用条件编译是为了防止include的循环引用问题
#ifndef __BASIC_ERROR_ID_H__
#define __BASIC_ERROR_ID_H__
/*
**基本标准库
*/
#include <stdio.h>
#include <stdlib.h>



/*
**宏定义错误编码
*/
#define Matrix_Addition_Error_001 "@ERROR(ID:001): Matrix Dimensions Wrong!\n\tDetails: The number of rows and columns of two matrices is not equal!\n"
#define Matrix_Muliply_Error_002 "@ERROR(ID:002): Matrix Dimensions Wrong!\n\tDetails: The dimension of the column about first matrix is not equal to the dimension of the row about second matrix!\n"
#define Matrix_Square_Error_003 "@ERROR(ID:003): Matrix Dimensions Wrong!\n\tDetails: The matrix is not a square array!\n"



//结束条件编译
#endif
/***************************************************************************************************************************************************/
/***************************************************************************************************************************************************/


