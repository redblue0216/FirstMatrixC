//模块名称: matrix_decomposition_operation
//设计：施华
//编码：施华
//简介：这是first_matrix_c矩阵库的矩阵分解运算模块，提供矩阵分解运算操作
//备注：无



/*
**功能介绍
**
**提供矩阵分解运算操作，包括
**（1）LU分解
**（2）Cholesky分解
**（3）QR分解
**（4）特征值分解
**（5）奇异值分解
**
*********************************************************
**
**技术介绍
**
**声明函数
*/


//使用条件编译是为了防止include的循环引用问题
#ifndef __MATRIX_DECOMPOSITION_OPERATION_H__
#define __MATRIX_DECOMPOSITION_OPERATION_H__
/*
**基本标准库
*/
#include <stdio.h>
#include <stdlib.h>


/*
**矩阵存储库
*/
#include "matrix_store.h"


/*
**矩阵分解运算相关操作函数声明
*/
// LU矩阵分解
LUMatrix MatrixDecompositionLU(Matrix a);
// Cholesky分解
CholeskyMatrix MatrixDecompositionCholesky(Matrix a);
// QR分解
QRMatrix MatrixDecompositionQR(Matrix a);



//结束条件编译
#endif
/***********************************************************************************/
/***********************************************************************************/


