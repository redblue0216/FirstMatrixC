//模块名称: matrix_basic_operation
//设计：施华
//编码：施华
//简介：这是first_matrix_c矩阵库的矩阵基本运算模块，提供矩阵基本运算操作函数声明
//备注：无



/*
**功能介绍
**
**提供矩阵基本运算操作，包括
**（1）矩阵加法
**（2）矩阵标量乘法
**（3）矩阵乘法
**（4）矩阵转置
**（5）矩阵求逆
**（6）矩阵行列式
**（7）矩阵范数
**
*********************************************************
**
**技术介绍
**
**声明函数
*/



//使用条件编译是为了防止include的循环引用问题
#ifndef __MATRIX_BASIC_OPERATION_H__
#define __MATRIX_BASIC_OPERATION_H__
/*
**基本标准库
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


/*
**矩阵存储库
*/
#include "matrix_store.h"


/*
**矩阵基本运算相关操作函数声明
*/
// 矩阵加法
Matrix MatrixAdd(Matrix a,Matrix b);
// 矩阵标量乘法
Matrix MatrixScalarMultiply(REAL scalar,Matrix a);
// 矩阵乘法
Matrix MatrixMultiply(Matrix a,Matrix b);
// 矩阵转置
Matrix MatrixTranspose(Matrix a);
// 矩阵代数余子式计算
Cofactor AlgebraicCofactor(Matrix a);
// 代数余子式转变为矩阵数组
MatrixSeries CofactorToMatrixSeries(Cofactor b);
// 矩阵行列式计算
REAL MatrixDeterminant(Matrix a);
// 矩阵余子式计算
Matrix MatrixCofactor(Matrix a,INDEX i,INDEX j);
// 伴随矩阵计算
AdjointMatrix MatrixAdjoint(Matrix a);
// 矩阵求逆计算
Matrix MatrixInverse(Matrix a);



//结束条件编译
#endif
/*******************************************************/
/*******************************************************/


