//模块名称: matrix_store
//设计：施华
//编码：施华
//简介：这是first_matrix_c矩阵库的矩阵存储模块，提供存储相关操作函数实现
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
**功能函数
*/



/*
**矩阵存储头文件
*/
#include "matrix_store.h"



/*
**矩阵存储相关操作函数
*/


// 创建矩阵结构
Matrix CreateMatrix(INDEX row,INDEX column)
/*
**函数功能：
**定义一个创建矩阵结构的函数
**参数：
**row (INDEX): 行索引
**column (INDEX): 列索引
**返回：
**matrix (Matrix): 已分配好内存的矩阵结构
*/
{   
    Matrix matrix;
    matrix.row = row;
    matrix.column = column;
    matrix.data = (REAL *)malloc(row * column * sizeof(REAL));  // 动态分配内存

    return matrix;
}


// 设置矩阵数据
Matrix SetMatrixData(Matrix matrix,const REAL *array)
/*
**函数功能：
**定义一个设置矩阵数据的函数
**参数：
**matrix (Matrix): 已分配好内存的矩阵结构
**array (REAL): 一维实数矩阵
**返回：
**无
*/
{
    INDEX i;
    for (i=0;i<matrix.row*matrix.column;i++)
    {
        matrix.data[i] = array[i];
    }

    // return matrix; 
}


// 转换矩阵索引
INDEX TransformerIndex(INDEX matrix_column,INDEX i,INDEX j)
/*
**函数功能：
**定义一个转换矩阵索引的函数
**参数：
**matrix_column (INDEX): 矩阵列数
**i (INDEX): 行索引
**j (INDEX): 列索引
**返回：
**index (INDEX): 一维矩阵中对应的索引
*/
{
    INDEX index;
    index = ((i - 1) * matrix_column + j) - 1;

    return index;
}


// 获取矩阵元素
REAL GetMatrixElement(Matrix matrix,INDEX i,INDEX j)
/*
**函数功能：
**定义一个获取矩阵元素的函数
**参数：
**matrix (Matrix): 已分配好内存的矩阵结构
**i (INDEX): 行索引
**j (INDEX): 列索引
**返回：
**element (REAL): 矩阵元素
*/
{
    REAL element;
    INDEX ij;
    ij = TransformerIndex(matrix.column,i,j);
    element = matrix.data[ij];

    return element;
}


// 设置矩阵归零
void SetMatrixZero(Matrix matrix)
/*
**函数功能：
**定义一个设置矩阵归零的函数
**参数：
**matrix (Matrix): 已分配好内存的矩阵结构
**返回：
**无
*/
{
    INDEX i;
    for(i = 0;i < matrix.row * matrix.column;i++)
    {
        matrix.data[i] = 0.0;
    }
}



/*******************************************************/
/*******************************************************/


