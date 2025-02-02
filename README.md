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
		+ **特征值分解**
		+ **奇异值分解**
	+ 矩阵变换运算
		+ 矩阵Gram-Schmidt变换
		+ 矩阵HouseHolder变换
		+ 矩阵Givens变换
	+ 矩阵特殊运算
		+ 矩阵Kronecker积
+ FirstMatrixC也基于ctypes提供了相应的python接口，以下是python接口使用示例：
```python
from firstmatrixc.utils import create_matrix,float_array_to_list,list_to_float_array
from firstmatrixc import lib



# 示例使用
if __name__ == "__main__":
    # 创建两个矩阵
    data1 = [1.0, 2.0, 3.0, 4.0]
    data2 = [5.0, 6.0, 7.0, 8.0]
    
    matrix1 = create_matrix(2, 2, data1)
    matrix2 = create_matrix(2, 2, data2)
    
    # 矩阵加法
    result_add = lib.MatrixAdd(matrix1, matrix2)
    result_add_data = float_array_to_list(result_add.data, result_add.row * result_add.column)
    print("Matrix Addition Result:", result_add_data)
    
    # 矩阵标量乘法
    scalar = 2.0
    result_scalar_mul = lib.MatrixScalarMultiply(scalar, matrix1)
    result_scalar_mul_data = float_array_to_list(result_scalar_mul.data, result_scalar_mul.row * result_scalar_mul.column)
    print("Matrix Scalar Multiplication Result:", result_scalar_mul_data)
    
    # 矩阵乘法
    result_mul = lib.MatrixMultiply(matrix1, matrix2)
    result_mul_data = float_array_to_list(result_mul.data, result_mul.row * result_mul.column)
    print("Matrix Multiplication Result:", result_mul_data)
    
    # 矩阵转置
    result_transpose = lib.MatrixTranspose(matrix1)
    result_transpose_data = float_array_to_list(result_transpose.data, result_transpose.row * result_transpose.column)
    print("Matrix Transpose Result:", result_transpose_data)
    
    # 矩阵行列式
    determinant = lib.MatrixDeterminant(matrix1)
    print("Matrix Determinant:", determinant)
    
    # LU 分解
    lu_result = lib.MatrixDecompositionLU(matrix1)
    print("LU Decomposition Result:")
    print("L Matrix:", float_array_to_list(lu_result.L.data, lu_result.L.row * lu_result.L.column))
    print("U Matrix:", float_array_to_list(lu_result.U.data, lu_result.U.row * lu_result.U.column))
    
    # Cholesky 分解
    cholesky_result = lib.MatrixDecompositionCholesky(matrix1)
    print("Cholesky Decomposition Result:")
    print("L Matrix:", float_array_to_list(cholesky_result.L.data, cholesky_result.L.row * cholesky_result.L.column))
    
    # QR 分解
    qr_result = lib.MatrixDecompositionQR(matrix1)
    print("QR Decomposition Result:")
    print("Q Matrix:", float_array_to_list(qr_result.Q.data, qr_result.Q.row * qr_result.Q.column))
    print("R Matrix:", float_array_to_list(qr_result.R.data, qr_result.R.row * qr_result.R.column))
    
    # SVD 分解
    svd_result = lib.MatrixDecompositionSVD(matrix1, 1e-6)
    print("SVD Decomposition Result:")
    print("U Matrix:", float_array_to_list(svd_result.U.data, svd_result.U.row * svd_result.U.column))
    print("V Matrix:", float_array_to_list(svd_result.V.data, svd_result.V.row * svd_result.V.column))
    
    # EVD 分解
    evd_result = lib.MatrixDecompositionEVD(matrix1, 1e-6)
    print("EVD Decomposition Result:")
    print("V Matrix:", float_array_to_list(evd_result.V.data, evd_result.V.row * evd_result.V.column))

```
