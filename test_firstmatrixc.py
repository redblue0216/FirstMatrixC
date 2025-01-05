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
