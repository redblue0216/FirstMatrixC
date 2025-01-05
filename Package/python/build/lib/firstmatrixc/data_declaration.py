import ctypes
from firstmatrixc import lib


# 定义 C 结构体和函数原型

# Matrix 结构体
class Matrix(ctypes.Structure):
    _fields_ = [
        ("row", ctypes.c_int),
        ("column", ctypes.c_int),
        ("data", ctypes.POINTER(ctypes.c_float))
    ]

# Cofactor 结构体
class Cofactor(ctypes.Structure):
    _fields_ = [
        ("order", ctypes.c_int),
        ("element", ctypes.POINTER(ctypes.c_float)),
        ("reverse_order_number", ctypes.POINTER(ctypes.c_float)),
        ("data", ctypes.POINTER(ctypes.c_float))
    ]

# MatrixSeries 结构体
class MatrixSeries(ctypes.Structure):
    _fields_ = [
        ("matrix_number", ctypes.c_int),
        ("order", ctypes.c_int),
        ("element", ctypes.POINTER(ctypes.c_float)),
        ("reverse_order_number", ctypes.POINTER(ctypes.c_float)),
        ("data", ctypes.POINTER(Matrix))
    ]

# AdjointMatrix 结构体
class AdjointMatrix(ctypes.Structure):
    _fields_ = [
        ("reverse_order_number", ctypes.POINTER(ctypes.c_float)),
        ("data", ctypes.POINTER(ctypes.c_float)),
        ("order", ctypes.c_int),
        ("matrix_determinant", ctypes.c_float)
    ]

# LUMatrix 结构体
class LUMatrix(ctypes.Structure):
    _fields_ = [
        ("L", Matrix),
        ("U", Matrix)
    ]

# CholeskyMatrix 结构体
class CholeskyMatrix(ctypes.Structure):
    _fields_ = [
        ("L", Matrix)
    ]

# QRMatrix 结构体
class QRMatrix(ctypes.Structure):
    _fields_ = [
        ("Q", Matrix),
        ("R", Matrix)
    ]

# SVDMatrix 结构体
class SVDMatrix(ctypes.Structure):
    _fields_ = [
        ("A", Matrix),
        ("U", Matrix),
        ("V", Matrix)
    ]

# EVDMatrix 结构体
class EVDMatrix(ctypes.Structure):
    _fields_ = [
        ("A", Matrix),
        ("V", Matrix)
    ]

# 定义函数原型

# 矩阵基本运算
lib.MatrixAdd.argtypes = [Matrix, Matrix]
lib.MatrixAdd.restype = Matrix

lib.MatrixScalarMultiply.argtypes = [ctypes.c_float, Matrix]
lib.MatrixScalarMultiply.restype = Matrix

lib.MatrixMultiply.argtypes = [Matrix, Matrix]
lib.MatrixMultiply.restype = Matrix

lib.MatrixTranspose.argtypes = [Matrix]
lib.MatrixTranspose.restype = Matrix

lib.AlgebraicCofactor.argtypes = [Matrix]
lib.AlgebraicCofactor.restype = Cofactor

lib.CofactorToMatrixSeries.argtypes = [Cofactor]
lib.CofactorToMatrixSeries.restype = MatrixSeries

lib.MatrixDeterminant.argtypes = [Matrix]
lib.MatrixDeterminant.restype = ctypes.c_float

lib.MatrixCofactor.argtypes = [Matrix, ctypes.c_int, ctypes.c_int]
lib.MatrixCofactor.restype = Matrix

lib.MatrixAdjoint.argtypes = [Matrix]
lib.MatrixAdjoint.restype = AdjointMatrix

lib.MatrixInverse.argtypes = [Matrix]
lib.MatrixInverse.restype = Matrix

# 矩阵分解运算
lib.MatrixDecompositionLU.argtypes = [Matrix]
lib.MatrixDecompositionLU.restype = LUMatrix

lib.MatrixDecompositionCholesky.argtypes = [Matrix]
lib.MatrixDecompositionCholesky.restype = CholeskyMatrix

lib.MatrixDecompositionQR.argtypes = [Matrix]
lib.MatrixDecompositionQR.restype = QRMatrix

lib.MatrixDecompositionSVD.argtypes = [Matrix, ctypes.c_float]
lib.MatrixDecompositionSVD.restype = SVDMatrix

lib.MatrixDecompositionEVD.argtypes = [Matrix, ctypes.c_float]
lib.MatrixDecompositionEVD.restype = EVDMatrix