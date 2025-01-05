import ctypes
from firstmatrixc.data_declaration import Matrix


# 辅助函数：将 Python 列表转换为 C 数组
def list_to_float_array(lst):
    return (ctypes.c_float * len(lst))(*lst)

# 辅助函数：将 C 数组转换为 Python 列表
def float_array_to_list(arr, length):
    return [arr[i] for i in range(length)]

# 创建 Matrix 对象
def create_matrix(rows, cols, data):
    matrix = Matrix()
    matrix.row = rows
    matrix.column = cols
    matrix.data = list_to_float_array(data)
    return matrix

