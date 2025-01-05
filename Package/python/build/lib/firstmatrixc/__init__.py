import ctypes
import os

# 加载共享库
_lib_path = os.path.join(os.path.dirname(__file__), "first_matrix_c.so")
lib = ctypes.CDLL(_lib_path)