# 相关编译指令
SOURCE_MATRIX_STORE := matrix_store.c
OBJS_MATRIX_STORE := matrix_store.o   
SOURCE_MATRIX_BASIC_OPERATION := matrix_basic_operation.c
OBJS_MATRIX_BASIC_OPERATION := matrix_basic_operation.o
TARGET := libFirstMatrixC.a


# 编译相关参数
CC := gcc
LIBS := 
LDFLAGS := -L
DEFINES :=
INCULDE :=
CFLAGS := 
CXXFLAGS := 


# 归档命令
AR := ar


# 编译操作
	cd ./src

$(OBJS_MATRIX_STORE):$(SOURCE_MATRIX_STORE)
	$(CC) $(SOURCE_MATRIX_STORE) -o $(OBJS_MATRIX_STORE) -c

$(OBJS_MATRIX_BASIC_OPERATION):$(SOURCE_MATRIX_BASIC_OPERATION)
	$(CC) $(SOURCE_MATRIX_BASIC_OPERATION) -o $(OBJS_MATRIX_BASIC_OPERATION) -c


# 链接操作


# 归档操作
$(AR) -rcs $(TARGET) $(OBJS_MATRIX_STORE) $(OBJS_MATRIX_BASIC_OPERATION)


# 清空缓存
	rm -rf *.o