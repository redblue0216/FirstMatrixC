# 定义编译器
CC = gcc

# 定义归档工具
AR = ar

# 定义编译选项
CFLAGS = -Wall

# 定义静态库的名称
LIB_NAME = libFirstMatrixC.a

# 定义源文件列表（假设所有源文件都在 src 目录下）
SOURCES = matrix_store.c matrix_basic_operation.c matrix_decomposition_operation.c
# 定义源文件的路径
SOURCES_PATH = ./src/

# 定义目标文件列表（.o 文件）
OBJECTS = $(patsubst %.c,$(SOURCES_PATH)%.o,$(SOURCES))

# 默认目标
all: $(LIB_NAME)

# 静态库规则
$(LIB_NAME): $(OBJECTS)
	$(AR) rcs $@ $^

# 从源文件生成目标文件的规则
$(SOURCES_PATH)%.o: $(SOURCES_PATH)%.c
	$(CC) $(CFLAGS) -c $< -o $@

# 清理规则
.PHONY: clean
clean:
	rm -f $(OBJECTS) $(LIB_NAME)