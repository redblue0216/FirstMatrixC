//模块名称: matrix_basic_operation
//设计：施华
//编码：施华
//简介：这是first_matrix_c矩阵库的矩阵基本运算模块，提供矩阵基本运算操作函数实现
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
**功能函数
*/



/*
**矩阵基本运算头文件
*/
#include "matrix_basic_operation.h"



// 矩阵加法
Matrix MatrixAdd(Matrix a,Matrix b)
/*
**函数功能：
**定义一个实现矩阵加法的函数
**参数：
**a (Matrix): a矩阵
**b (Matrix): b矩阵
**返回：
**c (Matrix): c矩阵
*/
{
    // 防御性编程-检查a,b矩阵维度是否一样
    if (a.row != b.row && a.column != b.column)
    {
        printf(Matrix_Addition_Error_001);
    }
    // 实现矩阵加法
    else
    {
        // step1: 创建一个新的矩阵
        Matrix c;
        c.row = a.row;
        c.column = a.column;
        c.data = (REAL *)malloc(a.row * a.column * sizeof(REAL));  // 动态分配内存
        // step2: 将对应坐标点元素相加后添加到矩阵存储Matrix的data中
        INDEX i,j;
        for (i=1;i<=a.row;i++)
        {   
            for (j=1;j<=a.column;j++)
            {
                INDEX c_index;
                c_index = TransformerIndex(a.column,i,j);
                c.data[c_index] = a.data[c_index] + b.data[c_index];
            }
        }
        return c;

    }
}


// 矩阵标量乘法
Matrix MatrixScalarMultiply(REAL scalar,Matrix a)
/*
**函数功能：
**定义一个矩阵标量乘法的函数
**参数：
**scalar (REAL): 标量
**a (Matrix): a矩阵
**返回：
**b (Matrix): b矩阵
*/
{
    // 防御性编程-忽略
    // 实现矩阵标量乘法
    // step1: 创建一个新的矩阵
    Matrix b;
    b.row = a.row;
    b.column = a.column;
    b.data = (REAL *)malloc(a.row * a.column * sizeof(REAL));  // 动态分配内存
    // step2: 将对应元素乘以标量后添加到矩阵存储Matrix的data中
    INDEX i,j;
    for (i=1;i<=a.row;i++)
    {
        for (j=1;j<=a.column;j++)
        {
            INDEX b_index;
            b_index = TransformerIndex(a.column,i,j);
            b.data[b_index] = a.data[b_index] * scalar;
        }
    }
    return b;
}


// 矩阵乘法
Matrix MatrixMultiply(Matrix a,Matrix b)
/*
**函数功能：
**定义一个矩阵乘法的函数
**参数：
**a (Matrix): a矩阵
**b (Matrix): b矩阵
**返回：
**c (Matrix): c矩阵
*/
{
    // 防御性编程-检查矩阵对应内部维度是否一致例如mn,nl
    if (a.column != b.row)
    {
        printf(Matrix_Muliply_Error_002);
    }
    // 实现矩阵乘法
    else{
        // step1: 创建一个新的矩阵，维度mn,nl->ml
        Matrix c;
        c.row = a.row;
        c.column = b.column;
        c.data = (REAL *)malloc(a.row * b.column * sizeof(REAL));  // 动态分配内存
        // step2: 根据矩阵乘法的维度变化ik,kj->ij进行对应行列元素相乘，组建新的矩阵
        INDEX i,j;
        for (i=1;i<=a.row;i++)
        {
            for (j=1;j<=b.column;j++)
            {   
                // 构建新矩阵角标
                INDEX c_index;
                c_index = TransformerIndex(b.column,i,j);  // 此处使用第二个矩阵的列作为转换窗口
                REAL c_element_ij = 0;
                INDEX k;
                // 根据角标n进行乘积求和计算
                for (k=1;k<=a.column;k++)
                {
                    // 根据角标获取矩阵元素
                    REAL a_ik;
                    REAL b_kj;  
                    a_ik = GetMatrixElement(a,i,k);
                    b_kj = GetMatrixElement(b,k,j);
                    c_element_ij += a_ik * b_kj;
                }
                // 将计算好的矩阵元素加载进新矩阵c
                c.data[c_index] = c_element_ij;
            }
        }
        return c;
    }
}


// 矩阵转置
Matrix MatrixTranspose(Matrix a)
/*
**函数功能：
**定义一个转置矩阵的函数
**参数：
**a (Matrix): a矩阵
**返回：
**b (Matrix): b矩阵
*/
{
    // 防御性编程-忽略
    // 实现矩阵转置
    // step1: 创建一个新的矩阵
    Matrix b;
    b.row = a.column;
    b.column = a.row;
    b.data = (REAL *)malloc(a.row * a.column * sizeof(REAL));  // 动态分配内存 
    // step2: 根据转置的维度变化(行列对调)，将对应元素加载进新建的矩阵中
    INDEX i,j;
    for (i=1;i<=a.row;i++)
    {
        for (j=1;j<=a.column;j++)
        {
            // 根据角标获取矩阵元素
            REAL b_ji;
            b_ji = GetMatrixElement(a,i,j);
            // 构建新矩阵角标
            INDEX b_index;
            b_index = TransformerIndex(b.column,j,i);
            //将获取的矩阵元素加载进新矩阵
            b.data[b_index] = b_ji;
        }
    }
    return b;
}


// 矩阵代数余子式计算
Cofactor AlgebraicCofactor(Matrix a)
/*
**函数功能：
**定义一个矩阵代数余子式计算的函数
**参数：
**a (Matrix): a矩阵
**返回：
**algebraic_cofactor (Cofactor): 代数余子式
*/
{
    // 防御性编程-检查矩阵是否为方阵
    if (a.row != a.column)
    {
        printf(Matrix_Square_Error_003);
    }
    // 实现代数余子式计算
    else
    {
        // step1: 根据矩阵a得到对应代数余子式中的行列式的阶
        Cofactor algebraic_cofactor;
        INTEGER algebraic_cofactor_order;
        algebraic_cofactor_order = a.row - 1;
        // step2: 根据阶数采取不同的计算方案。二阶行列式采取交叉乘积差，高阶直接使用一维数组存储元素
        if (algebraic_cofactor_order == 0)
        {
            // 确定余子式阶数
            algebraic_cofactor.order = 1;
            // 动态分配内存空间为逆序数、元素和余子式
            algebraic_cofactor.element = (REAL *)malloc(sizeof(REAL));
            algebraic_cofactor.reverse_order_number = (REAL *)malloc(sizeof(REAL));
            algebraic_cofactor.data = (REAL *)malloc(sizeof(REAL));
            // 加载数据到cofactor数据结构中
            algebraic_cofactor.data[0] = a.data[0];
            algebraic_cofactor.order = 1;
            algebraic_cofactor.element[0] = a.data[0];
            algebraic_cofactor.reverse_order_number[0] = 1;
            return algebraic_cofactor;
        }
        else if(algebraic_cofactor_order == 1)
        {
            // 确定余子式阶数
            algebraic_cofactor.order = 2;
            // 动态分配内存空间为逆序数、元素和余子式
            algebraic_cofactor.element = (REAL *)malloc(2 * sizeof(REAL));
            algebraic_cofactor.reverse_order_number = (REAL *)malloc(2 * sizeof(REAL));
            algebraic_cofactor.data = (REAL *)malloc(4 * sizeof(REAL));
            // 加载数据到cofactor数据结构中
            INDEX i,j;
            REAL a_ij;
            for (i=1;i<=2;i++)
            {
                for (j=1;j<=2;j++)
                {   
                    // 利用一维数组和元素角标整合代数余子式
                    a_ij = GetMatrixElement(a,i,j);
                    INDEX element_index;
                    element_index = TransformerIndex(a.column,i,j);
                    algebraic_cofactor.element[element_index] = a_ij;
                    algebraic_cofactor.reverse_order_number[element_index] = pow(-1,i+j);
                    INDEX data_index;
                    data_index = TransformerIndex(a.column,i,j);
                    algebraic_cofactor.data[data_index] = a_ij;
                }
            }
            return algebraic_cofactor;
        }
        else
        {   
            // 确定余子式阶数
            algebraic_cofactor.order = algebraic_cofactor_order;
            // 使用第一行元素作为代数余子式基础元素
            REAL a_1j;
            INDEX j;
            // 动态分配内存空间为逆序数、元素和余子式
            algebraic_cofactor.element = (REAL *)malloc(a.column * sizeof(REAL));
            algebraic_cofactor.reverse_order_number = (REAL *)malloc(a.column * sizeof(REAL));
            algebraic_cofactor.data = (REAL *)malloc(a.column * (a.column-1) * (a.row-1) * sizeof(REAL));
            // 加载数据到cofactor数据结构中
            INDEX data_index;  // 数据元素索引要放到所有循环外
            data_index = 0;
            for (j=1;j<=a.column;j++)
            {   
                // 利用一维数组和元素角标整合代数余子式
                a_1j = GetMatrixElement(a,1,j);
                INDEX element_index;
                element_index = TransformerIndex(a.column,1,j);
                algebraic_cofactor.element[element_index] = a_1j;
                algebraic_cofactor.reverse_order_number[element_index] = pow(-1,1+j);
                INDEX m,n;
                for (m=2;m<=a.row;m++)
                {   
                    for (n=1;n<=a.column;n++)
                    {   
                        if (n==j)
                        {
                            continue;
                        }
                        else
                        {
                            REAL a_mn;
                            a_mn = GetMatrixElement(a,m,n);
                            algebraic_cofactor.data[data_index] = a_mn;
                            // printf("======%d,%d======%f======%d\n",m,n,a_mn,data_index);
                            data_index += 1;   
                        }   
                    }
                }
            }
            return algebraic_cofactor;
        }
    }
}


// 代数余子式转变为矩阵数组
MatrixSeries CofactorToMatrixSeries(Cofactor b)
/*
**函数功能：
**定义一个将代数余子式Cofator数据结构转化为矩阵数组的函数
**参数：
**b (Cofactor): 代数余子式b
**返回：
**c (MatrixSeries): 代数余子式的矩阵数组
*/
{   
    // 创建一个矩阵序列
    MatrixSeries c;
    INTEGER b_matrix_number = b.order + 1;
    // 动态分配内存
    c.element = (REAL *)malloc((b.order+1) * sizeof(REAL));
    c.reverse_order_number = (REAL *)malloc((b.order+1) * sizeof(REAL));
    Matrix tmp_b_i;
    tmp_b_i.row = b.order;
    tmp_b_i.column = b.order;
    tmp_b_i.data = (REAL *)malloc(tmp_b_i.row * tmp_b_i.column * sizeof(REAL)); 
    c.data = (Matrix *)malloc((b.order+1) * sizeof(tmp_b_i));
    // 加载矩阵序列基础变量，包括元素、逆序数和矩阵个数
    c.element = b.element;
    c.reverse_order_number = b.reverse_order_number;    
    c.matrix_number = b_matrix_number;
    c.order = b.order;
    // 使用循环将矩阵数据加载进矩阵序列
    INDEX i;
    for (i=0;i<b_matrix_number;i++)
        {
            // 从代数余子式中获取对应矩阵角标
            Matrix b_i;
            b_i.row = b.order;
            b_i.column = b.order;
            b_i.data = (REAL *)malloc(b_i.row * b_i.column * sizeof(REAL));  // 动态分配内存
            INDEX j;
            for (j=0;j<b_i.row * b_i.column;j++)
            {
                INDEX ij;
                ij = i * b_i.row * b_i.column + j;
                // 根据矩阵角标获取矩阵b_i
                b_i.data[j] = b.data[ij];
            }
            // 将获取的矩阵加载进矩阵序列c中
            c.data[i] = b_i;
        }
    // 返回
    return c;
}


// 矩阵行列式计算
REAL MatrixDeterminant(Matrix a)
/*
**函数功能：
**定义一个矩阵行列式计算的函数
**参数：
**a (Matrix): 矩阵a
**返回：
**determiant_a (REAL): 行列式值
*/
{   
    // 防御性编程-检查矩阵是否为方阵
    if (a.row != a.column)
    {
        printf(Matrix_Square_Error_003);
    }
    // 实现代数余子式计算
    else
    {    
        // 使用递归技术实现Laplace行列式计算
        // step0: 声明a矩阵的行列式值
        REAL determinant_a;
        // step1: 先做一次代数余子式的计算，为了保证一、二阶和高阶的cofactor数据类型保持一致
        if (a.column == 1)
        {
            determinant_a = a.data[0];
        }
        else if (a.column == 2)
        {
            determinant_a = a.data[0] * a.data[3] - a.data[1] * a.data[2];
        }
        // step2: 使用递归技术进行高阶代数余子式的迭代计算，以得到最终矩阵行列式的值
        else
        {
            // 使用递归技术，在函数中循环自己调用自己---->>关键编程点
            INDEX j;
            // 使用初始零值加循环实现求和操作
            determinant_a = 0;
            for (j=0;j<a.column;j++)
            {   
                // 根据矩阵a获得代数余子式cofactor
                Cofactor cofator_c;
                cofator_c = AlgebraicCofactor(a);
                // 将代数余子式转化为矩阵序列
                MatrixSeries c_matrix_series;
                c_matrix_series = CofactorToMatrixSeries(cofator_c);
                // 使用递归技术实现Laplace公式计算
                determinant_a += c_matrix_series.reverse_order_number[j] * c_matrix_series.element[j] * MatrixDeterminant(c_matrix_series.data[j]);
            }
        }
        // 返回矩阵行列式的计算值
        return determinant_a;    
    }
}


// 矩阵余子式计算
Matrix MatrixCofactor(Matrix a,INDEX i,INDEX j)
/*
**函数功能：
**定义一个矩阵余子式计算的函数
**
**参数：
**a (Matrix): 矩阵
**i (INDEX): 矩阵行
**j (INDEX): 矩阵列
**返回：
**b (Matrix): 余子式矩阵b
*/
{
    // 创建一个余子式矩阵
    Matrix b;
    b.row = a.row - 1;
    b.column = a.column - 1;
    b.data = (REAL *)malloc(b.row * b.column * sizeof(REAL));
    // 在双循环中使用if语句跳过目标元素所在的行和列
    INDEX r,c;
    REAL a_rc;
    INDEX data_index;
    data_index = 0;
    for (r=1;r<=a.row;r++)
    {
        if (r == i)
        {
            continue;
        }
        else
        {
            for (c=1;c<=a.column;c++)
            {
                if (c == j)
                {
                    continue;
                }
                else
                {
                    // 根据给定坐标获取矩阵元素
                    a_rc = GetMatrixElement(a,r,c);
                    // 加载进矩阵b中
                    b.data[data_index] = a_rc;
                    // 矩阵b的一维data索引加一
                    data_index += 1;
                }
            }
        }
    }
    return b;
}


// 伴随矩阵计算
AdjointMatrix MatrixAdjoint(Matrix a)
/*
**函数功能：
**定义一个伴随矩阵计算函数  
**参数：
**a (Matrix): 矩阵a
**返回：
**b_matrix_matrix: 伴随矩阵(矩阵的矩阵)
*/
{   
    // 创建一个伴随矩阵
    AdjointMatrix c_adjoint_matrix;
    // 动态分配内存
    c_adjoint_matrix.reverse_order_number = (REAL *)malloc(a.column * a.row * sizeof(REAL));
    c_adjoint_matrix.data = (REAL *)malloc(a.column * a.row * sizeof(REAL));
    // 遍历矩阵a计算对应的余子式和逆序数并将余子式加载到矩阵的矩阵b_matrix_matirx中
    INDEX i,j;
    INDEX data_index;
    data_index = 0;
    for (i=1;i<=a.row;i++)
    {
        for (j=1;j<=a.column;j++)
        {
            // 根据角标获取对应的余子式
            Matrix a_cofactor;
            a_cofactor = MatrixCofactor(a,i,j);
            // 为伴随矩阵添加逆序数
            c_adjoint_matrix.reverse_order_number[data_index] = pow(-1,i+j);
            // 计算目标子余子式矩阵c_adjoint_matrix的行列式值，并加载进伴随矩阵矩阵c_adjoint_matrix
            REAL determinant_c_adjoint_matrix;
            determinant_c_adjoint_matrix = MatrixDeterminant(a_cofactor);
            c_adjoint_matrix.data[data_index] = determinant_c_adjoint_matrix;
            // 以下注释用于显示子余子式的矩阵元素
            // INDEX m;
            // printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%d,%d\n",i,j);
            // for (m=0;m<a_cofactor.column*a_cofactor.row;m++)
            // {
            //     printf("~~~~~~~~~~~~~%f\n",a_cofactor.data[m]);
            // }
            // 索引加一
            data_index += 1;
        }
    }
    // 为伴随矩阵添加矩阵行列式的值
    REAL determinant_a;
    determinant_a = MatrixDeterminant(a);
    c_adjoint_matrix.matrix_determinant = determinant_a;
    // 为伴随矩阵添加子余子式的阶数
    c_adjoint_matrix.order = a.row - 1;
    // 返回伴随矩阵
    return c_adjoint_matrix;
}


// 矩阵求逆计算
Matrix MatrixInverse(Matrix a)
/*
**函数功能：
**定义一个矩阵求逆计算函数
**参数：
**a (Matrix): 矩阵a
**返回：
**b (Matrix): 矩阵b
*/
{
    // 创建一个逆矩阵
    Matrix b;
    b.row = a.row;
    b.column = a.column;
    b.data = (REAL *)malloc(a.row * a.column * sizeof(REAL));
    // 利用伴随矩阵和行列式进行矩阵求逆计算
    AdjointMatrix adjoint_matrix_a;
    adjoint_matrix_a = MatrixAdjoint(a);
    REAL determinant_a;
    determinant_a = adjoint_matrix_a.matrix_determinant;
    // 使用循环计算逆矩阵的元素，并加载进逆矩阵b中
    INDEX i;
    for (i=0;i<a.row*a.column;i++)
    {
        //  使用伴随矩阵、行列式和逆序数进行逆矩阵元素计算
        REAL sign,data,inverse_matrix_element;
        sign = adjoint_matrix_a.reverse_order_number[i];
        data = adjoint_matrix_a.data[i];
        inverse_matrix_element = sign * data / determinant_a;
        b.data[i] = inverse_matrix_element;
    }
    // 返回逆矩阵
    return b;
}
/********************************************************************************************/
/********************************************************************************************/


