//模块名称: matrix_basic_operation
//设计：施华
//编码：施华
//简介：这是first_matrix_c矩阵库的矩阵分解运算模块，提供矩阵分解运算操作函数实现
//备注：无



/*
**功能介绍
**
**提供矩阵基本运算操作，包括
**（1）矩阵LU分解
**（2）矩阵Cholesky分解
**（3）矩阵QR分解
**（4）矩阵SVD分解
**（5）矩阵特征值分解
**
*********************************************************
**
**技术介绍
**
**功能函数
*/



/*
**矩阵分解运算头文件
*/
#include "matrix_basic_operation.h"



// LU分解
LUMatrix MatrixDecompositionLU(Matrix a)
/*
**函数功能：
**定义一个LU矩阵分解函数  
**参数：
**a (Matrix): 矩阵a
**返回：
**lu_matrix (LUMatrix): LU矩阵分解结果数据结构
*/
{
    // 声明索引变量和相关数据结构
    INTEGER n,i,j,k,nn;
    REAL temp1,temp2;
    LUMatrix lu_matrix;
    nn = a.column;
    REAL U_Matrix_element_nj;
    REAL L_Matrix_element_in;
    // 创建LU矩阵数据结构
    Matrix L_Matrix = CreateMatrix(nn,nn);
    L_Matrix.data = (REAL *)malloc(nn * nn * sizeof(REAL));
    Matrix U_Matrix = CreateMatrix(nn,nn);
    U_Matrix.data = (REAL *)malloc(nn * nn * sizeof(REAL));
    // 使用递归等式计算具体元素，开始矩阵LU分解运算
    for (i=1;i<=nn;i++)
        for (j=1;j<=nn;j++)
        {
            SetMatrixElement(L_Matrix,i,j,0);
            SetMatrixElement(U_Matrix,i,j,0);
        }
    for (n=1;n<=nn;n++)
    {
        for (j=n;j<=nn;j++)
        {
            temp1 = 0.0;
            for (k=1;k<=n-1;k++)
            {
                temp1 = temp1 + GetMatrixElement(L_Matrix,n,k) * GetMatrixElement(U_Matrix,k,j);
            }
            U_Matrix_element_nj = GetMatrixElement(a,n,j) - temp1;
            SetMatrixElement(U_Matrix,n,j,U_Matrix_element_nj);
        }
        for (i=n+1;i<=nn;i++)
        {
            temp2 = 0.0;
            for (k=1;k<=n-1;k++)
            {
                temp2 = temp2 + GetMatrixElement(L_Matrix,i,k) * GetMatrixElement(U_Matrix,k,n);
            }
            L_Matrix_element_in = (GetMatrixElement(a,i,n) - temp2) / GetMatrixElement(U_Matrix,n,n);
            SetMatrixElement(L_Matrix,i,n,L_Matrix_element_in);
        }
    }
    for (i=1;i<=nn;i++)
    {
        SetMatrixElement(L_Matrix,i,i,1.0);
    }
    lu_matrix.L = L_Matrix;
    lu_matrix.U = U_Matrix;
    
    return lu_matrix;
}


// Cholesky分解
CholeskyMatrix MatrixDecompositionCholesky(Matrix a)
/*
**函数功能：
**定义一个cholesky分解函数
**参数：
**a (Matrix): 矩阵a
**返回:
**cholesky_matrix (CholeskMatrix): cholesky矩阵数据结构
*/
{
    // 声明索引变量和相关数据结构
    INTEGER i,j,k,u,l,n;
    CholeskyMatrix cholesky_matrix;
    n = a.column;
    // 创建cholesky矩阵数据结构
    Matrix ac = CreateMatrix(n,n);
    ac.data = (REAL *)malloc(n * n * sizeof(REAL));
    // 复制a矩阵到ac矩阵中，避免修改原始矩阵
    for (i=0;i<=ac.row*ac.column;i++)
    {
        ac.data[i] = a.data[i];
    }
    // 使用递归等式计算具体元素，开始矩阵cholesky运算
    if ((ac.data[0] + 1.0 == 1.0) || (ac.data[0] < 0.0))
    {
        printf("Cholesky decomposition failed! \n");
        cholesky_matrix.L = ac;
        return cholesky_matrix;
    }
    ac.data[0] = sqrt(ac.data[0]);
    for (i=1;i<=n-1;i++)
    {
        u = i * n;
        ac.data[u] = ac.data[u] / ac.data[0];
    }
    for (j=1;j<=n-1;j++)
    {
        l = j * n + j;
        for (k=0;k<=j-1;k++)
        {
            u = j * n + k;
            ac.data[l] = ac.data[l] - ac.data[u] * ac.data[u];
        }
        if ((ac.data[l] + 1.0 == 1.0) || (ac.data[l] < 0.0))
        {
            printf("Cholesky decomposition failed! \n");
            cholesky_matrix.L = ac;
            return cholesky_matrix;
        }
        ac.data[l] = sqrt(ac.data[l]);
        for (i=j+1;i<=n-1;i++)
        {
            u = i * n + j;
            for (k=0;k<=j-1;k++)
            {
                ac.data[u] = ac.data[u] - ac.data[i*n+k] * ac.data[j*n+k];
            }
            ac.data[u] = ac.data[u] / ac.data[l];
        }
    }
    for (i=0;i<=n-2;i++)
    {
        for (j=i+1;j<=n-1;j++)
        {
            ac.data[i*n+j] = 0.0;
        }
    }
    cholesky_matrix.L = ac;
    return cholesky_matrix;
}


// QR分解
QRMatrix MatrixDecompositionQR(Matrix a)
/*
**函数功能：
**定义一个QR分解函数
**参数：
**a (Matrix): 矩阵a
**返回:
**qr_matrix (QRMatrix): QR矩阵数据结构
*/
{
    // 声明索引变量和相关数据结构
    INTEGER i,j,k,l,nn,p,jj,m,n;
    REAL u,alpha,w,t;
    QRMatrix qr_matrix;
    m = a.row;
    n = a.column;
    // 创建QR矩阵数据结构
    Matrix aq = CreateMatrix(m,n);
    aq.data = (REAL *)malloc(m * n * sizeof(REAL));
    Matrix q = CreateMatrix(m,m);
    q.data = (REAL *)malloc(m * m * sizeof(REAL));
    // 复制a矩阵到ac矩阵中，避免修改原始矩阵
    for (i=0;i<=aq.row*aq.column;i++)
    {
        aq.data[i] = a.data[i];
    }
    // 使用递归等式计算具体元素，开始矩阵QR运算
    if (m<n)
    {
        printf("QR decomposition failed! \n");
    }
    for (i=0;i<=m-1;i++)
    {
        for (j=0;j<=m-1;j++)
        {
            l = i * m + j;
            q.data[l] = 0.0;
            if (i==j)
            {
                q.data[l] = 1.0;
            }
        }
    }
    nn=n;
    if (m==n)
    {
        nn = m - 1;
    }
    for (k=0;k<=nn-1;k++)
    {
        u = 0.0;
        l = k * n + k;
        for (i=k;i<=m-1;i++)
        {
            w = fabs(aq.data[i*n+k]);
            if (w>u)
            {
                u=w;
            }
        }
        alpha = 0.0;
        for (i=k;i<=m-1;i++)
        {
            t = aq.data[i*n+k]/u;
            alpha = alpha + t * t;
        }
        if (aq.data[l]>0.0)
        {
            u = -u;
        }
        alpha = u * sqrt(alpha);
        if (fabs(alpha)+1.0==1.0)
        {
            printf("QR decomposition failed! \n");
        }
        u = sqrt(2.0 * alpha * (alpha - aq.data[l]));
        if ((u+1.0)!=1.0)
        {
            aq.data[l] = (aq.data[l] - alpha) / u;
            for (i=k+1;i<=m-1;i++)
            {
                p = i * n + k;
                aq.data[p] = aq.data[p] / u;
            }            
            for (j=0;j<=m-1;j++)
            {
                t = 0.0;
                for (jj=k;jj<=m-1;jj++)
                {
                    t = t + aq.data[jj*n+k]*q.data[jj*m+j];
                }
                for (i=k;i<=m-1;i++)
                {
                    p = i * m + j;
                    q.data[p] = q.data[p] - 2.0 * t * aq.data[i*n+k];
                }
            }
            for (j=k+1;j<=n-1;j++)
            {
                t=0.0;
                for (jj=k;jj<=m-1;jj++)
                {
                    t = t + aq.data[jj*n+k] * aq.data[jj*n+j];
                }
                for (i=k;i<=m-1;i++)
                {
                    p = i * n + j;
                    aq.data[p] = aq.data[p] - 2.0 * t * aq.data[i * n + k];
                }
            }
            aq.data[l] = alpha;
            for (i=k+1;i<=m-1;i++)
            {
                aq.data[i*n+k] = 0.0;
            }
        }
    }
    for (i=0;i<=m-2;i++)
    {
        for (j=i+1;j<=m-1;j++)
        {
            p = i * m + j;
            l = j * m + i;
            t = q.data[p];
            q.data[p] = q.data[l];
            q.data[l] = t;
        }
    }
    qr_matrix.Q = q;
    qr_matrix.R = aq;
    return qr_matrix;
}



/*********************************************************************************************/
/*********************************************************************************************/


