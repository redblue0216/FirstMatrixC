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


// SVD分解-ppp函数
VOID ppp(Matrix a,REAL *e,REAL *s,Matrix v,INTEGER m,INTEGER n)
/*
**函数功能：
**定义一个SVD辅助ppp函数
**参数：
**a (Matrix): 矩阵a
**e (Array): 数组e
**s (Array): 数组s
**v (Matrix): 矩阵v
**m (Int): 索引m
**n (Int): 索引n
**返回:
**无
*/
{
    // 声明索引变量和相关数据结构
    INTEGER i,j,p,q;
    REAL d;
    if (m>=n)
    {
        i=n;
    }
    else
    {
        i=m;
    }
    for (j=1;j<=i-1;j++)
    {
        a.data[(j-1)*n+j-1] = s[j-1];
        a.data[(j-1)*n+j] = e[j-1];
    }
    a.data[(i-1)*n+i-1] = s[i-1];
    if (m<n)
    {
        a.data[(i-1)*n+i] = e[i-1];
    }
    for (i=1;i<=n-1;i++)
    {
        for (j=i+1;j<=n;j++)
        {
            p=(i-1)*n+j-1;
            q=(j-1)*n+i-1;
            d = v.data[p];
            v.data[p] = v.data[q];
            v.data[q] = d;
        }
    }
    return;
}


// SVD分解-sss函数
VOID sss(REAL fg[2],REAL cs[2])
/*
**函数功能：
**定义一个SVD辅助ss函数
**参数：
**fg (Array): 数组fg
**cs (Array): 数组cs
**返回:
**无
*/
{
    // 声明索引变量和相关数据结构
    REAL r,d;
    if ((fabs(fg[0])+fabs(fg[1]))==0.0)
    {
        cs[0]=1.0;
        cs[1]=0.0;
        d=0.0;
    }
    else
    {
        d=sqrt(fg[0]*fg[0]+fg[1]*fg[1]);
        if (fabs(fg[0])>fabs(fg[1]))
        {
            d=fabs(d);
            if (fg[0]<0.0)
            {
                d=-d;
            }
        }
        if (fabs(fg[1]>=fabs(fg[0])))
        {
            d=fabs(d);
            if (fg[1]<0.0)
            {
                d=-d;
            }
        }
        cs[0]=fg[0]/d;
        cs[1]=fg[1]/d;
    }
    r=1.0;
    if (fabs(fg[0]>fabs(fg[1])))
    {
        r=cs[1];
    }
    else
    {
        if (cs[0]!=0.0)
        {
            r=1.0/cs[0];
        }
    }
    fg[0]=d;
    fg[1]=r;
    return;
}


// SVD分解
SVDMatrix MatrixDecompositionSVD(Matrix a,REAL eps)
/*
**函数功能：
**定义一个SVD矩阵分解函数
**参数：
**a (Matrix): 矩阵a
**eps (REAL): 精度eps
**返回:
**svd_matrix (SVDMatrix): SVD矩阵
*/
{ 
    // 声明索引变量和相关数据结构
    INTEGER i,j,k,l,it,ll,kk,ix,iy,mm,nn,iz,m1,ks;
    REAL d,dd,t,sm,sm1,em1,sk,ek,b,c,shh,fg[2],cs[2];
    REAL *s,*e,*w;
    REAL u_data_ij,v_data_ij;
    INTEGER m,n,ka;
    m = a.row;
    n = a.column;
    ka = MAX(m,n);
    Matrix u = CreateMatrix(m,m);
    u.data = (REAL *)malloc(m * m * sizeof(REAL));
    Matrix v = CreateMatrix(n,n);
    v.data = (REAL *)malloc(n * n * sizeof(REAL));
    // 创建缓存矩阵数据结构
    Matrix aa = CreateMatrix(a.row,a.column);
    aa.data = (REAL *)malloc(a.row * a.column * sizeof(REAL));
    SVDMatrix svd_matrix;
    // 用零值填充U,V矩阵
    SetMatrixZero(u);
    SetMatrixZero(v);
    // 复制a矩阵到ac矩阵中，避免修改原始矩阵
    for (i=0;i<=aa.row*aa.column;i++)
    {
        aa.data[i] = a.data[i];
    }
    // 重置SVD矩阵
    svd_matrix.A = aa;
    svd_matrix.U = u;
    svd_matrix.V = v;
    // 开始运算
    s=(REAL *)malloc(ka * sizeof(REAL));
    e=(REAL *)malloc(ka * sizeof(REAL));
    w=(REAL *)malloc(ka * sizeof(REAL));
    it=60; 
    k=n;
    if (m-1<n) 
        {
            k=m-1;
        }
    l=m;
    if (n-2<m) 
        {
            l=n-2;
        }
    if (l<0) 
        {
            l=0;
        }
    ll=k;
    if (l>k) 
        {
            ll=l;
        }
    if (ll>=1)
        { for (kk=1; kk<=ll; kk++)
            { if (kk<=k)
                { d=0.0;
                for (i=kk; i<=m; i++)
                    { ix=(i-1)*n+kk-1; d=d+aa.data[ix]*aa.data[ix];}
                s[kk-1]=sqrt(d);
                if (s[kk-1]!=0.0)
                    { ix=(kk-1)*n+kk-1;
                    if (aa.data[ix]!=0.0)
                        { s[kk-1]=fabs(s[kk-1]);
                        if (aa.data[ix]<0.0) 
                            {s[kk-1]=-s[kk-1];}
                        }
                    for (i=kk; i<=m; i++)
                        { iy=(i-1)*n+kk-1;
                        aa.data[iy]=aa.data[iy]/s[kk-1];
                        }
                    aa.data[ix]=1.0+aa.data[ix];
                    }
                s[kk-1]=-s[kk-1];
                }
            if (n>=kk+1)
                { for (j=kk+1; j<=n; j++)
                    { if ((kk<=k)&&(s[kk-1]!=0.0))
                        { d=0.0;
                        for (i=kk; i<=m; i++)
                            { ix=(i-1)*n+kk-1;
                            iy=(i-1)*n+j-1;
                            d=d+aa.data[ix]*aa.data[iy];
                            }
                        d=-d/aa.data[(kk-1)*n+kk-1];
                        for (i=kk; i<=m; i++)
                            { ix=(i-1)*n+j-1;
                            iy=(i-1)*n+kk-1;
                            aa.data[ix]=aa.data[ix]+d*aa.data[iy];
                            }
                        }
                    e[j-1]=aa.data[(kk-1)*n+j-1];
                    }
                }
            if (kk<=k)
                { for (i=kk; i<=m; i++)
                    { ix=(i-1)*m+kk-1; iy=(i-1)*n+kk-1;
                    u.data[ix]=aa.data[iy];
                    }
                }
            if (kk<=l)
                { d=0.0;
                for (i=kk+1; i<=n; i++)
                    {d=d+e[i-1]*e[i-1];}
                e[kk-1]=sqrt(d);
                if (e[kk-1]!=0.0)
                    { if (e[kk]!=0.0)
                        { e[kk-1]=fabs(e[kk-1]);
                        if (e[kk]<0.0) 
                            {e[kk-1]=-e[kk-1];}
                        }
                    for (i=kk+1; i<=n; i++)
                        {e[i-1]=e[i-1]/e[kk-1];}
                    e[kk]=1.0+e[kk];
                    }
                e[kk-1]=-e[kk-1];
                if ((kk+1<=m)&&(e[kk-1]!=0.0))
                    { 
                    for (i=kk+1; i<=m; i++) 
                        {w[i-1]=0.0;}
                    for (j=kk+1; j<=n; j++)
                        {for (i=kk+1; i<=m; i++)
                            {w[i-1]=w[i-1]+e[j-1]*aa.data[(i-1)*n+j-1];}}
                    for (j=kk+1; j<=n; j++)
                        {for (i=kk+1; i<=m; i++)
                        { ix=(i-1)*n+j-1;
                            aa.data[ix]=aa.data[ix]-w[i-1]*e[j-1]/e[kk];
                        }}
                    }
                for (i=kk+1; i<=n; i++)
                    {v.data[(i-1)*n+kk-1]=e[i-1];}
                }
            }
        }
    mm=n;
    if (m+1<n) 
        {
            mm=m+1;
        }
    if (k<n) 
        {
            s[k]=aa.data[k*n+k];
        }
    if (m<mm) 
        {
            s[mm-1]=0.0;
        }
    if (l+1<mm) 
        {
            e[l]=aa.data[l*n+mm-1];
        }
    e[mm-1]=0.0;
    nn=m;
    if (m>n) 
        {
            nn=n;
        }
    if (nn>=k+1)
        { for (j=k+1; j<=nn; j++)
            { for (i=1; i<=m; i++)
                {u.data[(i-1)*m+j-1]=0.0;}
            u.data[(j-1)*m+j-1]=1.0;
            }
        }
    if (k>=1)
        { for (ll=1; ll<=k; ll++)
            { kk=k-ll+1; iz=(kk-1)*m+kk-1;
            if (s[kk-1]!=0.0)
                { if (nn>=kk+1)
                    for (j=kk+1; j<=nn; j++)
                    { d=0.0;
                        for (i=kk; i<=m; i++)
                        { ix=(i-1)*m+kk-1;
                          iy=(i-1)*m+j-1;
                            d=d+u.data[ix]*u.data[iy]/u.data[iz];
                        }
                        d=-d;
                        for (i=kk; i<=m; i++)
                        { ix=(i-1)*m+j-1;
                            iy=(i-1)*m+kk-1;
                            u.data[ix]=u.data[ix]+d*u.data[iy];
                        }
                    }
                    for (i=kk; i<=m; i++)
                    { ix=(i-1)*m+kk-1; u.data[ix]=-u.data[ix];}
                    u.data[iz]=1.0+u.data[iz];
                    if (kk-1>=1)
                        {for (i=1; i<=kk-1; i++)
                            {u.data[(i-1)*m+kk-1]=0.0;}}
                }
            else
                { for (i=1; i<=m; i++)
                    {u.data[(i-1)*m+kk-1]=0.0;}
                u.data[(kk-1)*m+kk-1]=1.0;
                }
            }
        }
    for (ll=1; ll<=n; ll++)
        { kk=n-ll+1; iz=kk*n+kk-1;
        if ((kk<=l)&&(e[kk-1]!=0.0))
            { for (j=kk+1; j<=n; j++)
                { d=0.0;
                for (i=kk+1; i<=n; i++)
                    { ix=(i-1)*n+kk-1; iy=(i-1)*n+j-1;
                    d=d+v.data[ix]*v.data[iy]/v.data[iz];
                    }
                d=-d;
                for (i=kk+1; i<=n; i++)
                    { ix=(i-1)*n+j-1; iy=(i-1)*n+kk-1;
                    v.data[ix]=v.data[ix]+d*v.data[iy];
                    }
                }
            }
        for (i=1; i<=n; i++)
            {v.data[(i-1)*n+kk-1]=0.0;}
        v.data[iz-n]=1.0;
        }
    for (i=1; i<=m; i++)
        {for (j=1; j<=n; j++)
            {aa.data[(i-1)*n+j-1]=0.0;}}
    m1=mm; it=60;
    while (l==1)
        { if (mm==0)
            { ppp(aa,e,s,v,m,n);
            free(s); free(e); free(w); 
            return svd_matrix;
            }
        if (it==0)
            { ppp(aa,e,s,v,m,n);
            free(s); free(e); free(w);
            return svd_matrix;
            }
        kk=mm-1;
    while ((kk!=0)&&(fabs(e[kk-1])!=0.0))
            { d=fabs(s[kk-1])+fabs(s[kk]);
            dd=fabs(e[kk-1]);
            if (dd>eps*d) {kk=kk-1;}
            else {e[kk-1]=0.0;}
            }
        if (kk==mm-1)
            { kk=kk+1;
            if (s[kk-1]<0.0)
                { s[kk-1]=-s[kk-1];
                for (i=1; i<=n; i++)
                    { ix=(i-1)*n+kk-1; v.data[ix]=-v.data[ix];}
                }
            while ((kk!=m1)&&(s[kk-1]<s[kk]))
                { d=s[kk-1]; s[kk-1]=s[kk]; s[kk]=d;
                if (kk<n)
                    {for (i=1; i<=n; i++)
                    { ix=(i-1)*n+kk-1; iy=(i-1)*n+kk;
                        d=v.data[ix]; v.data[ix]=v.data[iy]; v.data[iy]=d;
                    }}
                if (kk<m)
                    {for (i=1; i<=m; i++)
                    { ix=(i-1)*m+kk-1; iy=(i-1)*m+kk;
                        d=u.data[ix]; u.data[ix]=u.data[iy]; u.data[iy]=d;
                    }}
                kk=kk+1;
                }
            it=60;
            mm=mm-1;
            }
        else
            { ks=mm;
            while ((ks>kk)&&(fabs(s[ks-1])!=0.0))
                { d=0.0;
                if (ks!=mm) {d=d+fabs(e[ks-1]);}
                if (ks!=kk+1) {d=d+fabs(e[ks-2]);}
                dd=fabs(s[ks-1]);
                if (dd>eps*d) {ks=ks-1;}
                else {s[ks-1]=0.0;}
                }
            if (ks==kk)
                { kk=kk+1;
                d=fabs(s[mm-1]);
                t=fabs(s[mm-2]);
                if (t>d) {d=t;}
                t=fabs(e[mm-2]);
                if (t>d) {d=t;}
                t=fabs(s[kk-1]);
                if (t>d) {d=t;}
                t=fabs(e[kk-1]);
                if (t>d) {d=t;}
                sm=s[mm-1]/d; sm1=s[mm-2]/d;
                em1=e[mm-2]/d;
                sk=s[kk-1]/d; ek=e[kk-1]/d;
                b=((sm1+sm)*(sm1-sm)+em1*em1)/2.0;
                c=sm*em1; c=c*c; shh=0.0;
                if ((b!=0.0)||(c!=0.0))
                    { shh=sqrt(b*b+c);
                    if (b<0.0) {shh=-shh;}
                    shh=c/(b+shh);
                    }
                fg[0]=(sk+sm)*(sk-sm)-shh;
                fg[1]=sk*ek;
                for (i=kk; i<=mm-1; i++)
                    { sss(fg,cs);
                    if (i!=kk) {e[i-2]=fg[0];}
                    fg[0]=cs[0]*s[i-1]+cs[1]*e[i-1];
                    e[i-1]=cs[0]*e[i-1]-cs[1]*s[i-1];
                    fg[1]=cs[1]*s[i];
                    s[i]=cs[0]*s[i];
                    if ((cs[0]!=1.0)||(cs[1]!=0.0))
                        {for (j=1; j<=n; j++)
                        { ix=(j-1)*n+i-1;
                            iy=(j-1)*n+i;
                            d=cs[0]*v.data[ix]+cs[1]*v.data[iy];
                            v.data[iy]=-cs[1]*v.data[ix]+cs[0]*v.data[iy];
                            v.data[ix]=d;
                        }}
                    sss(fg,cs);
                    s[i-1]=fg[0];
                    fg[0]=cs[0]*e[i-1]+cs[1]*s[i];
                    s[i]=-cs[1]*e[i-1]+cs[0]*s[i];
                    fg[1]=cs[1]*e[i];
                    e[i]=cs[0]*e[i];
                    if (i<m)
                        {if ((cs[0]!=1.0)||(cs[1]!=0.0))
                        for (j=1; j<=m; j++)
                            { ix=(j-1)*m+i-1;
                            iy=(j-1)*m+i;
                            d=cs[0]*u.data[ix]+cs[1]*u.data[iy];
                            u.data[iy]=-cs[1]*u.data[ix]+cs[0]*u.data[iy];
                            u.data[ix]=d;
                            }}
                    }
                e[mm-2]=fg[0];
                it=it-1;
                }
            else
                { if (ks==mm)
                    { kk=kk+1;
                    fg[1]=e[mm-2]; e[mm-2]=0.0;
                    for (ll=kk; ll<=mm-1; ll++)
                        { i=mm+kk-ll-1;
                        fg[0]=s[i-1];
                        sss(fg,cs);
                        s[i-1]=fg[0];
                        if (i!=kk)
                            { fg[1]=-cs[1]*e[i-2];
                            e[i-2]=cs[0]*e[i-2];
                            }
                        if ((cs[0]!=1.0)||(cs[1]!=0.0))
                            {for (j=1; j<=n; j++)
                            { ix=(j-1)*n+i-1;
                                iy=(j-1)*n+mm-1;
                                d=cs[0]*v.data[ix]+cs[1]*v.data[iy];
                                v.data[iy]=-cs[1]*v.data[ix]+cs[0]*v.data[iy];
                                v.data[ix]=d;
                            }}
                        }
                    }
                else
                    { kk=ks+1;
                    fg[1]=e[kk-2];
                    e[kk-2]=0.0;
                    for (i=kk; i<=mm; i++)
                        { fg[0]=s[i-1];
                        sss(fg,cs);
                        s[i-1]=fg[0];
                        fg[1]=-cs[1]*e[i-1];
                        e[i-1]=cs[0]*e[i-1];
                        if ((cs[0]!=1.0)||(cs[1]!=0.0))
                            {for (j=1; j<=m; j++)
                            { ix=(j-1)*m+i-1;
                                iy=(j-1)*m+kk-2;
                                d=cs[0]*u.data[ix]+cs[1]*u.data[iy];
                                u.data[iy]=-cs[1]*u.data[ix]+cs[0]*u.data[iy];
                                u.data[ix]=d;
                            }}
                        }
                    }
                }
            }
        }
    // // 修改符号
    // for (i=1;i<u.row;i++)
    // {
    //     for (j=1;j<u.column;j++)
    //     {
    //         u_data_ij = GetMatrixElement(u,i,j);
    //         SetMatrixElement(u,i,j,-u_data_ij);
    //     }
    // }
    // for (i=1;i<v.row;i++)
    // {
    //     for (j=1;j<v.column;j++)
    //     {
    //         v_data_ij = GetMatrixElement(v,i,j);
    //         SetMatrixElement(v,i,j,-v_data_ij);
    //     }
    // }
    // 重新配置svd_matrix数据
    svd_matrix.A = aa;
    svd_matrix.U = u;
    svd_matrix.V = v;
    return svd_matrix;
    }



/*********************************************************************************************/
/*********************************************************************************************/


