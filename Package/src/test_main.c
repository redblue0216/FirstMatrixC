// 载入依赖头文件
#include "first_matrix_c.h"



/*************************************/
/**           程序入口              **/
/*************************************/
void main(void)
{
    /*
    **矩阵存储模块测试
    */
   // 创建矩阵
    printf("===============================>>>>>> This is a CreateMatrix function test!\n");
    Matrix a = CreateMatrix(3,5);
    Matrix b = CreateMatrix(3,5);
    float array[] = {1,2,3,5,2,
                     4,5,6,4,3,
                     7,8,9,7,7};                    
    SetMatrixData(a,array);
    SetMatrixData(b,array);
    // 获取矩阵元素
    float element = GetMatrixElement(a,3,2);
    printf("The element of the matrix is %f\n",element);
    // SetMatrixZero(matrix);
    // float element_0 = GetMatrixElement(matrix,3,3);
    // printf("The element of the matrix is %f\n",element_0);
    int m;
    for (m=0;m<a.row*a.column;m++){
        printf("The matrix a element is %f\n",a.data[m]);
        printf("The matrix b element is %f\n",b.data[m]);
    }
    /*
    **矩阵基本运算测试
    */
    //矩阵加法
    printf("==============================>>>>>> This is a Matrixadd function test!\n");
    Matrix c = MatrixAdd(a,b);
    printf("The matrix row is %d\n",c.row);
    printf("The matrix column is %d\n",c.column);
    int i;
    for (i=0;i<a.row*a.column;i++){
        printf("The matrix c element is %f\n",c.data[i]);    
    }
    // 矩阵标量乘法
    printf("==============================>>>>>> This is a MatrixScalaMultiply function test!\n");
    int tmp_scalar = 10;
    Matrix bb = MatrixScalarMultiply(tmp_scalar,a);
    int j;
    for (j=0;j<a.row*a.column;j++){
        printf("The matrix bb element is %f\n",bb.data[j]);    
    }
    // 矩阵乘法
    printf("==============================>>>>>> This is a MatrixMultiply function test!\n");
    Matrix e = CreateMatrix(3,3);
    Matrix f = CreateMatrix(3,2);
    float array_e[] = {1,2,3,
                       4,5,6,
                       7,8,9};
    float array_f[] = {1,4,
                       2,5,
                       3,6};                                         
    SetMatrixData(e,array_e);
    SetMatrixData(f,array_f);
    Matrix cc = MatrixMultiply(e,f);
    printf("The matrix cc's row is %d\n",cc.row);
    printf("The matrix cc's column is %d\n",cc.column);
    int k;
    for (k=0;k<cc.row*cc.column;k++){
        printf("The matrix cc element is %f\n",cc.data[k]);    
    }
    // 矩阵转置
    printf("==============================>>>>>> This is a MatrixTranspose function test!\n");
    Matrix ft = MatrixTranspose(f);
    int n;
    printf("The matrix ft's row is %d\n",ft.row);
    printf("The matrix ft's column is %d\n",ft.column);
    for (n=0;n<ft.row*ft.column;n++){
        printf("The matrix ft element is %f\n",ft.data[n]);    
    }    
    // 矩阵行列式计算
    printf("==============================>>>>>> This is a AlgebraicCofactor function test!\n");
    Matrix g1 = CreateMatrix(1,1);
    float array_g1[] = {2};                                        
    SetMatrixData(g1,array_g1);
    Cofactor cofactor_g1;
    cofactor_g1 = AlgebraicCofactor(g1);
    printf("The g1 matrix cofactor is %f\n",cofactor_g1.data[0]);     
    Matrix g2 = CreateMatrix(2,2);
    float array_g2[] = {2,8,
                       7,9};                                        
    SetMatrixData(g2,array_g2);
    Cofactor cofactor_g2;
    cofactor_g2 = AlgebraicCofactor(g2);
    int l;
    for (l=0;l<=3;l++)
    {
        printf("The g2 matrix cofactor is %f\n",cofactor_g2.data[l]);
    }  
    Matrix g = CreateMatrix(4,4);
    float array_g[] = {1,2,3,4,
                       5,6,7,8,
                       9,10,11,12,
                       13,14,15,16};                                                              
    SetMatrixData(g,array_g);
    Cofactor cofactor_g;
    cofactor_g = AlgebraicCofactor(g);
    int h;
    for (h=0;h<g.column*(g.column-1)*(g.row-1);h++)
    {
        printf("&&&&&&,%d\n",h);
        printf("The g matrix cofactor is %f\n",cofactor_g.data[h]);
    }
    printf("The g matrix order is %d\n",cofactor_g.order);
    // 代数余子式转变为矩阵序列
    printf("==============================>>>>>> This is a CofactorToMatrixSeries function test!\n");
    MatrixSeries matrix_series_g;
    matrix_series_g = CofactorToMatrixSeries(cofactor_g);
    int matrix_num;
    printf("======= matrix number is %d\n",matrix_series_g.matrix_number);
    for (matrix_num=0;matrix_num<matrix_series_g.matrix_number;matrix_num++)
    {
        Matrix tmp_matrix;
        tmp_matrix = matrix_series_g.data[matrix_num];
        INDEX matrix_i;
        printf("==========================================NO.%d matrix\n",matrix_num);
        for (matrix_i=0;matrix_i<(matrix_series_g.order*matrix_series_g.order);matrix_i++)
        {
            printf("~~~~~---->> %f\n",tmp_matrix.data[matrix_i]);
        }
    }
    // 矩阵行列式计算
    printf("==============================>>>>>> This is a MatrixDeterminant function test!\n");
    REAL det_g1;
    det_g1 = MatrixDeterminant(g1);
    printf("The g2 matrix determinant is %f\n",det_g1);
    REAL det_g2;
    det_g2 = MatrixDeterminant(g2);
    printf("The g2 matrix determinant is %f\n",det_g2); 
    Matrix g3 = CreateMatrix(4,4);
    // det = 16
    // float array_g3[] = {2,-1,0,0,
    //                     0,3,1,0,
    //                     0,0,-2,4,
    //                     5,1,2,-3};  
    // det = -198        
    // float array_g3[] = {2,-1,1,3,
    //                     1,4,0,0,
    //                     -2,0,2,0,
    //                     5,0,0,-3};   
    // det = -130                         
    float array_g3[] = {5,0,0,-3,
                        0,3,1,0,
                        0,4,-2,0,
                        1,0,0,2};                                                        
    SetMatrixData(g3,array_g3); 
    REAL det_g3;
    det_g3 = MatrixDeterminant(g3);
    printf("The g3 matrix determinant is %f\n",det_g3); 
    Matrix g4 = CreateMatrix(4,3);
    float array_g4[] = {2,-1,0,
                        0,3,1,
                        0,0,-2,
                        5,1,2,};                                        
    SetMatrixData(g4,array_g4); 
    REAL det_g4;
    det_g4 = MatrixDeterminant(g4);
    printf("The g4 matrix determinant is %f\n",det_g4);
    // det = -12
    Matrix g5 = CreateMatrix(5,5);
    float array_g5[] = {3,-7,8,9,-6,
                        0,2,-5,7,3,
                        0,0,1,5,0,
                        0,0,2,4,-1,
                        0,0,0,-2,0};  
    SetMatrixData(g5,array_g5); 
    REAL det_g5;                        
    det_g5 = MatrixDeterminant(g5);
    printf("The g5 matrix determinant is %f-----dim:%d\n",det_g5,g5.column);  
    printf("==============================>>>>>> This is a MatrixCofactor function test!\n");
    Matrix cofactor_matrix_g3;
    cofactor_matrix_g3 = MatrixCofactor(g3,2,2);
    INDEX g3_i;
    // 5,0,-3,0,-2,0,1,0,2
    for (g3_i=0;g3_i<(cofactor_matrix_g3.row*cofactor_matrix_g3.column);g3_i++)
    {
        printf("The g3 cofactor matrix's element is %f\n",cofactor_matrix_g3.data[g3_i]);
    }
    printf("==============================>>>>>> This is a MatrixAdjoint function test!\n");
    AdjointMatrix adjoint_matirx_g5;
    adjoint_matirx_g5 = MatrixAdjoint(g5);
    printf("The g5's determinant is %f\n",adjoint_matirx_g5.matrix_determinant);
    printf("The g5's sign is %f\n",adjoint_matirx_g5.reverse_order_number[23]);
    printf("The g5's data is %f\n",adjoint_matirx_g5.data[23]);
    printf("==============================>>>>>> This is a MatrixAdjoint function test!\n");
    Matrix inverse_matrix_g5;
    inverse_matrix_g5 = MatrixInverse(g5);
    INDEX inverse_g5_i;
    printf("The g5's inverse matrix row is %d\n",inverse_matrix_g5.row);
    printf("The g5's inverse matrix column is %d\n",inverse_matrix_g5.column);
    for (inverse_g5_i=0;inverse_g5_i<g5.column*g5.row;inverse_g5_i++)
    {   
        printf("The g5's inverse matrix element is %f\n",inverse_matrix_g5.data[inverse_g5_i]);
    }
    /*
    **矩阵分解运算测试
    */
    // float e = addthree(2,3);
    // printf("The result is %f\n",e);
}