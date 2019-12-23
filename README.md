"# maLib" 
 Mathematical Analysis Lib for C Ver2.0
#使用说明
##复数函数
1. void Complex_conjugate(int n,struct Complex in[],struct Complex out[]); //求共轭
2. void Complex_abs(struct Complex a,double *abs);//求模
3. void Complex_plus(struct Complex a,struct Complex b,struct Complex *c);//求复数和
4. void Complex_sub(struct Complex a,struct Complex b,struct Complex *c);//求复数差
5. void Complex_multiply(struct Complex a,struct Complex b,struct Complex *c);//求复数积
##矩阵函数
1. void Matrix_plus(int **Matrix_a,int **Matrix_b,int **Matrix_sum,int row,int column);//求矩阵和
2. void Matrix_multiply(int **Matrix_a,int **Matrix_b,int **Matrix_product,int common,int row,int column);//求矩阵积
3. void Matrix_transpose(int **Matrix_a,int **Matrix_at, int row,int column); //求转置矩阵
##傅里叶变换函数
1. void DFT(struct Complex *x,struct Complex *X,int N);//离散傅里叶变换
2. void FFT(int N,struct Complex f[],struct Complex *out );//快速傅里叶变换
3. void IFFT(int N,struct Complex f[],struct Complex *out);//快速傅里叶逆变换
