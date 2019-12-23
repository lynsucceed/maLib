 /************************************************************************
    # Mathematical Analysis Lib	  by v2.0
    # Author : lynsucceed
    # function : Linear Algebra (Matrix Calculate) 
                 Fourier transform
			     Complex number Calculate 
    # https://github.com/lynsucceed
    # E-mail:lynsucceed@gmail.com
    # Copyright (c) 2017 lynsucceed
 **************************************************************************
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*****************************************************************************/
#include "malib.h"
#include "math.h"
#define pi 3.1415926535


void Complex_conjugate(int n,struct Complex in[],struct Complex out[])  
{  
  int i = 0;  
  for(i=0;i<n;i++)  
  {  
    out[i].real = in[i].real; 
	out[i].imaginary = -in[i].imaginary;  
  }   
}  

void Complex_abs(struct Complex a,double *abs)
{
  *abs = sqrt(pow(a.real,2)+pow(a.imaginary,2));
 }


void Complex_plus(struct Complex a,struct Complex b,struct Complex *c)
{
	c->real = a.real + b.real;
	c->imaginary = a.imaginary + b.imaginary;
}

void Complex_sub(struct Complex a,struct Complex b,struct Complex *c)  
{  
  c->real = a.real - b.real;  
  c->imaginary = a.imaginary - b.imaginary;   
}  
  

void Complex_multiply(struct Complex a,struct Complex b,struct Complex *c)
{
  c->real = a.real*b.real-a.imaginary*b.imaginary;
  c->imaginary = a.imaginary*b.real+a.real*b.imaginary;
 }


void Matrix_plus(int **Matrix_a,int **Matrix_b,int **Matrix_sum,int row,int column )
{
	int i,j;
	for(i=0;i<row;i++)
	{
	  for(j=0;j<column;j++)
	  {
	      //Matrix_sum[i][j]= Matrix_a[i][j]+Matrix_b[i][j];
		  *((int*)Matrix_sum+column*i+j)= *((int*)Matrix_a+column*i+j)+*((int*)Matrix_b+column*i+j);
	  }
	}
}


void Matrix_multiply(int **Matrix_a,int **Matrix_b,int **Matrix_product,int common,int row,int column)
{
   int i,j,k;
   int temp;
	for(i=0;i<row;i++)
	{
	  for(j=0;j<column;j++)
	  {
	    temp = 0;
		for(k=0;k<common;k++)
		{
          temp += *((int*)Matrix_a+common*i+k)*(*((int*)Matrix_b+column*k+j));//temp += Matrix_a[i][k]*Matrix_b[k][j]
	    }
		*((int*)Matrix_product+column*i+j)=temp;//Matrix_product[i][j]=temp
	 }
	}
}
	
void Matrix_transpose(int **Matrix_a,int **Matrix_at, int row,int column)  
{  
    int i, j;     
    for (i = 0; i < row; i++)  
    {  
        for(j = 0; j < column; j++)  
        {  
            *((int*)Matrix_at+column*j+i) =  *((int*)Matrix_a+column*i+j);
			//Matrix_at[j][i] =  Matrix_a[i][j]
        }  
    }  
}  	

void DFT(struct Complex *x,struct Complex *X,int N)
{
   int n,k;
   struct Complex angle;
   angle.real = cos(2*pi*k*n/N);
   angle.imaginary = -sin(2*pi*k*n/N);
   for(k=0;k<N;k++)
   {
    //X[k].real = 0;
    //X[k].imaginary = 0;
	for(n=0;n<N;n++)
	{
	  X[k].real += x[n].real*angle.real-x[n].imaginary*angle.imaginary;
	  X[k].imaginary += x[n].imaginary*angle.real+x[n].real*angle.imaginary;
	  }
   }
}


void Wn_i(int n,int i,struct Complex *Wn,char flag)  
{  
  Wn->real = cos(2*pi*i/n);  
  if(flag == 1)  
  Wn->imaginary = -sin(2*pi*i/n);  
  else if(flag == 0)  
  Wn->imaginary = -sin(2*pi*i/n);  
}  


void FFT(int N,struct Complex f[],struct Complex *out )  
{  
  struct Complex t,wn;//中间变量  
  int i,j,c,k,m,n,l,r,M;  
  int la,lb,lc;  
  /*----计算分解的级数M=log2(N)----*/  
  for(i=N,M=1;(i=i/2)!=1;M++);   
  /*----按照倒位序重新排列原信号----*/  
  for(i=1,j=N/2;i<=N-2;i++)  
  {  
    if(i<j)  
    {  
      t=f[j];  
      f[j]=f[i];  
      f[i]=t;  
    }  
    k=N/2;  
    while(k<=j)  
    {  
      j=j-k;  
      k=k/2;  
    }  
    j=j+k;  
  } 

 /*----把倒位序的原信号赋给输出数组----*/  
  for(c=0;c<N;c++) 
  {
	int num= 0;  
	out[num] = f[num];
  }

  /*----FFT算法----*/  
  for(m=1;m<=M;m++)  
  {  
    la = pow((double)2,m); //la=2^m代表第m级每个分组所含节点数       
    lb=la/2;    //lb代表第m级每个分组所含蝶形单元数  
                 //同时它也表示每个蝶形单元上下节点之间的距离  
    /*----蝶形运算----*/  
    for(l=1;l<=lb;l++)  
    {  
      r=(l-1)*pow((double)2,M-m);     
      for(n=l-1;n<N-1;n=n+la) //遍历每个分组，分组总数为N/la  
      {  
        lc=n+lb;  //n,lc分别代表一个蝶形单元的上、下节点编号       
        Wn_i(N,r,&wn,1);//wn=Wnr  
        Complex_multiply(out[lc],wn,&t);//t = f[lc] * wn复数运算  
        Complex_sub(out[n],t,&(out[lc]));//f[lc] = f[n] - f[lc] * Wnr  
        Complex_plus(out[n],t,&(out[n]));//f[n] = f[n] + f[lc] * Wnr  
      }  
    }  
  }  
} 

void IFFT(int N,struct Complex f[],struct Complex *out)  
{  
  int i=0;  
  Complex_conjugate(N,f,f);  
  FFT(N,f,out);  
  Complex_conjugate(N,out,out);  
  for(i=0;i<N;i++)  
  {  
    out[i].real = (out[i].real)/N; 
	out[i].imaginary = (out[i].imaginary)/N;   
  }  
}  
