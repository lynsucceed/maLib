 /************************************************************************
    # Mathematical Analysis Lib	  by v1.0
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
#define pi 3.14159

//struct Complexnum {
//     int real;
//	 int imaginary;
//	 }; 







double Complex_abs(struct Complex a)
{
  double abs;
  abs = sqrt(pow(a.real,2)+pow(a.imaginary,2));
  return abs;
 }

struct Complex Complex_plus(struct Complex a,struct Complex b)
{
    struct Complex sum;
	sum.real = a.real + b.real;
	sum.imaginary = a.imaginary + b.imaginary;
	return sum;
}

struct Complex Complex_multiply(struct Complex a,struct Complex b)
{
  struct Complex product;
  product.real = a.real*b.real-a.imaginary*b.imaginary;
  product.imaginary = a.imaginary*b.real+a.real*b.imaginary;
  return product;
 }


void Matrix_plus(int **Matrix_a,int **Matrix_b,int **Matrix_sum,int row,int column )
{
	int i,j;
	for(i=0;i<row;i++)
	{
	  for(j=0;j<column;j++)
	  {
	    Matrix_sum[i][j]= Matrix_a[i][j]+Matrix_b[i][j];
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
		  temp += Matrix_a[i][k]+Matrix_b[k][j];
	    }
		Matrix_product[i][j]=temp;
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
            Matrix_at[j][i] =  Matrix_a[i][j];
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

void FFT(struct Complex *x,struct Complex *X, struct Complex *xe,struct Complex *xo,struct Complex *X1,struct Complex *X2,int N)
{
   int r,k;
   struct Complex value;
   value.real = cos(2*pi*k/N);
   value.imaginary = -sin(2*pi*k/N);
   for(r = 0;r<N;r+=2)
   {
	   int i;
	   xe[i].real = x[r].real;
	   xe[i].imaginary = x[r].imaginary;
	   i++;
   }
    for(r = 1;r<N;r+=2)
   {
	   int j;
	   xo[j].real = x[r].real;
	   xo[j].imaginary = x[r].imaginary;
	   j++;
   }
   
    DFT(xe,X1,N/2);
   	DFT(xo,X2,N/2);
   
   	for(k=0;k<N;k++)
	{
		X[k].real =	X1[k].real +(X2[k].real*value.real-X2[k].imaginary*value.imaginary); 
		X[k].imaginary = X1[k].imaginary + (X2[k].imaginary*value.real+X2[k].real*value.imaginary);
	}
}
