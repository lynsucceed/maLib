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
#ifndef __MALIB_H
#define __MALIB_H


struct Complex {
     double real;
	 double imaginary;
	 }; 

// Complex number Calculate
void Complex_conjugate(int n,struct Complex in[],struct Complex out[]);
void Complex_abs(struct Complex a,double *abs);
void Complex_plus(struct Complex a,struct Complex b,struct Complex *c);
void Complex_sub(struct Complex a,struct Complex b,struct Complex *c);  
void Complex_multiply(struct Complex a,struct Complex b,struct Complex *c);
// Matrix Calculate
void Matrix_plus(int **Matrix_a,int **Matrix_b,int **Matrix_sum,int row,int column);
void Matrix_multiply(int **Matrix_a,int **Matrix_b,int **Matrix_product,int common,int row,int column);
void Matrix_transpose(int **Matrix_a,int **Matrix_at, int row,int column); 
// Fourier transform
void DFT(struct Complex *x,struct Complex *X,int N);
void FFT(int N,struct Complex f[],struct Complex *out );   
void IFFT(int N,struct Complex f[],struct Complex *out);

#endif 
