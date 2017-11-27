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
#ifndef __MALIB_H
#define __MALIB_H


struct Complex {
     double real;
	 double imaginary;
	 }; 

// Complex number Calculate
double Complex_abs(struct Complex a);
struct Complex Complex_plus(struct Complex a,struct Complex b);
struct Complex Complex_multiply(struct Complex a,struct Complex b);
// Matrix Calculate
void Matrix_plus(int **Matrix_a,int **Matrix_b,int **Matrix_sum,int row,int column);
void Matrix_multiply(int **Matrix_a,int **Matrix_b,int **Matrix_product,int common,int row,int column);
void Matrix_transpose(int **Matrix_a,int **Matrix_at, int row,int column); 
// Fourier transform
void DFT(struct Complex *x,struct Complex *X,int N);
void FFT(struct Complex *x,struct Complex *X, struct Complex *xe,struct Complex *xo,struct Complex *X1,struct Complex *X2,int N);
	
#endif 
