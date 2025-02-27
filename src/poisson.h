#include <iostream>
#include <fstream>
#include <complex>
#include <algorithm>
#include <string>
#include <time.h>
#include <iomanip>

#include <fftw3.h>
#include <mdb.h>

using namespace std;

void meshgrid_wavenumber(int N_x, int N_y, double x_domain, double y_domain, double **k_square);

void subtract_const(int N_x, int N_y, double **vec_c);

void poisson_solver(double **vec_rhs, double x_domain, double y_domain, int N_x, int N_y, double **u_fft);

void divide_fftw_vec(int N_x, int N_y, fftw_complex *fftwIn, double **vec_x);

