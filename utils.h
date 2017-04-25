#pragma once
#ifndef UTILS_H
#define UTILS_H

#include <fftw3.h>


double deg2rad (double);
double rad2deg (double);
void cvector_to_fftw(int, std::vector<std::complex<double>>, fftw_complex*);
void fftw_to_cvector(int, fftw_complex*, std::vector<std::complex<double>>&);
int fftshift(std::vector<std::complex<double>>&, int);
int fftshift(fftw_complex *, int);

#endif // UTILS_H