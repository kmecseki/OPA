#pragma once
#ifndef UTILS_H
#define UTILS_H

#include <fftw3.h>


double deg2rad (double);
double rad2deg (double);
void cvector_to_fftw(int, std::vector<std::complex<double>>, fftw_complex*);
void fftw_to_cvector(int, fftw_complex*, std::vector<std::complex<double>>&);
void fftshift(std::vector<std::complex<double>>&, int);
void fftshift(fftw_complex *, int);
double FindMax(const std::vector<double> &);
double get_FWHM(const std::vector<std::complex<double>> &, std::vector<double>);
void writeToFile(const char *ofname, std::vector<double> &data1, std::vector<std::complex<double>> &data2);
void writeToFile(const char *ofname, std::vector<double> &data1, std::vector<double> &data2);

#endif // UTILS_H