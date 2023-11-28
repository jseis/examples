/*
C++ implementation of a radix-2 Cooley-Tukey Fast Fourier Transform
and Inverse Fast Fourier Transform
Author: Jonah Eisen
*/

#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>

enum fft_direction_t {fwd = 1, inverse = -1};
const double PI = std::acos(-1);
const std::complex<double> i (0,1);

std::complex<double>* fft(const std::complex<double> * , const int);
std::complex<double>* ifft(const std::complex<double> * , const int);
std::complex<double>* fftHelper(const std::complex<double> *, const int, const int, const fft_direction_t);
void printFft(const std::string, const std::complex<double> *, const int);

int main()
{
    int n = 4;
    std::complex<double> * f = new std::complex<double>[n]{{1,0}, {1,0}, {1,0}, {1,0}};
    std::cout << "f:" << '\n';
    printFft("f", f, n);
    std::complex<double> * fhat = fft(f, n);
    std::cout << "fft(f):" << '\n';
    printFft("fhat", fhat, n);
    f = ifft(fhat, n);
    std::cout << "ifft(fft(f)):" << '\n';
    printFft("f", f, n);
    delete f;
    delete fhat;
    return 0;

}

std::complex<double>* fft(const std::complex<double> * f, const int n) {
    std::complex<double>* fhat = fftHelper(f, n, 1, fft_direction_t::fwd);
    //rescale by factor of 1/sqrt(n)
    //this definition ensures that f and fhat have the same l^2 norm
    for(int j = 0; j < n; j++){
        fhat[j] /= std::sqrt(n);
    }
    return fhat;
}

std::complex<double>* ifft(const std::complex<double> * fhat, const int n) {
    std::complex<double>* f = fftHelper(fhat, n, 1, fft_direction_t::inverse);
    //rescale by factor of 1/sqrt(n)
    //this definition ensures that f and fhat have the same l^2 norm
    for(int j = 0; j < n; j++){
        f[j] /= std::sqrt(n);
    }
    return f;
}

std::complex<double>* fftHelper(const std::complex<double> * f, const int r, const int stride, const fft_direction_t dir) {
    //the i'th level of recursive calls to fftHelper has 2^i calls
    //for 0 <= i <= logn, stride = 2^i and r = n/stride
    //assign the j_i'th index of the calls at level i according to j_i = \sum_{k = 0}^{i-1}c_k*2^k,
    //where c_k = 1 if the i'th FFT was taken on odd indices, 0 otherwise
    //then, the j_i'th call takes the FFT of the elements in the set
    //J_i = {x| x = j_i + 2is, 0<= x <= N-1}
    if(r == 1) {
        std::complex<double>* fhat = new std::complex<double>[1];
        fhat[0] = *f;
        return fhat;
    } else {
        std::complex<double> * fhat = new std::complex<double>[r];
        //take fft of even indices
        std::complex<double> * evenhat = fftHelper(f, r/2, 2*stride, dir);
        //take fft of odd indices
        std::complex<double> * oddhat = fftHelper(f + stride, r/2, 2*stride, dir);
        //use odd and even fft's to compute full fft
        std::complex<double> W(1,0);
        double theta_n = dir*(-2*PI*(1/((double)r)));
        std::complex<double> Wn = std::exp(theta_n*i);
        for(int j = 0; j < r/2; j++) {
            fhat[j] = evenhat[j] + W*oddhat[j];
            fhat[j+r/2] = evenhat[j] - W*oddhat[j];
            W *= Wn;
        }
        delete evenhat;
        delete oddhat;
        return fhat;
    }
}

void printFft(const std::string name, const std::complex<double> * func, const int n) {
    for(int j = 0; j < n; j++){
        std::cout << name << "[" << j << "]="  << func[j] << '\n';
    }
}

