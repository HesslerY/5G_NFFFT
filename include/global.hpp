#ifndef _GLOBAL_H_INCLUDED
#define _GLOBAL_H_INCLUDED

#include <string>
#include <vector>
#include <set>
#include <map>
#include <unordered_set>
#include <unordered_map>
#include <complex>
#include <assert.h>
#include <cmath>
#include <math.h>

#include <iostream>
#include <fstream>
#include <unistd.h>

typedef std::complex<double> Complexd;

inline void ERR(std::string err_msg) {
  std::cerr << err_msg << std::endl;
}

inline void MSG(std::string err_msg) {
  std::cout << err_msg << std::endl;
}

// spherical hankel function
inline std::complex<double> sph_hankel_2(unsigned n, double x){
    std::complex<double> result = std::sph_bessel(n,x) - std::complex<double> (0,1) * std::sph_neumann(n,x);
    return result;
}

inline double d_check(double val, double check_less = 1e-10){
    if(abs(val) < check_less){
        return 0;
    }else{
        return val;
    }
}

inline Complexd c_check(Complexd val, double check_less = 1e-10){
    if(abs(val.real()) < check_less){
        val.real(0);
    }
    if(abs(val.imag()) < check_less){
        val.imag(0);
    }
    return val;
}

#endif // _GLOBAL_H_INCLUDED

