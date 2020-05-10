#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <stdio.h>
#include <stdlib.h>

#include <Eigen/Dense>
#include "data_field.hpp"
#include "global.h"

int main(int argc, char** argv){
    std::cout << "koike made this file\n";

    std::cout << "argc=" << argc << std::endl;


    std::string filename_in = "../data/result.txt";
    data_field::field field;
    calcu_field::calcu koike;

    // Enear.read_file(filename_in);
    koike.neardata.read_file(filename_in);
    // koike.calcu_error(field,field);
    koike.set_matrix();


    
// sample matrix 
//     Matrix<std::complex<double> ,2,2> m1;
//     m1(0,0) = std::complex<double> (3,1);
//     m1(1,0) = std::complex<double> (2,1);
//     m1(0,1) = std::complex<double> (-1,0);
//     m1(1,1) = std::complex<double> (4,0);
//     std::cout << m1 << std::endl;
//     Matrix<std::complex<double> ,2,2> m2;
//     m2(0,0) = std::complex<double> (0,1);
//     m2(1,0) = std::complex<double> (0,1);
//     m2(0,1) = std::complex<double> (-1,0);
//     m2(1,1) = std::complex<double> (-1,0);
//     std::cout << "m1*m2 = \n" << m1*m2 << std::endl;

// // calcu singular value (M = U * singular matrix * V.adjoint)
//     JacobiSVD<Matrix< std::complex<double>,2,2 >> SVD(m1,ComputeFullU | ComputeFullV);
//     // std::cout << "this is singular value = \n" << SVD.singularValues().asDiagonal() << std::endl; 

//     std::cout << "this is U = \n" << SVD.matrixU() << std::endl; 
//     std::cout << "this is V = \n" << SVD.matrixV() << std::endl; 
//     std::cout << SVD.matrixU() * SVD.singularValues().asDiagonal() * SVD.matrixV().adjoint() << std::endl;

    return 0;
}