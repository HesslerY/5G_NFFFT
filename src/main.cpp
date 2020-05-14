#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <stdio.h>
#include <stdlib.h>

#include <Eigen/Dense>
#include "data_field.hpp"
#include "global.hpp"

int main(int argc, char** argv){

    // std::cout << "argc=" << argc << std::endl;

    std::string file_near = "../data/neardata.txt";
    std::string file_far_ref = "../data/fardata.txt";
    data_field::field field;
    calcu_field::calcu calcu1;

    calcu1.neardata.read_file(file_near);
    calcu1.far_ref.read_file(file_far_ref);
    calcu1.set_matrix();
    calcu1.calcu_ansbySVD();
    calcu1.print_info();

    std::cout << "finish all calculation" << std::endl;


    
// sample matrix 
    // Matrix<std::complex<double> ,2,2> m1;
    // m1(0,0) = std::complex<double> (3,1);
    // m1(1,0) = std::complex<double> (2,1);
    // m1(0,1) = std::complex<double> (-1,0);
    // m1(1,1) = std::complex<double> (4,0);

    // Matrix<std::complex<double>,2,2> m2;
    // m2 = m1.array().inverse();
    // std::cout << m2 <<std::endl;
    // std::cout << m1 << std::endl;
    // Matrix<std::complex<double> ,2,2> m2;
    // m2(0,0) = std::complex<double> (0,1);
    // m2(1,0) = std::complex<double> (0,1);
    // m2(0,1) = std::complex<double> (-1,0);
    // m2(1,1) = std::complex<double> (-1,0);
    // std::cout << "m1*m2 = \n" << m1*m2 << std::endl;

// calcu singular value (M = U * singular matrix * V.adjoint)
    // JacobiSVD<Matrix< std::complex<double>,2,2 >> SVD(m1,ComputeFullU | ComputeFullV);
    // // std::cout << "this is singular value = \n" << SVD.singularValues().asDiagonal() << std::endl; 

    // std::cout << "this is U = \n" << SVD.matrixU() << std::endl; 
    // std::cout << "this is V = \n" << SVD.matrixV() << std::endl; 
    // std::cout << SVD.matrixU() * SVD.singularValues().asDiagonal() * SVD.matrixV().adjoint() << std::endl;

// 球面での数値計算
    // const double pai = 3.1415926535;
    // int L = 1000;
    // std::vector<Matrix<double,3,1>> vec_k;
    // std::vector<double> vec_sin_theata;
    // double delta_theata = pai/(L+1);
    // double delta_phai = 2*pai/(2*L + 1);
    // double w_theate = pai/(L+1);
    // double w_phai = 2*pai/(2*L+1);
    // for(int n_theata = 0 ; n_theata < L ; n_theata++){
    //     std::cout << n_theata << ": " << "(" << delta_theata*(n_theata+1) << ")" ;
    //     std::cout << std::sin(delta_theata*(n_theata+1)) << std::endl;
    //     for(int n_phai = 0 ; n_phai < 2*L + 1 ; n_phai++){
    //         Matrix<double,3,1> k;
    //         double x = std::sin(delta_theata*(n_theata+1)) * std::cos(delta_phai * n_phai);
    //         double y = std::sin(delta_theata*(n_theata+1)) * std::sin(delta_phai * n_phai);
    //         double z = std::cos(delta_theata*(n_theata+1));
    //         k << x,y,z;
    //         vec_k.push_back(k);
    //         vec_sin_theata.push_back(std::sin(delta_theata*(n_theata+1)));

    //         // std::cout << "!!!!!!!!!! norm = " << k.norm() << "!!!!!!!!!!!!!" <<std::endl; 
    //         // std::cout << k << std::endl;
    //     }
    // }
    // double result = 0;
    // for(int i = 0 ; i < vec_k.size() ; i++){
    //     result += vec_sin_theata[i]*w_theate*w_phai;
    // }
    // std::cout << "result = " << result << std::endl;
    // std::cout << "ideal result = " << 4*pai << std::endl;

    return 0;
}