#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <stdio.h>
#include <stdlib.h>

// #include <Eigen/Dense>
#include "data_field.hpp"
#include "global.hpp"

int main(int argc, char** argv){

    // std::cout << "argc=" << argc << std::endl;
    // if(argc < 2){
    //     std::cout << "error (./field_trans nearfield.txt farfield.txt)" <<std::endl;
    // }

    std::string file_near = "../data/near_15_32.txt";
    std::string file_far_ref = "../data/far_50_16.txt";
    calcu_field::calcu calcu1;

    calcu1.neardata.read_file(file_near);
    calcu1.far_ref.read_file(file_far_ref);
    calcu1.start_calcu();
    calcu1.set_matrix();
    // calcu1.calcu_ansbySVD();
    calcu1.calcu_ansbyEigen();
    calcu1.calcu_fardata(calcu1.fardata,calcu1.far_ref);
    // calcu1.calcu_fardata2(calcu1.fardata,calcu1.far_ref);

    calcu1.print_info();
    calcu1.calcu_error(calcu1.fardata,calcu1.far_ref);

    std::cout << "finish all calculation" << std::endl;

    // ルジャンドル*ハンケル
    // double krm = 0.5;
    // double kdotr = 0.5;
    // Complexd sum = 0;
    // for(int i = 0 ; i < 50 ; i++){
    //     Complexd result =  (double)(2*i+1)*sph_hankel_2(i,krm)*std::legendre(i,kdotr)*std::pow(Complexd(0,-1),i);
    //     sum += result;
    //     std::cout << "(" << i << ") = " << result ;
    //     std::cout << " :(sum) = " << sum <<std::endl;
    //     // std::cout << i <<"hankel_2 = " << sph_hankel_2(i,krm) <<std::endl;;
    // }

    // MatrixXd mat(2,2);
    // mat << 1 ,3 ,4 ,5;
    // std::cout << mat * 3 << std::endl;
    // return 0;
}