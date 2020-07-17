#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <stdio.h>
#include <stdlib.h>

// #include <Eigen/Dense>
#include "data_field.hpp"
#include "calcu_field.hpp"
#include "global.hpp"
#include <cmath>

int main(int argc, char** argv){

    std::string file_near_ref;
    std::string file_far_ref;
    int val_L = 0;

    if(argc == 3){
        file_near_ref = argv[1];
        file_far_ref = argv[2];
    }else if(argc == 4){
        val_L = std::stoi(argv[1]);
        file_near_ref = argv[2];
        file_far_ref = argv[3];
    }else{
        ERR("error (./field_trans (int)L nearfield.txt farfield.txt)");
    }

    std::cout << "val_L = " << val_L << std::endl;


    std::cout.precision(10);


    calcu_field::calcu calcu1;
    calcu1.near_ref.read_file(file_near_ref);
    calcu1.far_ref.read_file(file_far_ref);
    // calcu1.near_ref.print_info();
    calcu1.start_calcu(val_L);
    calcu1.set_matrix(calcu1.A,calcu1.near_ref.Rxyz);
    // calcu1.calcu_ansbySVD();
    calcu1.calcu_ansbyEigen();
    calcu1.calcu_fardata(calcu1.fardata,calcu1.far_ref);
    // calcu1.calcu_fardata2(calcu1.fardata,calcu1.far_ref);
    // calcu1.calcu_fardata3(calcu1.fardata,calcu1.far_ref);
    // calcu1.calcu_fardata_U(calcu1.fardata,calcu1.far_ref);

    calcu1.print_info();
    std::string title = file_near_ref + " to " + file_far_ref;
    calcu1.calcu_error(calcu1.fardata,calcu1.far_ref,title);
    // calcu1.calcu_phase(calcu1.fardata,calcu1.far_ref);

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
    // // std::cout << mat.norm() << std::endl;
    // std::cout << mat.rowwise().norm() << std::endl;
    // std::cout << mat.array().log10() <<std::endl;
    // std::cout << mat.topRows(2) << std::endl;

    // Complexd temp(-2,-1e-10);
    // std::cout << c_check(temp) <<std::endl;

    return 0;
}
