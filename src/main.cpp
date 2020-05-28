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
#include <cmath>

int main(int argc, char** argv){

    // std::cout << "argc=" << argc << std::endl;
    // if(argc < 2){
    //     std::cout << "error (./field_trans nearfield.txt farfield.txt)" <<std::endl;
    // }

    std::string file_near = "../data/near_15_32_short.txt";
    std::string file_far_ref = "../data/far_50_32_short.txt";
    calcu_field::calcu calcu1;

    calcu1.neardata.read_file(file_near);
    calcu1.far_ref.read_file(file_far_ref);
    calcu1.start_calcu();
    calcu1.set_matrix(calcu1.A,calcu1.neardata.Rxyz);
    // calcu1.calcu_ansbySVD();
    calcu1.calcu_ansbyEigen();
    calcu1.calcu_fardata(calcu1.fardata,calcu1.far_ref);
    // calcu1.calcu_fardata2(calcu1.fardata,calcu1.far_ref);
    // calcu1.calcu_fardata_U(calcu1.fardata,calcu1.far_ref);

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

    //球面の数値計算 
    // double result = 0;
    // double pai = 3.141592;
    // int L = 15;
    // int p_theata = L;
    // int p_phai = 2*L;
    // double delta_theata = pai/(p_theata);
    // double delta_phai = 2*pai/(p_phai);
    // double w_theata = pai / (p_theata);
    // double w_phai = 2*pai/ (p_phai);
    // for(int i = 1 ; i <= p_theata ; i++){
    //     for(int j = 1 ; j <= p_phai ; j++){
    //         result += std::sin(delta_theata * i - delta_theata/2) * w_theata * w_phai;
    //         std::cout << std::sin(delta_theata * i - delta_theata/2) <<std::endl;
    //     }
    // }
    // std::cout << result << std::endl;

    //  check T_L
    // Complexd result = 0;
    // int i = 10;
    // result = std::pow(Complexd(0,-1),i) * (double)(2*i+1) * std::legendre(i,0.5) * sph_hankel_2(i,2.2);
    // std::cout << result << std::endl;
    return 0;
}
