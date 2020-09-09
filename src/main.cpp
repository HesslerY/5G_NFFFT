#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <stdio.h>
#include <stdlib.h>

// #include <Eigen/Dense>
#include "data_field.hpp"
#include "fiafta.hpp"
#include "global.hpp"
#include "byft.hpp"
#include <cmath>

int main(int argc, char** argv){

    std::string file_near_ref;
    std::string file_far_ref;
    
    bool flag_fiafta = false;
    bool flag_fourier = true;

    if(flag_fiafta){
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
        calcu_field::fiafta fiafta1;
        fiafta1.near_ref.read_file(file_near_ref);
        fiafta1.far_ref.read_file(file_far_ref);
        // fiafta1.near_ref.print_info();
        fiafta1.start_calcu(val_L);
        fiafta1.set_matrix(fiafta1.A,fiafta1.near_ref.Rxyz);
        // fiafta1.calcu_ansbySVD();
        fiafta1.calcu_ansbyEigen();
        fiafta1.calcu_fardata(fiafta1.fardata,fiafta1.far_ref);
        // fiafta1.calcu_fardata2(fiafta1.fardata,fiafta1.far_ref);
        // fiafta1.calcu_fardata3(fiafta1.fardata,fiafta1.far_ref);
        // fiafta1.calcu_fardata_U(fiafta1.fardata,fiafta1.far_ref);

        fiafta1.print_info();
        std::string title = file_near_ref + " to " + file_far_ref;
        fiafta1.calcu_error(fiafta1.fardata,fiafta1.far_ref,title);

        // std::cout << "finish all calculation" << std::endl;
        // fiafta1.savetxt_csv(fiafta1.A,"data_A",true);
        // fiafta1.savetxt_csv(fiafta1.U_mea,"data_Umea",true);
        // fiafta1.savetxt_csv(fiafta1.ans,"data_ans",true);
    }
    if(flag_fourier){
        if(argc != 3){
            ERR("error ./field_trans nearfield.txt");
        }
        file_near_ref = argv[1];
        file_far_ref = argv[2];

        calcu_field::byft sample;
        sample.near_ref.read_file(file_near_ref);
        sample.far_ref.read_file(file_far_ref);
        sample.calcu_pws();
    }



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
