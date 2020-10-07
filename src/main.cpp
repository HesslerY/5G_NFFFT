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
    
    bool flag_fiafta = true;
    bool flag_fourier = false;

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
        fiafta1.start_calcu(val_L);
        fiafta1.set_matrix(fiafta1.A,fiafta1.near_ref.Rxyz);
        // fiafta1.calcu_ansbySVD();
        fiafta1.calcu_ansbyEigen();
        // fiafta1.calcu_fardata(fiafta1.fardata,fiafta1.far_ref);
        fiafta1.calcu_fardata2(fiafta1.fardata,fiafta1.far_ref);
        // fiafta1.calcu_fardata3(fiafta1.fardata,fiafta1.far_ref);
        // fiafta1.calcu_fardata_U(fiafta1.fardata,fiafta1.far_ref);

        fiafta1.print_info();
        std::string title = file_near_ref + " to " + file_far_ref;
        fiafta1.calcu_error(fiafta1.fardata,fiafta1.far_ref,title);
        data_field::make_graph_xcut(fiafta1.fardata,fiafta1.far_ref,title);
        data_field::make_graph_ycut(fiafta1.fardata,fiafta1.far_ref,title);

        // std::cout << "finish all calculation" << std::endl;
        // fiafta1.savetxt_csv(fiafta1.A,"data_A",true);
        // fiafta1.savetxt_csv(fiafta1.U_mea,"data_Umea",true);
        // fiafta1.savetxt_csv(fiafta1.ans,"data_ans",true);
    }
    if(flag_fourier){
        if(argc != 3){
            ERR("error ./field_trans nearfield.txt farfield.txt");
        }
        file_near_ref = argv[1];
        file_far_ref = argv[2];

        calcu_field::byft sample;
        sample.near_ref.read_file(file_near_ref);
        sample.far_ref.read_file(file_far_ref);
        std::string title = file_near_ref + " to " + file_far_ref;
        sample.calcu_pws(title);
    }
    return 0;
}
