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
    data_field::field Enear;
    Enear.read_file(filename_in);


    
    // sample matrix 
    // MatrixXd m(2,2);
    // m(0,0) = 3;
    // m(1,0) = 2.5;
    // m(0,1) = -1;
    // m(1,1) = m(1,0) + m(0,1);
    // std::cout << m << std::endl;

    return 0;
}