#include "data_field.hpp"
#include "global.h"


namespace data_field{
    int field::read_file(std::string filename){
        std::ifstream inf(filename);
        if(!inf){
            std::cout << "error cannot open the file" << std::endl;
        }

    //line number 1,2,3
        std::string str;
        char delim = '=';
        int line = 0;
        while(line < 3){
            getline(inf, str);
            std::stringstream ss(str);
            std::string element;
            getline(ss,element,delim);
            getline(ss,element,delim);
            if(line == 0){
                z = std::stod(element);
            }else if(line == 1){
                n_sample = std::stoi(element);
            }
            line++;
        }

    //line number >3 (x,y,z,realEx,realEy,realEz,imgEx,imgEy,imgEz)
        Rxyz = Matrix<double,3,Dynamic>::Zero(3,n_sample*n_sample);
        Exyz = Matrix<double,3,Dynamic>::Zero(3,n_sample*n_sample);
        delim = ',';
        int column = 0;
        while(getline(inf, str)){
            std::stringstream ss(str);
            std::string element;
            int index = 0;

            while(getline(ss,element,delim)){
                if(index < 3){
                    Rxyz(index,column) = std::stod(element);
                }else if(index < 6){
                    Exyz(index-3,column) += std::complex<double> (std::stod(element),0);
                }else if(index < 9){
                    Exyz(index-6,column) += std::complex<double> (0,std::stod(element));
                }

                if(index >= 9){
                    ERR("error input file has more column than 9");
                }
                index++;
            }
            column++;
        }

        std::cout << Exyz.transpose() << std::endl;
        // std::cout << Rxyz.transpose() << std::endl;
        std::cout << "finish read_file" << std::endl;
        return 0;
    }

}