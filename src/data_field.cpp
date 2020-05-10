#include "data_field.hpp"
#include "global.h"


namespace data_field{
    int field::read_file(std::string filename){
        std::ifstream inf(filename);
        if(!inf){
            std::cout << "error cannot open the file" << std::endl;
            return 0;
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

        // std::cout << Exyz.transpose() << std::endl;
        // std::cout << Rxyz.transpose() << std::endl;
        std::cout << "finish read_file" << std::endl;
        return 1;
    }
}

namespace calcu_field{
    int calcu::calcu_error(data_field::field field_calcu, data_field::field ref){
        std::cout << "this is calcu_error()" << std::endl;
        return 0;
    }

    int calcu::set_matrix(){
        // for(int i = 0 ; i < 5 ; i++){
        //     std::cout << std::legendre(i,0.1) << std::endl;
        // }

        std::cout << "neardata.n_sample = " << neardata.n_sample << std::endl;
        std::vector<Matrix<double,3,1>> vec_k;
        double delta_theata = pai/L+1;
        double delta_phai = 2*pai/(2*L + 1);
        for(int n_theata = 0 ; n_theata < L ; n_theata++){
            for(int n_phai = 0 ; n_phai < 2*L + 1 ; n_phai++){
                Matrix<double,3,1> k;
                double x = std::sin(delta_theata*(n_theata+1)) * std::cos(delta_phai * n_phai);
                double y = std::sin(delta_theata*(n_theata+1)) * std::sin(delta_phai * n_phai);
                double z = std::cos(delta_theata*(n_theata+1));
                k << x,y,z;
                vec_k.push_back(k);
                std::cout << "!!!!!!!!!! norm = " << k.norm() << "!!!!!!!!!!!!!" <<std::endl; 
                std::cout << k << std::endl;
            }
        }
        Matrix<double,3,1> k;
        k << 3,0,4;
        Matrix<double,3,1> o;
        o << 2,1,3;        
        k.normalize();
        std::cout << k.array()*o.array() << std::endl;
        

        A.resize(neardata.n_sample*neardata.n_sample , 2*L*(2*L+1));
        std::cout << "A size = " << A.rows() << " * " << A.cols() << std::endl;
        for(int j = 0 ; j < A.cols() ; j++){
            for(int i = 0 ; i < A.rows() ; i++){
                A(i,j) = i+j;
            }
        }
        // std::cout << A << std::endl;
        std::cout << std::sin(pai/2) << std::endl;

        //wavenumber k
        // Matrix<std::complex<double>,3,1> k;
        // k << 1,0,0;

        // calcu_T(k,k);

        return 0;
    }

    std::complex<double> calcu::calcu_T(Matrix<double,3,1> k, Matrix<double,3,1> r_R){
        std::complex<double> result = 0;
        std::cout << k.dot(r_R) << std::endl;
        std::cout << "k.normalized =\n" <<k.normalized() << std::endl;
        for(int i = 0 ; i < L+1 ; i++){
            result += std::pow(std::complex<double> (0,-1),i) * (double)(2*i+1) * std::legendre(i,k.dot(r_R));
        }
        std::cout << result << std::endl;

        return result;
    }
}