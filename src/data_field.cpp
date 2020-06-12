#include "data_field.hpp"
#include "calcu_field.hpp"
#include "global.hpp"


namespace data_field{
    int field::read_file(std::string filename){
        std::ifstream inf(filename);
        if(!inf){
            std::cout << "error cannot open the file" << std::endl;
            return 0;
        }
    //line number 1,2,3,4
        std::string str;
        char delim = '=';
        int line = 0;
        while(line < 4){
            getline(inf, str);
            std::stringstream ss(str);
            std::string element;
            getline(ss,element,delim);
            getline(ss,element,delim);
            if(line == 0){
                z = std::stod(element);
            }else if(line == 1){
                freq = std::stod(element);
            }else if(line == 2){
                n_sample = std::stoi(element);
            }
            line++;
        }
    //line number > 4 (x,y,z,realEx,realEy,realEz,imgEx,imgEy,imgEz)
        Rxyz = Matrix<double,3,Dynamic>::Zero(3,n_sample);
        Exyz = Matrix<double,3,Dynamic>::Zero(3,n_sample);
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
        std::cout << "#finish read file = " << filename <<std::endl;

        return 1;
    }

    int field::calcu_polar(){
        Rpolar = Matrix<double,3,Dynamic>::Zero(3,n_sample);
        Epolar = Matrix<double,3,Dynamic>::Zero(3,n_sample);
        for(int i = 0 ; i < n_sample ; i++){
            double r_r = Rxyz.col(i).norm();
            double r_theata = std::acos( Rxyz(2,i)/Rxyz.col(i).norm() );
            double r_phai;
            if(Rxyz(1,i) == 0 && Rxyz(0,i) == 0){
                r_phai = 0;
            }else{
                r_phai = std::atan2(Rxyz(1,i) , Rxyz(0,i));
            }
            Rpolar(0,i) = r_r;
            Rpolar(1,i) = r_theata;
            Rpolar(2,i) = r_phai;
            Epolar(0,i) = Exyz(0,i)*std::sin(r_theata)*std::cos(r_phai) + Exyz(1,i)*std::sin(r_theata)*std::sin(r_phai) + Exyz(2,i)*std::cos(r_theata);
            Epolar(1,i) = Exyz(0,i)*std::cos(r_theata)*std::cos(r_phai) + Exyz(1,i)*std::cos(r_theata)*std::sin(r_phai) - Exyz(2,i)*std::sin(r_theata);
            Epolar(2,i) = -Exyz(0,i)*std::sin(r_phai) + Exyz(1,i)*std::cos(r_phai);
        }
        return 0;
    }

    int field::calcu_cart(){
        return 0;
    }

    int field::copy_data(const field& in_field){
        Rxyz = in_field.Rxyz;
        n_sample = in_field.n_sample;
        freq = in_field.freq;
        Exyz = Matrix<double,3,Dynamic>::Zero(3,n_sample);
        return 0;
    }
}
