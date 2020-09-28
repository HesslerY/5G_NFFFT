#include "data_field.hpp"

namespace data_field{
    int field::read_file(std::string filename){
        std::ifstream inf(filename);
        if(!inf){
            ERR("error cannot open the file :"+ filename);
            return 0;
        }
    //antenna info format
        std::string str;
        char delim = '=';
        while(getline(inf, str)){
            if(str == "#data"){
                break;
            }

            // antenna info =int or string
            std::stringstream ss(str);
            std::string element;
            int info_type = 0;
            while(getline(ss,element,delim)){
                if(info_type == 0){
                    if(element == "antenna name"){
                        info_type = 1;
                    }else if(element == "antenna size"){
                        info_type = 2;
                    }else if(element == "measurement surface"){
                        info_type = 3;
                    }else if(element == "freq"){
                        info_type = 4;
                    }else if(element == "sampling point"){
                        info_type = 5;
                    }else if(element == "sampling space"){
                        info_type = 6;
                    }else{
                        ERR("unknown file format : " + element);
                    }
                }else{
                    switch (info_type){
                        case 1:
                            antenna_name = element;
                            break;
                        case 2:
                            antenna_size = std::stod(element);
                            break;
                        case 3:
                            surface_type = 0; //surface type is not used yet
                            break;
                        case 4:
                            freq = std::stod(element);
                            break;
                        case 5:
                            n_sample = std::stoi(element);
                            break;
                        case 6:
                            space = std::stod(element);
                            break;
                        default:
                            ERR("reading file format is not correct (case default)\n");
                            break;
                    }
                }
            }
        }

        getline(inf,str);
    //data csv format (x,y,z,realEx,realEy,realEz,imgEx,imgEy,imgEz)
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
                    // Exyz(index-6,column) += 0;
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
        if(Rpolar.size() != 0 || Epolar.size() != 0 ){
            std::cout << "polar info is already calculated" << std::endl;
            return 0;
        }else{
            std::cout << "#calcu Epolar and Rpolar from xyz" << std::endl;
        }

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
        if(Exyz != Matrix<double,3,Dynamic>::Zero(3,n_sample)){
            std::cout << "Exyz info is already calculated" << std::endl;
            return 0;
        }else{
            std::cout << "#calcu Exyz from Epolar" << std::endl;
        }

        for(int i = 0 ; i < n_sample ; i++){
            double r_r = Rpolar(0,i);
            double r_theata = Rpolar(1,i);
            double r_phai = Rpolar(2,i);

            Exyz(0,i) = Epolar(0,i)*std::sin(r_theata)*std::cos(r_phai) + Epolar(1,i)*std::cos(r_theata)*std::cos(r_phai) - Epolar(2,i)*std::sin(r_phai);
            Exyz(1,i) = Epolar(0,i)*std::sin(r_theata)*std::sin(r_phai) + Epolar(1,i)*std::cos(r_theata)*std::sin(r_phai) + Epolar(2,i)*std::cos(r_phai);
            Exyz(2,i) = Epolar(0,i)*std::cos(r_theata) - Epolar(1,i)*std::sin(r_theata);
        }
        return 0;
    }

    int field::copy_data(const field& in_field){
        Rxyz = in_field.Rxyz;
        n_sample = in_field.n_sample;
        freq = in_field.freq;
        Exyz = Matrix<double,3,Dynamic>::Zero(3,n_sample);
        return 0;
    }

    int field::print_info(){
        std::cout << "========== file info ===========" <<std::endl;
        std::cout << "antenna name = " << antenna_name << std::endl;
        std::cout << "antenna size = " << antenna_size << std::endl;
        std::cout << "surface type = " << surface_type << std::endl;
        std::cout << "n_sample = " << n_sample << std::endl;
        std::cout << "space = " << space << std::endl;
        std::cout << "freq = " << freq << std::endl;

        std::cout << "data size = " << Rxyz.cols() <<std::endl;
        std::cout << "Exyz = " << Exyz.transpose() << std::endl;
        std::cout << "========== end info ===========" << std::endl;
    }
}
