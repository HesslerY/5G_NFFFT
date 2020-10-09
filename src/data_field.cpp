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

    int make_graph_xcut(data_field::field far_ref,data_field::field fardata,std::string title){
        // std::cout << "this is make grapah function" << std::endl;
        int n_point = (int)std::sqrt(fardata.n_sample);
        // std::cout << "n_point = " << n_point << std::endl;

        MatrixXd alldata_db = Matrix<double,2,Dynamic>::Zero(2,fardata.n_sample);
        double maxE = far_ref.Exyz.colwise().norm().maxCoeff();
        alldata_db.row(0) = (far_ref.Exyz.colwise().norm() / maxE).array().log10()*20;
        // alldata_db.row(1) = (fardata.Exyz.colwise().norm() / maxE).array().log10()*20;
        alldata_db.row(1) = (fardata.Exyz.colwise().norm() / fardata.Exyz.colwise().norm().maxCoeff()).array().log10()*20;
        
        MatrixXd data = Matrix<double,2,Dynamic>::Zero(2,n_point);
        Matrix<double,1,Dynamic> val_x = Matrix<double,1,Dynamic>::Zero(1,n_point);

        for(int i = 0 ; i < n_point ; i++){
            data(0,i) = alldata_db(0,n_point*i + n_point/2 - 1);
            data(1,i) = alldata_db(1,n_point*i + n_point/2 - 1);
            val_x(0,i) = fardata.Rxyz(0,n_point*i + n_point/2 - 1); // plot by x
        }

        std::vector<std::string> key_info{"amp_{ref}","amp_{cal}","amp_{error}","phase_{ref}","phase_{cal}"};
        std::replace(title.begin(),title.end(),'_','-');
        title += " (x cut)";
        std::vector<std::string> graph_info{title,"x[m]","relative mag[dB]","phase[{/Symbol \260}]"};
        plot_field_global(val_x,data,key_info,graph_info);
        return 0;
    }

    int make_graph_ycut(const data_field::field far_ref,const data_field::field fardata,std::string title){
        // std::cout << "this is make grapah function" << std::endl;
        int n_point = (int)std::sqrt(fardata.n_sample);
        // std::cout << "n_point = " << n_point << std::endl;

        MatrixXd alldata_db = Matrix<double,2,Dynamic>::Zero(2,fardata.n_sample);
        double maxE = far_ref.Exyz.colwise().norm().maxCoeff();
        alldata_db.row(0) = (far_ref.Exyz.colwise().norm() / maxE).array().log10()*20;
        // alldata_db.row(1) = (fardata.Exyz.colwise().norm() / maxE).array().log10()*20;
        alldata_db.row(1) = (fardata.Exyz.colwise().norm() / fardata.Exyz.colwise().norm().maxCoeff()).array().log10()*20;
        
        MatrixXd data = Matrix<double,2,Dynamic>::Zero(2,n_point);
        Matrix<double,1,Dynamic> val_x = Matrix<double,1,Dynamic>::Zero(1,n_point);

        for(int i = 0 ; i < n_point ; i++){
            data(0,i) = alldata_db(0,n_point*(n_point/2 - 1) + i);
            data(1,i) = alldata_db(1,n_point*(n_point/2 - 1) + i);
            val_x(0,i) = fardata.Rxyz(1,n_point*(n_point/2 - 1) + i); // plot by y
        }

        std::vector<std::string> key_info{"amp_{ref}","amp_{cal}","amp_{error}","phase_{ref}","phase_{cal}"};
        std::replace(title.begin(),title.end(),'_','-');
        title += " (ycut)";
        std::vector<std::string> graph_info{title,"y[m]","relative mag[dB]","phase[{/Symbol \260}]"};
        plot_field_global(val_x,data,key_info,graph_info);
        return 0;
    }

}
