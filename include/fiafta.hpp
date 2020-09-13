#ifndef _FIAFTA_HPP_INCLUDED
#define _FIAFTA_HPP_INCLUDED

#include "global.hpp"
#include "data_field.hpp"


namespace calcu_field{
    class fiafta{
    public:
        fiafta(){
            std::cout << "this is calcu field constructa"  << std::endl;
        }
        
    public:
        static constexpr double pai = 3.141592653589793;
        static constexpr double myu = 1.256637062121e-6; //myu_0 [N A^(-2)]
        static constexpr double eps = 8.854187812813e-12; // eps_0 [F m^(-1)]
        static constexpr double freq = 27e9;// 27Ghz;
        static constexpr double accur_SVD = 300;// singular value less than this is set to be 0 

    private:
        int L;
        int P_theata;
        int P_phai;
        int P;
        double w_theata;
        double w_phai;
        double k_0;
        Complexd coeff_A;


        // std::vector<double> weight_phai;
        // std::vector<double> weight_theata;
        std::vector<double> weight;
        std::vector<Matrix<double,3,1>> vec_k;
        std::vector<double> vec_sin; // vec of sin theata
        std::vector<Matrix<double,2,1>> vec_k_angle; // vec[i][0] = theata_i , vec[i][1] = phai_i
    
    public:
        Mat_XC A;// b = Ax
        Mat_XC ans;
        Mat_XC U_mea; //prove voltage at neardata points

    public:
        data_field::field near_ref;
        data_field::field fardata;
        data_field::field far_ref;

        int start_calcu(int val_L = 0);
        int set_matrix(Mat_XC& mat_cup,const Matrix<double,3,Dynamic>& Rxyz); //set matrix A and other value (P,k_0 ...etc)
        int calcu_ansbySVD();
        int calcu_fardata(data_field::field& field_calcu,const data_field::field& ref);
        int calcu_fardata2(data_field::field& field_calcu,const data_field::field& ref);
        int calcu_fardata3(data_field::field& field_calcu,const data_field::field& ref);
        int calcu_fardata_U(data_field::field& field_calcu,const data_field::field& ref);
        int calcu_error(const data_field::field& field_calcu,const data_field::field& ref,std::string title);
        int calcu_ansbyEigen();
        int print_info();
        int savetxt_csv(Mat_XC data, std::string filename, bool cflag);

        Complexd calcu_T(Matrix<double,3,1> , Matrix<double,3,1>);
        Complexd provepattern_theata(Matrix<double,3,1> , int);
        Complexd provepattern_phai(Matrix<double,3,1> , int);
    };
}

#endif // _CIRCUIT_H_INCLUDED