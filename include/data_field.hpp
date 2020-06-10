#ifndef _CIRCUIT_H_INCLUDED
#define _CIRCUIT_H_INCLUDED

#include "global.hpp"
#include <Eigen/Dense>
#include <cmath>

using namespace Eigen;

typedef Matrix<std::complex<double>,Dynamic,Dynamic> Mat_XC;
typedef std::complex<double> Complexd;

namespace data_field{
    class field{
    public:
        // Circuit() {};
        // virtual ~Circuit() { Clear(); };
        // void Clear() {
        // inputs.clear();
        // outputs.clear();
        // ffs.clear();
        // all_nodes_map.clear();
        // for (long i = 0; i < all_nodes.size(); i++) {
        //     delete all_nodes[i];
        // }
        // all_nodes.clear();
        // }
        // Circuit *GetDuplicate(std::string input_prefix, std::string output_prefix, std::string internal_prefix);

    public:
        double z;
        int n_sample;
        double freq;
        Matrix<double,3,Dynamic> Rxyz;
        Matrix<Complexd,3,Dynamic> Exyz;
        Matrix<double,3,Dynamic> Rpolar;
        Matrix<Complexd,3,Dynamic> Epolar;

        int read_file(std::string);
        int calcu_polar();
        int calcu_cart();
        int copy_data(const field& in_field);
    };
}

namespace calcu_field{
    class calcu{
    public:    
        calcu(){
            std::cout << "this is calcu field constructa"  << std::endl;
        }
        
    public:
        static constexpr double pai = 3.141592653589793;
        static constexpr double myu = 1.25663706212e-6; //myu_0 [N A^(-2)]
        static constexpr double eps = 8.8541878128e-12; // eps_0 [F m^(-1)]
        static constexpr double freq = 27e9;// 27Ghz;
        static constexpr double accur_SVD = 1;// singular value less than this is set to be 0 

    private:
        int L;
        int P_theata;
        int P_phai;
        int P;
        double w_theata;
        double w_phai;
        double k_0;
        Complexd coeff_A;

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

        int start_calcu();
        int set_matrix(Mat_XC& mat_cup,const Matrix<double,3,Dynamic>& Rxyz); //set matrix A and other value (P,k_0 ...etc)
        int calcu_ansbySVD();
        int calcu_fardata(data_field::field& field_calcu,const data_field::field& ref);
        int calcu_fardata2(data_field::field& field_calcu,const data_field::field& ref);
        int calcu_fardata_U(data_field::field& field_calcu,const data_field::field& ref);
        int calcu_error(const data_field::field& field_calcu,const  data_field::field& ref);
        int calcu_ansbyEigen();
        int print_info();
        int plot_field();//plot E field by gnuplot


        Complexd calcu_T(Matrix<double,3,1> , Matrix<double,3,1>);
        Complexd provepattern_theata(Matrix<double,3,1> , int);
        Complexd provepattern_phai(Matrix<double,3,1> , int);
    };
}

#endif // _CIRCUIT_H_INCLUDED
