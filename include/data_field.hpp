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

        int read_file(std::string);
    };
}

namespace calcu_field{
    class calcu{
    public:    
        calcu(){
            std::cout << "this is calcu field constructa"  << std::endl;
        }
        
    public:
        data_field::field neardata;
        data_field::field fardata;
        data_field::field far_ref;

    public:
        static constexpr double pai = 3.141592653589793;
        static constexpr double myu = 1.25663706212e-10; //myu_0 [N A^(-2)]
        static constexpr double eps = 8.8541878128e-12; // eps_0 [F m^(-1)]
        static constexpr double freq = 28e9;// 28Ghz;
        static constexpr double accur_SVD = 1e2;// singular value less than this is set to be 0 

    private:
        int L = 3;
        int P_theata;
        int P_phai;
        int P;
        double k_0;
        Mat_XC A;// b = Ax
        Mat_XC ans;
        Mat_XC U_mea; //prove voltage at neardata points

    public:
        int calcu_error(data_field::field field_calcu, data_field::field ref);
        int set_matrix(); //set matrix A and other value (P,k_0 ...etc)
        int calcu_ansbySVD();
        int calcu_fardata();
        int print_info();
        Complexd calcu_T(Matrix<double,3,1> , Matrix<double,3,1>);
    };
}

#endif // _CIRCUIT_H_INCLUDED
