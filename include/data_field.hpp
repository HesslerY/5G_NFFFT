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
        std::string antenna_name;
        double antenna_size;
        int surface_type; //surface type is not used yet
        int n_sample;
        double space;
        double freq;
        Matrix<double,3,Dynamic> Rxyz;
        Matrix<Complexd,3,Dynamic> Exyz;
        Matrix<double,3,Dynamic> Rpolar;
        Matrix<Complexd,3,Dynamic> Epolar;

        int read_file(std::string);
        int calcu_polar();
        int calcu_cart();
        int copy_data(const field& in_field);
        int print_info();
    };
}

#endif
