#ifndef _DATA_FIELD_HPP_INCLUDED
#define _DATA_FIELD_HPP_INCLUDED

#include "global.hpp"

#include <Eigen/Dense>
#include <cmath>

using namespace Eigen;

typedef Matrix<std::complex<double>,Dynamic,Dynamic> Mat_XC;


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
        Matrix<double,3,Dynamic> Rpolar; //(r,theta,phi)
        Matrix<Complexd,3,Dynamic> Epolar; //(Er,Etheta,Ephi)

        int read_file(std::string);
        int calcu_polar();
        int calcu_cart();
        int copy_data(const field& in_field);
        int print_info();
    };

    int make_graph_xcut(const data_field::field fardata,const data_field::field far_ref,std::string title);
    int make_graph_ycut(const data_field::field fardata,const data_field::field far_ref,std::string title);
}

long CountNumbersofTextLines(std::string filepath);

#endif
