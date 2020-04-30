#ifndef _CIRCUIT_H_INCLUDED
#define _CIRCUIT_H_INCLUDED

#include "global.h"
#include <Eigen/Dense>

using namespace Eigen;

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
        Matrix<double,3,Dynamic> Rxyz;
        Matrix<std::complex<double>,3,Dynamic> Exyz;

        int read_file(std::string);

        // // find a node by name, returns NULL if not found
        // Node *GetNode(std::string name) {
        // std::map<std::string, Node *>::iterator it = all_nodes_map.find(name);
        // if (it != all_nodes_map.end())
        //     return it->second;
        // return NULL;
        // }
    };
}

#endif // _CIRCUIT_H_INCLUDED
