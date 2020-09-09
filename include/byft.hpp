#ifndef _BYFT_HPP_INCLUDED
#define _BYFT_HPP_INCLUDED

#include "data_field.hpp"

namespace calcu_field{
    class byft{
        public:
        data_field::field near_ref;
        data_field::field fardata;
        data_field::field far_ref;

        private:
        static constexpr double pai = 3.141592653589793;
        static constexpr double myu = 1.256637062121e-6; //myu_0 [N A^(-2)]
        static constexpr double eps = 8.854187812813e-12; // eps_0 [F m^(-1)]
        static constexpr double freq = 27e9;// 27Ghz;

        double k;

        Mat_XC A_pws;

        public:
        int calcu_pws();
    };
}

#endif