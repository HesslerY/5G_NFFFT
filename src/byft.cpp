#include "data_field.hpp"
#include "global.hpp"
#include "byft.hpp"

namespace calcu_field{
    int byft::calcu_pws(){
        std::cout << "this is calcu_pws" << std::endl;

        near_ref.calcu_polar();
        far_ref.calcu_polar();
        near_ref.print_info();
        k = 2*pai*freq*std::sqrt(myu*eps);

        MatrixXd absEnear = near_ref.Exyz.colwise().norm();
        MatrixXd absEfar_ref = far_ref.Exyz.colwise().norm();
        MatrixXd absEfardata = MatrixXd::Zero(1,absEfar_ref.cols());
        std::cout << "absEnear = \n" << absEnear << std::endl;
        std::cout << "absEfar_ref = \n" << absEfar_ref << std::endl;

        A_pws = Mat_XC::Zero(1,absEnear.cols());
        
        for(int i = 0 ; i < A_pws.cols() ; i++){
            for(int j = 0 ; j < absEnear.cols() ; j++){
                double theta = near_ref.Rpolar(1,i);
                double phi = near_ref.Rpolar(2,i);
                double kx = k * sin(theta)*cos(phi);
                double ky = k * sin(theta)*sin(phi);
                A_pws(0,i) += absEnear(0,j) * exp(Complexd(0,-1) * (kx*near_ref.Rxyz(0,j) + ky*near_ref.Rxyz(0,j)) ) ;
            }
        }
        A_pws = A_pws / (A_pws.cols());

        std::cout << A_pws << std::endl;
        for(int i = 0; i < absEfar_ref.cols() ; i++){
            double r = far_ref.Rpolar(0,i);
            double kz = k * cos(near_ref.Rpolar(1,i));
            // absEfardata(0,i) = Complexd(0,-1) * exp(Complexd(0,-1)*k*r) / r * kz * A_pws(0,i);
            absEfardata(0,i) = kz * A_pws;
        }






        return 0;
    }
}