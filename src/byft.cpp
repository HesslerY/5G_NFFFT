#include "byft.hpp"

namespace calcu_field{
    int byft::calcu_pws(std::string title){
        std::cout << "this is calcu_pws" << std::endl;
        near_ref.calcu_polar();
        far_ref.calcu_polar();
        fardata.copy_data(far_ref);
        fardata.calcu_polar();
        // near_ref.print_info();
        k = 2*pai*freq*std::sqrt(myu*eps);

        A_pws = Mat_XC::Zero(3,near_ref.Exyz.cols());
        
        for(int i = 0 ; i < near_ref.Exyz.cols() ; i++){
            for(int j = 0 ; j < near_ref.Exyz.cols() ; j++){
                double theta = near_ref.Rpolar(1,i);
                double phi = near_ref.Rpolar(2,i);
                double kx = k * sin(theta)*cos(phi);
                double ky = k * sin(theta)*sin(phi);
                A_pws(0,i) += near_ref.Exyz(0,j) * exp(Complexd(0,-1) * (kx*near_ref.Rxyz(0,j) + ky*near_ref.Rxyz(1,j)) );
                A_pws(1,i) += near_ref.Exyz(1,j) * exp(Complexd(0,-1) * (kx*near_ref.Rxyz(0,j) + ky*near_ref.Rxyz(1,j)) );
                A_pws(2,i) += near_ref.Exyz(2,j) * exp(Complexd(0,-1) * (kx*near_ref.Rxyz(0,j) + ky*near_ref.Rxyz(1,j)) );
            }
        }
        A_pws = A_pws / near_ref.Exyz.cols() / (2*pai);
        std::cout << "A_pws = \n" << A_pws << std::endl;

        for(int i = 0; i < fardata.Exyz.cols() ; i++){
            double kz = k * cos(fardata.Rpolar(1,i));
            double r = fardata.Rpolar(0,i);
            fardata.Exyz(0,i) = Complexd(0,1) * exp(Complexd(0,-1) * k * r) / r * kz * A_pws(0,i);
            fardata.Exyz(1,i) = Complexd(0,1) * exp(Complexd(0,-1) * k * r) / r * kz * A_pws(1,i);
            fardata.Exyz(2,i) = Complexd(0,1) * exp(Complexd(0,-1) * k * r) / r * kz * A_pws(2,i);
        }
        fardata.Exyz = fardata.Exyz / (2*pai);

        std::cout << "fardata Exyz = \n" << fardata.Exyz << std::endl;
        std::cout << "far_ref Exyz = \n" << far_ref.Exyz << std::endl;

        Matrix<double,1,Dynamic> val_x = Matrix<double,1,Dynamic>::Zero(1,far_ref.n_sample);
        Matrix<double,2,Dynamic> mat_power = Matrix<double,2,Dynamic>::Zero(2,far_ref.n_sample);
        for(int i = 0 ; i < far_ref.n_sample ; i++){
            val_x(i) = i;
        }

        mat_power.row(0) = far_ref.Exyz.colwise().norm();
        mat_power.row(1) = fardata.Exyz.colwise().norm();
        Matrix<double,3,Dynamic> power_db = Matrix<double,3,Dynamic>::Zero(3,far_ref.n_sample);
        double max_power = mat_power.row(0).maxCoeff();

        // power_db.topRows(2) = 20*(mat_power.array() / max_power).log10();
        power_db.row(0) = 20*(mat_power.row(0).array() / mat_power.row(0).maxCoeff()).log10();
        power_db.row(1) = 20*(mat_power.row(1).array() / mat_power.row(1).maxCoeff()).log10();

        // power_db.row(2) = 20*(abs(mat_power.row(0).array() - mat_power.row(1).array())/max_power).log10();
        power_db.row(2) = 20*(abs(mat_power.row(0).array() - mat_power.row(1).array())*mat_power.row(0).array().inverse()).log10();
        // power_db.row(2) = 20*(error.colwise().norm().array()*mat_power.row(0).array().inverse()).log10();

        std::vector<std::string> key_info{"amp_{ref}","amp_{cal}","amp_{error}"};
        std::replace(title.begin(),title.end(),'_','-');
        title = "NFFFT by FT (" + title + ")";
        std::vector<std::string> graph_info{title,"sample points","relative mag[dB]",""};

        plot_field_global(val_x,power_db,key_info,graph_info);
        
        return 0;
    }
}