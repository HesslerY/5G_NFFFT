#include "data_field.hpp"
#include "global.hpp"


namespace data_field{
    int field::read_file(std::string filename){
        std::ifstream inf(filename);
        if(!inf){
            std::cout << "error cannot open the file" << std::endl;
            return 0;
        }
    //line number 1,2,3,4
        std::string str;
        char delim = '=';
        int line = 0;
        while(line < 4){
            getline(inf, str);
            std::stringstream ss(str);
            std::string element;
            getline(ss,element,delim);
            getline(ss,element,delim);
            if(line == 0){
                z = std::stod(element);
            }else if(line == 1){
                freq = std::stod(element);
            }else if(line == 2){
                n_sample = std::stoi(element);
            }
            line++;
        }
    //line number > 4 (x,y,z,realEx,realEy,realEz,imgEx,imgEy,imgEz)
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
        return 0;
    }

    int field::copy_data(const field& in_field){
        Rxyz = in_field.Rxyz;
        n_sample = in_field.n_sample;
        freq = in_field.freq;
        Exyz = Matrix<double,3,Dynamic>::Zero(3,n_sample);
        return 0;
    }
}

namespace calcu_field{
    int calcu::start_calcu(){
        // define constant value
        k_0 = 2*pai*freq*std::sqrt(myu*eps);
        int d_0 = 6;
        double d = 1.2e-2;
        L = (int) k_0*d + 1.8*std::pow(d_0,(double)2/3)*std::pow(k_0*d,(double(1/3)));
        L = 20;
        std::cout << "L = " << L << std::endl;
        P_theata = L;
        // P_phai = 2*L + 1;
        P_phai = 2*L;
        P = P_theata * P_phai;
        w_theata = pai/(P_theata);
        w_phai = 2*pai/(P_phai);

        // set measured U vector ()
        near_ref.calcu_polar();
        // U_mea = MatrixXd::Zero(near_ref.n_sample,1);
        // for(int i = 0 ; i < near_ref.n_sample ; i++){
        //     // U_mea(i) = near_ref.Exyz.col(i)(0);
        //     // U_mea(i) = near_ref.Exyz.col(i).sum();
        //     // U_mea(i) = near_ref.Exyz.col(i).norm();
        //     U_mea(i) = near_ref.Epolar(1,i) + near_ref.Epolar(2,i);
        // }
        U_mea = MatrixXd::Zero(near_ref.n_sample * 2 ,1);
        for(int i = 0 ; i < near_ref.n_sample * 2 ; i++){
            if(i < near_ref.n_sample){
                U_mea(i) = near_ref.Epolar(1,i);
            }else{
                U_mea(i) = near_ref.Epolar(2,i-near_ref.n_sample);
            }
        }

        // make wavenumber vec_k
        double delta_theata = pai/(P_theata);
        double delta_phai = 2*pai/(P_phai);
        for(int n_theata = 1 ; n_theata <= P_theata ; n_theata++){
            for(int n_phai = 1 ; n_phai <= P_phai ; n_phai++){
                Matrix<double,3,1> k;
                Matrix<double,2,1> angle;
                double theata = delta_theata*n_theata - delta_theata / 2;
                double phai = delta_phai * n_phai - delta_phai / 2;
                double x = std::sin(theata) * std::cos(phai);
                double y = std::sin(theata) * std::sin(phai);
                double z = std::cos(theata);
                k << x,y,z;
                angle << theata, phai;
                vec_k.push_back(k);
                vec_sin.push_back(std::sin(theata));
                vec_k_angle.push_back(angle);
            }
        }

        return 0;
    }

// set coupling matrix mat_cup
    int calcu::set_matrix(Mat_XC& mat_cup,const Matrix<double,3,Dynamic>& Rxyz){
        std::cout << "#set matrix mat_cup" << std::endl;
        int n_sample = Rxyz.cols();
        mat_cup.resize(n_sample*2 , 2*P);
        // A.resize(near_ref.n_sample , P);
        // std::cout << "A size = " << A.rows() << " * " << A.cols() << std::endl;
        for(int j = 0 ; j < mat_cup.cols() ; j++){
            if(j < P){//about theata
                for(int i = 0 ; i < mat_cup.rows() ; i++){
                    // first polarization (theta) 
                    if(i < n_sample){
                        std::complex<double> temp_T = calcu_T(vec_k[j],Rxyz.col(i));
                        std::complex<double> temp_W = provepattern_theata(vec_k[j],i);
                        mat_cup(i,j) = w_theata*w_phai*temp_T*temp_W*vec_sin[j];
                    }else if(n_sample <= i && i < 2*n_sample){
                    // second polarization (phi)
                        std::complex<double> temp_T = calcu_T(vec_k[j],Rxyz.col(i - n_sample));
                        std::complex<double> temp_W = provepattern_theata(vec_k[j],i - n_sample);
                        // mat_cup(i,j) = w_theata*w_phai*temp_T*temp_W*vec_sin[j];
                        mat_cup(i,j) = 0;
                    }else{
                        std::cout << "error index of mat_cup is wrong" << std::endl;
                    }
                }

            }else if( j >= P){//about phai
                for(int i = 0 ; i < mat_cup.rows() ; i++){
                    if(i < n_sample){
                        std::complex<double> temp_T = calcu_T(vec_k[j-P],Rxyz.col(i));
                        std::complex<double> temp_W = provepattern_phai(vec_k[j-P],i);
                        // mat_cup(i,j) = w_theata*w_phai*temp_T*temp_W*vec_sin[j-P];
                        mat_cup(i,j) = 0;
                    }else if(n_sample <= i && i < 2*n_sample){
                        std::complex<double> temp_T = calcu_T(vec_k[j-P],Rxyz.col(i - n_sample));
                        std::complex<double> temp_W = provepattern_phai(vec_k[j-P],i - n_sample);
                        mat_cup(i,j) = w_theata*w_phai*temp_T*temp_W*vec_sin[j-P];
                    }else{
                        std::cout << "error index of mat_cup is wrong" << std::endl;
                    }
               }
            }
        }
        coeff_A = Complexd(0,-1) * freq * myu / (double)2;
        // A = A * coeff_A;
        // std::cout << "A = " << A << std::endl;
        return 0;
    }

// args are not normalized
    Complexd calcu::calcu_T(Matrix<double,3,1> k, Matrix<double,3,1> r_R){
        Complexd result = 0;
        Complexd coeff_T = Complexd(0,-1) * k_0 / (4*pai); //coff_T = -j*k_0/(4*pai)
        // std:: cout << "k*r_R = " << k_0*r_R.norm() <<std::endl;
        // std::cout << "k.normalized =\n" <<k.normalized() << std::endl;
        // std::cout << "r_R.norm = " << r_R.norm() << std::endl;
        for(int i = 0 ; i < L+1 ; i++){
            // std::cout << "!!!!!!!!!!!!" << i << "!!!!!!!!!!!" << std::endl;
            // std::cout << "k = \n" << k << std::endl;
            // std::cout << "k.norm = " << k.norm() <<std::endl;
            // std::cout << "r_R.normalized = \n" << r_R.normalized() << std::endl;
            // std::cout << "k.norm = " << r_R.normalized().norm() <<std::endl;
            // std::cout << "dot = " << k.dot(r_R.normalized()) << std::endl;
            // std::cout << "legendre = " << std::legendre(i,k.normalized().dot(r_R.normalized())) <<std::endl;
            // std::cout << "hankel_2 = " << sph_hankel_2(i,k_0*r_R.norm()) <<std::endl;
            result += std::pow(Complexd(0,-1),i) * (double)(2*i+1) * std::legendre(i,k.normalized().dot(r_R.normalized())) * sph_hankel_2(i,k_0*r_R.norm());
        }
        // result = result * coeff_T;
        // std::cout << "result = " << result <<std::endl;

        return result;
    }
    
    int calcu::calcu_ansbySVD(){
        std::cout << "========== SVD info ==============" << std::endl;
        std::cout << "this is calcu_ansbySVD" <<std::endl;

    // calcu singular value (M = U * singular matrix * V.adjoint)
        JacobiSVD<Mat_XC> SVD(A,ComputeThinU | ComputeThinV);

        std::cout << "sigular value = \n" << SVD.singularValues() << std::endl;
        // std::cout << "sigular value = \n" << (MatrixXd) SVD.singularValues().asDiagonal() << std::endl;
        std::cout << "SV row = " << SVD.singularValues().rows() << std::endl;
        std::cout << "SV col = " << SVD.singularValues().cols() << std::endl;
        // std::cout << "this is U = \n" << SVD.matrixU() << std::endl; 
        std::cout << "U.row = " << SVD.matrixU().rows() << std::endl; 
        std::cout << "U.col = " << SVD.matrixU().cols() << std::endl; 
        // std::cout << "this is V = \n" << SVD.matrixV() << std::endl; 
        std::cout << "V.row = " << SVD.matrixV().rows() << std::endl; 
        std::cout << "V.col = " << SVD.matrixV().cols() << std::endl; 

        // std::cout << "A = " << A << std::endl;
        // std::cout << "A by SVD = \n" << SVD.matrixU() * (MatrixXd)SVD.singularValues().asDiagonal() * SVD.matrixV().adjoint()  << std::endl;        
        // std::cout << "error = \n" << SVD.matrixU() * (MatrixXd)SVD.singularValues().asDiagonal() * SVD.matrixV().adjoint() - A << std::endl;

        // MatrixXd singular_diag = (MatrixXd)SVD.singularValues().asDiagonal().array();

    // inverse singular value more than accur_SVD is set to be 0
        MatrixXd inv_sing_trun = SVD.singularValues().array().inverse().matrix();
        std::cout << "inv_sing = \n" << inv_sing_trun << std::endl;
        for(int i = 0 ; i < inv_sing_trun.size() ; i++){
            // if(inv_sing_trun(i) >= accur_SVD){
            //     inv_sing_trun(i) = 0;
            // }
            // if( i > 100){
            //     inv_sing_trun(i) = 0;
            // }
            ;
        }
        std::cout << "inv_sing_trun = \n" << inv_sing_trun << std::endl;

    // set diagonal matrix
        MatrixXd diag_inv_sing_trun = inv_sing_trun.asDiagonal();
        // std::cout << "diag_inv_sing_trun = \n" << diag_inv_sing_trun << std::endl;

    // calcu ans matrix
        ans = SVD.matrixV() * diag_inv_sing_trun * SVD.matrixU().adjoint() * U_mea;
        std::cout << "ans = " << ans.transpose() << std::endl;
        std::cout << "error U_mea - y*ans = " << (U_mea - A*ans).transpose() <<std::endl;
        std::cout << "ans.row = " << ans.rows() << std::endl;
        std::cout << "ans.col = " << ans.cols() << std::endl;
        std::cout << "========== SVD info end ==============" << std::endl;

        return 0;
    }

    int calcu::calcu_ansbyEigen(){
        Mat_XC x;
    // ans by LU
        // FullPivLU<Mat_XC> lu(A);
        // x = lu.solve(U_mea);

    // ans by SVD of Eigen
        BDCSVD<Mat_XC> SVD(A,ComputeThinU | ComputeThinV);
        x = SVD.solve(U_mea);
        ans = x;

    // ans by QR decomposition
        // Mat_XC y;
        // y = A.colPivHouseholderQr().solve(U_mea);
        // y = A.fullPivHouseholderQr().solve(U_mea);
        // std::cout << "y.size() = " << y.size() <<std::endl;
        // std::cout << "y = "<< y.transpose() << std::endl;
        // ans = y;

        return 0;
    }

// calculation of prove pattern weight
// ideal prove always 1 (for all direction)
    Complexd calcu::provepattern_theata(Matrix<double,3,1> k , int index_row){
        return 1;
    }
    Complexd calcu::provepattern_phai(Matrix<double,3,1> k , int index_row){
        return 1;
    }

// calcu fardata from ans (U_mea = A*ans) angle by r_m
    int calcu::calcu_fardata(data_field::field& field_calcu,const data_field::field& ref){
        Mat_XC Epolar(3,ref.n_sample);//E_r , E_theata , E_phai @field_calcu
        field_calcu.copy_data(ref);

        // 各r_iに対して計算 T*Dj
        for(int i = 0 ; i < field_calcu.Rxyz.cols() ; i++){
            Complexd E_theata = 0;
            Complexd E_phai = 0;
            for(int j = 0 ; j < A.cols() ; j++){
                if (j < P){
                    E_theata += w_theata*w_phai*calcu_T(vec_k[j],field_calcu.Rxyz.col(i))*vec_sin[j]*ans(j);
                }else if(j >= P){
                    E_phai += w_theata*w_phai*calcu_T(vec_k[j-P],field_calcu.Rxyz.col(i))*vec_sin[j-P]*ans(j);
                }
            }
            // cos(theata) = z/r : tan(phai) = y/x
            double r_theata = std::acos( field_calcu.Rxyz(2,i)/field_calcu.Rxyz.col(i).norm() );
            double r_phai;
            if(field_calcu.Rxyz(1,i) == 0 && field_calcu.Rxyz(0,i) == 0){
                r_phai = 0;
            }else{
                r_phai = std::atan2(field_calcu.Rxyz(1,i) , field_calcu.Rxyz(0,i));
            }

            field_calcu.Exyz(0,i) = E_theata*std::cos(r_theata)*std::cos(r_phai) - E_phai*std::sin(r_phai);
            field_calcu.Exyz(1,i) = E_theata*std::cos(r_theata)*std::sin(r_phai) + E_phai*std::cos(r_phai);
            field_calcu.Exyz(2,i) = -E_theata*std::sin(r_theata);
            // std::cout << "theata,phai = " << r_theata << " " << r_phai << std::endl;
            // std::cout << "E_theata,E_phai = " << E_theata << " " << E_phai << std::endl;

        }
        // field_calcu.Exyz = field_calcu.Exyz * coeff_A;
        std::cout << "finish calculating fardata" <<std::endl;
        // std::cout << field_calcu.Exyz <<std::endl;
        return 0;
    }

// calcu fardata from ans (U_mea = A*ans) angle by k
    int calcu::calcu_fardata2(data_field::field& field_calcu,const data_field::field& ref){
        Mat_XC Epolar(3,ref.n_sample);//E_r , E_theata , E_phai @field_calcu
        field_calcu.copy_data(ref);

        // 各r_iに対して計算 T*Dj
        for(int i = 0 ; i < field_calcu.Rxyz.cols() ; i++){
            for(int j = 0 ; j < ans.rows() ; j++){
                Complexd E_theata = 0;
                Complexd E_phai = 0;
                if (j < P){
                    E_theata = w_theata*w_phai*calcu_T(vec_k[j],field_calcu.Rxyz.col(i))*vec_sin[j]*ans(j);
                    field_calcu.Exyz(0,i) += E_theata*std::cos(vec_k_angle[j](0))*std::cos(vec_k_angle[j](1)) ;
                    field_calcu.Exyz(1,i) += E_theata*std::cos(vec_k_angle[j](0))*std::sin(vec_k_angle[j](1)) ;
                    field_calcu.Exyz(2,i) += -E_theata*std::sin(vec_k_angle[j](0));
                }else if(j >= P){
                    E_phai = w_theata*w_phai*calcu_T(vec_k[j-P],field_calcu.Rxyz.col(i))*vec_sin[j-P]*ans(j);
                    field_calcu.Exyz(0,i) +=  - E_phai*std::sin(vec_k_angle[j-P](1));
                    field_calcu.Exyz(1,i) +=  E_phai*std::cos(vec_k_angle[j-P](1));
                    field_calcu.Exyz(2,i) += 0;
                }
            }
        }
        // field_calcu.Exyz = field_calcu.Exyz * coeff_A;
        return 0;
    }

    int calcu::calcu_fardata_U(data_field::field& field_calcu,const data_field::field& ref){
        Mat_XC A_far;
        field_calcu.copy_data(ref);
        set_matrix(A_far , field_calcu.Rxyz);
        Mat_XC result = MatrixXd::Zero(ref.n_sample,1);
        result = A_far * ans;
        std::cout << "fardata_U.cols() = " << result.cols() << std::endl;
        std::cout << "fardata_U.rows() = " << result.rows() << std::endl;
        std::cout << "fardata_U = \n" << result.transpose() << std::endl;
        return 0;
    }

    int calcu::calcu_error(const data_field::field& field_calcu,const data_field::field& ref){
        std::cout << "======= calcu error between two field  ========" << std::endl;
        Matrix<Complexd,3,Dynamic> error = Matrix<double,3,Dynamic>::Zero(3,ref.n_sample);

        if(field_calcu.Rxyz != ref.Rxyz || field_calcu.Exyz.size() != ref.Exyz.size() ){
            std::cout << "error(calcu_error) field between two fardata should the same" << std::endl;
            std::cout << "======= calcu error between two field  end ========" << std::endl;
            return -1;
        }else{
            error = field_calcu.Exyz - ref.Exyz;
            double sum_error = 0;
            double sum_ref = 0;
            for(int i = 0 ; i < ref.Exyz.cols() ; i++){
                sum_error += error.col(i).norm();
                sum_ref += ref.Exyz.col(i).norm();
            }

            std::cout << "two fields Rxyz and Exyz.size() match" <<std::endl;
            std::cout << "far_ref Exyz = \n" << ref.Exyz << std::endl;
            std::cout << "fardata Exyz = \n" << field_calcu.Exyz << std::endl;
            std::cout << "error (fardata - far_ref) = \n" << error << std::endl; 
            std::cout << "sum_error = " << sum_error << std::endl;
            std::cout << "sum_ref = " << sum_ref << std::endl;
            std::cout << "error (%) = " << sum_error / sum_ref * 100 << std::endl;

            Matrix<double,2,Dynamic> mat_power = Matrix<double,2,Dynamic>::Zero(2,ref.n_sample);
            for(int i = 0 ; i < ref.n_sample ; i++){
                mat_power(0,i) = ref.Exyz.col(i).norm();
                mat_power(1,i) = field_calcu.Exyz.col(i).norm();
            }
            std::cout << "ref and calcu power = \n" << mat_power.transpose() <<std::endl;
            std::cout << "======= calcu error between two field  end ========" << std::endl;
        }
        return 0;
    }

// print info 
    int calcu::print_info(){
        std::cout << "========= info calcu=========" << std::endl;
        std::cout << "sampling points = " << near_ref.n_sample << std::endl;
        std::cout << "freq = " << freq << std::endl;
        std::cout << "k_0 = " << k_0 << std::endl;
        std::cout << "L = " << L << std::endl;
        std::cout << "P_theata = " << P_theata << std::endl;
        std::cout << "P_phai = " << P_phai << std::endl;
        
        // std::cout << "A.row = " << A.rows() << std::endl;
        // std::cout << "A.col = " << A.cols() << std::endl;
        // std::cout << "A.col(0) = " << A.col(0).transpose() <<std::endl;
        // std::cout << "A.col(465) = " << A.col(465).transpose() <<std::endl;
        // std::cout << "A = " << A << std::endl;

        std::cout << "========= end =========" << std::endl;

        std::cout << "========= info ans =========" << std::endl;
        std::cout << "ans.row = " << ans.rows() << std::endl;
        std::cout << "ans.col = " << ans.cols() << std::endl;
        std::cout << "ans = " << ans.transpose() << std::endl;
        std::cout << "U_mea = " << U_mea.transpose() << std::endl;
        std::cout << " A*ans = " << (A*ans).transpose() << std::endl;
        std::cout << "error U_mea - A*ans = " << (U_mea - A*ans).transpose() <<std::endl;
        std::cout << "========= end info ans =========" << std::endl;

        std::cout << "========= info debug =========" << std::endl;
        far_ref.calcu_polar();
        // std::cout << "far_Rxyz = \n" << far_ref.Rxyz.col(527) << std::endl;
        // std::cout << "far_Exyz = \n" << far_ref.Exyz.col(527) << std::endl;
        // std::cout << "far_Rpolar = \n" << far_ref.Rpolar.col(527) << std::endl;
        // std::cout << "far_Epolar = \n" << far_ref.Epolar.col(527) << std::endl;
        std::cout << "far_Rpolar = \n" << far_ref.Rpolar << std::endl;
        std::cout << "far_Epolar = \n" << far_ref.Epolar << std::endl;
        near_ref.calcu_polar();
        std::cout << "near_ref_Rpolar = \n" << near_ref.Rpolar << std::endl;
        std::cout << "near_ref_Epolar = \n" << near_ref.Epolar << std::endl;
        fardata.calcu_polar();
        std::cout << "fardata_Rpolar = \n" << fardata.Rpolar << std::endl;
        std::cout << "fardata_Epolar = \n" << fardata.Epolar << std::endl;
        std::cout << "========= end info debug =========" << std::endl;
        return 0;
    }

    int calcu::plot_field(){
        return 0;
    }
}