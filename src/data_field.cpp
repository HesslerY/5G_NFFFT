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
        std::cout << Exyz.transpose() << std::endl;
        // std::cout << Rxyz.transpose() << std::endl;
        std::cout << "finish read_file" << std::endl;
        return 1;
    }
}

namespace calcu_field{
    int calcu::calcu_error(data_field::field field_calcu, data_field::field ref){
        std::cout << "this is calcu_error()" << std::endl;
        return 0;
    }

    int calcu::set_matrix(){
        k_0 = 2*pai*freq*std::sqrt(myu*eps);
        P_theata = L;
        P_phai = 2*L + 1;
        P = P_theata * P_phai;
        std::vector<Matrix<double,3,1>> vec_k;
        std::vector<double> vec_sin; // vec of sin theata
    //recheck 分割 L or L+1 and n_theata + 1 or not 
        double delta_theata = pai/P_theata;
        double delta_phai = 2*pai/P_phai;

        double w_theata = pai/P_theata;
        double w_phai = 2*pai/P_phai;
        for(int n_theata = 0 ; n_theata < L+1 ; n_theata++){
            for(int n_phai = 0 ; n_phai < 2*L+1 ; n_phai++){
                Matrix<double,3,1> k;
                double x = std::sin(delta_theata*(n_theata+1)) * std::cos(delta_phai * n_phai);
                double y = std::sin(delta_theata*(n_theata+1)) * std::sin(delta_phai * n_phai);
                double z = std::cos(delta_theata*(n_theata+1));
                k << x,y,z;
                vec_k.push_back(k);
                vec_sin.push_back(std::sin(delta_theata*(n_theata+1)));
                // std::cout << "!!!!!!!!!! norm = " << k.norm() << "!!!!!!!!!!!!!" <<std::endl; 
                // std::cout << k << std::endl;
            }
        }

        A.resize(neardata.n_sample , 2*P);
        // A.resize(neardata.n_sample , P);
        // std::cout << "A size = " << A.rows() << " * " << A.cols() << std::endl;
        for(int j = 0 ; j < A.cols() ; j++){
            if(j < P){
                for(int i = 0 ; i < A.rows() ; i++){
                    std::complex<double> temp = calcu_T(vec_k[j],neardata.Rxyz.col(i));
                    A(i,j) = w_theata*w_phai*temp*vec_sin[j];
                }
            }else if( j >= P){
                for(int i = 0 ; i < A.rows() ; i++){
                    std::complex<double> temp = calcu_T(vec_k[j-P],neardata.Rxyz.col(i));
                    A(i,j) = w_theata*w_phai*temp*vec_sin[j-P];
                }
            }
        }
        // std::cout << A.transpose() << std::endl;

        return 0;
    }
    
    int calcu::calcu_ansbySVD(){
        std::cout << "========== SVD info ==============" << std::endl;
        std::cout << "this is calcu_ansbySVD" <<std::endl;

    // calcu singular value (M = U * singular matrix * V.adjoint)
        JacobiSVD<Mat_XC> SVD(A,ComputeThinU | ComputeThinV);

        // std::cout << "sigular value = \n" << SVD.singularValues() << std::endl;
        // std::cout << "sigular value = \n" << (MatrixXd) SVD.singularValues().asDiagonal() << std::endl;
        std::cout << "SV row = " << SVD.singularValues().rows() << std::endl;
        std::cout << "SV col = " << SVD.singularValues().cols() << std::endl;
        // std::cout << "this is U = \n" << SVD.matrixU() << std::endl; 
        std::cout << "U.row = " << SVD.matrixU().rows() << std::endl; 
        std::cout << "U.col = " << SVD.matrixU().cols() << std::endl; 
        // std::cout << "this is V = \n" << SVD.matrixV() << std::endl; 
        std::cout << "V.row = " << SVD.matrixV().rows() << std::endl; 
        std::cout << "V.col = " << SVD.matrixV().cols() << std::endl; 
        
        // std::cout << "error = \n" << SVD.matrixU() * (MatrixXd)SVD.singularValues().asDiagonal() * SVD.matrixV().adjoint() - A << std::endl;

        // MatrixXd singular_diag = (MatrixXd)SVD.singularValues().asDiagonal().array();

    // inverse singular value more than accur_SVD is set to be 0
        MatrixXd inv_sing_trun = SVD.singularValues().array().inverse().matrix();
        std::cout << "sing = \n" << inv_sing_trun << std::endl;
        for(int i = 0 ; i < inv_sing_trun.size() ; i++){
            if(inv_sing_trun(i) >= accur_SVD){
                inv_sing_trun(i) = 0;
            }
        }
        std::cout << "inv_sing_trun = \n" << inv_sing_trun << std::endl;

    // set diagonal matrix
        MatrixXd diag_inv_sing_trun = inv_sing_trun.asDiagonal();
        // std::cout << "diag_inv_sing_trun = \n" << diag_inv_sing_trun << std::endl;

    // set measured U vector ()
        U_mea = MatrixXd::Ones(neardata.n_sample,1);
        for(int i = 0 ; i < neardata.n_sample ; i++){
            U_mea(i) = neardata.Exyz.col(i).sum();
        }
        std::cout << "U_mea = \n" << U_mea.transpose() << std::endl;

    // calcu ans matrix
        ans = SVD.matrixV() * diag_inv_sing_trun * SVD.matrixU().adjoint() * U_mea;
        // std::cout << "ans = \n" << ans << std::endl;
        // std::cout << "ans.row = " << ans.rows() << std::endl;
        // std::cout << "ans.col = " << ans.cols() << std::endl;
        std::cout << "========== SVD info end ==============" << std::endl;

        return 0;
    }

// args are not normalized
    Complexd calcu::calcu_T(Matrix<double,3,1> k, Matrix<double,3,1> r_R){
        std::complex<double> result = 0;
        // std::cout << "k.normalized =\n" <<k.normalized() << std::endl;
        // std::cout << "r_R.norm = " << r_R.norm() << std::endl;
        for(int i = 0 ; i < L+1 ; i++){
            // std::cout << "!!!!!!!!!!!!" << i << "!!!!!!!!!!!" << std::endl;
            // std::cout << "k = \n" << k << std::endl;
            // std::cout << "k.norm = " << k.norm() <<std::endl;
            // std::cout << "r_R.normalized = \n" << r_R.normalized() << std::endl;
            // std::cout << "k.norm = " << r_R.normalized().norm() <<std::endl;
            // std::cout << "dot = " << k.dot(r_R.normalized()) << std::endl;
            result += std::pow(Complexd(0,-1),i) * (double)(2*i+1) * std::legendre(i,k.normalized().dot(r_R.normalized())) * sph_hankel_2(i,k_0*r_R.norm());
        }
        return result;
    }

// calcu fardata from ansbySVD
    int calcu::calcu_fardata(){
        return 0;
    }

// print info 
    int calcu::print_info(){
        std::cout << "========= info calcu=========" << std::endl;
        std::cout << "sampling points = " << neardata.n_sample << std::endl;
        std::cout << "freq = " << freq << std::endl;
        std::cout << "k_0 = " << k_0 << std::endl;
        std::cout << "L = " << L << std::endl;
        std::cout << "P_theata = " << P_theata << std::endl;
        std::cout << "P_phai = " << P_phai << std::endl;
        
        std::cout << "A.row = " << A.rows() << std::endl;
        std::cout << "A.col = " << A.cols() << std::endl;
        // std::cout << "A = " << A << std::endl;

        std::cout << "========= end =========" << std::endl;

        std::cout << "========= info ans =========" << std::endl;
        std::cout << "ans.row = " << ans.rows() << std::endl;
        std::cout << "ans.col = " << ans.cols() << std::endl;
        std::cout << "ans = \n" << ans << std::endl;
        std::cout << "========= end info ans=========" << std::endl;

        return 0;

    }
}