#ifndef _CALCU_FIELD_HPP_INCLUDED
#define _CALCU_FIELD_HPP_INCLUDED

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
        int clacu_fardata3(data_field::field& field_calcu,const data_field::field& ref);
        int calcu_fardata_U(data_field::field& field_calcu,const data_field::field& ref);
        int calcu_error(const data_field::field& field_calcu,const data_field::field& ref,std::string title);
        int calcu_ansbyEigen();
        int print_info();
        int plot_field(const MatrixXd& val_x,const MatrixXd& val_y , std::string title = "graph title");//plot x y and graph title by gnuplot


        Complexd calcu_T(Matrix<double,3,1> , Matrix<double,3,1>);
        Complexd provepattern_theata(Matrix<double,3,1> , int);
        Complexd provepattern_phai(Matrix<double,3,1> , int);
    };
}

#endif // _CIRCUIT_H_INCLUDED