#include "global.hpp"

int plot_field_global(const MatrixXd& val_x, const MatrixXd& val_y, std::vector<std::string> key_info, std::vector<std::string> graph_info){
    // std::cout << val_y.transpose() << std::endl;
    if(val_x.cols() != val_y.cols()){
        ERR("error : plot data should have the same cols");
    }

    if(graph_info.size() == 0){
        graph_info.resize(4);
        graph_info[0] = "sample title";
        graph_info[1] = "xlabel";
        graph_info[2] = "y1label";
        graph_info[3] = "y2label";
    }
    if(key_info.size() == 0 ){
        key_info.resize(val_y.size());
        key_info[0] = "ref";
        for(int i = 1 ; i < key_info.size() ; i++){
            key_info[i] = "data" + std::to_string(i);
        }
    }

//  make datafile for gnuplot
    FILE* dataf;
    const char* filename = "datafile";
    dataf = fopen(filename,"w");
    fprintf(dataf,"#x");
    for(int i = 0; i < key_info.size(); i++){
        fprintf(dataf," %s",key_info[i].c_str());
    }
    fprintf(dataf,"\n");
    for(int i = 0 ; i < val_y.cols() ; i++){
        fprintf(dataf,"%f ",val_x(i));
        for(int j = 0 ; j < val_y.rows() ; j++){
            fprintf(dataf,"%f ",val_y(j,i));
        }
        fprintf(dataf,"\n");
    }
    fclose(dataf);

// gnuplot make graph
    FILE* gp;
    gp = popen("gnuplot -persist","w");
    fprintf(gp,"set grid \n");
    fprintf(gp,"set title \"%s\"\n", graph_info[0].c_str());
    fprintf(gp,"set xlabel \"%s\"\n",graph_info[1].c_str());
    fprintf(gp,"set ylabel \"%s\"\n",graph_info[2].c_str());
    fprintf(gp,"set y2label \"%s\"\n",graph_info[3].c_str());
    // fprintf(gp,"%s\n",cmd.str().c_str());

    for(int i = 0 ; i < val_y.rows() ; i++){
        if(i == 0){
            fprintf(gp,"plot \"%s\" using 1:2 with linespoints title \"%s\"\n",filename,key_info[0].c_str());
        }else{
            fprintf(gp,"replot \"%s\" using 1:%d with linespoints title \"%s\"\n",filename,i+2,key_info[i].c_str());
        }
    }
    pclose(gp);

    // std::cout << cmd.str() << std::endl;
    std::cout << "finish gnuplot" << std::endl;
    return 0;
}

int plot_field_twoaxis_global(const MatrixXd& val_x, const MatrixXd& val_y, std::vector<int> axis, std::vector<std::string> key_info, std::vector<std::string> graph_info){
    // std::cout << val_y.transpose() << std::endl;
    if(val_x.cols() != val_y.cols() || val_y.rows() != axis.size()){
        ERR("error : plot data should have the same cols");
    }
    // error 処理　axis,keyinfo,graphinfo size
    
    if(graph_info.size() == 0){
        graph_info.resize(4);
        graph_info[0] = "sample title";
        graph_info[1] = "xlabel";
        graph_info[2] = "y1label";
        graph_info[3] = "y2label";
    }
    if(key_info.size() == 0 ){
        key_info.resize(val_y.size());
        key_info[0] = "ref";
        for(int i = 1 ; i < key_info.size() ; i++){
            key_info[i] = "data" + std::to_string(i);
        }
    }

//  make datafile for gnuplot
    FILE* dataf;
    const char* filename = "datafile";
    dataf = fopen(filename,"w");
    fprintf(dataf,"# x");
    for(int i = 0; i < key_info.size(); i++){
        fprintf(dataf," %s",key_info[i].c_str());
    }
    fprintf(dataf,"\n");

    for(int i = 0 ; i < val_y.cols() ; i++){
        fprintf(dataf,"%f ",val_x(i));
        for(int j = 0 ; j < val_y.rows() ; j++){
            fprintf(dataf,"%f ",val_y(j,i));
        }
        fprintf(dataf,"\n");
    }
    fclose(dataf);

// gnuplot make graph
    FILE* gp;
    gp = popen("gnuplot -persist","w");
    fprintf(gp,"set grid \n");
    fprintf(gp,"set title \"%s\"\n", graph_info[0].c_str());
    fprintf(gp,"set xlabel \"%s\"\n",graph_info[1].c_str());
    fprintf(gp,"set ylabel \"%s\"\n",graph_info[2].c_str());
    fprintf(gp,"set y2label \"%s\"\n",graph_info[3].c_str());
    fprintf(gp,"set y2tics\n");

    for(int i = 0 ; i < val_y.rows() ; i++){
        int y_axis = axis[i];
        if(y_axis != 1 && y_axis != 2){
            ERR("gnuplot error: axis should be 1 or 2");
        }

        if(i == 0){
            fprintf(gp,"plot \"%s\" using 1:2 axis x1y%d with linespoints title \"%s\"\n",filename,y_axis,key_info[0].c_str());
        }else{
            fprintf(gp,"replot \"%s\" using 1:%d axis x1y%d with linespoints title \"%s\"\n",filename,i+2,y_axis,key_info[i].c_str());
        }
    }
    pclose(gp);

    // std::cout << cmd.str() << std::endl;
    std::cout << "finish gnuplot" << std::endl;
    return 0;
}