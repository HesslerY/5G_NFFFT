#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cassert>
#include <vector>
#include <random>

using namespace std;

int main(int argc, char** argv){
    std::cout << "koike made this file\n";


    std::string filename_in = "data/result.txt";
    std::ifstream inf(filename_in);
    if(!inf){
        std::cout << "error cannot open the file" << std::endl;
    }
    std::string str;
    while(getline(inf, str)){
        // std::cout << str << std::endl;
        std::stringstream ss(str);
        std::string element;
        getline(ss,element,'=');
        getline(ss,element,'=');
        std::cout << "z = " << element << std::endl;
    }
}