// Cstruct2.cpp
#include <iostream>
#include <vector>
#include <cmath>
#include "numerical.hpp"
#include <fstream>
#include <string>
#include <iomanip>
#include <map>
#include <filesystem>
#include <chrono>
#include "yaml-cpp/yaml.h"
#include <prettyprint.hpp>
// # include "practice.h"
using namespace std;
fs::path OUTPUT_FILE_PATH, LOG_PATH;
ofstream OUTPUT_FILE, LOG_FILE;

int main(int argc, char* argv[])
{
    if (argc < 2){
        printf("must specify config file in command line argument!\n");
        exit(1);
    }
    string line = "-----------------------------";

    printf("load confing file from %s\n", argv[1]);
    cout << line <<"\n";
    YAML::Node config = YAML::LoadFile(argv[1]);
    cout << config << endl;
    cout << line <<"\n";

    int space_dim = config["dim"].as<int>();
    float t_f;
    double t_d;
    string type_name = config["type"].as<string>();

    cout << "ok" << "\n";
    if (space_dim==1){
        if (type_name == "float"){
            pde(config, argv, t_f);
        }else if (type_name == "double"){
            pde(config, argv, t_d);
        }else{
            printf("config argument \'type_name\' must be in (\'float\', \'double\')!\n");
            exit(1);
        }
    }else if(space_dim==2){
        if (type_name == "float"){
            pde_2d(config, argv, t_f);
        }else if (type_name == "double"){
            pde_2d(config, argv, t_d);
        }else{
            printf("config argument \'type_name\' must be in (\'float\', \'double\')!\n");
            exit(1);
        }
    }else{
        printf("config argument \'type_name\' must be in (\'float\', \'double\')!\n");
        exit(1);
    }
    return 0;
}