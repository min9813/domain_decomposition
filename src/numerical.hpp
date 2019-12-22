#ifndef _NUMERICAL_MINTEI
#define _NUMERICAL_MINTEI
#include <iostream>
#include <vector>
#include <cmath>
#include <utility>
#include <map>
#include <filesystem> // make directory
#include <fstream>
#include "yaml-cpp/yaml.h"
using namespace std;
namespace fs = filesystem;

tuple<vector<vector<float>>, vector<int>> lu_partition(vector<vector<float>> coef_matrix);
vector<vector<float>> cholesky_factorize(vector<vector<float>> coef_matrix);
vector<vector<float>> exact_solution_1d(vector<vector<float>> *field, map<string, pair<float, float>> spatio_bound, map<string, float> delta, float coef);
vector<vector<vector<float>>> exact_solution_2d(vector<vector<vector<float>>> *field, map<string, pair<float, float>> spatio_bound, map<string, float> delta, float coef);
string float_to_string(float value, int precision);
int pde_2d_main(YAML::Node config, char* argv[], float type_value);
int pde(YAML::Node config, char* argv[], float type_value);
int pde_2d(YAML::Node config, char* argv[], float type_value);

extern fs::path OUTPUT_FILE_PATH, LOG_PATH;
extern ofstream OUTPUT_FILE, LOG_FILE;

// vector<vector<float>> lu_partition(vector<vector<float>> coef_matrix);

// class DomainPDE
// {
//     public:
//         void explicit_solve(int time);
//         void implicit_solve(int time);
// }

#endif
