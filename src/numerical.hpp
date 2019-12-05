#ifndef _NUMERICAL_MINTEI
#define _NUMERICAL_MINTEI
#include <iostream>
#include <vector>
#include <cmath>
#include <utility>
using namespace std;

tuple<vector<vector<float>>, vector<int>> lu_partition(vector<vector<float>> coef_matrix);
vector<vector<float>> cholesky_factorize(vector<vector<float>> coef_matrix);
// vector<vector<float>> lu_partition(vector<vector<float>> coef_matrix);

// class DomainPDE
// {
//     public:
//         void explicit_solve(int time);
//         void implicit_solve(int time);
// }

#endif
