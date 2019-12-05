// Cstruct2.cpp
#include <iostream>
#include <queue>
#include <vector>
#include <prettyprint.hpp>
// # include "practice.h"
using namespace std;

// void visualize(int num){
//     printf("Now: %d \n", num);
// }

// void add_class(Keisan k){
//     visualize(2);    
//     int c;
//     c = k.add();
//     printf("c = %d", c);
//     visualize(3);
// }



int main()
{
    int n_rows = 5;
    int n_cols = 5;
    vector<vector<int> > coef_matrix(n_rows,vector<int>(n_cols, 0));
    coef_matrix[0].assign(3,1);
    cout << coef_matrix << endl;


}