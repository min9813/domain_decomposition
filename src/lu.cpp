#include <iostream>
#include <vector>
#include <cmath>
#include <utility>
using namespace std;

tuple<vector<vector<float>>, vector<int> > lu_partition(vector<vector<float>> coef_matrix){
// vector<vector<float>> lu_partition(vector<vector<float>> coef_matrix)

    int n_rows = coef_matrix.size();
    int n_cols = coef_matrix[0].size();
    if (n_rows != n_cols)
    {
        cout << "expect columns = rows, actual " << n_cols << " = " << n_rows << endl;
        exit(1);
    }
    int i, j, k, tmp_pivot_idx;
    bool need_pivot = false;
    float tmp_c, pivot_value, abs_cand, omega;
    vector<int> pivot_indexes(n_rows, 0);
    for (int i = 0; i < n_rows; i++)
    {
        pivot_indexes[i] = i;
    };
    // start LU, with row's pivot exchange.
    for (i = 0; i < n_rows-1; i++)
    {
        pivot_value = abs(coef_matrix[i][i]);
        for(tmp_pivot_idx=i,k=i+1;k<n_rows;k++){
            abs_cand = abs(coef_matrix[k][i]);
            if(abs_cand>pivot_value){
                pivot_value = abs_cand;
                tmp_pivot_idx = k;
            }
        }
        swap(coef_matrix[i], coef_matrix[tmp_pivot_idx]);
        swap(pivot_indexes[i], pivot_indexes[tmp_pivot_idx]);
        omega = 1/ coef_matrix[i][i];

        for(j=i+1;j<n_rows;j++){
            coef_matrix[j][i] *= omega; // calc L
            for(k=i+1;k<n_cols;k++){
                coef_matrix[j][k] -= coef_matrix[i][k] * coef_matrix[j][i];
            }
        }
    }
    return forward_as_tuple(coef_matrix, pivot_indexes);
    // return coef_matrix;

}

vector<vector <float>> cholesky_factorize(vector<vector<float>> coef_matrix){
    int n_rows = coef_matrix.size();
    int n_cols = coef_matrix[0].size();
    // vector<vector<float> > cholesky_factorized_vec(n_rows, vector<float>(n_rows,0));
    if (n_rows != n_cols)
    {
        cout << "expect columns = rows, actual " << n_cols << " = " << n_rows << endl;
        exit(1);
    }

    int i,j,k;
    float omega;

    for(i=0;i<n_rows;i++){
        omega = sqrt(coef_matrix[i][i]);
        cout << omega << endl;
        for(j=i+1;j<n_rows;j++){
            coef_matrix[j][i] /= omega;
            for(k=i+1;k<=j;k++){
                coef_matrix[j][k] -= coef_matrix[j][i] * coef_matrix[k][i];
            }
        }
        coef_matrix[i][i] = omega;
    }

    return coef_matrix;
}
