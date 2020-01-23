#include <iostream>
#include <vector>
#include <cmath>
#include <utility>
#include <prettyprint.hpp>
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

    int count = 0;
    for(i=0;i<n_rows;i++){
        omega = sqrt(coef_matrix[i][i]);
        // cout << omega << endl;
        for(j=i+1;j<n_rows;j++){
            coef_matrix[j][i] /= omega;
            for(k=i+1;k<=j;k++){
                coef_matrix[j][k] -= coef_matrix[j][i] * coef_matrix[k][i];
                count ++ ;
            }
        }
        coef_matrix[i][i] = omega;
    }
    // printf("count: %d %d", count, n_rows);
    // exit(0);

    return coef_matrix;
}

vector<vector <float>> sparse_lu_factorize(vector<vector<float>> coef_matrix){
    int n_rows = coef_matrix.size();
    int n_cols = coef_matrix[0].size();
    int center_idx = n_cols / 2;
    if (n_rows < n_cols)
    {
        cout << "expect rows > cols, actual " << n_rows << " <= " << n_cols << endl;
        exit(1);
    }

    int i,j,k, col_index, dense_row_index, max_col_index, tmp;
    float omega;

    for(i=0;i<n_rows;i++){
        // for (j = 0; j < n_rows; j++)
        // {
            // cout << "y=" << j << coef_matrix[j] << endl;
        // }
        omega = sqrt(coef_matrix[i][center_idx]);
        max_col_index = min(n_cols, center_idx + n_rows - i);
        // max_col_index = min(max_col_index, center_idx+i);
        // cout << "i:"<<i<<" max_col:" << max_col_index << endl;
        for(j=center_idx+1;j<max_col_index;j++){
            // cout << "j=" <<j;
            coef_matrix[i][j] /= omega;
            for(k=center_idx+1;k<=j;k++){
                // cout << " k=" << k;
                tmp = k - center_idx;
                dense_row_index = i + tmp;
                // cout << " (" << dense_row_index << "," << j-k << ") " ;
                coef_matrix[dense_row_index][j-tmp] -= coef_matrix[i][j] * coef_matrix[i][k];
                // count ++ ;
            }
            // cout << "\n";
        }
        coef_matrix[i][center_idx] = omega;
    }

    return coef_matrix;
}


vector<vector <float>> sparse_cholesky_factorize(vector<vector<float>> coef_matrix){
    int n_rows = coef_matrix.size();
    int n_cols = coef_matrix[0].size();
    // vector<vector<float> > cholesky_factorized_vec(n_rows, vector<float>(n_rows,0));
    // if (n_rows < n_cols)
    // {
    //     cout << "expect rows > cols, actual " << n_rows << " <= " << n_cols << endl;
    //     exit(1);
    // }

    int i,j,k, col_index, dense_row_index, max_col_index;
    float omega;

    // int count = 0;
    for(i=0;i<n_rows;i++){
        // for (j = 0; j < n_rows; j++)
        // {
            // cout << "y=" << j << coef_matrix[j] << endl;
        // }
        omega = sqrt(coef_matrix[i][0]);
        max_col_index = min(n_cols, n_rows - i);
        // cout << "i:"<<i<<" max_col:" << max_col_index << endl;
        for(j=1;j<max_col_index;j++){
            // cout << "j=" <<j;
            coef_matrix[i][j] /= omega;
            for(k=1;k<=j;k++){
                // cout << " k=" << k;
                dense_row_index = i + k;
                // cout << " (" << dense_row_index << "," << j-k << ") " ;
                coef_matrix[dense_row_index][j-k] -= coef_matrix[i][j] * coef_matrix[i][k];
                // count ++ ;
            }
            // cout << "\n";
        }
        coef_matrix[i][0] = omega;
    }
    // printf("count: %d %d", count, n_rows);
    // exit(0);

    return coef_matrix;
}
