#include <iostream>
#include <vector>
#include <cmath>
#include <prettyprint.hpp>
#include "numerical.hpp"
using namespace std;

vector<vector<float> > explicit_solve(vector<vector<float>> *field, int now_step, float calc_time_step, float x_delta, float t_delta, pair<int, int> domain_col, float coef = 1.)
{
    int i, j, k;
    int calc_step = int(calc_time_step / t_delta);
    int start_col, end_col;
    if (domain_col.first == 0)
        start_col = 1;
    else
        start_col = domain_col.first;
    cout << "ok exp1" << endl;
    cout << "end col:" << domain_col.second << endl;
    cout << "start col:" << start_col << endl;

    float ratio = coef * t_delta / x_delta;
    // initialize
    printf("%lx %lx\n", (long)field, (long)&(*field)[0]);

    // solve
    for (i = now_step; i < calc_step + now_step; i++)
    {
        cout << "ok exp2 " << now_step <<" i=  " << i<< endl;
        for (j = start_col; j < domain_col.second; j++)
        {
            (*field)[i][j] = ratio * (*field)[i - 1][j - 1] + (1 - 2 * ratio) * (*field)[i - 1][j] + ratio * (*field)[i - 1][j + 1];
        }
    }

    cout << "end exp3" << endl;

    return *field;
}

vector<vector<float >> implicit_solve(vector<vector<float>> *field, int now_step, float calc_time_step, float x_delta, float t_delta, pair<int, int> domain_col, float coef = 1.)
{
    int i, j, k, coef_k, coef_j, field_k, field_j;
    int n_rows = (*field).size();
    int start_col, end_col;
    if (domain_col.first == 0)
        start_col = 1;
    else
        start_col = domain_col.first;
    end_col = domain_col.second;
    int f2c_idx = 0;
    int coef_n_cols = end_col - start_col;
    int calc_step = int(calc_time_step / t_delta);

    float ratio = coef * t_delta / x_delta;
    float substitute_sum;

    vector<float> coef_list{-ratio, 1 + 2 * ratio, -ratio};
    cout << "n_cols  " << (*field)[0].size() << endl;
    cout << "n_rows  " << n_rows << endl;
    // making coef matrix
    vector<vector<float>> coef_matrix(coef_n_cols, vector<float>(coef_n_cols, 0));

    vector<int> pivot_indexes;

    for (i = 0; i < coef_n_cols; i++)
    {
        if (i == 0)
        {
            coef_matrix[0][0] = 1 + 2 * ratio;
            coef_matrix[0][1] = -ratio;
        }
        else if (i == coef_n_cols-1)
        {
            coef_matrix[i][i] = 1 + 2 * ratio;
            coef_matrix[i][i - 1] = -ratio;
        }
        else
        {
            for (j = i - 1; j < i + 2; j++)
            {
                coef_matrix[i][j] = coef_list[1 + j - i];
            }
        }
    }

    // printf("ok1");

    cout << "coef matrix:" << coef_n_cols << "," << coef_n_cols << endl; 
    cout << "---------- coef mat initial ------------" << endl;
    cout << coef_matrix << endl;

    std::tie(coef_matrix, pivot_indexes) = lu_partition(coef_matrix);
    // coef_matrix = cholesky_factorize(coef_matrix);
    cout << "i = " << coef_matrix << endl;
    cout << "pv = " << pivot_indexes << endl;
    // exit(0);

    // start time step
    for (i = now_step; i <= now_step+calc_time_step; i++)
    {

        cout << "############### forward #############" << endl;
        cout << "imp start fsub i = " << i << endl;
        cout << "end col = " << end_col << " start col " << start_col << endl;

        cout << (*field)[i-1] <<endl;

        (*field)[i - 1][start_col] += (*field)[i - 1][start_col - 1] * ratio;
        (*field)[i - 1][end_col-1] += (*field)[i - 1][end_col] * ratio;
        // start forward substitution
        // cout << "i = " << i << endl;
        cout << "imp ok 1" << endl;

        cout << (*field)[i-1] <<endl;

        for (j = start_col; j < end_col; j++)
        {
            coef_j = j - start_col;
            field_j = pivot_indexes[coef_j] + f2c_idx + start_col;
            substitute_sum = 0.;
            cout << "j = " << j <<" coef_j " << coef_j << endl;

            for (k = start_col; k < j; k++)
            {
                field_k = pivot_indexes[coef_j] + f2c_idx + start_col;
                coef_k = k - start_col;
                cout << "j = " << j << " pivot index = " << pivot_indexes[coef_k] + f2c_idx << endl;
                substitute_sum += coef_matrix[coef_j][coef_k] * (*field)[i - 1][field_k];
            }
            cout << "substitute sum " << substitute_sum << endl;
            cout << "j = " << j << " pivot index after = " << pivot_indexes[coef_j] + f2c_idx << endl;
            cout << "field : " << (*field)[i - 1][field_j] << " pivot " << pivot_indexes[coef_j] << endl;
            cout << "bunshi : " << (*field)[i - 1][field_j] - substitute_sum << endl;
            // devide by coef_matrix[coef_j][coef_j] only once, backward, xor forward.
            (*field)[i][j + f2c_idx] = (*field)[i - 1][field_j] - substitute_sum ; 
        }
        cout << "imp ok 2" << endl;

        (*field)[i - 1][start_col] -= (*field)[i - 1][start_col - 1] * ratio;
        (*field)[i - 1][end_col-1] -= (*field)[i - 1][end_col] * ratio;

        cout << (*field)[i-1] <<endl;
        cout << (*field)[i] <<endl;

        cout << "imp ok 3" << endl;

        cout << "############### backward #############" << endl;


        // start backward substitution
        for (j = end_col-1; j >= start_col; j--)
        {
            // cout << "j = " << j << endl;
            coef_j = j - start_col;
            cout << "j = " << j <<" coef_j " << coef_j << endl;

            substitute_sum = 0.;
            for (k = end_col-1; k > j; k--)
            {
                coef_k = k - start_col;
                substitute_sum += coef_matrix[coef_j][coef_k] * (*field)[i][k + f2c_idx];
            }
            cout << "substitute sum " << substitute_sum << endl;
            cout << "j = " << j << " pivot index after = " << pivot_indexes[coef_j] + f2c_idx << endl;
            cout << "field : " << (*field)[i][j + f2c_idx] << " pivot " << pivot_indexes[coef_j] << endl;
            cout << "bunshi : " << (*field)[i][j + f2c_idx] - substitute_sum << endl;
            cout << "bunbo :" << coef_matrix[coef_j][coef_j] << endl;
            (*field)[i][j + f2c_idx] = ((*field)[i][j + f2c_idx] - substitute_sum) / coef_matrix[coef_j][coef_j];
        }
    }

    for (i=0;i<10;i++){
        cout << (*field)[i] << endl;

    }

    return *field;

}

vector<vector <float>> solve(vector<vector<float>> *field, vector<float> x_bound, float x_delta, float t_delta, vector<float> imp_boundary, vector<float> exp_boundary, float coef = 1., float exp_time_delta=-1, float imp_time_delta=-1)
{
    int now_step, step, j, k;
    int n_rows = (*field).size();
    int n_cols = (*field)[0].size();
    pair<int, int> imp_cols, exp_cols;
    int col0 = (int)((imp_boundary[0] - x_bound[0]) / x_delta);
    int col1 = (int)((imp_boundary[1] - x_bound[0]) / x_delta);
    imp_cols = make_pair(col0, col1);
    if (exp_time_delta < 0){
        exp_time_delta = t_delta;
    }
    if (imp_time_delta < 0){
        imp_time_delta = t_delta;
    }

    cout << "ok" << endl;

    col0 = (int)((exp_boundary[0] - x_bound[0]) / x_delta);
    col1 = (int)((exp_boundary[1] - x_bound[0]) / x_delta);
    exp_cols = make_pair(col0, col1);
    cout << "ok1" << endl;

    for (now_step = 1; now_step < n_rows; now_step++)
    {
        // explicit_solve(field, now_step, exp_time_delta, x_delta, t_delta, exp_cols, coef);

        implicit_solve(field, now_step, imp_time_delta, x_delta, t_delta, imp_cols, coef);
        break;
    }

    return *field;
}

vector<vector<float>> initialize(vector<vector<float>> *field, float initial_value, vector<float> x_bound){
    int i,j;
    int n_cols = (*field)[0].size();
    int n_rows = (*field).size();
    for(i=0;i<n_rows;i++){
        for(j=0;j<n_cols;j++){
            if(i==0 & j>0 & j<n_cols-1){
                (*field)[0][j] = initial_value;
            }else if(j==0){
                (*field)[i][j] = x_bound[0];
            }else if(j==n_cols-1){
                (*field)[i][j] = x_bound[1];
            }
        }
    }

    cout << "########## initial ###########" << endl;
    cout << (*field) << endl;

    return *field;
}

int main()
{

    // 初期値、境界値（x=0, 1, t=0 の3つ）、tiem delta、x_delta を引数にすると良さそう
    int n_cols = 11;         // 列数
    float initial_value = 1, coef=1; // 初期値
    float exp_time_delta = -1, imp_time_delta = -1;
    float t_delta, x_delta;
    int i, j, n_rows;
    vector<float> x_value = {0.0, 0.0};    // boundary value
    vector<float> x_bound = {0.0, 1.0};
    vector<float> t_bound = {0.0, 1.0};
    vector<float> imp_boundary = {0,1.0};
    vector<float> exp_boundary = {0.0,0.7};

    float x_length = x_bound[1] - x_bound[0];
    float time_length = t_bound[1] - t_bound[0];
    x_delta = x_length / float(n_cols-1);
    // t_delta = 0.1;
    t_delta = 0.01;
    // printf("x_delta = %f", x_delta);

    n_rows = (int)((t_bound[1] - t_bound[0]) / t_delta);
    vector<vector<float>> vec(n_rows, vector<float>(n_cols, 0));

    vec = initialize(&vec, initial_value, x_value);
    // vector<vector<float> > vec;

    vec = solve(&vec, x_bound, x_delta, t_delta, imp_boundary, exp_boundary, coef, exp_time_delta, imp_time_delta);
    // vec = implicit_solve(time_length, t_delta, x_length, x_delta, initial_value, x_value);

    for (int i = 0; i < vec.size(); i++)
    {
        printf("t=%f\t", i * t_delta);
        for (int j = 0; j < vec[0].size(); j++)
        {
            printf("%f ", vec[i][j]);
        }
        printf("\n");
    }
    return 0;
}

void print(vector<vector<float>> *vec)
{
    int n_rows = (*vec).size();
    int n_cols = (*vec)[0].size();
}