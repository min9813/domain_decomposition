#include <iostream>
#include <vector>
#include <cmath>
#include <prettyprint.hpp>
#include <matplotlibcpp.h>
#include "numerical.hpp"
#include <fstream>
#include <string>
#include <iomanip>
#include <map>
using namespace std;
// namespace plt = matplotlibcpp;

vector<vector<vector<float>>> explicit_solve(vector<vector<vector<float>>> *field, int now_step, float calc_time_step, map<string, float> delta, map<int, pair<int,int>> domain_col, float coef = 1.)
{
    int t, i, j;
    int calc_step = int(calc_time_step / delta["t"]);
    int start_xid, end_xid, start_yid, end_yid;
    float x_diff, ratio_x, ratio_y, y_diff;
    if (domain_col[0].first == 0)
        start_xid = 1;
    else
        start_xid = domain_col[0].first;
    if(domain_col[1].first == 0)
        start_yid = 1;
    else
        start_yid = domain_col[1].first;
    // cout << "ok exp1" << endl;
    // cout << "end col:" << domain_col.second << endl;
    // cout << "start col:" << start_col << endl;

    ratio_x = coef * delta["t"] / delta["x"];
    ratio_y = coef * delta["t"] / delta["y"];

    // initialize
    // printf("%lx %lx\n", (long)field, (long)&(*field)[0]);

    // solve
    for (t = now_step; t < calc_step + now_step; t++)
    {
        // cout << "ok exp2 " << now_step <<" i=  " << i<< endl;
        for (i = start_xid; i < domain_col[0].second; i++)
        {
            for(j = start_yid; j < domain_col[1].second;j++){
                x_diff = (*field)[t-1][i+1][j] + (*field)[t-1][i-1][j] - 2 * (*field)[t-1][i][j];
                y_diff = (*field)[t-1][i][j+1] + (*field)[t-1][i][j-1] - 2 * (*field)[t-1][i][j];
                (*field)[t][i][j] = (*field)[t-1][i][j] + ratio_x * x_diff + ratio_y * y_diff;
            }
        }
    }

    // cout << "end exp3" << endl;

    return *field;
}

vector<vector<vector<float>>> implicit_solve(vector<vector<vector<float>>> *field, int now_step, float calc_time_step, map<string, float> delta, map<int, pair<int,int>> domain_col, float coef = 1.)
{
    int t, i, j, k, index, x_idx, y_idx, coef_idx, field_xid, field_yid, coef_k;
    int calc_step = int(calc_time_step / delta["t"]);
    int start_xid, end_xid, start_yid, end_yid;
    float x_diff, ratio_x, ratio_y, y_diff;
    if (domain_col[0].first == 0)
        start_xid = 1;
    else
        start_xid = domain_col[0].first;
    if(domain_col[1].first == 0)
        start_yid = 1;
    else
        start_yid = domain_col[1].first;
    // cout << "ok exp1" << endl;
    // cout << "end col:" << domain_col.second << endl;
    // cout << "start col:" << start_col << endl;

    ratio_x = coef * delta["t"] / delta["x"];
    ratio_y = coef * delta["t"] / delta["y"];

    end_xid = domain_col[0].second;
    end_yid = domain_col[1].second;

    int f2c_idx = 0;
    int x_idx_num, y_idx_num, coef_n_cols;
    x_idx_num = end_xid - start_xid;
    y_idx_num = end_yid - start_yid;
    coef_n_cols = x_idx_num * y_idx_num;

    float substitute_sum;

    vector<float> coef_list{-ratio_y, -ratio_x, 1 + 2 * ratio_x + 2*ratio_y, -ratio_x, -ratio_y};
    // vector<float> x0_coef_list{-ratio_y, -ratio_x, 1 + 2 * ratio_x + 2*ratio_y, -ratio_x, -ratio_y};

    // cout << "n_cols  " << (*field)[0].size() << endl;
    // cout << "n_rows  " << n_rows << endl;
    // making coef matrix
    vector<vector<float>> coef_matrix(coef_n_cols, vector<float>(coef_n_cols, 0));
    vector<int> index_list{-x_idx_num, -1, 0, 1, x_idx_num};

    // vector<int> pivot_indexes;

    for (i = 0; i < coef_n_cols; i++)
    {
        // if (i == 0)
        // {
        //     coef_matrix[i][i] = coef_list[2];
        //     coef_matrix[i][i+1] = coef_list[3];
        //     coef_matrix[i][i+x_idx_num] = coef_list[4];

        // }else if(i / x_idx_num == 0){
        //     coef_matrix[i][i] = coef_list[2];
        //     coef_matrix[i][i-1] = coef_list[1];
        //     coef_matrix[i][i+1] = coef_list[3];
        //     coef_matrix[i][i+x_idx_num] = coef_list[4];
        // }else if(i % x_idx_num == 1){
        //     coef_matrix[i][i] = coef_list[2];
        //     coef_matrix[i][i-x_idx_num] = coef_list[0];
        //     coef_matrix[i][i+1] = coef_list[3];
        //     coef_matrix[i][i+x_idx_num] = coef_list[4];
        // }else if(i / x_idx_num == (y_idx_num - 1)){
        //     coef_matrix[i][i] = coef_list[2];
        //     coef_matrix[i][i-x_idx_num] = coef_list[0];
        //     coef_matrix[i][i-1] = coef_list[1];
        //     coef_matrix[i][i+1] = coef_list[3];
        // }else if(i % x_idx_num == 0){
        //     coef_matrix[i][i] = coef_list[2];
        //     coef_matrix[i][i-x_idx_num] = coef_list[0];
        //     coef_matrix[i][i-1] = coef_list[1];
        //     coef_matrix[i][i+x_idx_num] = coef_list[4];
        // }
        for (j = 0; j < index_list.size(); j++)
        {
            index = i + index_list[j];
            if(index < 0 || index >= coef_n_cols) continue;
            if((index / x_idx_num != i / x_idx_num) && (index % x_idx_num != i % x_idx_num)){
                continue;
            }
            coef_matrix[i][i+index_list[j]] = coef_list[j];
        }
    }

    // printf("ok1");

    // cout << "coef matrix:" << coef_n_cols << "," << coef_n_cols << endl; 
    // cout << "---------- coef mat initial ------------" << endl;
    // cout << coef_matrix << endl;

    // std::tie(coef_matrix, pivot_indexes) = lu_partition(coef_matrix);
    coef_matrix = cholesky_factorize(coef_matrix);
    // cout << "i = " << coef_matrix << endl;
    // cout << "pv = " << pivot_indexes << endl;
    // exit(0);

    // start time step
    for (i = now_step; i <= now_step+calc_time_step; i++)
    {

        // cout << "############### forward #############" << endl;
        // cout << "imp start fsub i = " << i << endl;
        // cout << "end col = " << end_col << " start col " << start_col << endl;

        // cout << (*field)[i-1] <<endl;
        for(x_idx=start_xid;x_idx<end_xid;x_idx++){
            if(x_idx == start_xid){
                (*field)[i - 1][start_xid][start_yid] += (*field)[i - 1][start_xid - 1][start_yid] * coef_list[1];
                (*field)[i - 1][start_xid][start_yid] += (*field)[i - 1][start_xid][start_yid-x_idx_num] * coef_list[0];
                (*field)[i - 1][start_xid][end_yid-1] += (*field)[i - 1][start_xid - 1][end_yid-1] * coef_list[3];
                (*field)[i - 1][start_xid][end_yid-1] += (*field)[i - 1][start_xid][end_yid] * coef_list[4];
            }else if(x_idx == end_xid-1){
                (*field)[i - 1][x_idx][start_yid] += (*field)[i - 1][x_idx+1][start_yid] * coef_list[1];
                (*field)[i - 1][x_idx][start_yid] += (*field)[i - 1][x_idx][start_yid-1] * coef_list[0];
                (*field)[i - 1][x_idx][end_yid-1] += (*field)[i - 1][x_idx+1][end_yid-1] * coef_list[3];
                (*field)[i - 1][x_idx][end_yid-1] += (*field)[i - 1][x_idx][end_yid] * coef_list[4];
                
            }else{
                (*field)[i-1][x_idx][start_yid] += (*field)[i-1][x_idx][start_yid-1] * coef_list[0];
                (*field)[i-1][x_idx][end_yid-1] += (*field)[i-1][x_idx][end_yid] * coef_list[4];
            }
        }
        for(y_idx=start_yid+1;y_idx<end_yid-1;y_idx++){
            (*field)[i-1][start_xid][y_idx] = (*field)[i-1][start_xid-1][y_idx] * coef_list[1];
            (*field)[i-1][end_xid-1][y_idx] = (*field)[i-1][end_xid][y_idx] * coef_list[3];
        }

        // start forward substitution
        // cout << "i = " << i << endl;
        // cout << "imp ok 1" << endl;

        // cout << (*field)[i-1] <<endl;

        for (x_idx = start_xid; x_idx < end_xid; x_idx++)
        {
            for(y_idx = start_yid; y_idx < end_yid; y_idx++){
                coef_idx = (x_idx - start_xid) + (y_idx - start_yid) * x_idx_num;
                field_xid = j + f2c_idx;
                substitute_sum = 0.;
                }
            
        }
        // cout << "imp ok 2" << endl;

        (*field)[i - 1][start_xid] -= (*field)[i - 1][start_xid - 1] * ratio;
        (*field)[i - 1][end_xid-1] -= (*field)[i - 1][end_xid] * ratio;

        // cout << (*field)[i-1] <<endl;
        // cout << (*field)[i] <<endl;

        // cout << "imp ok 3" << endl;

        // cout << "############### backward #############" << endl;


        // start backward substitution
        for (j = end_xid-1; j >= start_xid; j--)
        {
            // cout << "j = " << j << endl;
            coef_j = j - start_xid;
            // cout << "j = " << j <<" coef_j " << coef_j << endl;

            substitute_sum = 0.;
            for (k = end_xid-1; k > j; k--)
            {
                coef_k = k - start_xid;
                // coef_k = end_col-coef_k;
                // cout << "coef_j " << coef_j << " coef_k " << coef_k << " field k :" << k << endl;
                // cout << "coef value: " << coef_matrix[coef_k][coef_j] << " field valud : " << (*field)[i][k + f2c_idx] << endl;
                // cout << " product ; " << coef_matrix[coef_k][coef_j] * (*field)[i][k + f2c_idx] << endl;
                substitute_sum += coef_matrix[coef_k][coef_j] * (*field)[i][k + f2c_idx];
            }
            // cout << "substitute sum " << substitute_sum << endl;
            // cout << "j = " << j << " pivot index after = " << pivot_indexes[coef_j] + f2c_idx << endl;
            // cout << "j = " << j << endl;
            // cout << "field : " << (*field)[i][j + f2c_idx] << endl;
            // cout << "bunshi : " << (*field)[i][j + f2c_idx] - substitute_sum << endl;
            // cout << "bunbo :" << coef_matrix[coef_j][coef_j] << endl;
            (*field)[i][j + f2c_idx] = ((*field)[i][j + f2c_idx] - substitute_sum) / coef_matrix[coef_j][coef_j];
        }
    }

    // for (i=0;i<10;i++){
        // cout << (*field)[i] << endl;

    // }

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

    // cout << "ok" << endl;

    col0 = (int)((exp_boundary[0] - x_bound[0]) / x_delta);
    col1 = (int)((exp_boundary[1] - x_bound[0]) / x_delta);
    exp_cols = make_pair(col0, col1);
    // cout << "ok1" << endl;

    for (now_step = 1; now_step < n_rows; now_step++)
    {
        explicit_solve(field, now_step, exp_time_delta, x_delta, t_delta, exp_cols, coef);

        implicit_solve(field, now_step, imp_time_delta, x_delta, t_delta, imp_cols, coef);
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

    // cout << "########## initial ###########" << endl;
    // cout << (*field) << endl;

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
    vector<float> t_bound = {0.0, 2.0};
    vector<float> imp_boundary = {0.7,1.0};
    vector<float> exp_boundary = {0.0,0.7};

    string filename = "pde_output.txt";
    ofstream writing_file;

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

    float t_now, x_now;

    writing_file.open(filename, ios::out);

    cout << "writing " << filename << "..." << endl;

    for (i = 0; i < vec.size(); i++)
    {
        t_now = t_bound[0] + i * t_delta;
        if(i==0){
            printf("time,");
            writing_file << "time,";
            for(j=0;j < n_cols; j++){
                x_now = x_bound[0] + x_delta * j;
                cout << fixed << setprecision(3) << x_now;
                writing_file << "x=" <<fixed <<setprecision(3)<< x_now;
                if (j<n_cols-1){
                    cout << "\t";
                    writing_file << ",";
                }
            }
            writing_file << endl;
            cout << endl;
        }
        printf("t=%f\t", t_now);
        writing_file << t_now << ",";

        for (j = 0; j < vec[0].size(); j++)
        {
            if(j < vec[0].size()-1){
                printf("%f\t", vec[i][j]);
                writing_file << vec[i][j] << ",";
            }else{
                printf("%f", vec[i][j]);
                writing_file << vec[i][j];
            }
            
        }
        writing_file << endl;
        printf("\n");
    }

    writing_file.close();
    return 0;
}

void print(vector<vector<float>> *vec)
{
    int n_rows = (*vec).size();
    int n_cols = (*vec)[0].size();
}