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
#include <filesystem>
#include <chrono>
#include "yaml-cpp/yaml.h"
using namespace std;
namespace fs = filesystem;
// namespace plt = matplotlibcpp;
static vector<vector<float>> field;
static double choleskey_elapsed_time = 0, gauss_time = 0;
static bool is_sep = false;

void explicit_solve(int now_step, float calc_time_step, map<string, float> delta, map<string, pair<float, float>> domain_col, float coef = 1., pair<float, float> neighb_coef={1.0, 1.0})
{
    int t, i, j;
    int calc_step = int(calc_time_step / delta["t"]);
    int start_xid, end_xid;
    float x_diff, ratio_x;
    if (domain_col["x"].first == 0)
        start_xid = 1;
    else
        start_xid = domain_col["x"].first;
    // printf("ok exp1");
    // cout << "ok exp1" << endl;
    // cout << "end col:" << domain_col.second << endl;
    // cout << "start col:" << start_col << endl;
    float ratio = coef * delta["t"] / (delta["x"]*delta["x"]);
    float ratio_small = 2 * (coef * neighb_coef.first) / (coef + neighb_coef.first) * delta["t"] / (delta["x"]*delta["x"]);
    float ratio_large = 2 * (coef * neighb_coef.second) / (coef + neighb_coef.second) * delta["t"] / (delta["x"]*delta["x"]);
    if(ratio_large >= 1./2.){
        ratio_large = ratio;
    }
    if(ratio_small >= 1./2.){
        ratio_small = ratio;
    }
    if(ratio>=1./2.){
        printf("must satisfy von-neuman's stable conditions. actual ratio = %f", ratio);
        exit(1);
    }
    // printf("coef=%f, coef_large=%f, coef_small=%f\n", ratio, ratio_large, ratio_small);
    // cout << "exp ratio: " << ratio << endl;
    // initialize
    // printf("%lx %lx\n", (long)field, (long)&field[0]);

    // solve
    for (i = now_step; i < calc_step + now_step; i++)
    {
        // cout << "ok exp2 " << now_step <<" i=  " << i<< endl;
        if(start_xid < domain_col["x"].second){
            j = start_xid;
            field[i][j] = ratio_small * field[i - 1][j - 1] + (1 - 2 * ratio) * field[i - 1][j] + ratio * field[i - 1][j + 1];

            for (j = start_xid+1; j < domain_col["x"].second-1; j++)
            {
                field[i][j] = ratio * field[i - 1][j - 1] + (1 - 2 * ratio) * field[i - 1][j] + ratio * field[i - 1][j + 1];
            }
            if(j < domain_col["x"].second){
                field[i][j] = ratio * field[i - 1][j - 1] + (1 - 2 * ratio) * field[i - 1][j] + ratio_large * field[i - 1][j + 1];
            }

        }
        
    }
    // cout << field[now_step-1] <<endl;
    // cout << field[now_step] <<endl;

    // cout << "end exp3" << endl;

}

void implicit_solve(int now_step, float calc_time_step, map<string, float> delta, map<string, pair<float, float>> domain_col, float coef = 1., pair<float, float> neighb_coef={1.0, 1.0})
{
    int i, j, k, coef_k, coef_j, field_k, field_j;
    int n_rows = field.size();
    int start_col, end_col;
    if (domain_col["x"].first == 0)
        start_col = 1;
    else
        start_col = domain_col["x"].first;
    end_col = domain_col["x"].second;
    int f2c_idx = 0;
    int coef_n_cols = end_col - start_col;
    int calc_step = int(calc_time_step / delta["t"]);

    float ratio = coef * delta["t"] / (delta["x"]*delta["x"]);
    float ratio_small = 2 * (neighb_coef.first * coef) / (neighb_coef.first + coef) * delta["t"] / (delta["x"] * delta["x"]);
    float ratio_large = 2 * (neighb_coef.second * coef) / (neighb_coef.second + coef) * delta["t"] / (delta["x"] * delta["x"]);
    
    float substitute_sum;

    vector<float> coef_list{-ratio, 1 + 2 * ratio, -ratio};
    vector<float> coef_list_small{-ratio_small, 1 + 2 * ratio, -ratio};
    vector<float> coef_list_large{-ratio, 1 + 2 * ratio, -ratio_large};

    // cout << "n_cols  " << field[0].size() << endl;
    // cout << "n_rows  " << n_rows << endl;
    // making coef matrix
    vector<vector<float>> coef_matrix(coef_n_cols, vector<float>(coef_n_cols, 0));

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

        // cout << field[i-1] <<endl;
        // cout << "-------- prev field org----------" << endl;
        // cout << field[i-1] <<endl;
        field[i - 1][start_col] -= field[i - 1][start_col - 1] * coef_list[0];
        field[i - 1][end_col-1] -= field[i - 1][end_col] * coef_list[2];
        // start forward substitution
        // cout << "i = " << i << endl;
        // cout << "imp ok 1" << endl;
        
        // cout << field[i-1] <<endl;
        // cout << "-------- prev field after add ----------" << endl;
        // cout << field[i-1] <<endl;

        for (j = start_col; j < end_col; j++)
        {
            coef_j = j - start_col;
            field_j = j + f2c_idx;
            substitute_sum = 0.;
            // cout << "j = " << j <<" coef_j " << coef_j << endl;

            for (k = start_col; k < j; k++)
            {
                field_k = k + f2c_idx;
                coef_k = k - start_col;
                // cout << "j = " << j << " pivot index = " << pivot_indexes[coef_k] + f2c_idx << endl;
                // cout << "coef j = " << coef_j << " coef k " << coef_k << endl;
                // in forwarding, use the variable you want to solve. in this case, field[i], not [i-1]
                substitute_sum += coef_matrix[coef_j][coef_k] * field[i][field_k];
            }
            // cout << "substitute sum " << substitute_sum << endl;
            // cout << "j = " << j << " pivot index after = " << pivot_indexes[coef_j] + f2c_idx << endl;
            // cout << "field : " << field[i - 1][field_j] << " pivot " << pivot_indexes[coef_j] << endl;
            // cout << "bunshi : " << field[i - 1][field_j] - substitute_sum << endl;
            // cout << "bunbo :" << coef_matrix[coef_j][coef_j] << endl;

            // devide by coef_matrix[coef_j][coef_j] only once, backward, xor forward.
            field[i][field_j] = (field[i - 1][field_j] - substitute_sum) / coef_matrix[coef_j][coef_j]; 
        }
        // cout << "imp ok 2" << endl;

        field[i - 1][start_col] += field[i - 1][start_col - 1] * coef_list[0];
        field[i - 1][end_col-1] += field[i - 1][end_col] * coef_list[2];

        // cout << "-------- now field ----------" << endl;
        // cout << field[i] <<endl;

        // cout << "imp ok 3" << endl;

        // cout << "############### backward #############" << endl;


        // start backward substitution
        for (j = end_col-1; j >= start_col; j--)
        {
            // cout << "j = " << j << endl;
            coef_j = j - start_col;
            // cout << "j = " << j <<" coef_j " << coef_j << endl;

            substitute_sum = 0.;
            for (k = end_col-1; k > j; k--)
            {
                coef_k = k - start_col;
                // coef_k = end_col-coef_k;
                // cout << "coef_j " << coef_j << " coef_k " << coef_k << " field k :" << k << endl;
                // cout << "coef value: " << coef_matrix[coef_k][coef_j] << " field valud : " << field[i][k + f2c_idx] << endl;
                // cout << " product ; " << coef_matrix[coef_k][coef_j] * field[i][k + f2c_idx] << endl;
                substitute_sum += coef_matrix[coef_k][coef_j] * field[i][k + f2c_idx];
            }
            // cout << "substitute sum " << substitute_sum << endl;
            // cout << "j = " << j << " pivot index after = " << pivot_indexes[coef_j] + f2c_idx << endl;
            // cout << "j = " << j << endl;
            // cout << "field : " << field[i][j + f2c_idx] << endl;
            // cout << "bunshi : " << field[i][j + f2c_idx] - substitute_sum << endl;
            // cout << "bunbo :" << coef_matrix[coef_j][coef_j] << endl;
            field[i][j + f2c_idx] = (field[i][j + f2c_idx] - substitute_sum) / coef_matrix[coef_j][coef_j];
        }
    }

    // for (i=0;i<10;i++){
    // cout << field[now_step] << endl;

    // }
}

void sparse_implicit_solve(int now_step, float calc_time_step, map<string, float> delta, map<string, pair<float, float>> domain_col, float coef = 1., pair<float, float>neighb_coef={1.0, 1.0})
{
    int i, j, k, coef_k, coef_j, field_k, field_j;
    int n_rows = field.size();
    int start_col, end_col;
    if (domain_col["x"].first == 0)
        start_col = 1;
    else
        start_col = domain_col["x"].first;
    end_col = domain_col["x"].second;
    int f2c_idx = 0;
    int coef_n_rows = end_col - start_col;
    int coef_n_cols = 2;
    int calc_step = int(calc_time_step / delta["t"]);

    float ratio = coef * delta["t"] / (delta["x"]*delta["x"]);
    float ratio_small = 2 * (coef * neighb_coef.first) / (coef + neighb_coef.first) * delta["t"] / (delta["x"]*delta["x"]);
    float ratio_large = 2 * (coef * neighb_coef.second) / (coef + neighb_coef.second) * delta["t"] / (delta["x"]*delta["x"]);
    float substitute_sum;
    
    vector<float> coef_list_main{1 + 2 * ratio, -ratio};
    vector<float> coef_list_small{1 + 2 * ratio, -ratio_small};
    vector<float> coef_list_large{1 + 2 * ratio, -ratio_large};

    // cout << "n_cols  " << field[0].size() << endl;
    // cout << "n_rows  " << n_rows << endl;
    // making coef matrix
    vector<float> coef_list;
    vector<vector<float>> coef_matrix(coef_n_rows, vector<float>(coef_n_cols, 0));
    // cout << "coef list:" << coef_list_main <<endl;
    for (i = 0; i < coef_n_rows; i++)
    {
        if (i == 0)
        {
            coef_matrix[0][0] = 1 + 2 * ratio;
            coef_matrix[0][1] = -ratio;
        }
        else if (i == coef_n_rows-1)
        {
            coef_matrix[i][0] = 1 + 2 * ratio;
        }
        else
        {
            for (j = 0; j < coef_list_main.size(); j++)
            {
                coef_matrix[i][j] = coef_list_main[j];
            }
        }
    }

    // printf("ok1");

    // cout << "coef matrix:" << coef_n_cols << "," << coef_n_cols << endl; 
    // cout << "---------- coef mat initial ------------" << endl;
    // cout << coef_matrix << endl;
    // cout << "------- before csolve ---------" << endl;
    // for (i = 0; i < coef_n_rows; i++)
    // {
        // cout << "y=" << i << coef_matrix[i] << endl;
    // }
    // exit(0);
    // std::tie(coef_matrix, pivot_indexes) = lu_partition(coef_matrix);
    std::chrono::system_clock::time_point start, end;
    float elapsed = 0.;
    start = std::chrono::system_clock::now();

    coef_matrix = sparse_cholesky_factorize(coef_matrix);
    
    end = std::chrono::system_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0;
    choleskey_elapsed_time += elapsed;

    // cout << "i = " << coef_matrix << endl;
    // cout << "pv = " << pivot_indexes << endl;
    // exit(0);
    // cout << "------- after csolve ---------" << endl;
    // for (i = 0; i < coef_n_rows; i++)
    // {
    //     cout << "y=" << i << coef_matrix[i] << endl;
    // }
    int max_coef_row_index, add_index;
    // start time step
    start = std::chrono::system_clock::now();

    for (i = now_step; i <= now_step+calc_time_step; i++)
    {

        // cout << "############### forward #############" << endl;
        // cout << "imp start fsub i = " << i << endl;
        // cout << "end col = " << end_col << " start col " << start_col << endl;

        // cout << field[i-1] <<endl;
        // cout << "-------- prev field org----------" << endl;
        // cout << field[i-1] <<endl;
        field[i - 1][start_col] -= field[i - 1][start_col - 1] * coef_list_small[1];
        field[i - 1][end_col-1] -= field[i - 1][end_col] * coef_list_large[1];
        // start forward substitution
        // cout << "i = " << i << endl;
        // cout << "imp ok 1" << endl;
        
        // cout << field[i-1] <<endl;
        // cout << "-------- prev field after add ----------" << endl;
        // cout << field[i-1] <<endl;

        for (j = start_col; j < end_col; j++)
        {
            coef_j = j - start_col;
            field_j = j + f2c_idx;
            substitute_sum = 0.;
            max_coef_row_index = min(coef_j, coef_n_cols-1);
            // cout << "max coef row index:" << max_coef_row_index << endl;
            if (max_coef_row_index > 0){
                substitute_sum += coef_matrix[coef_j-1][1] * field[i][field_j-1];
            }
            field[i][field_j] = (field[i-1][field_j] - substitute_sum) / coef_matrix[coef_j][0];
            // cout << "j = " << j <<" coef_j " << coef_j <<" field j " << field_j << " sub sum: " << substitute_sum << " " << field[i][field_j]<< endl;

            // for (k = 0; k < max_coef_row_index; k++)
            // {
            //     tmp_field_id = k + add_index;
            //     // printf("k=%d coef=(%d, %d) tmp_y=%d tmp_x = %d\n", k, add_index+k, max_gauss_coef_index-k, tmp_field_yid, tmp_field_xid);

            //     substitute_sum += coef_matrix[tmp_field_id][max_gauss_coef_index-k] * field[i][tmp_field_yid][tmp_field_xid];
            // }
            // cout << "substitute sum " << substitute_sum << endl;
            // cout << "j = " << j << " pivot index after = " << pivot_indexes[coef_j] + f2c_idx << endl;
            // cout << "field : " << field[i - 1][field_j] << " pivot " << pivot_indexes[coef_j] << endl;
            // cout << "bunshi : " << field[i - 1][field_j] - substitute_sum << endl;
            // cout << "bunbo :" << coef_matrix[coef_j][coef_j] << endl;

            // devide by coef_matrix[coef_j][coef_j] only once, backward, xor forward.
        }

        // cout << "imp ok 2" << endl;

        field[i - 1][start_col] += field[i - 1][start_col - 1] * coef_list_small[1];
        field[i - 1][end_col-1] += field[i - 1][end_col] * coef_list_large[1];

        // cout << "-------- now field before bacjward ----------" << endl;
        // cout << field[i] <<endl;

        // cout << "imp ok 3" << endl;

        // cout << "############### backward #############" << endl;


        // start backward substitution
        for (j = end_col-1; j >= start_col; j--)
        {
            // cout << "j = " << j << endl;
            coef_j = j - start_col;
            field_j = j + f2c_idx;
            max_coef_row_index = min(coef_n_rows-coef_j-1, coef_n_cols-1);

            // cout << "j = " << j <<" coef_j " << coef_j << endl;

            substitute_sum = 0.;
            if (max_coef_row_index > 0){
                substitute_sum += coef_matrix[coef_j][1] * field[i][field_j+1];
            }
            // cout << "j = " << j <<" coef_j " << coef_j <<" field j " << field_j <<"max:" <<max_coef_row_index << " sub sum: " << substitute_sum << " " << field[i][field_j]<< endl;

            // cout << "substitute sum " << substitute_sum << endl;
            // cout << "j = " << j << " pivot index after = " << pivot_indexes[coef_j] + f2c_idx << endl;
            // cout << "j = " << j << endl;
            // cout << "field : " << field[i][j + f2c_idx] << endl;
            // cout << "bunshi : " << field[i][j + f2c_idx] - substitute_sum << endl;
            // cout << "bunbo :" << coef_matrix[coef_j][coef_j] << endl;
            field[i][field_j] = (field[i][field_j] - substitute_sum) / coef_matrix[coef_j][0];
        }
    }
    end = std::chrono::system_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0;
    gauss_time += elapsed;
}

void explicit_solve_whole(int now_step, float calc_time_step, map<string, float> delta, vector<map<string, pair<float, float>>> domain_col){
    int t, i, j, k;
    int calc_step = int(calc_time_step / delta["t"]);
    int start_xid, end_xid;
    float x_diff, ratio, ratio_tx, coef, coef_2, ratio_neighb;
    
    // printf("ok exp1");
    // cout << "ok exp1" << endl;
    // cout << "end col:" << domain_col.second << endl;
    // cout << "start col:" << start_col << endl;

    // cout << "exp ratio: " << ratio << endl;
    // cout << "explicit " << domain_col << endl;
    // initialize
    // printf("%lx %lx\n", (long)field, (long)&field[0]);
    ratio_tx = delta["t"] / (delta["x"]*delta["x"]);

    // solve
    for (i = now_step; i < calc_step + now_step; i++)
    {
        for (j = 0; j < domain_col.size(); j++)
        {
            
            if (domain_col[j]["is_exp"].first == 0){
                printf("only for explicit method.\n");
                exit(1);
            }
            coef = domain_col[j]["coef"].first;

            if (domain_col[j]["x"].first == 0)
                start_xid = 1;
            else
                start_xid = domain_col[j]["x"].first;
            // if(domain_col[j]["t_delta"].first > 0){
            //     delta["t"] = domain_col[j]["t_delta"].first;
            // }
            ratio = ratio_tx * coef;
            if (ratio >= 1./2.){
                printf("must satisfy von-neuman's stable conditions. actual ratio = %f", ratio);
                exit(1);
            }
            if(j>0){
                coef_2 = domain_col[j-1]["coef"].first;
                ratio_neighb = ratio_tx * 2 * (coef_2 * coef) / (coef_2 + coef);
            }else{
                ratio_neighb = ratio;
            }
            k = start_xid;
            field[i][k] = ratio_neighb * field[i - 1][k - 1] + (1 - 2 * ratio) * field[i - 1][k] + ratio * field[i - 1][k + 1];

            for (k = start_xid+1; k < domain_col[j]["x"].second-1; k++)
            {
                field[i][k] = ratio * field[i - 1][k - 1] + (1 - 2 * ratio) * field[i - 1][k] + ratio * field[i - 1][k + 1];
            }
            if(j<domain_col.size()-1){
                coef_2 = domain_col[j+1]["coef"].first;
                ratio_neighb = ratio_tx * 2 * (coef_2 * coef) / (coef_2 + coef);
            }else{
                ratio_neighb = ratio;
            }
            field[i][k] = ratio * field[i - 1][k - 1] + (1 - 2 * ratio) * field[i - 1][k] + ratio_neighb * field[i - 1][k + 1];
            
        }
    }
}

void implicit_solve_whole(int now_step, float calc_time_step, map<string, float> delta, vector<map<string, pair<float, float>>> domain_col){
    int i, j, k, coef_k, coef_j, field_k, field_j;
    int n_rows = field.size();
    int start_col, end_col;
    if (domain_col[0]["x"].first == 0)
        start_col = 1;
    else
        start_col = domain_col[0]["x"].first;
    end_col = domain_col[domain_col.size()-1]["x"].second;
    int f2c_idx = 0;
    int coef_n_rows = end_col - start_col;
    int coef_n_cols = 2;
    int calc_step = int(calc_time_step / delta["t"]);
    bool small_edge, large_edge;

    float ratio_tx = delta["t"] / (delta["x"]*delta["x"]), coef1, coef2;
    float substitute_sum, ratio, large_edge_ratio, small_edge_ratio;

    vector<float> coef_list{1 + 2 * ratio, -ratio};
    // making coef matrix
    vector<vector<float>> coef_matrix(coef_n_rows, vector<float>(coef_n_cols, 0));

    vector<int> pivot_indexes;

    k = 0;
    for (i = 0; i < coef_n_rows; i++)
    {
        small_edge = (i+start_col) == domain_col[k]["x"].second;
        large_edge = (i+start_col) == (domain_col[k]["x"].second - 1);
        large_edge &= k < domain_col.size() - 1;
        if(small_edge){
            k ++;
        }else if(large_edge){
            coef1 = domain_col[k]["coef"].first;
            coef2 = domain_col[k+1]["coef"].first;
            // cout << "large " << i << endl;
            large_edge_ratio = 2. * (coef1 * coef2) / (coef1 + coef2);
            large_edge_ratio = large_edge_ratio * ratio_tx;
            // cout << "large ratio:" << large_edge_ratio;
        }
        ratio = domain_col[k]["coef"].first * ratio_tx;
        // cout << "ratio " << ratio << endl;
        if (i == 0)
        {
            coef_matrix[i][0] = 1 + 2 * ratio;
            coef_matrix[i][1] = - ratio;
        }
        else if (i == coef_n_rows-1)
        {
            coef_matrix[i][0] = 1 + 2 * ratio;
        }else if(large_edge){
            coef_matrix[i][0] = 1 + 2 * ratio;
            coef_matrix[i][1] = - large_edge_ratio;
        }
        else
        {
            coef_matrix[i][0] = 1 + 2 * ratio;
            coef_matrix[i][1] = - ratio;
        }
    }

    // cout << coef_matrix << endl;
    // exit(0);

    // printf("ok1");

    // cout << "coef matrix:" << coef_n_cols << "," << coef_n_cols << endl; 
    // cout << "---------- coef mat initial ------------" << endl;
    // cout << "------- before csolve ---------" << endl;
    // cout << coef_matrix <<endl;
    // std::tie(coef_matrix, pivot_indexes) = lu_partition(coef_matrix);
    std::chrono::system_clock::time_point start, end;
    float elapsed = 0.;
    start = std::chrono::system_clock::now();

    coef_matrix = sparse_cholesky_factorize(coef_matrix);
    
    end = std::chrono::system_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0;
    choleskey_elapsed_time += elapsed;
    // cout << "------- after csolve ---------" << endl;
    // cout << coef_matrix <<endl;
    // cout << "pv = " << pivot_indexes << endl;
    // exit(0);
    int max_coef_row_index, add_index;
    // start time step
    start = std::chrono::system_clock::now();

    for (i = now_step; i <= now_step+calc_time_step; i++)
    {

        // cout << "############### forward #############" << endl;
        // cout << "imp start fsub i = " << i << endl;
        // cout << "end col = " << end_col << " start col " << start_col << endl;

        // cout << field[i-1] <<endl;
        // cout << "-------- prev field org----------" << endl;
        // cout << field[i-1] <<endl;
        field[i - 1][start_col] -= field[i - 1][start_col - 1] * coef_list[1];
        field[i - 1][end_col-1] -= field[i - 1][end_col] * coef_list[1];
        // start forward substitution
        // cout << "i = " << i << endl;
        // cout << "imp ok 1" << endl;
        
        // cout << field[i-1] <<endl;
        // cout << "-------- prev field after add ----------" << endl;
        // cout << field[i-1] <<endl;

        for (j = start_col; j < end_col; j++)
        {
            coef_j = j - start_col;
            field_j = j + f2c_idx;
            substitute_sum = 0.;
            max_coef_row_index = min(coef_j, coef_n_cols-1);
            // cout << "max coef row index:" << max_coef_row_index << endl;
            if (max_coef_row_index > 0){
                substitute_sum += coef_matrix[coef_j-1][1] * field[i][field_j-1];
            }
            field[i][field_j] = (field[i-1][field_j] - substitute_sum) / coef_matrix[coef_j][0];
            // cout << "j = " << j <<" coef_j " << coef_j <<" field j " << field_j << " sub sum: " << substitute_sum << " " << field[i][field_j]<< endl;

            // for (k = 0; k < max_coef_row_index; k++)
            // {
            //     tmp_field_id = k + add_index;
            //     // printf("k=%d coef=(%d, %d) tmp_y=%d tmp_x = %d\n", k, add_index+k, max_gauss_coef_index-k, tmp_field_yid, tmp_field_xid);

            //     substitute_sum += coef_matrix[tmp_field_id][max_gauss_coef_index-k] * field[i][tmp_field_yid][tmp_field_xid];
            // }
            // cout << "substitute sum " << substitute_sum << endl;
            // cout << "j = " << j << " pivot index after = " << pivot_indexes[coef_j] + f2c_idx << endl;
            // cout << "field : " << field[i - 1][field_j] << " pivot " << pivot_indexes[coef_j] << endl;
            // cout << "bunshi : " << field[i - 1][field_j] - substitute_sum << endl;
            // cout << "bunbo :" << coef_matrix[coef_j][coef_j] << endl;

            // devide by coef_matrix[coef_j][coef_j] only once, backward, xor forward.
        }

        // cout << "imp ok 2" << endl;

        field[i - 1][start_col] += field[i - 1][start_col - 1] * coef_list[1];
        field[i - 1][end_col-1] += field[i - 1][end_col] * coef_list[1];

        // cout << "-------- now field before bacjward ----------" << endl;
        // cout << field[i] <<endl;

        // cout << "imp ok 3" << endl;

        // cout << "############### backward #############" << endl;


        // start backward substitution
        for (j = end_col-1; j >= start_col; j--)
        {
            // cout << "j = " << j << endl;
            coef_j = j - start_col;
            field_j = j + f2c_idx;
            max_coef_row_index = min(coef_n_rows-coef_j-1, coef_n_cols-1);

            // cout << "j = " << j <<" coef_j " << coef_j << endl;

            substitute_sum = 0.;
            if (max_coef_row_index > 0){
                substitute_sum += coef_matrix[coef_j][1] * field[i][field_j+1];
            }
            // cout << "j = " << j <<" coef_j " << coef_j <<" field j " << field_j <<"max:" <<max_coef_row_index << " sub sum: " << substitute_sum << " " << field[i][field_j]<< endl;

            // cout << "substitute sum " << substitute_sum << endl;
            // cout << "j = " << j << " pivot index after = " << pivot_indexes[coef_j] + f2c_idx << endl;
            // cout << "j = " << j << endl;
            // cout << "field : " << field[i][j + f2c_idx] << endl;
            // cout << "bunshi : " << field[i][j + f2c_idx] - substitute_sum << endl;
            // cout << "bunbo :" << coef_matrix[coef_j][coef_j] << endl;
            field[i][field_j] = (field[i][field_j] - substitute_sum) / coef_matrix[coef_j][0];
        }
    }
    end = std::chrono::system_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0;
    gauss_time += elapsed;
    // cout << field[now_step] << endl;
    // exit(0);
    // }
}

void solve(map<string, pair<float, float>> spatio_bound, map<string, float> delta, vector<map<string, pair<float, float>>> domain_boundary)
{
    int now_step, step, j, k, tmp_col1, tmp_col2;
    int n_rows = field.size();
    int n_cols = field[0].size();
    map<string, map<string, float>> delta_info;
    float coef;

    // printf("okok\n");
    vector<map<string, pair<float, float>>> domain_col(domain_boundary.size());
    
    bool is_same = !is_sep;
    int is_exp = domain_boundary[0]["is_exp"].first;
    for (j = 0; j < domain_boundary.size(); j++)
    {
        map<string, pair<float, float>> tmp_domain_col;
        for (auto itr = domain_boundary[j].begin(); itr != domain_boundary[j].end(); itr++)
        {
            if (itr->first == "x"){
                tmp_col1 = (itr->second.first - spatio_bound[itr->first].first) / (delta[itr->first]);
                tmp_col2 = (itr->second.second - spatio_bound[itr->first].first) / (delta[itr->first]);
                tmp_domain_col[itr->first] = make_pair(tmp_col1, tmp_col2);
            }else if(itr->first == "is_exp"){
                is_same &= is_exp == (itr->second).first;
                tmp_domain_col[itr->first] = itr->second;
            }
            else{
                tmp_domain_col[itr->first] = itr->second;                
            }
        }
        domain_col[j] = tmp_domain_col;
    }

    // cout << domain_col << endl;
    // cout << delta << endl;

    delta_info["exp"] = delta;
    delta_info["imp"] = delta;

    // printf("okok1\n");

    std::chrono::system_clock::time_point  start, end;
    float elapsed = 0., coef_small, coef_large;
    pair<float, float> neighb_coef;
    vector<float> domain_elapsed(domain_col.size(), 0.);
    // vector<int> counter(domain_col.size(), 0.);
    // cout << "is sep:" << is_sep <<endl;
    // cout << "is same:" << is_same <<endl;
    for (now_step = 1; now_step < n_rows; now_step++)
    {
        // printf("now step %d\n", now_step);
        if (is_same){
            start = std::chrono::system_clock::now();
            if(is_exp){
                explicit_solve_whole(now_step, delta_info["exp"]["t"], delta_info["exp"], domain_col);
            }else{
                implicit_solve_whole(now_step, delta_info["exp"]["t"], delta_info["imp"], domain_col);
            }
            end = std::chrono::system_clock::now();
            elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end-start).count();
            domain_elapsed[0] = elapsed;
        }else{
            for (j = 0; j < domain_col.size(); j++)
            {
                // cout << "j="<<j<<endl;
                coef = domain_col[j]["coef"].first;
                if(j>0){
                    coef_small = domain_col[j-1]["coef"].first;
                }else{
                    coef_small = coef;
                }
                if(j<domain_col.size()-1){
                    coef_large = domain_col[j+1]["coef"].first;
                }else{
                    coef_large = coef;
                }
                neighb_coef = make_pair(coef_small, coef_large);
                if (domain_col[j]["is_exp"].first==1){
                    if(domain_col[j]["t_delta"].first > 0){
                        delta_info["exp"]["t"] = domain_col[j]["t_delta"].first;
                    }
                    start = std::chrono::system_clock::now();
                    explicit_solve(now_step, delta_info["exp"]["t"], delta_info["exp"], domain_col[j], coef, neighb_coef);
                    end = std::chrono::system_clock::now();
                    // elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end-start).count();
                }else{
                    if(domain_col[j]["t_delta"].first > 0){
                        delta_info["imp"]["t"] = domain_col[j]["t_delta"].first;
                    }
                    start = std::chrono::system_clock::now();
                    sparse_implicit_solve(now_step, delta_info["imp"]["t"], delta_info["imp"], domain_col[j], coef, neighb_coef);
                    end = std::chrono::system_clock::now();
                    elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end-start).count();
                }
                domain_elapsed[j] += elapsed / (n_rows-1);
            }
        }
        
        // exit(0);
        // printf(" finish !\n");
    }

    string tmp, tmp2;
    LOG_FILE.open(LOG_PATH, ios::app);
    LOG_FILE << "--- mean time | total (" << n_rows-1 << " step) time ---\n";
    for(j=0;j<domain_elapsed.size();j++){
        // tmp = float_to_string(domain_elapsed[j], 2);
        LOG_FILE << domain_elapsed[j] << "micro s | ";
        // tmp = float_to_string(domain_elapsed[j]*(n_rows-1), 2);
        LOG_FILE << domain_elapsed[j]*(n_rows-1) <<"micro s" << "\n";
    }
    LOG_FILE << "chol: " << choleskey_elapsed_time / (float)(n_rows - 1) << "ms | ";
    LOG_FILE << choleskey_elapsed_time << "ms"
             << "\n";
    LOG_FILE << "gauss: " << gauss_time / (float)(n_rows - 1) << "ms | ";
    LOG_FILE << gauss_time << "ms"
             << "\n";
    LOG_FILE.close();
}

void initialize(float initial_value, map<string, pair<float, float>> bound_value){
    int i,j;
    int n_rows = field.size();
    int n_cols = field[0].size();
    if(n_cols==0 & n_rows==0){
        printf("must initialize field vector first !\n");
        exit(1);
    }

    for(i=0;i<n_rows;i++){
        for(j=0;j<n_cols;j++){
            if(i==0 & j>0 & j<n_cols-1){
                field[0][j] = initial_value;
            }else if(j==0){
                field[i][j] = bound_value["x"].first;
            }else if(j==n_cols-1){
                field[i][j] = bound_value["x"].second;
            }
        }
        // cout << field[i] << endl;
        // cout << bound_value << endl;
    }

    // cout << "########## initial ###########" << endl;
    // cout << field << endl;
}

int pde(YAML::Node config, char *argv[], float type_value)
{
    cout << "start pde" << "\n";
    // 初期値、境界値（x=0, 1, t=0 の3つ）、tiem delta、x_delta を引数にすると良さそう
    int n_cols, n_rows,i,j, try_time;
    bool is_exact_solution = (config["is_exact"].as<int>()>0);
    try_time = config["try_time"].as<int>();

    try{
        is_sep = config["is_sep"].as<int>() == 1;
    }catch(YAML::TypedBadConversion<int>){
        is_sep = false;
    }
    string line = "-----------------------------";


    OUTPUT_FILE_PATH = argv[2];
    fs::path save_dir = OUTPUT_FILE_PATH.parent_path();
    if (!fs::exists(save_dir)) fs::create_directories(save_dir);
    cout << "save result to " << save_dir << "\n";

    LOG_PATH = OUTPUT_FILE_PATH.string<char>() + ".log";
    LOG_FILE.open(LOG_PATH, ios::out);
    LOG_FILE << "save result to directory: " << save_dir << "\n";
    LOG_FILE << "load config config file from " << argv[1] << "\n";
    LOG_FILE << line << "\n";
    LOG_FILE << config << "\n";
    LOG_FILE << line << "\n";
    LOG_FILE.close();

    float tmp;
    map<string, float> delta = config["delta"].as<map<string, float>>();

    auto spatio_bound = config["space"].as<map<string, pair<float,float>>>();
    auto t_bound = config["time"].as<pair<float, float>>();


    auto bound_value = config["condition"]["bound_value"].as<map<string, pair<float, float>>>();
    float initial_value = config["condition"]["initial_value"]["t"].as<float>();

    vector<map<string, pair<float, float>>> domain_list;
    cout << "solver num:" << config["solver"].size() << endl;
    for(std::size_t i=0;i<config["solver"].size();i++){
        auto tmp_domain_info = config["solver"][i]["boundary"].as<map<string, pair<float, float>>>();
        
        tmp = config["solver"][i]["coef"].as<float>();
        tmp_domain_info["coef"] = make_pair(tmp, tmp);
        
        tmp = (float)(config["solver"][i]["type"].as<string>()=="exp");
        tmp_domain_info["is_exp"] = make_pair(tmp, tmp);
        
        tmp = config["solver"][i]["t_delta"].as<float>();
        tmp_domain_info["t_delta"] = make_pair(tmp, tmp);
        
        domain_list.push_back(tmp_domain_info);
    }
    float x_length = spatio_bound["x"].second - spatio_bound["x"].first;
    float time_length = t_bound.second - t_bound.first;
    n_cols = (int)(x_length / delta["x"]) + 1;
    n_rows = (int)(time_length / delta["t"]) + 1;
    if(try_time==1){
        vector<vector<float>> tmp_vec(n_rows, vector<float>(n_cols, 0));
        field = tmp_vec;
        tmp_vec.clear();

        if (!is_exact_solution){
            cout << "initialize" << endl;
            initialize(initial_value, bound_value);
        }

        // vector<vector<float> > vec;
        printf("start solving ... \n");
        auto start = std::chrono::system_clock::now();
        if (is_exact_solution) exact_solution_1d(&field, spatio_bound, delta, domain_list[0]["coef"].first);
        else solve(spatio_bound, delta, domain_list);
        auto end = std::chrono::system_clock::now();
        double elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end-start).count(); //処理に要した時間をミリ秒に変換
        string e = "end in " + float_to_string(elapsed / 1000, 4) + "ms";
        cout << e << "\n";
        // vec = implicit_solve(time_length, t_delta, x_length, x_delta, initial_value, x_value);
        LOG_FILE.open(LOG_PATH, ios::app);
        LOG_FILE << e << "\n";
        LOG_FILE.close();
    }else{
        double total_elapsed = 0, elapsed;
        vector<double> elapsed_list;
        for(i=0; i<try_time;i++){
            vector<vector<float>> tmp_vec(n_rows, vector<float>(n_cols, 0));
            field = tmp_vec;
            tmp_vec.clear();

            if (!is_exact_solution){
                // cout << "initialize" << endl;
                initialize(initial_value, bound_value);
            }

            // vector<vector<float> > vec;
            // printf("start solving ... \n");
            auto start = std::chrono::system_clock::now();
            if (is_exact_solution) exact_solution_1d(&field, spatio_bound, delta, domain_list[0]["coef"].first);
            else solve(spatio_bound, delta, domain_list);
            auto end = std::chrono::system_clock::now();
            elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end-start).count(); //処理に要した時間をミリ秒に変換
            elapsed_list.push_back(elapsed);
            total_elapsed += elapsed / (double)try_time;
            if(i % 5 == 0){
                printf("finish %d / %d \n", i, try_time);
            }
        }
        string e = "end in " + float_to_string(total_elapsed / 1000 * try_time, 4) + "ms";
        cout << e << "\n";
        cout << "detail:" << "\n";
        LOG_FILE.open(LOG_PATH, ios::app);
        LOG_FILE << e << "\n";
        for(i=0;i<elapsed_list.size();i++){
            e = "trial " + to_string(i+1) + " : " + to_string(elapsed_list[i]/1000) + " ms";
            cout << e << "\n";
            LOG_FILE << e << "\n";
        }
        e = "mean " + to_string(i+1) + " : " + to_string(total_elapsed/1000) + " ms";
        cout << e << "\n";
        LOG_FILE << e << "\n";
        // vec = implicit_solve(time_length, t_delta, x_length, x_delta, initial_value, x_value);
        LOG_FILE.close();
    }


    

    OUTPUT_FILE.open(OUTPUT_FILE_PATH, ios::out);
    float t_now, x_now;

    cout << "writing " << OUTPUT_FILE_PATH << "..." << endl;

    for (i = 0; i < field.size(); i++)
    {
        t_now = t_bound.first + i * delta["t"];
        if(i==0){
            // printf("time,");
            OUTPUT_FILE << "time,";
            for(j=0;j < n_cols; j++){
                x_now = spatio_bound["x"].first + delta["x"] * j;
                // cout << fixed << setprecision(3) << x_now;
                OUTPUT_FILE << "x=" <<fixed <<setprecision(5)<< x_now;
                if (j<n_cols-1){
                    // cout << "\t";
                    OUTPUT_FILE << ",";
                }
            }
            OUTPUT_FILE << endl;
            // cout << endl;
        }
        // printf("t=%f\t", t_now);
        OUTPUT_FILE << fixed <<setprecision(5) << t_now << ",";

        for (j = 0; j < field[0].size(); j++)
        {
            if(j < field[0].size()-1){
                // printf("%f\t", vec[i][j]);
                OUTPUT_FILE << field[i][j] << ",";
            }else{
                // printf("%f", vec[i][j]);
                OUTPUT_FILE << field[i][j];
            }
            
        }
        OUTPUT_FILE << endl;
        // printf("\n");
    }

    OUTPUT_FILE.close();
    return 0;
}
