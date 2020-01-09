#include <iostream>
#include <vector>
#include <cmath>
#include <prettyprint.hpp>
#include "numerical.hpp"
#include <fstream> // file io
#include <string>
#include <sstream> // digit adjust
#include <iomanip> //exit?
#include <map>
#include <chrono>     // time calculation
#include <filesystem> // make directory
#include "yaml-cpp/yaml.h"
using namespace std;
namespace fs = filesystem;
// namespace plt = matplotlibcpp;

// template <typename float>
// int pde_2d(YAML::Node config, char* argv[], float type_value);
static vector<vector<vector<float>>> field;
float choleskey_elapsed_time = 0, gauss_time = 0;
bool is_sparse_lu = true;

inline void preprocess_marginal_field(int time_index, int start_xid, int start_yid, int end_xid, int end_yid, vector<float> coef_list)
{
    int x_idx, y_idx;
    for (x_idx = start_xid; x_idx < end_xid; x_idx++)
    {
        if (x_idx == start_xid)
        {
            field[time_index - 1][start_yid][start_xid] -= field[time_index - 1][start_yid][start_xid - 1] * coef_list[1];
            field[time_index - 1][start_yid][start_xid] -= field[time_index - 1][start_yid - 1][start_xid] * coef_list[0];
            field[time_index - 1][end_yid - 1][start_xid] -= field[time_index - 1][end_yid - 1][start_xid - 1] * coef_list[1];
            field[time_index - 1][end_yid - 1][start_xid] -= field[time_index - 1][end_yid][start_xid] * coef_list[4];
        }
        else if (x_idx == end_xid - 1)
        {
            field[time_index - 1][start_yid][x_idx] -= field[time_index - 1][start_yid][x_idx + 1] * coef_list[3];
            field[time_index - 1][start_yid][x_idx] -= field[time_index - 1][start_yid - 1][x_idx] * coef_list[0];
            field[time_index - 1][end_yid - 1][x_idx] -= field[time_index - 1][end_yid - 1][x_idx + 1] * coef_list[3];
            field[time_index - 1][end_yid - 1][x_idx] -= field[time_index - 1][end_yid][x_idx] * coef_list[4];
        }
        else
        {
            field[time_index - 1][start_yid][x_idx] -= field[time_index - 1][start_yid - 1][x_idx] * coef_list[0];
            field[time_index - 1][end_yid - 1][x_idx] -= field[time_index - 1][end_yid][x_idx] * coef_list[4];
        }
    }
    for (y_idx = start_yid + 1; y_idx < end_yid - 1; y_idx++)
    {
        field[time_index - 1][y_idx][start_xid] -= field[time_index - 1][y_idx][start_xid - 1] * coef_list[1];
        field[time_index - 1][y_idx][end_xid - 1] -= field[time_index - 1][y_idx][end_xid] * coef_list[3];
    }
}

inline void postprocess_marginal_field(int time_index, int start_xid, int start_yid, int end_xid, int end_yid, vector<float> coef_list)
{
    int x_idx, y_idx;
    for (x_idx = start_xid; x_idx < end_xid; x_idx++)
    {
        if (x_idx == start_xid)
        {
            field[time_index - 1][start_yid][start_xid] += field[time_index - 1][start_yid][start_xid - 1] * coef_list[1];
            field[time_index - 1][start_yid][start_xid] += field[time_index - 1][start_yid - 1][start_xid] * coef_list[0];
            field[time_index - 1][end_yid - 1][start_xid] += field[time_index - 1][end_yid - 1][start_xid - 1] * coef_list[1];
            field[time_index - 1][end_yid - 1][start_xid] += field[time_index - 1][end_yid][start_xid] * coef_list[4];
        }
        else if (x_idx == end_xid - 1)
        {
            field[time_index - 1][start_yid][x_idx] += field[time_index - 1][start_yid][x_idx + 1] * coef_list[3];
            field[time_index - 1][start_yid][x_idx] += field[time_index - 1][start_yid - 1][x_idx] * coef_list[0];
            field[time_index - 1][end_yid - 1][x_idx] += field[time_index - 1][end_yid - 1][x_idx + 1] * coef_list[3];
            field[time_index - 1][end_yid - 1][x_idx] += field[time_index - 1][end_yid][x_idx] * coef_list[4];
        }
        else
        {
            field[time_index - 1][start_yid][x_idx] += field[time_index - 1][start_yid - 1][x_idx] * coef_list[0];
            field[time_index - 1][end_yid - 1][x_idx] += field[time_index - 1][end_yid][x_idx] * coef_list[4];
        }
    }
    for (y_idx = start_yid + 1; y_idx < end_yid - 1; y_idx++)
    {
        field[time_index - 1][y_idx][start_xid] += field[time_index - 1][y_idx][start_xid - 1] * coef_list[1];
        field[time_index - 1][y_idx][end_xid - 1] += field[time_index - 1][y_idx][end_xid] * coef_list[3];
    }
}

void explicit_solve_2d(int now_step, float calc_time_step, map<string, float> delta, map<string, pair<float, float>> domain_col, float coef = 1.)
{
    int t, i, j;
    // int calc_step = int(calc_time_step / delta["t"]);
    int calc_step = 1;
    int start_xid, end_xid, start_yid, end_yid;
    int n_cols1 = field[0].size();
    float x_diff, ratio_x, ratio_y, y_diff;
    if (domain_col["x"].first == 0)
        start_xid = 1;
    else
        start_xid = domain_col["x"].first;
    if (domain_col["y"].first == 0)
        start_yid = 1;
    else
        start_yid = domain_col["y"].first;
    // cout << "ok exp1" << endl;
    // cout << "---------- exp -----------" << endl;
    // cout << "start y:" << domain_col["y"] << endl;
    // cout << "start x:" << domain_col["x"] << endl;
    // cout << "calc step:" << calc_step << " calc time step: " << calc_time_step << endl;

    ratio_x = coef * delta["t"] / (delta["x"] * delta["x"]);
    ratio_y = coef * delta["t"] / (delta["y"] * delta["y"]);
    // cout << "ratio x, y = " << ratio_x << "," << ratio_y << "\n";
    // exit(0);
    // cout << "explicit" << endl; 
    // cout << domain_col << endl;
    // initialize
    // printf("%lx %lx\n", (long)field, (long)&field[0]);
    // int count = 0;
    // solve
    // cout << "field shape:" << field.size() << ","<<field[0].size() << "," <<field[0][0].size() << " now step:" << now_step << endl;
    // for(int i=0;i<n_cols1;i++){
    // cout << "y=" << i << " " << field[now_step][i] << endl;
    // }
    for (t = now_step; t < calc_step + now_step; t++)
    {
        // cout << "ok exp2 " << now_step <<" i=  " << i<< endl;
        for (i = start_yid; i < domain_col["y"].second; i++)
        {
            for (j = start_xid; j < domain_col["x"].second; j++)
            {
                x_diff = field[t - 1][i][j + 1] + field[t - 1][i][j - 1] - 2 * field[t - 1][i][j];
                y_diff = field[t - 1][i + 1][j] + field[t - 1][i - 1][j] - 2 * field[t - 1][i][j];
                field[t][i][j] = field[t - 1][i][j] + ratio_x * x_diff + ratio_y * y_diff;
                // count ++;
            }
        }
        // cout << "t:" << t << " count:" <<count << endl;
    }

    // cout << "end exp3" << endl;
}

map<string, pair<int, int>> calc_column_status(map<string, pair<float, float>> domain_col){
    int start_xid, start_yid, end_xid, end_yid;
    map<string, pair<int, int>> result;
    if (domain_col["x"].first == 0)
        start_xid = 1;
    else
        start_xid = domain_col["x"].first;
    if (domain_col["y"].first == 0)
        start_yid = 1;
    else
        start_yid = domain_col["y"].first;

    end_xid = domain_col["x"].second;
    end_yid = domain_col["y"].second;

    result["x"] = make_pair(start_xid, end_xid);
    result["y"] = make_pair(start_yid, end_yid);


    return result;
}

void implicit_solve_2d(int now_step, float calc_time_step, map<string, float> delta, map<string, pair<float, float>> domain_col, float coef = 1.)
{
    int t, i, j, k, index, x_idx, y_idx, coef_idx, field_xid, field_yid, coef_k;
    int calc_step = int(calc_time_step / delta["t"]);
    int start_xid, end_xid, start_yid, end_yid, tmp_field_xid, tmp_field_yid;
    float x_diff, ratio_x, ratio_y, y_diff;
    map<string, pair<int, int>> col_status;

    col_status = calc_column_status(domain_col);
    start_xid = col_status["x"].first;
    end_xid = col_status["x"].second;

    start_yid = col_status["y"].first;
    end_yid = col_status["y"].second;


    ratio_x = coef * delta["t"] / (delta["x"] * delta["x"]);
    ratio_y = coef * delta["t"] / (delta["y"] * delta["y"]);

    int f2c_idx = 0;
    int x_idx_num, y_idx_num, coef_n_cols;
    x_idx_num = end_xid - start_xid;
    y_idx_num = end_yid - start_yid;
    coef_n_cols = x_idx_num * y_idx_num;

    float substitute_sum;

    vector<float> coef_list{-ratio_y, -ratio_x, 1 + 2 * ratio_x + 2 * ratio_y, -ratio_x, -ratio_y};
    vector<vector<float>> coef_matrix(coef_n_cols, vector<float>(coef_n_cols, 0));
    vector<int> index_list{-x_idx_num, -1, 0, 1, x_idx_num};

    for (i = 0; i < coef_n_cols; i++)
    {
        for (j = 0; j < index_list.size(); j++)
        {
            index = i + index_list[j];
            if (index < 0 || index >= coef_n_cols)
                continue;
            if ((index / x_idx_num != i / x_idx_num) && (index % x_idx_num != i % x_idx_num))
            {
                // 例えば(2,0) と (1, x_idx_num-1) が関連づけられちゃうのを防ぐ。
                continue;
            }
            coef_matrix[i][index] = coef_list[j];
        }
    }
    // cout << "------- before csolve ---------" << endl;
    // for (i = 0; i < coef_n_cols; i++)
    // {
    //     cout << "y=" << i << coef_matrix[i] << endl;
    // }
    std::chrono::system_clock::time_point start, end;
    float elapsed = 0.;
    start = std::chrono::system_clock::now();
    coef_matrix = cholesky_factorize(coef_matrix);
    end = std::chrono::system_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0;
    choleskey_elapsed_time += elapsed;
    // cout << "------- after csolve ---------" << endl;
    // for (i = 0; i < coef_n_cols; i++)
    // {
    //     cout << "y=" << i << coef_matrix[i] <<setprecision(6) << endl;
    // }

    // start time step
    start = std::chrono::system_clock::now();
    for (i = now_step; i <= now_step + calc_time_step; i++)
    {

        preprocess_marginal_field(i, start_xid, start_yid, end_xid, end_yid, coef_list);
        for (y_idx = start_yid; y_idx < end_yid; y_idx++)
        {
            for (x_idx = start_xid; x_idx < end_xid; x_idx++)
            {
                coef_idx = (y_idx - start_yid) * x_idx_num + (x_idx - start_xid);
                substitute_sum = 0.;
                for (k = 0; k < coef_idx; k++)
                {
                    tmp_field_xid = k % x_idx_num + start_xid;
                    tmp_field_yid = k / x_idx_num + start_yid;
                    substitute_sum += coef_matrix[coef_idx][k] * field[i][tmp_field_yid][tmp_field_xid];
                }
                field[i][y_idx][x_idx] = (field[i - 1][y_idx][x_idx] - substitute_sum) / coef_matrix[coef_idx][coef_idx];
            }
        }

        postprocess_marginal_field(i, start_xid, start_yid, end_xid, end_yid, coef_list);
        // printf("---------- after forward ----------\n");
        // for (int jj = 0; jj < field[now_step].size(); jj++)
        // {
        //     cout << "y=" << jj << " " << field[i][jj] << endl;
        // }

        // start backward substitution
        for (y_idx = end_yid - 1; y_idx >= start_yid; y_idx--)
        {
            for (x_idx = end_xid - 1; x_idx >= start_xid; x_idx--)
            {
                coef_idx = (y_idx - start_yid) * x_idx_num + (x_idx - start_xid);
                substitute_sum = 0.;
                for (k = coef_n_cols - 1; k > coef_idx; k--)
                {
                    tmp_field_xid = k % x_idx_num + start_xid;
                    tmp_field_yid = k / x_idx_num + start_yid;
                    substitute_sum += coef_matrix[k][coef_idx] * field[i][tmp_field_yid][tmp_field_xid];
                }
                field[i][y_idx][x_idx] = (field[i][y_idx][x_idx] - substitute_sum) / coef_matrix[coef_idx][coef_idx];
            }
        }
        // printf("---------- after backward ------------ \n");
        // for (int jj = 0; jj < field[now_step].size(); jj++)
        // {
        //     cout << "y=" << jj << " " << field[i][jj] << endl;
        // }


    }
    end = std::chrono::system_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0;
    gauss_time += elapsed;
}

void sparse_symmetric_implicit_solve_2d_xy(int now_step, float calc_time_step, map<string, float> delta, map<string, pair<float, float>> domain_col, float coef = 1.)
{
    int t, i, j, k, index, x_idx, y_idx, coef_idx, field_xid, field_yid, coef_k, dense_index, max_gauss_coef_index, add_index, coef_row_base;
    int calc_step = int(calc_time_step / delta["t"]);
    int start_xid, end_xid, start_yid, end_yid, tmp_field_xid, tmp_field_yid;
    float x_diff, ratio_x, ratio_y, y_diff;
    map<string, pair<int, int>> col_status;
    col_status = calc_column_status(domain_col);
    start_xid = col_status["x"].first;
    end_xid = col_status["x"].second;

    start_yid = col_status["y"].first;
    end_yid = col_status["y"].second;

    ratio_x = coef * delta["t"] / (delta["x"] * delta["x"]);
    ratio_y = coef * delta["t"] / (delta["y"] * delta["y"]);

    int f2c_idx = 0;
    int x_idx_num, y_idx_num, coef_n_rows, coef_n_cols;
    x_idx_num = end_xid - start_xid;
    y_idx_num = end_yid - start_yid;
    if (y_idx_num < x_idx_num)
    {
        printf("the number of x columns must be smaller than the number oh y columns!\n");
        exit(0);
    }
    coef_n_rows = x_idx_num * y_idx_num;
    coef_n_cols = x_idx_num + 1;

    float substitute_sum;

    vector<float> coef_list{-ratio_y, -ratio_x, 1 + 2 * ratio_x + 2 * ratio_y, -ratio_x, -ratio_y};
    vector<float> coef_sparse_list{1 + 2 * ratio_x + 2 * ratio_y, -ratio_x, -ratio_y};
    vector<vector<float>> coef_matrix(coef_n_rows, vector<float>(coef_n_cols, 0));
    vector<int> index_list{0, 1, x_idx_num};

    for (i = 0; i < coef_n_rows; i++)
    {
        for (j = 0; j < index_list.size(); j++)
        {
            index = index_list[j];
            dense_index = index + i;
            if (dense_index >= coef_n_rows)
            {
                continue;
            }
            if(i / x_idx_num != dense_index / x_idx_num && i % x_idx_num != dense_index % x_idx_num){
                // prevent making relation between (2,0) and (1, x_idx_num-1)
                // printf("i=%d, j=%d\n", i,dense_index);
                continue;
            }
            coef_matrix[i][index] = coef_sparse_list[j];
        }
    }

    // cout << "------- before csolve ---------" << endl;
    // for (i = 0; i < coef_n_rows; i++)
    // {
    //     cout << "y=" << i << coef_matrix[i] << endl;
    // }
    // exit(0);

    std::chrono::system_clock::time_point start, end;
    float elapsed = 0.;
    start = std::chrono::system_clock::now();

    coef_matrix = sparse_cholesky_factorize(coef_matrix);
    // cout << "------- after csolve ---------" << endl;
    // for (i = 0; i < coef_n_rows; i++)
    // {
    //     cout << "y=" << i << coef_matrix[i] << endl;
    // }
    end = std::chrono::system_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0;
    choleskey_elapsed_time += elapsed;

    // start time step
    start = std::chrono::system_clock::now();
    for (i = now_step; i <= now_step + calc_time_step; i++)
    {

        preprocess_marginal_field(i, start_xid, start_yid, end_xid, end_yid, coef_list);

        for (y_idx = start_yid; y_idx < end_yid; y_idx++)
        {
            for (x_idx = start_xid; x_idx < end_xid; x_idx++)
            {
                coef_idx = (y_idx - start_yid) * x_idx_num + (x_idx - start_xid);
                substitute_sum = 0.;
                max_gauss_coef_index = min(coef_idx, coef_n_cols-1);
                add_index = max(0, coef_idx - max_gauss_coef_index);
                // printf("-------------\n");
                // printf("y=%d x=%d coef_idx= %d max_idx=%d \n", y_idx, x_idx, coef_idx, max_gauss_coef_index);
                for (k = 0; k < max_gauss_coef_index; k++)
                {
                    tmp_field_xid = (k + add_index) % x_idx_num + start_xid;
                    tmp_field_yid = (k + add_index) / x_idx_num + start_yid;
                    // printf("k=%d coef=(%d, %d) tmp_y=%d tmp_x = %d\n", k, add_index+k, max_gauss_coef_index-k, tmp_field_yid, tmp_field_xid);

                    substitute_sum += coef_matrix[add_index+k][max_gauss_coef_index-k] * field[i][tmp_field_yid][tmp_field_xid];
                }
                field[i][y_idx][x_idx] = (field[i - 1][y_idx][x_idx] - substitute_sum) / coef_matrix[coef_idx][0];
            }
        }

        postprocess_marginal_field(i, start_xid, start_yid, end_xid, end_yid, coef_list);

        // printf("---------- after forward ------------ \n");
        // for (int jj = 0; jj < field[now_step].size(); jj++)
        // {
        //     cout << "y=" << jj << " " << field[i][jj] << endl;
        // }



        // start backward substitution
        for (y_idx = end_yid - 1; y_idx >= start_yid; y_idx--)
        {
            for (x_idx = end_xid - 1; x_idx >= start_xid; x_idx--)
            {
                coef_idx = (y_idx - start_yid) * x_idx_num + (x_idx - start_xid);
                substitute_sum = 0.;
                max_gauss_coef_index = min(coef_n_rows-coef_idx, coef_n_cols);
                // printf("y=%d x=%d coef_idx= %d max_idx=%d\n", y_idx, x_idx, coef_idx, max_gauss_coef_index);
                for (k = max_gauss_coef_index - 1; k > 0; k--)
                {
                    tmp_field_xid = (k + coef_idx) % x_idx_num + start_xid;
                    tmp_field_yid = (k + coef_idx) / x_idx_num + start_yid;
                    substitute_sum += coef_matrix[coef_idx][k] * field[i][tmp_field_yid][tmp_field_xid];
                }
                field[i][y_idx][x_idx] = (field[i][y_idx][x_idx] - substitute_sum) / coef_matrix[coef_idx][0];
            }
        }
        // printf("---------- after backward ------------ \n");
        // for (int jj = 0; jj < field[now_step].size(); jj++)
        // {
        //     cout << "y=" << jj << " " << field[i][jj] << endl;
    }
    end = std::chrono::system_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0;
    gauss_time += elapsed;
}

void sparse_symmetric_implicit_solve_2d_yx(int now_step, float calc_time_step, map<string, float> delta, map<string, pair<float, float>> domain_col, float coef = 1.)
{
    int t, i, j, k, index, x_idx, y_idx, coef_idx, field_xid, field_yid, coef_k, dense_index, max_gauss_coef_index, add_index, coef_row_base;
    int calc_step = int(calc_time_step / delta["t"]);
    int start_xid, end_xid, start_yid, end_yid, tmp_field_xid, tmp_field_yid;
    float x_diff, ratio_x, ratio_y, y_diff;
    map<string, pair<int, int>> col_status;
    col_status = calc_column_status(domain_col);
    start_xid = col_status["x"].first;
    end_xid = col_status["x"].second;

    start_yid = col_status["y"].first;
    end_yid = col_status["y"].second;

    ratio_x = coef * delta["t"] / (delta["x"] * delta["x"]);
    ratio_y = coef * delta["t"] / (delta["y"] * delta["y"]);

    int f2c_idx = 0;
    int x_idx_num, y_idx_num, coef_n_rows, coef_n_cols;
    x_idx_num = end_xid - start_xid;
    y_idx_num = end_yid - start_yid;
    if (x_idx_num < y_idx_num)
    {
        printf("the number of y columns must be smaller than the number oh x columns!\n");
        exit(0);
    }
    coef_n_rows = x_idx_num * y_idx_num;
    coef_n_cols = y_idx_num + 1;

    float substitute_sum;

    vector<float> coef_list{-ratio_y, -ratio_x, 1 + 2 * ratio_x + 2 * ratio_y, -ratio_x, -ratio_y};
    vector<float> coef_sparse_list{1 + 2 * ratio_x + 2 * ratio_y, -ratio_y, -ratio_x};
    vector<vector<float>> coef_matrix(coef_n_rows, vector<float>(coef_n_cols, 0));
    vector<int> index_list{0, 1, y_idx_num};

    for (i = 0; i < coef_n_rows; i++)
    {
        for (j = 0; j < index_list.size(); j++)
        {
            index = index_list[j];
            dense_index = index + i;
            if (dense_index >= coef_n_rows)
            {
                continue;
            }
            if(i / y_idx_num != dense_index / y_idx_num && i % y_idx_num != dense_index % y_idx_num){
                // prevent making relation between (2,0) and (1, x_idx_num-1)
                // printf("i=%d, j=%d\n", i,dense_index);
                continue;
            }
            coef_matrix[i][index] = coef_sparse_list[j];
        }
    }

    // cout << "------- before csolve ---------" << endl;
    // for (i = 0; i < coef_n_rows; i++)
    // {
    //     cout << "y=" << i << coef_matrix[i] << endl;
    // }
    // exit(0);

    std::chrono::system_clock::time_point start, end;
    float elapsed = 0.;
    start = std::chrono::system_clock::now();

    coef_matrix = sparse_cholesky_factorize(coef_matrix);
    // cout << "------- after csolve ---------" << endl;
    // for (i = 0; i < coef_n_rows; i++)
    // {
    //     cout << "y=" << i << coef_matrix[i] << endl;
    // }
    end = std::chrono::system_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0;
    choleskey_elapsed_time += elapsed;

    // start time step
    start = std::chrono::system_clock::now();
    for (i = now_step; i <= now_step + calc_time_step; i++)
    {

        preprocess_marginal_field(i, start_xid, start_yid, end_xid, end_yid, coef_list);

        for (x_idx = start_xid; x_idx < end_xid; x_idx++)
        {
            for (y_idx = start_yid; y_idx < end_yid; y_idx++)
            {
                coef_idx = (x_idx - start_xid) * y_idx_num + (y_idx - start_yid);
                substitute_sum = 0.;
                max_gauss_coef_index = min(coef_idx, coef_n_cols-1);
                add_index = max(0, coef_idx - max_gauss_coef_index);
                // printf("-------------\n");
                // printf("y=%d x=%d coef_idx= %d max_idx=%d \n", y_idx, x_idx, coef_idx, max_gauss_coef_index);
                for (k = 0; k < max_gauss_coef_index; k++)
                {
                    tmp_field_yid = (k + add_index) % y_idx_num + start_yid;
                    tmp_field_xid = (k + add_index) / y_idx_num + start_xid;
                    // printf("k=%d coef=(%d, %d) tmp_y=%d tmp_x = %d\n", k, add_index+k, max_gauss_coef_index-k, tmp_field_yid, tmp_field_xid);

                    substitute_sum += coef_matrix[add_index+k][max_gauss_coef_index-k] * field[i][tmp_field_yid][tmp_field_xid];
                }
                field[i][y_idx][x_idx] = (field[i - 1][y_idx][x_idx] - substitute_sum) / coef_matrix[coef_idx][0];
            }
        }

        postprocess_marginal_field(i, start_xid, start_yid, end_xid, end_yid, coef_list);

        // printf("---------- after forward ------------ \n");
        // for (int jj = 0; jj < field[now_step].size(); jj++)
        // {
        //     cout << "y=" << jj << " " << field[i][jj] << endl;
        // }



        // start backward substitution
        for (x_idx = end_xid - 1; x_idx >= start_xid; x_idx--)
        {
            for (y_idx = end_yid - 1; y_idx >= start_yid; y_idx--)
            {
                coef_idx = (x_idx - start_xid) * y_idx_num + (y_idx - start_yid);
                substitute_sum = 0.;
                max_gauss_coef_index = min(coef_n_rows-coef_idx, coef_n_cols);
                // printf("y=%d x=%d coef_idx= %d max_idx=%d\n", y_idx, x_idx, coef_idx, max_gauss_coef_index);
                for (k = max_gauss_coef_index - 1; k > 0; k--)
                {
                    tmp_field_yid = (k + coef_idx) % y_idx_num + start_yid;
                    tmp_field_xid = (k + coef_idx) / y_idx_num + start_xid;
                    substitute_sum += coef_matrix[coef_idx][k] * field[i][tmp_field_yid][tmp_field_xid];
                }
                field[i][y_idx][x_idx] = (field[i][y_idx][x_idx] - substitute_sum) / coef_matrix[coef_idx][0];
            }
        }
        // printf("---------- after backward ------------ \n");
        // for (int jj = 0; jj < field[now_step].size(); jj++)
        // {
        //     cout << "y=" << jj << " " << field[i][jj] << endl;
        // }


    }
    end = std::chrono::system_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0;
    gauss_time += elapsed;
}

vector<vector<vector<float>>> solve_2d(map<string, pair<float, float>> spatio_bound, map<string, float> delta, vector<map<string, pair<float, float>>> domain_boundary)
{
    int now_step, step, j, k, tmp_col1, tmp_col2;
    float coef;
    int n_rows = (field).size();
    int n_cols1 = (field)[0].size();
    int n_cols2 = (field)[0][0].size();
    float ini = 0;
    map<string, map<string, float>> delta_info;
    map<string, pair<int, int>> col_status;

    vector<map<string, pair<float, float>>> domain_col(domain_boundary.size());
    for (j = 0; j < domain_boundary.size(); j++)
    {
        map<string, pair<float, float>> tmp_domain_col;
        for (auto itr = domain_boundary[j].begin(); itr != domain_boundary[j].end(); itr++)
        {
            if (itr->first == "x" || itr->first == "y")
            {
                tmp_col1 = (float)((itr->second.first - spatio_bound[itr->first].first) / (delta[itr->first]));
                tmp_col2 = (float)((itr->second.second - spatio_bound[itr->first].first) / (delta[itr->first]));
                tmp_domain_col[itr->first] = make_pair(tmp_col1, tmp_col2);
            }
            else
            {
                tmp_domain_col[itr->first] = itr->second;
            }
        }
        domain_col[j] = tmp_domain_col;
    }

    delta_info["exp"] = delta;
    delta_info["imp"] = delta;
    // cout << "-------------- start ------------" <<endl;
    // for(int i=0;i<n_cols1;i++){
    // cout << "y=" << i << " " << field[0][i] << endl;
    // }
    // cout << "ok1" << endl;
    std::chrono::system_clock::time_point start, end;
    float elapsed = 0.;
    vector<float> domain_elapsed(domain_col.size(), 0.);
    // vector<int> counter(domain_col.size(), 0.);
    for (now_step = 1; now_step < n_rows; now_step++)
    {
        // cout << "-------------- start step = " << now_step << " ------------" <<endl;
        // for(int i=0;i<n_cols1;i++){
        // cout << "y=" << i << " " << field[now_step][i] << endl;
        // }
        for (j = 0; j < domain_col.size(); j++)
        {
            coef = domain_col[j]["coef"].first;
            // cout << "coef = " << coef << " now step:" << now_step <<"\n";
            // cout << domain_col[j] << "\n";
            if (domain_col[j]["is_exp"].first == 1)
            {
                if (domain_col[j]["t_delta"].first > 0)
                {
                    delta_info["exp"]["t"] = domain_col[j]["t_delta"].first;
                }
                start = std::chrono::system_clock::now();
                explicit_solve_2d(now_step, delta_info["exp"]["t"], delta_info["exp"], domain_col[j], coef);
                end = std::chrono::system_clock::now();
                elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0;
            }
            else
            {
                if (domain_col[j]["t_delta"].first > 0)
                {
                    delta_info["imp"]["t"] = domain_col[j]["t_delta"].first;
                }

                start = std::chrono::system_clock::now();
                // cout << "imp" << endl;
                if (is_sparse_lu)
                {
                    col_status = calc_column_status(domain_col[j]);
                    if (col_status["x"].second - col_status["x"].first <= col_status["y"].second - col_status["y"].first){
                        sparse_symmetric_implicit_solve_2d_xy(now_step, delta_info["imp"]["t"], delta_info["imp"], domain_col[j], coef);
                    }else{
                        sparse_symmetric_implicit_solve_2d_yx(now_step, delta_info["imp"]["t"], delta_info["imp"], domain_col[j], coef);
                    }
                }
                else
                {
                    implicit_solve_2d(now_step, delta_info["imp"]["t"], delta_info["imp"], domain_col[j], coef);
                }
                end = std::chrono::system_clock::now();
                elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0;
                // cout << "-------------- "<< j << " imp ------------" <<endl;
                // for(int i=0;i<n_cols1;i++){
                // cout << "y=" << i << " " << field[now_step][i] << endl;
                // }
            }
            // cout << "elapsed:" << elapsed << endl;
            domain_elapsed[j] += elapsed / (n_rows - 1);
        }

        // if (now_step>50){
        // exit(1);
        // }

        // cout << "-------------- all ------------" <<endl;
        // for (int i = 0; i < n_cols1; i++)
        // {
        //     cout << "y=" << i << " " << field[now_step][i] << endl;
        // }
        // if (now_step > 10)
        // {
        //     exit(2);
        // }
    }

    string tmp, tmp2;
    LOG_FILE.open(LOG_PATH, ios::app);
    LOG_FILE << "--- mean time | total (" << n_rows - 1 << " step) time ---\n";
    for (j = 0; j < domain_elapsed.size(); j++)
    {
        tmp = float_to_string(domain_elapsed[j], 2);
        LOG_FILE << tmp << "ms | ";
        tmp = float_to_string(domain_elapsed[j] * (n_rows - 1), 2);
        LOG_FILE << tmp << "ms"
                 << "\n";
    }
    tmp = float_to_string(choleskey_elapsed_time / (float)(n_rows - 1), 2);
    LOG_FILE << "chol: " << tmp << "ms | ";
    tmp = float_to_string(choleskey_elapsed_time, 2);
    LOG_FILE << tmp << "ms"
             << "\n";
    tmp = float_to_string(gauss_time / (float)(n_rows - 1), 2);
    LOG_FILE << "gauss: " << tmp << "ms | ";
    tmp = float_to_string(gauss_time, 2);
    LOG_FILE << tmp << "ms"
             << "\n";
    LOG_FILE.close();
    return field;
}

void initialize_2d(float initial_value, map<string, pair<float, float>> bound_value)
{
    int i, j, k;
    int n_rows = field.size();
    int n_cols = field[0].size();
    int n_cols2 = field[0][0].size();

    for (i = 0; i < n_rows; i++)
    {
        for (j = 0; j < n_cols; j++)
        {
            for (k = 0; k < n_cols2; k++)
            {
                if (i == 0 & j > 0 & j<n_cols - 1 & k> 0 & k < n_cols2 - 1)
                {
                    field[0][j][k] = initial_value;
                }
                else if (k == 0)
                {
                    field[i][j][k] = bound_value["x"].first;
                }
                else if (k == n_cols2 - 1)
                {
                    field[i][j][k] = bound_value["x"].second;
                }
                else if (j == 0)
                {
                    field[i][j][k] = bound_value["y"].first;
                }
                else if (j == n_cols - 1)
                {
                    field[i][j][k] = bound_value["y"].second;
                }
            }
        }
    }

    // cout << "########## initial ###########" << endl;
    // cout << field << endl;

    // return field;
}

int pde_2d(YAML::Node config, char *argv[], float type_value)
{
    // if (argc < 2){
    //     printf("must specify config file in command line argument!\n");
    //     exit(1);
    // }
    int i, j, k, n_rows, n_cols_y, n_cols_x;
    bool is_exact_solution;

    // OUTPUT_FILE.open(OUTPUT_FILE_PATH, ios::out);
    // OUTPUT_FILE.close();

    // printf("load confing file from %s\n", argv[1]);
    // cout << line <<"\n";
    // YAML::Node config = YAML::LoadFile(argv[1]);
    // cout << config << endl;
    // cout << line <<"\n";
    string line = "-----------------------------";

    is_exact_solution = (config["is_exact"].as<int>() > 0);

    OUTPUT_FILE_PATH = argv[2];
    fs::path save_dir = OUTPUT_FILE_PATH.parent_path();
    if (!fs::exists(save_dir))
        fs::create_directories(save_dir);
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

    auto spatio_bound = config["space"].as<map<string, pair<float, float>>>();
    auto t_bound = config["time"].as<pair<float, float>>();

    auto bound_value = config["condition"]["bound_value"].as<map<string, pair<float, float>>>();
    float initial_value = config["condition"]["initial_value"]["t"].as<float>();

    vector<map<string, pair<float, float>>> domain_list;
    for (std::size_t i = 0; i < config["solver"].size(); i++)
    {
        auto tmp_domain_info = config["solver"][i]["boundary"].as<map<string, pair<float, float>>>();

        tmp = config["solver"][i]["coef"].as<float>();
        tmp_domain_info["coef"] = make_pair(tmp, tmp);

        tmp = (float)(config["solver"][i]["type"].as<string>() == "exp");
        tmp_domain_info["is_exp"] = make_pair(tmp, tmp);

        tmp = config["solver"][i]["t_delta"].as<float>();
        tmp_domain_info["t_delta"] = make_pair(tmp, tmp);

        domain_list.push_back(tmp_domain_info);
    }

    is_sparse_lu = config["is_sparse_lu"].as<bool>();

    float x_length = spatio_bound["x"].second - spatio_bound["x"].first;
    float y_length = spatio_bound["y"].second - spatio_bound["y"].first;
    float time_length = t_bound.second - t_bound.first;

    n_cols_x = (int)(x_length / delta["x"]) + 1;
    n_cols_y = (int)(y_length / delta["y"]) + 1;
    n_rows = (int)(time_length / delta["t"]) + 1;

    map<int, int> shape;
    shape[0] = n_rows;
    shape[1] = n_cols_y;
    shape[2] = n_cols_x;

    vector<vector<vector<float>>> tmp_vec(n_rows, vector<vector<float>>(n_cols_y, vector<float>(n_cols_x, 0)));
    field = tmp_vec;
    tmp_vec.clear();

    // cout << "ok pde 2 " <<endl;

    vector<vector<vector<float>>> vec;

    // cout << "ok pde 3 " <<endl;

    if (!is_exact_solution)
    {
        cout << "initialize" << endl;
        initialize_2d(initial_value, bound_value);
    }
    // cout << "ok pde 4 " <<endl;

    printf("start solving ... \n");
    auto start = std::chrono::system_clock::now();
    if (is_exact_solution)
        vec = exact_solution_2d(&field, spatio_bound, delta, domain_list[0]["coef"].first);
    else
        vec = solve_2d(spatio_bound, delta, domain_list);
    auto end = std::chrono::system_clock::now();
    double elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / 1000.0; //処理に要した時間をミリ秒に変換
    string e = "end in " + float_to_string(elapsed, 2) + "s";
    cout << e << "\n";

    LOG_FILE.open(LOG_PATH, ios::app);
    LOG_FILE << e << "\n";
    LOG_FILE.close();

    OUTPUT_FILE.open(OUTPUT_FILE_PATH, ios::out);
    std::cout << "writing " << OUTPUT_FILE_PATH << "..." << std::endl;

    float t_now, x_now, y_now;
    for (i = 0; i < vec.size(); i++)
    {
        t_now = t_bound.first + i * delta["t"];
        if (i == 0)
        {
            // printf("time,");
            OUTPUT_FILE << "time,";
            for (j = 0; j < n_cols_y; j++)
            {
                y_now = spatio_bound["y"].first + delta["y"] * j;
                for (k = 0; k < n_cols_x; k++)
                {
                    x_now = spatio_bound["x"].first + delta["x"] * k;
                    OUTPUT_FILE << "(x y)=" << fixed << setprecision(3) << x_now << " " << y_now;
                    if (j < n_cols_y - 1 || k < n_cols_x - 1)
                    {
                        OUTPUT_FILE << ",";
                    }
                }
                // cout << fixed << setprecision(3) << x_now;
            }
            OUTPUT_FILE << endl;
        }
        OUTPUT_FILE << t_now << ",";

        for (j = 0; j < n_cols_y; j++)
        {
            for (k = 0; k < n_cols_x; k++)
            {
                if (j < n_cols_y - 1 || k < n_cols_x - 1)
                {
                    // printf("%f\t", vec[i][j]);
                    OUTPUT_FILE << vec[i][j][k] << ",";
                }
                else
                {
                    // printf("%f", vec[i][j]);
                    OUTPUT_FILE << vec[i][j][k];
                }
            }
        }
        OUTPUT_FILE << endl;
    }

    OUTPUT_FILE.close();
    return 0;
}
