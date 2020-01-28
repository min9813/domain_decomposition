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
#include <thread>
using namespace std;
namespace fs = filesystem;
// namespace plt = matplotlibcpp;

// template <typename float>
// int pde_2d(YAML::Node config, char* argv[], float type_value);
static vector<vector<vector<float>>> field;
static float choleskey_elapsed_time = 0, gauss_time = 0;
static bool is_sparse_lu = true;
static bool is_sep = false;

inline float postprocess_marginal_field_mp(int i, int x_idx, int y_idx, int start_xid, int start_yid, int end_xid, int end_yid, vector<float> coef_list_large, vector<float> coef_list_small){
    float prev_value;
    if(x_idx==start_xid){
        if(y_idx==start_yid){
            prev_value = field[i-1][y_idx][x_idx] - field[i-1][y_idx][x_idx-1] * coef_list_small[1];
            prev_value = prev_value - field[i-1][y_idx][x_idx-1] * coef_list_small[0];
        }else if(y_idx==(end_yid-1)){
            prev_value = field[i-1][y_idx][x_idx] - field[i-1][y_idx][x_idx-1] * coef_list_small[1];
            prev_value = prev_value - field[i-1][end_yid][x_idx] * coef_list_small[4];
        }else{
            prev_value = field[i-1][y_idx][x_idx] - field[i-1][y_idx][x_idx-1] * coef_list_small[1];
        }
    }else if(x_idx==(end_xid-1)){
        if(y_idx==start_yid){
            prev_value = field[i-1][y_idx][x_idx] - field[i-1][y_idx][x_idx+1] * coef_list_large[3];
            prev_value = prev_value - field[i-1][y_idx-1][x_idx] * coef_list_large[0];
        }else if(y_idx==(end_yid-1)){
            prev_value = field[i-1][y_idx][x_idx] - field[i-1][y_idx][x_idx+1] * coef_list_large[3];
            prev_value = prev_value - field[i-1][end_yid][x_idx] * coef_list_large[4];
        }else{
            prev_value = field[i-1][y_idx][x_idx] - field[i-1][y_idx][end_xid] * coef_list_large[3];
        }
    }else{
        if(y_idx==start_yid){
            prev_value = field[i-1][y_idx][x_idx] - field[i-1][y_idx-1][x_idx] * coef_list_small[0];
        }else if(y_idx==(end_yid-1)){
            prev_value = field[i-1][y_idx][x_idx] - field[i-1][end_yid][x_idx] * coef_list_large[4];
        }else{
            prev_value = field[i-1][y_idx][x_idx];
        }
    }
    return prev_value;
}


void explicit_solve_2d_mp(int now_step, float calc_time_step, map<string, float> delta, map<string, pair<float, float>> domain_col, float coef = 1., pair<float, float> neighb_coef={1., 1.})
{
    int t, i, j;
    // int calc_step = int(calc_time_step / delta["t"]);
    int calc_step = 1;
    int start_xid, end_xid, start_yid, end_yid;
    int n_cols1 = field[0].size();
    float x_diff, ratio_x, ratio_y, y_diff, ratio_x_small, ratio_x_large;
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
    ratio_x_small = 2 * neighb_coef.first * coef / (neighb_coef.first + coef) * delta["t"] / (delta["x"] * delta["x"]);
    ratio_x_large = 2 * neighb_coef.second * coef / (neighb_coef.second + coef) * delta["t"] / (delta["x"] * delta["x"]);

    ratio_y = coef * delta["t"] / (delta["y"] * delta["y"]);
    // cout << "exp r_x:" << ratio_x << " r_y:" << ratio_y << endl;

    float ratio = ratio_x + ratio_y;
    if(ratio >= 1./2.){
        printf("must satisfy von-neuman's stable conditions. actual ratio = %f", ratio);
        exit(1);
    }
    ratio = ratio_x_small + ratio_y;
    if(ratio >= 1./2.){
        ratio_x_small = ratio_x;
    }
    ratio = ratio_x_large + ratio_y;
    if(ratio >= 1./2.){
        ratio_x_large = ratio_x;
    }

    // cout << "coef:" << coef << " neighb coef:" << neighb_coef <<endl;
    // exit(0);
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
    if(start_xid<domain_col["x"].second){
        for (t = now_step; t < calc_step + now_step; t++)
        {
        // cout << "ok exp2 " << now_step <<" i=  " << i<< endl;
            for (i = start_yid; i < domain_col["y"].second; i++)
            {
                j = start_xid;
                x_diff = ratio_x * field[t - 1][i][j + 1] + ratio_x_small * field[t - 1][i][j - 1] - ratio_x * 2 * field[t - 1][i][j];
                y_diff = field[t - 1][i + 1][j] + field[t - 1][i - 1][j] - 2 * field[t - 1][i][j];
                field[t][i][j] = field[t - 1][i][j] + x_diff + ratio_y * y_diff;

                for (j = start_xid+1; j < domain_col["x"].second-1; j++)
                {
                    x_diff = field[t - 1][i][j + 1] + field[t - 1][i][j - 1] - 2 * field[t - 1][i][j];
                    y_diff = field[t - 1][i + 1][j] + field[t - 1][i - 1][j] - 2 * field[t - 1][i][j];
                    field[t][i][j] = field[t - 1][i][j] + ratio_x * x_diff + ratio_y * y_diff;
                // count ++;
                }
            // if(j!=domain_col["x"].second-1){
                // exit(1);
            // }
                if(j<domain_col["x"].second){
                    x_diff = ratio_x_large * field[t - 1][i][j + 1] + ratio_x * field[t - 1][i][j - 1] - ratio_x * 2 * field[t - 1][i][j];
                    y_diff = field[t - 1][i + 1][j] + field[t - 1][i - 1][j] - 2 * field[t - 1][i][j];
                    field[t][i][j] = field[t - 1][i][j] + x_diff + ratio_y * y_diff;
                }
                
            }
        // cout << "t:" << t << " count:" <<count << endl;
        
        }
    }

    // printf("ratio_x %f ratio_x_large %f ratio_small %f\n", ratio_x, ratio_x_large, ratio_x_small);
   
    // exit(0);

    // cout << "end exp3" << endl;
}

map<string, pair<int, int>> calc_column_status_mp(map<string, pair<float, float>> domain_col){
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


void sparse_symmetric_implicit_solve_2d_xy_mp(int now_step, float calc_time_step, map<string, float> delta, map<string, pair<float, float>> domain_col, float coef = 1., pair<float, float> neighb_coef={1., 1.})
{
    int t, i, j, k, index, x_idx, y_idx, coef_idx, field_xid, field_yid, coef_k, dense_index, max_gauss_coef_index, add_index, coef_row_base;
    int calc_step = int(calc_time_step / delta["t"]);
    int start_xid, end_xid, start_yid, end_yid, tmp_field_xid, tmp_field_yid;
    float x_diff, ratio_x, ratio_y, y_diff, ratio_x_small, ratio_x_large;
    map<string, pair<int, int>> col_status;
    col_status = calc_column_status_mp(domain_col);
    start_xid = col_status["x"].first;
    end_xid = col_status["x"].second;

    start_yid = col_status["y"].first;
    end_yid = col_status["y"].second;

    ratio_x = coef * delta["t"] / (delta["x"] * delta["x"]);
    ratio_y = coef * delta["t"] / (delta["y"] * delta["y"]);
    ratio_x_small = 2 * (neighb_coef.first * coef) / (neighb_coef.first + coef) * delta["t"] / (delta["x"] * delta["x"]);
    ratio_x_large = 2 * (neighb_coef.second * coef) / (neighb_coef.second + coef) * delta["t"] / (delta["x"] * delta["x"]);
    // cout << "imp r_x:" << ratio_x << " r_y:" << ratio_y << endl;

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
    // printf("imp ratio_x %f ratio_x_large %f ratio_small %f\n", ratio_x, ratio_x_large, ratio_x_small);


    vector<float> coef_list{-ratio_y, -ratio_x, 1 + 2 * ratio_x + 2 * ratio_y, -ratio_x, -ratio_y};
    vector<float> coef_list_small{-ratio_y, -ratio_x_small, 1 + 2 * ratio_x + 2 * ratio_y, -ratio_x, -ratio_y};
    vector<float> coef_list_large{-ratio_y, -ratio_x, 1 + 2 * ratio_x + 2 * ratio_y, -ratio_x_large, -ratio_y};


    vector<float> coef_sparse_list_main{1 + 2 * ratio_x + 2 * ratio_y, -ratio_x, -ratio_y};
    vector<float> coef_sparse_list_small{1 + 2 * ratio_x + 2 * ratio_y, -ratio_x_small, -ratio_y};
    vector<float> coef_sparse_list_large{1 + 2 * ratio_x + 2 * ratio_y, -ratio_x_large, -ratio_y};

    vector<float> coef_sparse_list;
    vector<vector<float>> coef_matrix(coef_n_rows, vector<float>(coef_n_cols, 0));
    vector<int> index_list{0, 1, x_idx_num};

    // cout << "n rows:" << coef_n_rows <<" n cols:" << coef_n_cols<<  endl;
    // exit(0);
    for (i = 0; i < coef_n_rows; i++)
    {
        x_idx = i % x_idx_num;
        y_idx = i / x_idx_num;
        if(x_idx + 1 == x_idx_num){
            coef_sparse_list = coef_sparse_list_large;
        }
        // else if(x_idx == 0){
            // coef_sparse_list = coef_sparse_list_main;
        // }
        else{
            coef_sparse_list = coef_sparse_list_main;
        }
        for (j = 0; j < index_list.size(); j++)
        {
            index = index_list[j];
            dense_index = index + i;
            // cout << "i=" << i <<" j= "<<j<<" index="<<index<<endl;

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
    // cout << " 1 n rows:" << coef_n_rows <<" n cols:" << coef_n_cols<<  endl;

    // cout << "------- before csolve ---------" << endl;
    // for (i = 0; i < coef_n_rows; i++)
    // {
        // cout << "y=" << i << coef_matrix[i] << endl;
    // }
    // exit(0);

    std::chrono::system_clock::time_point start, end;
    float elapsed = 0., prev_value;
    start = std::chrono::system_clock::now();

    coef_matrix = sparse_cholesky_factorize(coef_matrix);
    // cout << "------- after csolve ---------" << endl;
    // for (i = 0; i < coef_n_rows; i++)
    // {
        // cout << "y=" << i << coef_matrix[i] << endl;
    // }
    // cout << " 2 n rows:" << coef_n_rows <<" n cols:" << coef_n_cols<<  endl;

    end = std::chrono::system_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0;
    choleskey_elapsed_time += elapsed;

    // start time step
    start = std::chrono::system_clock::now();
    for (i = now_step; i < now_step + calc_step; i++)
    {


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

                prev_value = postprocess_marginal_field_mp(i, x_idx, y_idx, start_xid, start_yid, end_xid, end_yid, coef_list_large, coef_list_small);
                

                field[i][y_idx][x_idx] = (prev_value - substitute_sum) / coef_matrix[coef_idx][0];
            }
        }
        // cout << " 3 n rows:" << coef_n_rows <<" n cols:" << coef_n_cols<<  endl;


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
        // cout << " 5 n rows:" << coef_n_rows <<" n cols:" << coef_n_cols<< " i :" << i << " now step:"<<now_step <<" calc:"<< calc_step << endl;


    }
    // cout << " 6 n rows:" << coef_n_rows <<" n cols:" << coef_n_cols<<  endl;

    end = std::chrono::system_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0;
    gauss_time += elapsed;
}

void solve_2d_mp(map<string, pair<float, float>> spatio_bound, map<string, float> delta, vector<map<string, pair<float, float>>> domain_boundary)
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
    bool is_same = !is_sep;
    int is_exp = domain_boundary[0]["is_exp"].first;
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
            }else if(itr->first == "is_exp"){
                is_same &= is_exp == (itr->second).first;
                tmp_domain_col[itr->first] = itr->second;
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
    float elapsed = 0., coef_large, coef_small;
    pair<float, float> neighb_coef;
    // vector<int> counter(domain_col.size(), 0.);
    for (now_step = 1; now_step < n_rows; now_step++)
    {
        // cout << "-------------- start step = " << now_step << " ------------" <<endl;
        // for(int i=0;i<n_cols1;i++){
        // cout << "y=" << i << " " << field[now_step][i] << endl;
        // }
        // cout << "is same:" <<is_same <<endl;
        start = std::chrono::system_clock::now();

        vector<thread> thread_list(domain_col.size());
        for (j = 0; j < domain_col.size(); j++)
        {
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
            if (domain_col[j]["is_exp"].first == 1)
            {
                if (domain_col[j]["t_delta"].first > 0)
                {
                    delta_info["exp"]["t"] = domain_col[j]["t_delta"].first;
                }
                // cout << "j="<<j<<" neighb coef:" << neighb_coef <<endl;
                thread_list[j] = thread(explicit_solve_2d_mp, now_step, delta_info["exp"]["t"], delta_info["exp"], domain_col[j], coef, neighb_coef);
                // explicit_solve_2d_mp(now_step, delta_info["exp"]["t"], delta_info["exp"], domain_col[j], coef, neighb_coef);
            
            }
            else
            {
                if (domain_col[j]["t_delta"].first > 0)
                {
                    delta_info["imp"]["t"] = domain_col[j]["t_delta"].first;
                }

                // start = std::chrono::system_cloc::now();
                col_status = calc_column_status_mp(domain_col[j]);
                if (col_status["x"].second - col_status["x"].first <= col_status["y"].second - col_status["y"].first){
                    thread_list[j] = thread(sparse_symmetric_implicit_solve_2d_xy_mp, now_step, delta_info["exp"]["t"], delta_info["exp"], domain_col[j], coef, neighb_coef);

                    // sparse_symmetric_implicit_solve_2d_xy_mp(now_step, delta_info["imp"]["t"], delta_info["imp"], domain_col[j], coef, neighb_coef);
                }else{
                    printf("NOT IMPLEMENTED for y columns <= x ");
                    exit(1);
                }

            }
            // cout << "elapsed:" << elapsed << endl;
        }
        for(thread &t: thread_list){
            t.join();
        }
        end = std::chrono::system_clock::now();
        elapsed += std::chrono::duration_cast<std::chrono::microseconds>(end-start).count();
        
        // exit(1);

        if(now_step % 1000 == 0){
            printf("now step / total : [%d/%d]\n", now_step, n_rows-1);
        }
        // cout << "now step / total :"<<"["<< now_step <<"/"<<n_rows-1 <<"]"<<" %\r";

        // if (now_step>50){
        // exit(1);
        // }

        // cout << "-------------- all ------------" <<endl;
        // for (int i = 0; i < n_cols1; i++)
        // {
            // cout << "y=" << i << " " << field[now_step][i] << endl;
        // }
        // if (now_step > 10)
        // {
            // exit(1);
        // }
        // cout << "now/ total = " << now_step << "/" << n_rows << endl;
    }

    string tmp, tmp2;
    LOG_FILE.open(LOG_PATH, ios::app);
    LOG_FILE << "--- mean time | total (" << n_rows - 1 << " step) time ---\n";
    tmp = float_to_string(elapsed / (float)(n_rows - 1), 2);
    tmp2 = float_to_string(elapsed, 2);
    LOG_FILE << tmp << "ms | " << tmp2 <<" ms" << endl;
    LOG_FILE.close();
    // cout << "finish" << endl;
}

void initialize_2d_mp(float initial_value, map<string, pair<float, float>> bound_value)
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

int pde_2d_mp(YAML::Node config, char *argv[], float type_value)
{
    // if (argc < 2){
    //     printf("must specify config file in command line argument!\n");
    //     exit(1);
    // }
    int i, j, k, n_rows, n_cols_y, n_cols_x, try_time;
    bool is_exact_solution;
    try_time = config["try_time"].as<int>();
    try{
        is_sep = config["is_sep"].as<int>() == 1;
    }catch(YAML::TypedBadConversion<int>){
        is_sep = false;
    }
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



    // cout << "ok pde 2 " <<endl;

    // cout << "ok pde 3 " <<endl;

    if(try_time == 1){
        vector<vector<vector<float>>> tmp_vec(n_rows, vector<vector<float>>(n_cols_y, vector<float>(n_cols_x, 0)));
        field = tmp_vec;
        tmp_vec.clear();
        if (!is_exact_solution)
        {
            cout << "initialize" << endl;
            initialize_2d_mp(initial_value, bound_value);
        }
        // cout << "ok pde 4 " <<endl;

        printf("start solving ... \n");
        auto start = std::chrono::system_clock::now();
        if (is_exact_solution)
            exact_solution_2d(&field, spatio_bound, delta, domain_list[0]["coef"].first);
        else
            solve_2d_mp(spatio_bound, delta, domain_list);
        auto end = std::chrono::system_clock::now();
        double elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / 1000.0; //処理に要した時間をミリ秒に変換
        string e = "end in " + float_to_string(elapsed, 2) + "s";
        cout << e << "\n";

        LOG_FILE.open(LOG_PATH, ios::app);
        LOG_FILE << e << "\n";
        LOG_FILE.close();
    }else{
        double total_elapsed = 0, elapsed;
        vector<double> elapsed_list;
        for(i=0; i<try_time;i++){
            vector<vector<vector<float>>> tmp_vec(n_rows, vector<vector<float>>(n_cols_y, vector<float>(n_cols_x, 0)));
            field = tmp_vec;
            tmp_vec.clear();

            if (!is_exact_solution){
                // cout << "initialize" << endl;
                initialize_2d_mp(initial_value, bound_value);
            }

            // vector<vector<float> > vec;
            // printf("start solving ... \n");
            auto start = std::chrono::system_clock::now();
            if (is_exact_solution) exact_solution_2d(&field, spatio_bound, delta, domain_list[0]["coef"].first);
            else solve_2d_mp(spatio_bound, delta, domain_list);
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
    std::cout << "writing " << OUTPUT_FILE_PATH << "..." << std::endl;

    float t_now, x_now, y_now;
    bool is_ok;
    for (i = 0; i < field.size(); i++)
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
                    OUTPUT_FILE << "(x y)=" << fixed << setprecision(6) << x_now << " " << y_now;
                    if (j < n_cols_y - 1 || k < n_cols_x - 1)
                    {
                        OUTPUT_FILE << ",";
                    }
                }
                // cout << fixed << setprecision(3) << x_now;
            }
            OUTPUT_FILE << endl;
        }
        if(n_rows > 100000){
            is_ok = i % 100 == 0;
        }else{
            is_ok = true;
        }
        if (is_ok){
            OUTPUT_FILE << t_now << ",";

            for (j = 0; j < n_cols_y; j++)
            {
                for (k = 0; k < n_cols_x; k++)
                {
                    if (j < n_cols_y - 1 || k < n_cols_x - 1)
                    {
                        // printf("%f\t", vec[i][j]);
                        OUTPUT_FILE << field[i][j][k] << ",";
                    }
                    else
                    {
                        // printf("%f", vec[i][j]);
                        OUTPUT_FILE << field[i][j][k];
                    }
                }
            }
            OUTPUT_FILE << endl;
            if(i % 10000 == 0){
                printf("now / total = %d / %d\n", i, (int)field.size());
            }
        }
        
    }

    OUTPUT_FILE.close();
    return 0;
}
