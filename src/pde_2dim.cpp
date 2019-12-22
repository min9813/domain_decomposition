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
#include <chrono> // time calculation
#include <filesystem> // make directory
#include "yaml-cpp/yaml.h"
using namespace std;
namespace fs = filesystem;
// namespace plt = matplotlibcpp;

// template <typename float>
// int pde_2d(YAML::Node config, char* argv[], float type_value);

vector<vector<vector<float>>> explicit_solve(vector<vector<vector<float>>> *field, int now_step, float calc_time_step, map<string, float> delta, map<string, pair<int, int>> domain_col, float coef = 1.)
{
    int t, i, j;
    int calc_step = int(calc_time_step / delta["t"]);
    int start_xid, end_xid, start_yid, end_yid;
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
    // cout << "end col:" << domain_col.second << endl;
    // cout << "start col:" << start_col << endl;

    ratio_x = coef * delta["t"] / (delta["x"]*delta["x"]);
    ratio_y = coef * delta["t"] / (delta["y"]*delta["y"]);
    // cout << "explicit" << endl;
    // cout << domain_col << endl;
    // initialize
    // printf("%lx %lx\n", (long)field, (long)&(*field)[0]);

    // solve
    for (t = now_step; t < calc_step + now_step; t++)
    {
        // cout << "ok exp2 " << now_step <<" i=  " << i<< endl;
        for (i = start_yid; i < domain_col["y"].second; i++)
        {
            for (j = start_xid; j < domain_col["x"].second; j++)
            {
                x_diff = (*field)[t - 1][i][j + 1] + (*field)[t - 1][i][j - 1] - 2 * (*field)[t - 1][i][j];
                y_diff = (*field)[t - 1][i + 1][j] + (*field)[t - 1][i - 1][j] - 2 * (*field)[t - 1][i][j];
                (*field)[t][i][j] = (*field)[t - 1][i][j] + ratio_x * x_diff + ratio_y * y_diff;
            }
        }
    }

    // cout << "end exp3" << endl;

    return *field;
}

vector<vector<vector<float>>> implicit_solve(vector<vector<vector<float>>> *field, int now_step, float calc_time_step, map<string, float> delta, map<string, pair<int, int>> domain_col, float coef = 1.)
{
    int t, i, j, k, index, x_idx, y_idx, coef_idx, field_xid, field_yid, coef_k;
    int calc_step = int(calc_time_step / delta["t"]);
    int start_xid, end_xid, start_yid, end_yid, tmp_field_xid, tmp_field_yid;
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
    // cout << "end col:" << domain_col.second << endl;
    // cout << "start col:" << start_col << endl;

    ratio_x = coef * delta["t"] / (delta["x"]*delta["x"]);
    ratio_y = coef * delta["t"] / (delta["y"]*delta["y"]);
    // cout << "ratio" << endl;
    // cout << ratio_x << ratio_y << delta << coef << endl;

    end_xid = domain_col["x"].second;
    end_yid = domain_col["y"].second;

    int f2c_idx = 0;
    int x_idx_num, y_idx_num, coef_n_cols;
    x_idx_num = end_xid - start_xid;
    y_idx_num = end_yid - start_yid;
    coef_n_cols = x_idx_num * y_idx_num;

    float substitute_sum;

    vector<float> coef_list{-ratio_y, -ratio_x, 1 + 2 * ratio_x + 2 * ratio_y, -ratio_x, -ratio_y};
    // vector<float> x0_coef_list{-ratio_y, -ratio_x, 1 + 2 * ratio_x + 2*ratio_y, -ratio_x, -ratio_y};

    // cout << "n_cols  " << (*field)[0].size() << endl;
    // cout << "n_rows  " << n_rows << endl;
    // making coef matrix
    vector<vector<float>> coef_matrix(coef_n_cols, vector<float>(coef_n_cols, 0));
    vector<int> index_list{-x_idx_num, -1, 0, 1, x_idx_num};

    // cout << "implicit" << endl;
    // cout << domain_col << endl;
    // vector<int> pivot_indexes;

    for (i = 0; i < coef_n_cols; i++)
    {
        for (j = 0; j < index_list.size(); j++)
        {
            index = i + index_list[j];
            if (index < 0 || index >= coef_n_cols)
                continue;
            if ((index / x_idx_num != i / x_idx_num) && (index % x_idx_num != i % x_idx_num))
            {
                continue;
            }
            coef_matrix[i][i + index_list[j]] = coef_list[j];
        }
    }
    // cout << coef_list << endl;
    // cout << "-------------" << "coef matrix ini" << "------------" << endl;
    // for (int i=0;i<coef_matrix.size();i++){
        // cout << coef_matrix[i] << endl;
    // }
    // exit(1);
    // cout << coef_matrix << endl;

    // printf("ok1");

    // cout << "coef matrix:" << coef_n_cols << "," << coef_n_cols << endl;
    // cout << "---------- coef mat initial ------------" << endl;
    // cout << coef_matrix << endl;

    // std::tie(coef_matrix, pivot_indexes) = lu_partition(coef_matrix);
    coef_matrix = cholesky_factorize(coef_matrix);
    // cout << "-------------" << "factorized coef matrix" << "------------" << endl;
    // for (int i=0;i<coef_matrix.size();i++){
        // cout << fixed << setprecision(4) << coef_matrix[i] << endl;
    // }
    // exit(1);
    // cout << "i = " << coef_matrix << endl;
    // cout << "pv = " << pivot_indexes << endl;
    // exit(0);
    // cout << "------- field prev step ---------" << endl;
    // for(int i=0;i<(*field)[now_step].size();i++){
        // cout << i << (*field)[now_step-1][i] << endl;
    // }

    // start time step
    for (i = now_step; i <= now_step + calc_time_step; i++)
    {

        // cout << "############### forward #############" << endl;
        // cout << "imp start fsub i = " << i << endl;
        // cout << "end col = " << end_col << " start col " << start_col << endl;

        // cout << (*field)[i-1] <<endl;
        for (x_idx = start_xid; x_idx < end_xid; x_idx++)
        {
            if (x_idx == start_xid)
            {
                (*field)[i - 1][start_yid][start_xid] -= (*field)[i - 1][start_yid][start_xid - 1] * coef_list[1];
                (*field)[i - 1][start_yid][start_xid] -= (*field)[i - 1][start_yid - 1][start_xid] * coef_list[0];
                (*field)[i - 1][end_yid - 1][start_xid] -= (*field)[i - 1][end_yid - 1][start_xid - 1] * coef_list[1];
                (*field)[i - 1][end_yid - 1][start_xid] -= (*field)[i - 1][end_yid][start_xid] * coef_list[4];
            }
            else if (x_idx == end_xid - 1)
            {
                (*field)[i - 1][start_yid][x_idx] -= (*field)[i - 1][start_yid][x_idx + 1] * coef_list[3];
                (*field)[i - 1][start_yid][x_idx] -= (*field)[i - 1][start_yid - 1][x_idx] * coef_list[0];
                (*field)[i - 1][end_yid - 1][x_idx] -= (*field)[i - 1][end_yid - 1][x_idx + 1] * coef_list[3];
                (*field)[i - 1][end_yid - 1][x_idx] -= (*field)[i - 1][end_yid][x_idx] * coef_list[4];
            }
            else
            {
                (*field)[i - 1][start_yid][x_idx] -= (*field)[i - 1][start_yid - 1][x_idx] * coef_list[0];
                (*field)[i - 1][end_yid - 1][x_idx] -= (*field)[i - 1][end_yid][x_idx] * coef_list[4];
            }
        }
        for (y_idx = start_yid + 1; y_idx < end_yid - 1; y_idx++)
        {
            (*field)[i - 1][y_idx][start_xid] -= (*field)[i - 1][y_idx][start_xid - 1] * coef_list[1];
            (*field)[i - 1][y_idx][end_xid - 1] -= (*field)[i - 1][y_idx][end_xid] * coef_list[3];
        }

        // cout << "------- field prev step after add ---------" << endl;
        // for(int jj=0;jj<(*field)[now_step].size();jj++){
            // cout << jj << (*field)[now_step-1][jj] << endl;
        // }

        // cout << "------- field now step after add ---------" << endl;
        // for(int jj=0;jj<(*field)[now_step].size();jj++){
            // cout << jj << (*field)[now_step][jj] << endl;
        // }

        // start forward substitution
        // cout << "i = " << i << endl;
        // cout << "imp ok 1" << endl;

        // cout << (*field)[i-1] <<endl;

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
                    substitute_sum += coef_matrix[coef_idx][k] * (*field)[i][tmp_field_yid][tmp_field_xid];
                    // cout << k << " fx " << tmp_field_xid<<" fy " << tmp_field_yid ;
                    // cout << " sub sum: " << substitute_sum << " coef mat " <<coef_matrix[coef_idx][k] << " field "<< (*field)[i][tmp_field_yid][tmp_field_xid]<<endl;
                }
                // cout << y_idx <<" " << x_idx <<" " <<coef_idx<<" substitute sum = " << substitute_sum << endl;
                (*field)[i][y_idx][x_idx] = ((*field)[i - 1][y_idx][x_idx] - substitute_sum) / coef_matrix[coef_idx][coef_idx];
                // cout << "*** field matrix ***" << endl;
                // cout << y_idx<<(*field)[i][y_idx] <<endl;
            }
        }
        // cout << "------------- after forward ----------------" << endl;
        // for(int jj=0;jj<(*field)[now_step].size();jj++){
            // cout << jj<< fixed << setprecision(4)  << (*field)[now_step][jj] << endl;
        // }

        // cout << "imp ok 2" << endl;

        for (x_idx = start_xid; x_idx < end_xid; x_idx++)
        {
            if (x_idx == start_xid)
            {
                (*field)[i - 1][start_yid][start_xid] += (*field)[i - 1][start_yid][start_xid - 1] * coef_list[1];
                (*field)[i - 1][start_yid][start_xid] += (*field)[i - 1][start_yid - 1][start_xid] * coef_list[0];
                (*field)[i - 1][end_yid - 1][start_xid] += (*field)[i - 1][end_yid - 1][start_xid - 1] * coef_list[3];
                (*field)[i - 1][end_yid - 1][start_xid] += (*field)[i - 1][end_yid][start_xid] * coef_list[4];
            }
            else if (x_idx == end_xid - 1)
            {
                (*field)[i - 1][start_yid][x_idx] += (*field)[i - 1][start_yid][x_idx + 1] * coef_list[1];
                (*field)[i - 1][start_yid][x_idx] += (*field)[i - 1][start_yid - 1][x_idx] * coef_list[0];
                (*field)[i - 1][end_yid - 1][x_idx] += (*field)[i - 1][end_yid - 1][x_idx + 1] * coef_list[3];
                (*field)[i - 1][end_yid - 1][x_idx] += (*field)[i - 1][end_yid][x_idx] * coef_list[4];
            }
            else
            {
                (*field)[i - 1][start_yid][x_idx] += (*field)[i - 1][start_yid - 1][x_idx] * coef_list[0];
                (*field)[i - 1][end_yid - 1][x_idx] += (*field)[i - 1][end_yid][x_idx] * coef_list[4];
            }
        }
        for (y_idx = start_yid + 1; y_idx < end_yid - 1; y_idx++)
        {
            (*field)[i - 1][y_idx][start_xid] += (*field)[i - 1][y_idx][start_xid - 1] * coef_list[1];
            (*field)[i - 1][y_idx][end_xid - 1] += (*field)[i - 1][y_idx][end_xid] * coef_list[3];
        }

        // cout << (*field)[i-1] <<endl;
        // cout << (*field)[i] <<endl;

        // cout << "imp ok 3" << endl;

        // cout << "############### backward #############" << endl;

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
                    substitute_sum += coef_matrix[k][coef_idx] * (*field)[i][tmp_field_yid][tmp_field_xid];
                    // cout << k << " fx " << tmp_field_xid<<" fy " << tmp_field_yid ;
                    // cout << " sub sum: " << substitute_sum << " coef mat " <<coef_matrix[coef_idx][k] << " field "<< (*field)[i][tmp_field_yid][tmp_field_xid]<<endl;

                }
                (*field)[i][y_idx][x_idx] = ((*field)[i][y_idx][x_idx] - substitute_sum) / coef_matrix[coef_idx][coef_idx];

                // cout << y_idx <<" " << x_idx <<" " <<coef_idx<<" substitute sum = " << substitute_sum << endl;
                // cout << "*** field matrix ***" << endl;
                // cout << y_idx<<(*field)[i][y_idx] <<endl;
            }
        }
        // cout << "------------- after backward ----------------" << endl;
        // for(int jj=0;jj<(*field)[now_step].size();jj++){
            // cout << jj<< fixed << setprecision(4)  << (*field)[now_step][jj] << endl;
        // }
    }

    // exit(1);

    // for (i=0;i<10;i++){
    // cout << (*field)[i] << endl;

    // }

    return *field;
}

vector<vector<vector<float>>> solve(vector<vector<vector<float>>> *field, map<string, pair<float, float>> spatio_bound, map<string, float> delta, vector<map<string, pair<float, float>>> domain_boundary)
{
    int now_step, step, j, k, tmp_col1, tmp_col2;
    float coef;
    int n_rows = (*field).size();
    int n_cols1 = (*field)[0].size();
    int n_cols2 = (*field)[0][0].size();
    map<string, map<string, float>> delta_info;

    vector<map<string, pair<int, int>>> domain_col(domain_boundary.size());
    for (j = 0; j < domain_boundary.size(); j++)
    {
        map<string, pair<int, int>> tmp_domain_col;
        for (auto itr = domain_boundary[j].begin(); itr != domain_boundary[j].end(); itr++)
        {
            if (itr->first == "x" || itr->first == "y"){
                tmp_col1 = (int)((itr->second.first - spatio_bound[itr->first].first) / (delta[itr->first]));
                tmp_col2 = (int)((itr->second.second - spatio_bound[itr->first].first) / (delta[itr->first]));
                tmp_domain_col[itr->first] = make_pair(tmp_col1, tmp_col2);
            }else{
                tmp_domain_col[itr->first] = itr->second;                
            }
        }
        domain_col[j] = tmp_domain_col;
    }

    delta_info["exp"] = delta;
    delta_info["imp"] = delta;

    // cout << "ok1" << endl;
    std::chrono::system_clock::time_point  start, end;
    float elapsed = 0.;
    vector<float> domain_elapsed(domain_col.size(), 0.);
    // vector<int> counter(domain_col.size(), 0.);
    for (now_step = 1; now_step < n_rows; now_step++)
    {
        for (j = 0; j < domain_col.size(); j++)
        {
            coef = domain_col[j]["coef"].first;

            if (domain_col[j]["is_exp"].first==1){
                if(domain_col[j]["t_delta"].first > 0){
                    delta_info["exp"]["t"] = domain_col[j]["t_delta"].first;
                }
                start = std::chrono::system_clock::now();
                explicit_solve(field, now_step, delta_info["exp"]["t"], delta_info["exp"], domain_col[j], coef);
                end = std::chrono::system_clock::now();
                elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();
                // cout << "-------------- "<< j << " exp ------------" <<endl;
                // for(int i=0;i<n_cols1;i++){
                    // cout << "y=" << i << " " << (*field)[now_step][i] << endl;
                // }
            }else{
                if(domain_col[j]["t_delta"].first > 0){
                    delta_info["imp"]["t"] = domain_col[j]["t_delta"].first;
                }
                
                start = std::chrono::system_clock::now();
                implicit_solve(field, now_step, delta_info["imp"]["t"], delta_info["imp"], domain_col[j], coef);
                end = std::chrono::system_clock::now();
                elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();
                // cout << "-------------- "<< j << " imp ------------" <<endl;
                // for(int i=0;i<n_cols1;i++){
                    // cout << "y=" << i << " " << (*field)[now_step][i] << endl;
                // }
            }
            domain_elapsed[j] += elapsed / (n_rows-1);
        }

        // cout << "-------------- all ------------" <<endl;
        // for(int i=0;i<n_cols1;i++){
            // cout << "y=" << i << " " << (*field)[now_step][i] << endl;
        // }
        // if(now_step > 1){
            // exit(1);
        // }
    }

    string tmp, tmp2;
    LOG_FILE.open(LOG_PATH, ios::app);
    LOG_FILE << "--- mean time | total (" << n_rows-1 << " step) time ---\n";
    for(j=0;j<domain_elapsed.size();j++){
        tmp = float_to_string(domain_elapsed[j], 2);
        LOG_FILE << tmp << "ms | ";
        tmp = float_to_string(domain_elapsed[j]*(n_rows-1), 2);
        LOG_FILE << tmp <<"ms" << "\n";
    }
    LOG_FILE.close();
    return *field;
}

vector<vector<vector <float>>> initialize(vector<vector<vector<float>>> *field, float initial_value, map<string, pair<float, float>> bound_value)
{
    int i, j, k;
    int n_rows = (*field).size();
    int n_cols = (*field)[0].size();
    int n_cols2 = (*field)[0][0].size();

    for (i = 0; i < n_rows; i++)
    {
        for (j = 0; j < n_cols; j++)
        {
            for (k = 0; k < n_cols2; k++)
            {
                if (i == 0 & j > 0 & j<n_cols - 1 & k> 0 & k < n_cols2 - 1)
                {
                    (*field)[0][j][k] = initial_value;
                }
                else if (k == 0)
                {
                    (*field)[i][j][k] = bound_value["x"].first;
                }
                else if (k == n_cols2 - 1)
                {
                    (*field)[i][j][k] = bound_value["x"].second;
                }
                else if (j == 0)
                {
                    (*field)[i][j][k] = bound_value["y"].first;
                }
                else if (j == n_cols - 1)
                {
                    (*field)[i][j][k] = bound_value["y"].second;
                }
            }
        }
    }

    // cout << "########## initial ###########" << endl;
    // cout << (*field) << endl;

    return *field;
}

int pde_2d(YAML::Node config, char* argv[], float type_value)
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

    is_exact_solution = (config["is_exact"].as<int>()>0);


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

    auto spatio_bound = config["space"].as<map<string, pair<float ,float>>>();
    auto t_bound = config["time"].as<pair<float, float>>();


    auto bound_value = config["condition"]["bound_value"].as<map<string, pair<float, float>>>();
    float initial_value = config["condition"]["initial_value"]["t"].as<float>();

    vector<map<string, pair<float, float>>> domain_list;
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
    float y_length = spatio_bound["y"].second - spatio_bound["y"].first;
    float time_length = t_bound.second - t_bound.first;

    n_cols_x = (int)(x_length / delta["x"]) + 1;
    n_cols_y = (int)(y_length / delta["y"]) + 1;
    n_rows = (int)(time_length / delta["t"]) + 1;

    cout << "ok pde 2 " <<endl;

    vector<vector<vector <float>>> vec(n_rows, vector<vector<float>>(n_cols_y, vector<float>(n_cols_x, 0)));
    cout << "ok pde 3 " <<endl;
    
    if (!is_exact_solution){
        cout << "initialize" << endl;
        vec = initialize(&vec, initial_value, bound_value);
    }
    cout << "ok pde 4 " <<endl;

    printf("start solving ... \n");
    auto start = std::chrono::system_clock::now();
    if (is_exact_solution) vec = exact_solution_2d(&vec, spatio_bound, delta, domain_list[0]["coef"].first);
    else vec = solve(&vec, spatio_bound, delta, domain_list);
    auto end = std::chrono::system_clock::now();
    double elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count() / 1000.0; //処理に要した時間をミリ秒に変換
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
                for(k=0; k<n_cols_x;k++){
                    x_now = spatio_bound["x"].first + delta["x"] * k;
                    OUTPUT_FILE << "(x y)=" << fixed << setprecision(3) << x_now << " " << y_now;
                    if (j < n_cols_y - 1 || k < n_cols_x - 1){
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
            for(k = 0; k < n_cols_x; k++){
                if (j < n_cols_y - 1 || k < n_cols_x-1)
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
