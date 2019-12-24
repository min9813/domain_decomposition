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

void explicit_solve(int now_step, float calc_time_step, map<string, float> delta, map<string, pair<int, int>> domain_col, float coef = 1.)
{
    int t, i, j;
    int calc_step = int(calc_time_step / delta["t"]);
    int start_xid, end_xid;
    float x_diff, ratio_x;
    if (domain_col["x"].first == 0)
        start_xid = 1;
    else
        start_xid = domain_col["x"].first;
    // cout << "ok exp1" << endl;
    // cout << "end col:" << domain_col.second << endl;
    // cout << "start col:" << start_col << endl;

    float ratio = coef * delta["t"] / (delta["x"]*delta["x"]);

    cout << "exp ratio: " << ratio << endl;
    // cout << "explicit " << domain_col << endl;
    // initialize
    // printf("%lx %lx\n", (long)field, (long)&field[0]);

    // solve
    for (i = now_step; i < calc_step + now_step; i++)
    {
        // cout << "ok exp2 " << now_step <<" i=  " << i<< endl;
        for (j = start_xid; j < domain_col["x"].second; j++)
        {
            field[i][j] = ratio * field[i - 1][j - 1] + (1 - 2 * ratio) * field[i - 1][j] + ratio * field[i - 1][j + 1];
        }
    }

    // cout << "end exp3" << endl;

}

void implicit_solve(int now_step, float calc_time_step, map<string, float> delta, map<string, pair<int, int>> domain_col, float coef = 1.)
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
    float substitute_sum;

    vector<float> coef_list{-ratio, 1 + 2 * ratio, -ratio};
    // cout << "n_cols  " << field[0].size() << endl;
    // cout << "n_rows  " << n_rows << endl;
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

vector<vector<float>> solve(map<string, pair<float, float>> spatio_bound, map<string, float> delta, vector<map<string, pair<float, float>>> domain_boundary)
{
    int now_step, step, j, k, tmp_col1, tmp_col2;
    int n_rows = field.size();
    int n_cols = field[0].size();
    map<string, map<string, float>> delta_info;
    float coef;

    vector<map<string, pair<int, int>>> domain_col(domain_boundary.size());
    for (j = 0; j < domain_boundary.size(); j++)
    {
        map<string, pair<int, int>> tmp_domain_col;
        for (auto itr = domain_boundary[j].begin(); itr != domain_boundary[j].end(); itr++)
        {
            if (itr->first == "x"){
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
                explicit_solve(now_step, delta_info["exp"]["t"], delta_info["exp"], domain_col[j], coef);
                end = std::chrono::system_clock::now();
                elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();
                // cout << "-------------- "<< j << " exp ------------" <<endl;
                // for(int i=0;i<n_cols1;i++){
                    // cout << "y=" << i << " " << field[now_step][i] << endl;
                // }
            }else{
                if(domain_col[j]["t_delta"].first > 0){
                    delta_info["imp"]["t"] = domain_col[j]["t_delta"].first;
                }
                
                start = std::chrono::system_clock::now();
                implicit_solve(now_step, delta_info["imp"]["t"], delta_info["imp"], domain_col[j], coef);
                end = std::chrono::system_clock::now();
                elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();
                // cout << "-------------- "<< j << " imp ------------" <<endl;
                // for(int i=0;i<n_cols1;i++){
                    // cout << "y=" << i << " " << field[now_step][i] << endl;
                // }
            }
            domain_elapsed[j] += elapsed / (n_rows-1);
        }
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

    return field;
}

void initialize(float initial_value, map<string, pair<float, float>> bound_value){
    int i,j;
    int n_cols = field.size();
    int n_rows = field[0].size();
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
    }

    // cout << "########## initial ###########" << endl;
    // cout << field << endl;
}

int pde(YAML::Node config, char *argv[], float type_value)
{
    cout << "start pde" << "\n";
    // 初期値、境界値（x=0, 1, t=0 の3つ）、tiem delta、x_delta を引数にすると良さそう
    int n_cols, n_rows,i,j;
    bool is_exact_solution = (config["is_exact"].as<int>()>0);
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

    vector<vector<float>> tmp_vec(n_rows, vector<float>(n_cols, 0));
    field = tmp_vec;
    tmp_vec.clear();

    vector<vector<float>> vec;
    if (!is_exact_solution){
        cout << "initialize" << endl;
        initialize(initial_value, bound_value);
    }

    // vector<vector<float> > vec;
    printf("start solving ... \n");
    auto start = std::chrono::system_clock::now();
    if (is_exact_solution) vec = exact_solution_1d(&field, spatio_bound, delta, domain_list[0]["coef"].first);
    else vec = solve(spatio_bound, delta, domain_list);
    auto end = std::chrono::system_clock::now();
    double elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count() / 1000.0; //処理に要した時間をミリ秒に変換
    string e = "end in " + float_to_string(elapsed, 2) + "s";
    cout << e << "\n";
    // vec = implicit_solve(time_length, t_delta, x_length, x_delta, initial_value, x_value);
    LOG_FILE.open(LOG_PATH, ios::app);
    LOG_FILE << e << "\n";
    LOG_FILE.close();

    OUTPUT_FILE.open(OUTPUT_FILE_PATH, ios::out);
    std::cout << "writing " << OUTPUT_FILE_PATH << "..." << std::endl;
    float t_now, x_now;

    cout << "writing " << OUTPUT_FILE_PATH << "..." << endl;

    for (i = 0; i < vec.size(); i++)
    {
        t_now = t_bound.first + i * delta["t"];
        if(i==0){
            printf("time,");
            OUTPUT_FILE << "time,";
            for(j=0;j < n_cols; j++){
                x_now = spatio_bound["x"].first + delta["x"] * j;
                cout << fixed << setprecision(3) << x_now;
                OUTPUT_FILE << "x=" <<fixed <<setprecision(3)<< x_now;
                if (j<n_cols-1){
                    cout << "\t";
                    OUTPUT_FILE << ",";
                }
            }
            OUTPUT_FILE << endl;
            cout << endl;
        }
        printf("t=%f\t", t_now);
        OUTPUT_FILE << t_now << ",";

        for (j = 0; j < vec[0].size(); j++)
        {
            if(j < vec[0].size()-1){
                printf("%f\t", vec[i][j]);
                OUTPUT_FILE << vec[i][j] << ",";
            }else{
                printf("%f", vec[i][j]);
                OUTPUT_FILE << vec[i][j];
            }
            
        }
        OUTPUT_FILE << endl;
        printf("\n");
    }

    OUTPUT_FILE.close();
    return 0;
}
