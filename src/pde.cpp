# include <iostream>
# include <vector>
# include <cmath>
# include <prettyprint.hpp>
# include "numerical.hpp"
using namespace std;

void explicit_solve(vector<vector<float> > *vec, float delta_x, float delta_t){
    int i,j;
    int n_rows = (*vec).size();
    int n_cols = (*vec)[0].size();
    float ratio = delta_t / pow(delta_x, 2);
    // initialize
    (*vec)[0].assign(n_cols, 1);
    for(i=0;i<n_rows;i++){
        (*vec)[i][0] = 0;
        (*vec)[i][n_cols-1] = 0;
    }
    printf("%lx %lx", (long)vec, (long)&(*vec)[0]);

    // solve    
    for(i=1;i<n_rows;i++){
        for(j=1;j<n_cols-1;j++){
            (*vec)[i][j] = ratio * (*vec)[i-1][j-1] + (1-2*ratio) * (*vec)[i-1][j] + ratio*(*vec)[i-1][j+1];
        }
    }
}

vector<vector<float> > implicit_solve(int time_length, float delta_time, float x_length, float delta_x, float t_ini, vector<float> x_ini, float coef=1.){
    int i,j,k;
    int n_rows = (int)((float)time_length / delta_time);
    int n_cols = (int)((float)x_length / delta_x);
    int f2c_idx = 1;
    int coef_n_cols = n_cols - 2;

    
    float ratio = coef * delta_time / delta_x;
    float substitute_sum;
    
    vector<float> coef_list{-ratio, 1+2*ratio, -ratio};
    cout << "n_rows  " << n_rows << " " << time_length<<" "<<delta_time<< endl;
    cout << "n_cols  " << n_cols << " " << x_length<<" "<<delta_x<< endl;
    cout << "x_ini " << x_ini << endl;
    // making matrix
    vector<vector<float> > field(n_rows, vector<float>(n_cols, 0));


    field[0].assign(n_cols, t_ini);

    for (i=0;i<n_rows;i++){
        field[i][0] = x_ini[0];
        field[i][n_cols-1] = x_ini[1];
    }
    cout << field[0] << endl;


    // making coef matrix
    vector<vector<float> > coef_matrix(coef_n_cols,vector<float>(coef_n_cols,0));

    vector<int> pivot_indexes;

    for(i=0;i<coef_n_cols;i++){
        if(i==0){
            coef_matrix[0][0] = 1 + 2 * ratio;
            coef_matrix[0][1] = - ratio;
        }else if (i==coef_n_cols){
            coef_matrix[coef_n_cols][coef_n_cols] = 1 + 2 * ratio;
            coef_matrix[coef_n_cols][coef_n_cols-1] = - ratio;
        }else{
            for(j=i-1;j<i+2;j++){
                coef_matrix[i][j] = coef_list[1+j-i];
            }
        }
        
    }

    // printf("ok1");

    std::tie(coef_matrix, pivot_indexes) = lu_partition(coef_matrix);
    cout << "i = " << coef_matrix << endl;
    cout << "pv = " << pivot_indexes << endl;



    // start time step
    for(i=1;i<n_rows;i++){
        // start forward substitution
        // cout << "i = " << i << endl;
        for(j=0;j<coef_n_cols;j++){
            // cout << "j = " << j << endl;

            substitute_sum = 0.;
            for(k=0;k<j;k++){
                substitute_sum += coef_matrix[j][k] * field[i-1][pivot_indexes[k]+f2c_idx];
            }
            field[i][j+f2c_idx] = (field[i-1][pivot_indexes[k]+f2c_idx] - substitute_sum) / coef_matrix[j][j];
        }
        // start backward substitution
        for(j=coef_n_cols-1;j>=0;j--){
            // cout << "j = " << j << endl;

            substitute_sum = 0.;
            for(k=coef_n_cols-1;k>j;k--){
                substitute_sum += coef_matrix[j][k] * field[i][k+f2c_idx];
            }
            field[i][j+f2c_idx] = (field[i][k+f2c_idx] - substitute_sum) / coef_matrix[j][j];
        }
    }

    
    return field;
}


int main(){

    // 初期値、境界値（x=0, 1, t=0 の3つ）、tiem delta、x_delta を引数にすると良さそう
    int n_cols = 10;     // 列数
    float initial_value = 1;      // 初期値
    float delta_time, delta_x;
    int i,j, n_rows;
    vector<float> x_value = {0.0, 0.0}; // boundary value
    vector<float> x_boundary = {0.0, 1.0}; // boundary value
    vector<float> t_boundary = {0.0, 1.0};

    float x_length = x_boundary[1] - x_boundary[0];
    float time_length = t_boundary[1] - t_boundary[0];
    delta_x = x_length/float(n_cols);
    // delta_time = 0.1;
    delta_time = 0.001;
    // printf("delta_x = %f", delta_x);

    n_rows = (int)((t_boundary[1] - t_boundary[0]) / delta_time);
    vector<vector<float> > vec(n_rows, vector<float>(n_cols));
    // vector<vector<float> > vec;

    explicit_solve(&vec, delta_x, delta_time);
    // vec = implicit_solve(time_length, delta_time, x_length, delta_x, initial_value, x_value);

    for(int i=0;i<vec.size();i++){
        printf("t=%f\t", i*delta_time);
        for(int j=0;j<vec[0].size();j++){
            printf("%f ", vec[i][j]);
        }
        printf("\n");
    }
    return 0;
}

void print(vector<vector<float>> *vec){
    int n_rows = (*vec).size();
    int n_cols = (*vec)[0].size();

}