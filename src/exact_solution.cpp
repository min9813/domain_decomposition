#include <vector>
#include <cmath>
#include <map>
#include <sstream> // digit adjust
#include <iomanip> //exit?
#include<iostream>
# include <prettyprint.hpp>
using namespace std;

const float MIN_EPS = 1e-5;
const float MAX_BOUND = 20.; 
const float ALPHA_MAX = 100.; 
const float GAMMA_MAX = 100000.; 



inline float cos_hypo(float value);
inline vector<float> cos_hypo(vector<float> values);
inline float sin_hypo(float value);
inline vector<float> sin_hypo(vector<float> values);
inline float hypo_div_stable(float value1, float value2);
inline float exp_sin_frac(float exp_value, float radian_x, float radian_y, float alpha, float beta);
inline float exp_sin_1d(float exp_value, float radian_x, float alpha);
float point_value(float x, float t, float coef, float x_bound);
float point_value_2d(float x, float y, float t, float coef, float x_bound, float y_bound);

vector<vector<float>> exact_solution_1d(vector<vector<float>> *field, map<string, pair<float, float>> spatio_bound, map<string, float> delta, float coef){
    int t_idx, x_idx;
    int n_rows = (*field).size();
    int n_cols = (*field)[0].size();
    float t_now, x_now, tmp;
    for (t_idx = 0; t_idx < n_rows; t_idx++)
    {
        t_now = t_idx * delta["t"];
        for (x_idx = 0; x_idx < n_cols; x_idx++)
        {
            if(x_idx == 0 || x_idx==n_cols-1){
                (*field)[t_idx][x_idx] = 0;
            }else if(t_idx==0){
                (*field)[t_idx][x_idx] = 1;
            }else{
                x_now = x_idx * delta["x"] + spatio_bound["x"].first;
                tmp = point_value(x_now, t_now, coef, spatio_bound["x"].second);
                (*field)[t_idx][x_idx] = tmp;
            }
            
        }
    }
    return (*field);
}

vector<vector<vector<float>>> exact_solution_2d(vector<vector<vector<float>>> *field, map<string, pair<float, float>> spatio_bound, map<string, float> delta, float coef)
{
    // cout << "sol 2d" << endl;
    int t_idx, x_idx, y_idx;
    int n_rows = (*field).size();
    int n_cols1 = (*field)[0].size();
    int n_cols2 = (*field)[0][0].size();
    float t_now, x_now, y_now, tmp;
    for (t_idx = 0; t_idx < n_rows; t_idx++)
    {
        t_now = t_idx * delta["t"];
        for (y_idx = 0; y_idx < n_cols1; y_idx++)
        {
            y_now = y_idx * delta["y"] + spatio_bound["y"].first;
            for (x_idx = 0; x_idx < n_cols2; x_idx++)
            {
                if(x_idx == 0 || x_idx==n_cols2-1 || y_idx==0 || y_idx == n_cols1 - 1){
                    (*field)[t_idx][y_idx][x_idx] = 0;
                }else if(t_idx==0){
                    (*field)[t_idx][y_idx][x_idx] = 1;
                }else{
                    x_now = x_idx * delta["x"] + spatio_bound["x"].first;
                    tmp = point_value_2d(x_now, y_now, t_now, coef, spatio_bound["x"].second, spatio_bound["y"].second);
                    (*field)[t_idx][y_idx][x_idx] = tmp;
                }
                
            }
        }
        // if(t_idx>0){
            // exit(1);
        // }
    }

    return *field;
}

inline float cos_hypo(float value)
{
    return (exp(value) + exp(-value)) / 2.;
}

inline vector<float> cos_hypo(vector<float> values)
{
    vector<float> new_values(values.size());
    for (int i = 0; i < values.size(); i++)
    {
        new_values[i] = (exp(values[i]) + exp(-values[i])) / 2.;
    }

    return new_values;
}

inline float sin_hypo(float value)
{
    return (exp(value) - exp(-value)) / 2.;
}

inline vector<float> sin_hypo(vector<float> values)
{
    vector<float> new_values(values.size());
    for (int i = 0; i < values.size(); i++)
    {
        new_values[i] = (exp(values[i]) - exp(-values[i])) / 2.;
    }

    return new_values;
}

inline float hypo_div_stable(float value1, float value2)
{
    return exp(value1 - value2);
}

inline float exp_sin_frac(float exp_value, float radian_x, float radian_y, float alpha, float beta)
{
    float bunshi = exp(exp_value) * sin(radian_x) * sin(radian_y);
    float bunbo = alpha * beta;

    return bunshi / bunbo;
}

inline float exp_sin_1d(float exp_value, float radian_x, float alpha)
{
    float bunshi = exp(exp_value) * sin(radian_x);

    return bunshi / alpha;
}

float point_value(float x, float t, float coef, float x_bound){
    int i;
    float s=0, alpha=0, radian_x= - M_PI * x, exp_value, tmp;
    cout << "----------------------------" << "\n";
    for(i=1;i<200;i+=2){
        alpha = i * M_PI/x_bound;
        if(t>0){
            exp_value = - coef * alpha * t;
            exp_value *= alpha;
        }else{
            exp_value = 0;
        }
        
        if(abs(exp_value) > MAX_BOUND) break;
        radian_x += 2 * M_PI * x;
        if (radian_x > 2*M_PI) radian_x -= 2 * M_PI;
        tmp = exp_sin_1d(exp_value, radian_x, alpha);
        s += 4 * tmp;
        // cout << "t="<<t <<" x="<<x<<" i="<<i<<" alpha="<<alpha<<" exp_value="<<exp_value<<" radx="<<radian_x<<" tmp="<<tmp << " sum="<<s<<"\n";
        if(1./alpha < MIN_EPS) break;
    }

    return s;

}

float point_value_2d(float x, float y, float t, float coef, float x_bound, float y_bound)
{
    int i, j, index;
    float s = 0, tmp, exp_value, alpha_tmp, alpha_square_tmp, radian_x_tmp, beta, beta_square;
    float radian_y = -M_PI * y, gamma;
    vector<float> alpha, alpha_square, radian_x;
    // cout << "exact calculate" << endl;
    for (i = 1;i<100; i += 2)
    {
        beta = i * M_PI / y_bound;
        beta_square = pow(beta, 2);
        radian_y += 2 * M_PI * y;
        if (radian_y > 2*M_PI) radian_y -= 2 * M_PI;
        radian_x_tmp = - M_PI * x;
        for (j = 1;j<100; j += 2)
        {
            index = (j-1)/ 2;
            if (index >= alpha.size())
            {
                alpha_tmp = j * M_PI / x_bound;
                alpha_square_tmp = pow(alpha_tmp, 2);
                radian_x_tmp += 2*M_PI * x;
                if (radian_x_tmp > 2*M_PI) radian_x_tmp -= 2 * M_PI;
                alpha.push_back(alpha_tmp);
                alpha_square.push_back(alpha_square_tmp);
                radian_x.push_back(radian_x_tmp);
            }
            else
            {
                alpha_tmp = alpha[index];
                alpha_square_tmp = alpha_square[index];
                radian_x_tmp = radian_x[index];
            }
            if(t>0){
                gamma = pow(alpha_tmp, 2) + beta_square;
                if(gamma > GAMMA_MAX) break;     
                exp_value = - coef * gamma * t;
            }else{
                exp_value = 0;
            }
            if(abs(exp_value) > MAX_BOUND) break;

            tmp = exp_sin_frac(exp_value, radian_x_tmp, radian_y, alpha_tmp, beta);
            s += tmp;
            // cout << "t="<<t <<" x="<<x<<" y="<<y <<" i="<<i<<" j="<<j<<" alpha="<<alpha_tmp<<" beta="<<beta<< " alpha*beta="<<1./(alpha_tmp*beta)<<" exp_value="<<exp_value<<" radx="<<radian_x_tmp<<" rady="<<radian_y<<" tmp="<<tmp << " sum="<<s<<"\n";

            if(1./(alpha_tmp*beta) < MIN_EPS) break;

        }
    }

    s *= 16. / (x_bound * y_bound);
    // cout << " s = " <<s <<"\n";
    // cout << "--------------------------------"<<endl;

    return s;
}

string float_to_string(float value, int precision){
    stringstream ss;
    ss << fixed << setprecision(precision) << value;
    string str = ss.str();
    return str;
}
