#include <thread>
#include <future>
#include <cmath>
#include <chrono>

int count_;

double ThreadA(int a, int b)
{
    double sum = 0;
    int i;
    for(i=a; i<b; ++i){
        sum += std::exp(((double)i)/((double)b));
    }
    return sum;
}


int main()
{
    std::chrono::system_clock::time_point start, end;
    // double result, result2;
    double r;
    std::future<double> result, result2;

    start = std::chrono::system_clock::now();

    // 何かの処理
    // std::future<double> result = std::async(std::launch::async, [] { return ThreadA(0,1000000); });
    // std::future<double> result2 = std::async(std::launch::async, [] { return ThreadA(0,1000000); });
    for(int i=0;i<100;i++){
        // result = ThreadA(0, 1000000);
        // result2 = ThreadA(0, 1000000);
        // result += result2;
        result = std::async(std::launch::async, [] { return ThreadA(0,1000000); });
        result2 = std::async(std::launch::async, [] { return ThreadA(0,1000000); });
        r = result.get() + result2.get();
    }
        

        // double r = result.get();
        // r += result2.get();

    end = std::chrono::system_clock::now();

    double time = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0);
    printf("time %lf[ms]\n", time/100.);

    

    printf("count_ : %lf, \n", r );


    return 0;
}