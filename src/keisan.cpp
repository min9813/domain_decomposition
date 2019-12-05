# include "practice.h"
# include <iostream>

Keisan::Keisan(int _a, int _b){
    a = _a;
    b = _b;
    printf("call constructor \n");
}

Keisan::Keisan(Keisan & rother){
    a = rother.a;
    b = rother.b;

    printf("call constructor copy\n");
}

Keisan::~Keisan(){
    printf("call destructor\n");
}


int Keisan::add(){
    return a + b;
}