#ifndef _KEISAN_H_
#define _KEISAN_H_
 
class Keisan{
public:
    Keisan(int _a, int _b);
    Keisan(Keisan & rother);
    ~Keisan();
    int a;
    int b;
    int add();
};

# endif