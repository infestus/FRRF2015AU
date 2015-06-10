#include <iostream>
#include "cstring"
#include "stdio.h"
#include "string"
#include "stdint.h"
using namespace std; 

int main(){
    int64_t bytes = 1000000000;
    int64_t i = 0; 
    int64_t* arr = new int64_t[bytes]; 
    memset(arr, 0, sizeof(int64_t)*(bytes));
    delete arr;  
}