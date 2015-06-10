#include "nussinov.hpp"
#include "RNAFolding.hpp"
#include "experiments.hpp"
#include "string"
using namespace std; 

int main(int argc, char * argv[]){
    Experiments * exp = new Experiments();
    if (argc == 5){
        string seq = argv[1]; //Arg 1: Sequence
        string out = argv[2]; //Arg 2: Output file name
        bool burn_in = false; 
        string burn = argv[3]; //Arg 3: Burn in binary
        string param = argv[4]; //Arg 4: Greater than 0, static q value. Less than 0: Dynamic with b = 5
        if (burn == "True" || burn == "true"){
            burn_in = true; 
        }
        exp->single_n(seq, out, burn_in, stoi(param));
    } else if (argc == 2){
        string param = argv[1];
        if (param == "Test" || param == "test"){
            exp->test_basic("./sequences/testseqs_big.fasta");
        } else {
            exp->experiment();
        }
    } else {
        exp->experiment();
    }

    return 0; 
}