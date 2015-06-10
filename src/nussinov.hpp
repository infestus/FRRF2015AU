#pragma once
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>

using namespace std;

class Nussinov
{
    ifstream rnafile;
    int n; 
    vector<int> * mat;
    vector<int> back;
    string seq; 
    uint score_keeper;

public:
    int fold_score(string sequence){
        return nussinov(sequence); 
    }

    int nussinov(string sequence) 
    {
        seq = sequence; 
        n = sequence.length();
        mat = new vector<int>(n*n, -1);
        int s = score(0, n-1);
        delete mat; 
        return s;
    };

    int score(int i, int j)
    {
        if (j-i < 2) {
            return 0;
        }
        if (mat->at(i*n+j) == -1) {
            int best = 0;
            int source = 0; 
            for (int k = i+1; k < j-1; k++){ 
                int temp = score(i, k)+score(k+1, j); //Divided step
                if (temp > best) {
                    best = temp;
                }
            }
            if (permitted_match(seq[i], seq[j]) && (score(i+1, j-1)+1 > best)) {
                best = score(i+1, j-1)+1; //Matching
            }
            if (score(i+1, j) > best) {
                best = score(i+1, j); // Unpaired I
            }
            if (score(i, j-1) > best) {
                best = score(i, j-1); // Unpaired J
            }
            mat->at(i*n+j) = best;
        }
        return mat->at(i*n+j);
    }

    bool permitted_match(char a, char b) {
        string t = string (1, a);
        t += b;
        bool res = (t == "AU" || t == "CG" || t == "GU" || t == "UA" || t == "GC" || t == "UG");
        return res;
    }

};