#pragma once
#include <string>
#include <iterator>
#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <cmath>
#include <map>
#include <sstream>
#include <boost/dynamic_bitset.hpp>
#include <cstdint>
#include <climits>
using namespace std;

class RNAFold
{
private:
    int n; 
    vector<int> * mat;
    string seq; 
    int score_match;
    int q_sq = 0; 
    int q = 0; 
    int * R; 

public:
    int fold_score(string seq, int size){
        this->seq = seq; 
        n = seq.length();
        mat = new vector<int>(n*n, 0);
        russian(seq, size);
        int result = mat->at(0*n+(n-1));
        delete mat; 
        return result; 
    }

    int fold_score_dynamic(string seq, int base = 5){
        this->seq = seq; 
        n = seq.length(); 
        mat = new vector<int>(n*n, 0);
        int group_size = round(logbn(n, base));
        russian(seq, group_size);
        int result = mat->at(0*n+(n-1));
        delete mat; 
        return result; 
    }

    void russian(string sequence, int size){
        q = size;
        q_sq = pow(2, q-1); 
        int max_q = static_cast<int>(ceil(n/q))+1; 
        int width, height, depth; 
        width = n; 
        height = max_q; 
        depth = q_sq; 
        uint64_t R_size = width * height * depth;
        if (R_size > INT_MAX){
            cerr << "Q-value: (" << q << ") too large. R-tabel would increase max_int in size" << endl;
            exit(1);
        }
        R = new int[width * height * depth];
        vector<boost::dynamic_bitset<> > vg_store(max_q, boost::dynamic_bitset<>(0) );
        for (int j = 2; j < n; j++) {
            // Independant
            for (int i = 0; i < (j-1); i++) {
                mat->at(i*n+j) = max(mat->at((i+1)*n+(j-1))+permitted_match(i, j), mat->at(i*n+(j-1))); // Rules A, and B    
            }
            // Dependant
            for (int i = j-2; i >= 0; i--){
                mat->at(i*n+j) = max(mat->at((i+1)*n+j), mat->at(i*n+j)); // Rules C is missing from pseudo-code. Added here. 
                int k_prime = -1;
                int k_prime_val = -1;
                int g; 
                for(int z = j-1; z > i+1; ){
                    int k_prime_temp = -1;
                    int k_prime_temp_val = -1;
                    int z_mod = z%q;
                    if (z_mod==q-1 && z-(i+1)>=q-1) {
                        g = floor(z/q);
                        k_prime_temp = R[i * height * depth + g * depth + vg_store[g].to_ulong()];
                        k_prime_temp_val = mat->at(i*n+(k_prime_temp-1)) + mat->at(k_prime_temp*n+j); 
                        z = z-q; 
                    } else if (z - next_group(z, z_mod) < q && i<next_group(z, z_mod)) {
                        int top_score = -1;
                        int temp = 0; 
                        for (int aux_z = z; aux_z > next_group(z, z_mod); aux_z--) {
                            temp = mat->at(i*n+(aux_z-1)) + mat->at(aux_z*n+j); 
                            if (top_score < temp) {
                                top_score = temp;
                                k_prime_temp = aux_z; 
                                k_prime_temp_val = temp; 
                            }
                        }
                        z = next_group(z, z_mod); 
                    } else if (z-(i+1) < (q-1)) {
                        int top_score = -1;
                        int temp = 0; 
                        for (int aux_z = z; aux_z > i+1; aux_z--) {
                            temp = mat->at(i*n+(aux_z-1)) + mat->at(aux_z*n+j); 
                            if (top_score < temp) {
                                top_score = temp; 
                                k_prime_temp = aux_z;
                                k_prime_temp_val = temp; 
                            }
                        }
                        z = i;
                    }
                    if (k_prime_temp > 0){
                        // k_prime_temp_val = mat->at(i*n+(k_prime_temp-1)) + mat->at(k_prime_temp*n+j); 
                        if (k_prime > 0){
                            // k_prime_val = mat->at(i*n+(k_prime-1)) + mat->at(k_prime*n+j);  
                            if (k_prime_temp_val > k_prime_val){
                                k_prime = k_prime_temp; 
                                k_prime_val = k_prime_temp_val; 
                            }
                        } else {
                            k_prime = k_prime_temp; 
                            k_prime_val = k_prime_temp_val; 
                        }
                    }
                }
                if (k_prime > 0) { 
                    mat->at(i*n+j) = max(mat->at(i*n+j), mat->at(i*n+(k_prime-1)) + mat->at(k_prime*n+j));
                }
                if ((i-1)% q == q-1) {
                    g = floor(i/q);
                    boost::dynamic_bitset<> vg(q-1, 0);
                    vg = encode(j, g, q);
                    vg_store[g] = vg; 
                }
            }
            // Table
            if (j%q == q-1){ //Cgroup j/q is complete. 
                int g = floor(j/q); // G is a Rgroup 
                vector<int> V_prime;
                int foo = 0; 
                for (int vs = 0; vs<q_sq; vs++) {
                    boost::dynamic_bitset<> v(q-1, vs);  
                    V_prime = decode(v, q);
                    for(int i = 0; i < j-1; i++){
                        // g = floor(i/q);
                        int r_max_val = -1;
                        int r_max_idx = -1; 
                        foo = 0;
                        for(int k = (j-(q-1)); k <= j; k++){
                            if (k != 0 && r_max_val < mat->at(i*n+(k-1)) + V_prime[foo]){
                                r_max_val = mat->at(i*n+(k-1)) + V_prime[foo]; 
                                r_max_idx = k; 
                            }
                            foo++; 
                        }
                        R[i * height * depth + g * depth + vs] = r_max_idx;
                    }
                }
            }
        }
        delete[] R; 
    };

    vector<int> decode(boost::dynamic_bitset<> v, int q){
        vector<int> res(q);
        res[0] = 0; 
        for (int i = 1; i<res.size(); i++) {
            res[i] = (res[i-1]+v[i-1]);
        }
        reverse(res.begin(), res.end()); 
        return res; 
    };

    boost::dynamic_bitset<> encode(int j, int g, int q){
        boost::dynamic_bitset<> v(q-1);
        int v_idx = 0;
        int val = 0; 
        for(int i = min(n-2, (g * q)+(q-2)); i>=g*q; i--){
            val = mat->at(i*n+j) - mat->at((i+1)*n+j);
            v[v_idx] = val; 
            v_idx++; 
        }
        return v; 
    };

    int permitted_match(int i, int j) {
        string t = string (1, seq[i]);
        t += seq[j];
        int res = (t == "AU" || t == "CG" || t == "GU" || t == "UA" || t == "GC" || t == "UG") ? 1 : 0;
        return res;
    }

    int next_group(int i, int z_mod) {
        return (i-(z_mod)-1);
    }

    double logbn(double x, double b){
        return log10(x)/log10(b); 
    }
};  