#pragma once
#include <stdio.h>
#include <string>
#include <cstring>
#include <iostream>
#include <fstream>
#include "nussinov.hpp"
#include "RNAFolding.hpp"
#include "sys/types.h"
#include <vector>
#include <cmath>
#include <algorithm> 
#include <chrono>
#include <ctime>
#include <ratio>

using namespace std;

class Experiments
{
public: 
    bool burn_in = false; 

    void experiment(){
        // linear_n("uneven2k.fasta", "uneven2k.dat", 1, 9, false, 2);
        // linear_n("unique2k.fasta", "all_unique2k(2).dat", 1, 9, false, 2);
        // linear_n("unique6k.fasta", "all_unique6k(2).dat", 1, 7, false);
        // linear_n("unique35h.fasta", "all_unique35h(2).dat", 1, 7, false, 3);
        // dynamic_n("unique6k_extended.fasta", "dynamic_unique6k.dat");
        linear_n("unique6k_extended.fasta", "linear_unique6k.dat");
    };

    void linear_n(string filename, string output, int max_iter = 1, int max_q = 9, bool v = false, int min_q = 2){
        filename = "./sequences/" + filename; 
        ofstream out_stream;
        out_stream.open ("./data/" + output);
        vector< pair<string, string> > sequences = readFasta(filename);
        int burn_in_iter = 0;
        if (burn_in){
            burn_in_iter = warm_up();
        }
        int max_n = sequences.back().second.length();
        out_stream << "n;q;ms;burn_in_iter:" << burn_in_iter << ";max_q:" << max_q << ";min_q:" << min_q << ";iter:" << max_iter << ";max_n:" << max_n << endl; 
        if (v){
            verify_all(filename, max_q, min_q);
        }
        cout << "Running on setting: " << filename << ", using " << max_iter << "iterations ";
        cout << "with q ranging from " << min_q << " to " << max_q; 
        for (int j = 0; j < sequences.size(); j++) {
            for (int i = 0; i < max_iter; i++) {
                for (int q = min_q; q <= max_q; q++){
                    cout << "Running: [q:" << q << ",iter:" << i << ",seqId:" << j << ",seqlen:" << sequences[j].second.length() << "]" << endl;
                    out_stream << sequences[j].second.length() << ";" << q << ";"; 
                    out_stream << timeRussian(sequences[j].second, q) << endl; 
                }
                cout << "Running: [q: 0" << ",iter:" << i << ",seqId:" << j << ",seqlen:" << sequences[j].second.length() << "]" << endl;
                out_stream <<  sequences[j].second.length() << ";0;" << timeNussinov(sequences[j].second) << endl;  
                cout << "Running: [q:-1" << ",iter:" << i << ",seqId:" << j << ",seqlen:" << sequences[j].second.length() << "]" << endl;
                out_stream <<  sequences[j].second.length() << ";-1;" << timeDynamic(sequences[j].second) << endl;  
            }
        }
        out_stream.close(); 
    }

    void dynamic_n(string filename, string output, int max_iter = 1, int max_b = 8) {
        const int min_b = 3;
        filename = "./sequences/" + filename; 
        ofstream out_stream;
        out_stream.open ("./data/" + output);
        vector< pair<string, string> > sequences = readFasta(filename);
        int burn_in_iter = 0;
        if (burn_in){
            burn_in_iter = warm_up();
        }
        int max_n = sequences.back().second.length();
        out_stream << "n;q;ms;burn_in_iter:" << burn_in_iter << ";max_b:" << max_b << ";min_b:" << min_b << ";iter:" << max_iter << ";max_n:" << max_n << endl; 
        cout << "Running on setting: " << filename << ", using " << max_iter << "iterations ";
        cout << "with q ranging from " << min_b << " to " << max_b; 
        for (int j = 0; j < sequences.size(); j++) {
            for (int i = 0; i < max_iter; i++) {
                for (int b = min_b; b <= max_b; b++){
                    cout << "Running: [b:" << b << ",iter:" << i << ",seqId:" << j << ",seqlen:" << sequences[j].second.length() << "]" << endl;
                    out_stream <<  sequences[j].second.length() << ";" << b << ";" << timeDynamic(sequences[j].second, b) << endl;  
                }
                cout << "Running: [nussinov" << ",iter:" << i << ",seqId:" << j << ",seqlen:" << sequences[j].second.length() << "]" << endl;
                out_stream <<  sequences[j].second.length() << ";0;" << timeNussinov(sequences[j].second) << endl;  
            }
        }
        out_stream.close(); 
    }

    void single_n(string sequence, string output, bool do_warm_up = true, int q = 5){
        ofstream out_stream;
        out_stream.open ("./data/" + output, ofstream::out | ofstream::app);
        int burn_in_iter = 0;
        if (do_warm_up){
            burn_in_iter = warm_up();
        }
        if (q == 0){
            cout << "Running: [q:" << q << ",seqlen:" << sequence.length() << "]" << endl;
            out_stream <<  sequence.length() << ";0; burn_in:" << do_warm_up << ";" << timeNussinov(sequence) << endl;  
        } else if (q == -1){
            cout << "Running: [q:" << q << ",seqlen:" << sequence.length() << "]" << endl;
            out_stream <<  sequence.length() << ";-1; burn_in:" << do_warm_up << ";" << timeDynamic(sequence) << endl;  
        } else {
            cout << "Running: [q:" << q << ",seqlen:" << sequence.length() << "]" << endl;
            out_stream << sequence.length() << ";" << q << "; burn_in:" << do_warm_up << ";" << timeRussian(sequence, q) << endl; 
        }
        out_stream.close();    
    }

    int warm_up(int iter = 2){
        vector< pair<string, string> > sequences = readFasta("./sequences/unique2k.fasta");
        string seq = sequences.back().second; 
        RNAFold * russian = new RNAFold(); 
        Nussinov * nussinov = new Nussinov(); 
        for (int i = 0; i < iter; i++){
            nussinov->fold_score(seq);
            russian->fold_score(seq, 5); 
        }
        return iter;
    }

    void verify_all(string filename, int max_q, int min_q) {
        cerr << "Starting verification..." << endl; 
        vector< pair<string, string> > sequences = readFasta(filename);
        RNAFold * russian = new RNAFold();
        Nussinov * nussinov = new Nussinov();

        for (int j = 0; j < sequences.size(); j++) {
            if (sequences[j].second.length() > 999){
                int nus_cpp = nussinov->fold_score(sequences[j].second);
                for (int q = min_q; q <= max_q; q++){
                    int rus_cpp = russian->fold_score(sequences[j].second, q);
                    if (nus_cpp != rus_cpp){
                        cerr << "Values did not match" << endl;
                        cerr << "[FRM;NUS]" << endl;
                        cerr << "[" << rus_cpp << ";" << nus_cpp << "]" << endl;
                    }
                }
            } else {
                int nus_py = python_result(sequences[j].second);
                int nus_cpp = nussinov->fold_score(sequences[j].second);
                for (int q = min_q; q <= max_q; q++){
                    int rus_cpp = russian->fold_score(sequences[j].second, q);
                    if ((nus_cpp != nus_py) || (nus_py != rus_cpp)){
                        cerr << "Values did not match" << endl;
                        cerr << "[FRM;NUS;NPY]" << endl;
                        cerr << "[" << rus_cpp << ";" << nus_cpp << ";" << nus_py << "]" << endl;
                    }
                }
            }
        }
        cerr << "Verification complete " << endl; 
    }

    double timeRussian(string sequence, int q){
        int russian_result = 0; 
        RNAFold * russian = new RNAFold(); 
        chrono::high_resolution_clock::time_point start = chrono::high_resolution_clock::now();
        russian_result = russian->fold_score(sequence, q);
        chrono::high_resolution_clock::time_point end = chrono::high_resolution_clock::now();
        delete russian; 
        chrono::duration<double> diff = chrono::duration_cast< chrono::duration <double> > (end - start);
        return diff.count(); 
    }

    double timeNussinov(string sequence){
        int nussinov_result = 0; 
        Nussinov * nussinov = new Nussinov(); 
        chrono::high_resolution_clock::time_point start = chrono::high_resolution_clock::now();
        nussinov_result = nussinov->fold_score(sequence);
        chrono::high_resolution_clock::time_point end = chrono::high_resolution_clock::now();
        delete nussinov; 
        chrono::duration<double> diff = chrono::duration_cast< chrono::duration <double> > (end - start);
        return diff.count(); 
    }

    double timeDynamic(string sequence, int b = 5){
        int russian_result = 0; 
        RNAFold * russian = new RNAFold(); 
        chrono::high_resolution_clock::time_point start = chrono::high_resolution_clock::now();
        russian_result = russian->fold_score_dynamic(sequence, b);
        chrono::high_resolution_clock::time_point end = chrono::high_resolution_clock::now();
        delete russian; 
        chrono::duration<double> diff = chrono::duration_cast< chrono::duration <double> > (end - start);
        return diff.count();    
    }

    void time_fasta(string filename, string output){
        ofstream out_stream;
        filename = "./sequences/" + filename; 
        output = "./data/" + output;
        out_stream.open (output);
        vector< pair<string, string> > sequences = readFasta(filename);
        const int max_iter = 5; 
        const int max_q = 10;
        const int min_q = 2; 
        double average_nussinov_time = 0; 
        out_stream << "n;Nussinov;"; 
        for (int q = min_q; q<=max_q; q++){
            out_stream << "Russian_" << q << ";"; 
        }
        out_stream << endl; 
        vector<double> * average_russian_time = new vector<double>(max_q+1, 0); 
        for (int seqId = 0; seqId < sequences.size(); seqId++){
            for (int iteration = 0; iteration <= max_iter; iteration++){
                for (int q = min_q; q<=max_q; q++){
                    average_russian_time->at(q) += timeRussian(sequences[seqId].second, q)/max_iter; 
                }
                average_nussinov_time += timeNussinov(sequences[seqId].second)/max_iter;
            }
            out_stream << sequences[seqId].second.length() << ";" << average_nussinov_time << ";"; 
            for (int q = min_q; q<=max_q; q++){
                out_stream << average_russian_time->at(q) << ";";
            }
            out_stream << endl;
        }
        delete average_russian_time; 
        out_stream.close();
    };

    void test_basic(string foo = "./sequences/testseqs.fasta", int q = 5){
        cerr << "test_basic(" << foo << "," << q << ");" << endl; 
        vector< pair<string, string> > sequences;
        sequences = readFasta(foo); 
        for (int seqId = 0; seqId < sequences.size(); seqId++){
            string sequence = sequences[seqId].second; 
            Nussinov * nussi = new Nussinov();
            RNAFold * russian = new RNAFold(); 
            int nussinov_result = 0; 
            int russian_result = 0;
            nussinov_result = nussi->fold_score(sequence);
            russian_result = russian->fold_score(sequence, q);
            if (nussinov_result != russian_result || (nussinov_result != python_result(sequence))){
                cerr << "Values doesn't match on sequence: " << seqId << ": " << sequences[seqId].first << endl; 
                cerr << "Nussinov: " << nussinov_result << ", Russians: " << russian_result << endl; 
                ppython_result(sequence); 
                cout << endl; 
                exit(1); 
            }
            delete nussi; 
            delete russian; 
            cout << "Sequence: " << sequences[seqId].first << endl << sequence << "  completed" << endl; 
        }
    };



    void test_q_scaling(){
        vector< pair<string, string> > sequences = readFasta("./sequences/testseqs.fasta");
        for (int seqId = 0; seqId < sequences.size(); seqId++){
            string sequence = sequences[seqId].second; 
            RNAFold * russian = new RNAFold(); 
            int python_res = python_result(sequence); 
            int log_3 = (int) (log(sequence.length())/log(3)); 
            for (int l = 2; l <= 10; l++){
                if (python_res != russian->fold_score(sequence, l)){
                    cerr << "Values doesn't match on sequence: " << sequence << endl; 
                    cerr << "Python: " << python_res << ", Russians: " << russian->fold_score(sequence, l) << "; with q = " << l << endl; 
                    exit(1); 
                } else {
                    cout << "Sequence: " << sequences[seqId].first << endl << sequence << "  completed" << endl; 
                }
            }
        }
    }

    void test(string sequence, int q = 5){
        Nussinov * nussi = new Nussinov();
        RNAFold * russian = new RNAFold(); 
        int nussinov_result = 0; 
        int russian_result = 0; 
        nussinov_result = nussi->fold_score(sequence);
        russian_result = russian->fold_score(sequence, q);
        if (nussinov_result != russian_result || (nussinov_result != python_result(sequence))){
                cerr << "Values doesn't match on sequence: " << sequence << endl; 
                cerr << "Nussinov: " << nussinov_result << ", Russians: " << russian_result << endl; 
                ppython_result(sequence); 
                exit(1); 
        }
    }

    vector< pair<string, string> > readFasta(string filename){
        int numberOfSequences = 0; 
        int seqId = -1; 
        ifstream ifs(filename);
        string temp; 
        while(getline(ifs, temp)){
            if (temp[0] == '>'){
                numberOfSequences++; 
            }
        }
        ifs.clear(); 
        ifs.seekg(0, ifs.beg);
        vector< pair<string, string> > sequences(numberOfSequences, make_pair("", "")); 
        while(getline(ifs, temp)) {
            if (temp[0] == '>') {
                seqId++;
                sequences[seqId] = make_pair(temp.substr(1), "");
            } else {
                upperCase(temp); 
                sequences[seqId].second += temp; 
            }
        }
        return sequences;    
    }

    int python_result(string sequence){
        try {
            string command = "./src/nussinov_single.py " + sequence; 
            FILE * python_response = popen(command.c_str(), "r"); 
            if (!python_response){
                cerr << "Python program didn't cuz errors"; 
            }

            char buffer[256];
            char *python_line = fgets(buffer, sizeof(buffer), python_response); 
            pclose(python_response); 
            string resp(buffer); 
            std::string::size_type sz; 
            int foo = stoi(resp);
            return foo; 
        } catch (...) {
            cerr << "couldn't run python program" << endl; 
        }
    };

    void ppython_result(string sequence){
        try {
            string command = "./src/nussinov_single.py " + sequence; 
            FILE * python_response = popen(command.c_str(), "r"); 
            if (!python_response){
                cerr << "Python program didn't cuz errors"; 
            }

            char buffer[256];
            char *python_line = fgets(buffer, sizeof(buffer), python_response); 
            pclose(python_response); 
            string resp(buffer); 
            std::string::size_type sz; 
            int foo = stoi(resp);
            cerr << "Correct answer is: "; 
            cerr << foo << endl; 
        } catch (...) {
            cerr << "couldn't run python program" << endl; 
        }; 
    }

    void upperCase(string &s){
        transform(s.begin(), s.end(), s.begin(), ::toupper); 
    }
};
