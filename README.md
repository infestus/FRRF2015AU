# FRRF2015AU
Sub-qubic Nussinov RNA secondary structure prediction algorithm using Four Russians technique

Compile using you preferred C++11 compiler. I.e. with Clang: 
clang++ -I "path_to_boost" -L "path_to_boost"/stage/lib -std=c++11 -O3 -DNDEBUG -w src/main.cpp -o bin/RNA;

Running the program:
Simply running the executable will run the function "experiments" in experiments.hpp.
4 arguments can also be given for manual control. 
RNA sequence output burn_in_binary q_scaling

sequence: The sequence to be folded. 
output: the filename the program will output to. Just put a filename, all files automatically goes in data-folder.
burn_in_binary: 1 or 0, for whether to use burn-ins or not. 
q_scaling: the q-value to be used by the program. Putting 0 will use the Nussinov algorith, putting -1 will select dynamic scaling using log_5(sequence length) as q. 

- Frederik 
