#include <iostream>
#include <fstream>
#include <vector>
#include <string.h>
#include <sstream>
#include <map>

#include "Sequence.h"
#include <eigen-3.4.0/Eigen/Dense>
#include <eigen-3.4.0/Eigen/Eigenvalues>

using namespace std;

Alignment seqs_read(string file_name) {

    Alignment align;

    vector<Sequence> seqs;
    string line;
    string seq;
    string id;
    bool first_line = true;

    //map to count all the nucleotide freqs
    std::map <std::string , int> nucl_dict;
    std::string nucl;
	nucl_dict["A"] = 0;
	nucl_dict["C"] = 1;
	nucl_dict["G"] = 2;
	nucl_dict["T"] = 3;

    ifstream seq_file(file_name);
    if (!seq_file) {
        cout <<"Unable to open file";
        exit(1);
    }

    while(getline(seq_file, line)) {
        if (first_line) {
            first_line = false;
        }
        else {
            //initialise nucleotide frequencies
            Eigen::Vector4d nucl_freqs(0.0,0.0,0.0,0.0);

            istringstream s(line);
            s >> id >> seq;

             for(int i=0; i<seq.length();i++) {
                nucl = seq[i]; //type conversion?? char to string
                nucl_freqs(nucl_dict[nucl]) += 1;
            }

            Sequence id(seq, nucl_freqs);
            seqs.push_back(id);

            //add up the global freqs (for the whole alignment)
            for (int i=0; i<4; i++) {
                align.global_freqs(i) += nucl_freqs(i);
            }
        }
    }
    //set up the seq vector
    align.sequences = seqs;
    return align;

}
