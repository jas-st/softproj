#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <string.h>

#include "Sequence.h"
#include <eigen-3.4.0/Eigen/Dense>
#include <eigen-3.4.0/Eigen/Eigenvalues>

Alignment seqs_read_bio(std::string file_name) {

    Alignment align;

    std::vector<Sequence> seqs;
    std::string line;
    std::string seq;
    std::string id;
    std::string new_seq;
    int counter = 0;


    //map to count all the nucleotide freqs
    std::map <std::string , int> nucl_dict;
    std::string nucl;
	nucl_dict["A"] = 0;
	nucl_dict["C"] = 1;
	nucl_dict["G"] = 2;
	nucl_dict["T"] = 3;


    std::ifstream seq_file(file_name);
    if (!seq_file) {
        std::cout <<"Unable to open file" << std::endl;
        exit(1);
    }

    while(getline(seq_file, line)) {
        if (counter < 5) {
            counter += 1;
        }
        else if(strcmp(line.c_str(), ";") == 0 | line.empty()) {
            break;
        }
        else {
            //nucleotide frequencies
            Eigen::Vector4d nucl_freqs(0.0,0.0,0.0,0.0);

            std::istringstream s(line);
            s >> id >> seq;
            new_seq = id;

            for(int i=0; i<seq.length();i++) {
                if (seq[i] != '-') {
                    nucl = seq[i]; //type conversion?? char to string
                    nucl_freqs(nucl_dict[nucl]) += 1;
                }
            }



            Sequence new_seq(id, seq, nucl_freqs);
            seqs.push_back(new_seq);


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

Alignment seqs_read(std::string file_name) {

    Alignment align;

    std::vector<Sequence> seqs;
    std::string line;
    std::string seq;
    std::string id;
    bool first_line = true;

    //map to count all the nucleotide freqs
    std::map <std::string , int> nucl_dict;
    std::string nucl;
	nucl_dict["A"] = 0;
	nucl_dict["C"] = 1;
	nucl_dict["G"] = 2;
	nucl_dict["T"] = 3;

    std::ifstream seq_file(file_name);
    if (!seq_file) {
        std::cout <<"Unable to open file";
        exit(1);
    }

    while(getline(seq_file, line)) {
        if (first_line) {
            first_line = false;
        }
        else {
            //initialise nucleotide frequencies
            Eigen::Vector4d nucl_freqs(0.0,0.0,0.0,0.0);

            std::istringstream s(line);
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
