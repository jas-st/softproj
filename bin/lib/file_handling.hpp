#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <string.h>
#include <set>

#include "Sequence.hpp"
#include <eigen-3.4.0/Eigen/Dense>
#include <eigen-3.4.0/Eigen/Eigenvalues>

//function to calculate absolute nucleotide frequences of a sequence
Eigen::Vector4d get_nucleotide_frequencies(std::string seq) {
    //map to count all the nucleotide freqs
    std::map <std::string , int> nucl_dict;
    //nucleotide frequencies vector
    Eigen::Vector4d nucl_freqs(0.0,0.0,0.0,0.0);
    std::string nucl;

	nucl_dict["A"] = 0;
	nucl_dict["C"] = 1;
	nucl_dict["G"] = 2;
	nucl_dict["T"] = 3;

    for(int i=0; i<seq.length();i++) {
        if (seq[i] != '-') {
            nucl = seq[i]; //type conversion?? char to string
            nucl_freqs(nucl_dict[nucl]) += 1;
        }
    }

    return nucl_freqs;
}

//function to find illegal characters in a string
bool check_illegal_char(char str_char) {
    //returns true if the char is illegal so it can be replaced
    if (!(isdigit(str_char) || isalpha(str_char) || (str_char == '.') || (str_char == '_') || (str_char == '-'))) {
        return true;
    }
    else { return false; }
}

//main function to read in the sequences
Alignment seqs_read(std::string file_name, std::string extension) {

    Alignment align;

    //declare variables
    std::vector<Sequence> seqs;
    std::string line;
    std::string seq;
    std::string id;
    Eigen::Vector4d nucl_freqs;
    std::string new_seq;
    int counter_phy = 0;
    int counter_nex = 0;

    //exit if file is not found
    std::ifstream seq_file(file_name);
    if (!seq_file) {
        exit(1);
    }

    while(getline(seq_file, line)) {
        //end the reading if a blank line or end is read
        std::size_t first = line.find_first_not_of( " \f\n\r\t\v;" );
        std::size_t end = line.find("end;");

        //stop reading if blank line/ or end/ or ;
        if((first == line.npos) || (end != line.npos)) {
                break;
        }

        if (strcmp(extension.c_str(), "nex") == 0 ) {
            if (counter_nex < 6) {
                counter_nex += 1;
            }
        } else {
            if (counter_phy < 2) {
                counter_phy += 1;
            }
        }

        if ((counter_nex == 6) | (counter_phy == 2)) {
            std::istringstream s(line);
            s >> id >> seq;

            //remove '' or "" from begining of the identifier
            id.erase(std::remove(id.begin(), id.end(), '\''), id.end());
            id.erase(std::remove(id.begin(), id.end(), '\"'), id.end());

            //remove illegal characters (only allow alpha/num/.-_)
            replace_if(id.begin(), id.end(), check_illegal_char, '_');

            //get nucleotide frequencies
            nucl_freqs = get_nucleotide_frequencies(seq);

            new_seq = id;
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
