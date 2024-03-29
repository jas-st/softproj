#pragma once
#include <string>
#include <vector>

#include <eigen-3.4.0/Eigen/Dense>
#include <eigen-3.4.0/Eigen/Eigenvalues>

class Sequence {
    public:
        std::string id;
        std::string seq;
        double length;

        //nucleotide frequencies for the specific sequence {A,C,G,T}
        Eigen::Vector4d nucl_freqs = Eigen::Vector4d::Zero();

        //nucleotide frequncies only looking at the alignment {A,C,G,T}
        Eigen::Vector4d nucl_freqs_al = Eigen::Vector4d::Zero();

        //constructor for initialization
        Sequence(std::string seq_id, std::string seq_str, Eigen::Vector4d seq_freq) {
            id = seq_id;
            seq = seq_str;
            length = seq_freq.sum();
            nucl_freqs = seq_freq;
        }

};

class Alignment {
    public:
        //the vector with all sequences
        std::vector<Sequence> sequences;

        //get nucleotide frequencies for the whole alignment {A,C,G,T}
        Eigen::Vector4d global_freqs = Eigen::Vector4d::Zero();

        //default constructor
        Alignment() {
        }
};
