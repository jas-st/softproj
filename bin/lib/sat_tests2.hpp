#pragma once

#include<iostream>
#include<string>
#include <eigen-3.4.0/Eigen/Dense>
#include <eigen-3.4.0/Eigen/Eigenvalues>

#include "Sequence.hpp"

double sat_test_cas1(Eigen::Matrix4d H, int n, Eigen::Vector4d freqs) {

    double stat = -1;
    double div;
    double crit_val;
    double nucl_sum = freqs.sum();

    // get relative frequences
    freqs = freqs/nucl_sum;

    for (int i = 0; i<4; i++) {
        div = n*freqs(i);
        stat += H(i,i)/div;
    }


    return stat;
}

double sat_test_cas2(Sequence Seq1, Sequence Seq2, Eigen::Matrix4d H, int n) {

    double freqs[4] = {};
    double stat = -1;
    double div;

    //relative frequences for the two sequences
    for (int i=0; i<4; i++) {
        freqs[i] = (Seq1.nucl_freqs_al[i]+Seq2.nucl_freqs_al[i])/(2*n);
    }

    for (int i = 0; i<4; i++) {
        div = n*freqs[i];
        stat += H(i,i)/div;
    }

    return stat;

}


double chi_test(Sequence Seq1, Sequence Seq2, Eigen::Matrix4d H, int n) {

    Eigen::Vector4d seq1_freqs_al = Seq1.nucl_freqs_al;
    Eigen::Vector4d seq2_freqs_al = Seq2.nucl_freqs_al;

    double stat = 0;
    double div;


    for(int i = 0; i<4; i++) {
        for(int j = 0; j<4; j++) {

            div = std::pow(H(i,j)-(seq1_freqs_al(i)*seq2_freqs_al(j)/n),2);
            stat += div/(seq1_freqs_al(i)*seq2_freqs_al(j)/n);
        }

    }

    return stat;

}
