#pragma once

#include<iostream>
#include<string>
#include <eigen-3.4.0/Eigen/Dense>
#include <eigen-3.4.0/Eigen/Eigenvalues>

#include "Sequence.h"

double sat_test_cas1(Eigen::Matrix4d H, int n, Eigen::Vector4d freqs) {
    //z_val for one sided confidence interval for alpha=0.05
    double z_val = 1.65;

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

    crit_val = sqrt(3/n) * z_val;

    return stat;
}

double sat_test_cas2(Sequence Seq1, Sequence Seq2, Eigen::Matrix4d H, int n) {
    //z_val for one sided confidence interval for alpha=0.05
    double z_val = 1.65;

    double freqs[4] = {};
    double stat = -1;
    double div;
    double crit_val;

    //relative frequences for the two sequences
    for (int i=0; i<4; i++) {
        freqs[i] = (Seq1.nucl_freqs[i]+Seq2.nucl_freqs[i])/(Seq1.length + Seq2.length);
    }

    for (int i = 0; i<4; i++) {
        div = n*freqs[i];
        stat += H(i,i)/div;
    }

    crit_val = sqrt(3/n) * z_val;

    return stat;

}

double chi_test(Sequence Seq1, Sequence Seq2, Eigen::Matrix4d H, int n) {

    Eigen::Vector4d seq1_freqs = Seq1.nucl_freqs/Seq1.length;
    Eigen::Vector4d seq2_freqs = Seq2.nucl_freqs/Seq2.length;

    double stat = 0;
    double div;
    double crit_val = 19.023;


    for(int i = 0; i<4; i++) {
        for(int j = 0; j<4; j++) {
            div = std::pow(H(i,j)-n*seq1_freqs(i)*seq2_freqs(j),2);
            stat += div/(n*seq1_freqs(i)*seq2_freqs(j));

        }

    }

    return stat;

}

double chi_test_draft(Sequence Seq1, Sequence Seq2, Eigen::Matrix4d H, int n) {

    Eigen::Vector4d seq1_freqs_al = Seq1.nucl_freqs_al;
    Eigen::Vector4d seq2_freqs_al = Seq2.nucl_freqs_al;

    double stat = 0;
    double div;
    double crit_val = 19.023;


    for(int i = 0; i<4; i++) {
        for(int j = 0; j<4; j++) {

            div = std::pow(H(i,j)-(seq1_freqs_al(i)*seq2_freqs_al(j)/n),2);
            stat += div/(seq1_freqs_al(i)*seq2_freqs_al(j)/n);
        }

    }

    return stat;
