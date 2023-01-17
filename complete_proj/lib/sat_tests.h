#pragma once

#include<iostream>
#include<string>
#include <eigen-3.4.0/Eigen/Dense>
#include <eigen-3.4.0/Eigen/Eigenvalues>

#include "Sequence.h"

std::pair<int,double> sat_test_cas1(Eigen::Matrix4d H, int n, Eigen::Vector4d freqs) {
    //z_val for one sided confidence interval for alpha=0.05
    double z_val = 1.65;

    double stat = -1;
    double div;
    double crit_val;
    double nucl_sum = freqs.sum();

    std::pair<int,double> cas1(0,0.0);

    // get relative frequences
    freqs = freqs/nucl_sum;

    for (int i = 0; i<4; i++) {
        div = n*freqs(i);
        stat += H(i,i)/div;
    }

    crit_val = sqrt(3/n) * z_val;

    cas1.second = stat;

    if (stat > crit_val) {
        return cas1;
    } else {
        cas1.first = 1;
        return cas1;
    }
}

std::pair<int,double> sat_test_cas2(Sequence Seq1, Sequence Seq2, Eigen::Matrix4d H, int n) {
    //z_val for one sided confidence interval for alpha=0.05
    double z_val = 1.65;

    double freqs[4] = {};
    double stat = -1;
    double div;
    double crit_val;

    std::pair<int,double> cas2(0,0.0);

    //relative frequences for the two sequences
    for (int i=0; i<4; i++) {
        freqs[i] = (Seq1.nucl_freqs[i]+Seq2.nucl_freqs[i])/(Seq1.length + Seq2.length);
    }

    for (int i = 0; i<4; i++) {
        div = n*freqs[i];
        stat += H(i,i)/div;
    }

    crit_val = sqrt(3/n) * z_val;

    cas2.second = stat;

    if (stat > crit_val) {
        return cas2;
    } else {
        cas2.first = 1;
        return cas2;
    }
}

std::pair<int, double> chi_test(Sequence Seq1, Sequence Seq2, Eigen::Matrix4d H, int n) {
    //Eigen::Vector4d seq1_freqs = Seq1.nucl_freqs2/Seq1.nucl_freqs2.sum();
    //Eigen::Vector4d seq2_freqs = Seq2.nucl_freqs2/Seq2.nucl_freqs2.sum();

    Eigen::Vector4d seq1_freqs = Seq1.nucl_freqs/Seq1.length;
    Eigen::Vector4d seq2_freqs = Seq2.nucl_freqs/Seq2.length;

    double stat = 0;
    double div;
    double crit_val = 19.023;

    std::pair<int,double> chi(0,0.0);

    for(int i = 0; i<4; i++) {
        for(int j = 0; j<4; j++) {
            div = std::pow(H(i,j)-n*seq1_freqs(i)*seq2_freqs(j),2);
            stat += div/(n*seq1_freqs(i)*seq2_freqs(j));

        }

    }

    chi.second = stat;

    if (stat > crit_val) {
        return chi;
    } else {
        chi.first = 1;
        return chi;
    }

}
