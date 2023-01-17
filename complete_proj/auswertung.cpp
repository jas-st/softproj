#include "file_bio.h"
#include "sat_tests2.h"

#include "../lib/div_mat.h"
#include "../lib/stats.h"

using namespace std;


int main(int argc, char **argv) {
    string file_name;
    Alignment alignment;
    string extension;

    //read in flag argument -F for the file name with sequences
    int i = 1;
    while(i < argc) {
        if (strcmp( argv[i], "-F") == 0){
            file_name = argv[i+1];
        }
        i++;
    }

    extension = file_name.substr(file_name.length()-3);

    if (strcmp(extension.c_str(), "nex") == 0 ) {
        alignment = seqs_read_bio(file_name);
    } else {
        alignment = seqs_read(file_name);
    }

    //vector containing all of the sequences
    //idea - to return an object that will store the seq map and the freq array

    vector<Sequence> sequences = alignment.sequences;
    Eigen::Vector4d global_freqs = alignment.global_freqs;


    Matrix4d H;  //diversity matrix

    Eigen::Matrix<double, 6, 6> B;
    Eigen::Matrix<double, 6, 1> m;
    Eigen::Matrix<double, 3, 3> V;
    Eigen::Matrix<double, 3, 1> D;
    Eigen::Matrix<double, 4, 4> P_hat;
    Eigen::Matrix<double, 4, 1> d;
    Eigen::Matrix<double, 4, 4> V2;
    double stat;
    double var;

    int size = sequences.size();  //number of sequences
    int n; // fix it so that it counts only the seq length? without gap ? i guess
    double stat_b;
    double stat_s;
    double stat_q;
    double int_sym;
    int res_b;
    int res_s;
    int res_q;
    int int_sum_q;

    // print column names
    std::cout << "Sequences compared" << '\t' <<"Sat test Cassius 1" << '\t' << "Value" << '\t'
    << "Sat test Cassius 2" << '\t' << "Value" << '\t' << "Chi test" << '\t' << "Value"<< '\t'
    << "Bowker_test" << '\t' << "Value" << '\t' << "Stuart_test" << '\t' << "Value"<< '\t'
    << "Internal Symmetry" << '\t' << "Value"<< '\t' << "Proposed_test" << '\t' << "Value" << endl;

    for (int i=0; i<size; i++) {
        for (int j=i+1; j<size; j++) {
            Sequence seq1 = sequences[i];
            Sequence seq2 = sequences[j];

            std::cout << "("<< seq1.id << ", " << seq2.id << ")" << "\t";

            H = div_mat(seq1.seq, seq2.seq);
            n = H.sum(); //alignment length

            if (n == 0) {
                std::cout << "invalid matrix" << endl;
                continue;
            }

            //Saturation tests
            pair<int, double> cas1 = sat_test_cas1(H, n, alignment.global_freqs);
            pair<int, double> cas2 = sat_test_cas2(seq1, seq2, H, n);
            pair<int, double> chi = chi_test(seq1, seq2, H, n);


            std::cout << cas1.first << '\t' << cas1.second << '\t';
            std::cout << cas2.first << '\t' << cas2.second << '\t';
            std::cout << chi.first << '\t' << chi.second << '\t';

            //Bowker Test
            m = get_m(H);
            B = get_B(H);
            stat_b = bowker_stat(m,B);


            if (1.237 < stat_b && stat_b < 14.449) {
                res_b = 1;
            } else {
                res_b = 0;
            }

            std::cout << res_b << '\t' << stat_b << '\t';

            // Stuart test for marginal symmetry
            if (n != 0) {
                V = get_V(B);
                D = get_D(m);
                stat_s = stuart_stat(D,V);
            } else {
                stat_s = -1;
            }

            if (0.216 < stat_s && stat_s < 9.348) {
                res_s = 1;
            } else {
                res_s = 0;
            }

            std::cout << res_s << '\t' << stat_s << '\t';

            // Internal symmetry
            int_sym = stat_b - stat_s;

            if (0.216 < int_sym && int_sym < 9.348) {
                int_sum_q = 1;
            } else {
                int_sum_q = 0;
            }

            std::cout << int_sum_q << '\t' << int_sym << '\t';


            // Proposed test for quasy
            if (n != 0) {
            P_hat =H/n;

            d = get_d(H);

            var = 4; //sum of variances

            stat = 0;


            for(int i=0; i<4; i++){
                for(int j=i; j<4; j++){
                    if (i == j){
                        V2(i,j) = get_var(P_hat, n, i+1);
                    }else{
                        V2(i,j) = get_covar(P_hat, n, i+1, j+1);
                        V2(j,i) = V2(i,j);
                    }
                }
            }

            for(int i=0; i<4; i++){
                stat += d(i,0)/sqrt(V2(i,i));


                for(int j = 0; j<4; j++){
                    if(i != j){
                        var += 1/sqrt(V2(i,i)) * 1/sqrt(V2(j,j)) * V2(i,j);
                    }
                }
            }

            stat_q = stat/sqrt(var);
            } else {
                stat_q = -2.0;
            }

            if (abs(stat_q) <= 1.96) {
                res_q = 1;
            } else {
                res_q = 0;
            }

            std::cout << res_q << '\t' << stat_q << std::endl;


        }


    }
}
