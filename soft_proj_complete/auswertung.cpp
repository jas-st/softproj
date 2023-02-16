#include "lib/file_bio.h"
#include "lib/sat_tests2.h"

#include "lib/div_mat.h"
#include "lib/stats.h"

using namespace std;


int main(int argc, char **argv) {
    string file_name;
    Alignment alignment;
    string extension;
    bool sat_tests = false;

    //read in flag argument -F for the file name with sequences
    //read in flag argument -s for saturation tests, default - false
    int i = 1;
    while(i < argc) {
        if (strcmp( argv[i], "-F") == 0){
            file_name = argv[i+1];
        }
        if (strcmp( argv[i], "-s") == 0){
            sat_tests = argv[i+1];
        }
        i++;
    }

    extension = file_name.substr(file_name.length()-3);


    //vector containing all of the sequences
    //idea - to return an object that will store the seq map and the freq array

    alignment = seqs_read(file_name, extension);
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

    // print column names
    std::cout << "Sequences compared" << ',';
    if (sat_tests) {
        std::cout << "Sat_test_Cassius1" << ',' << "Sat_test_Cassius2"
                  << ',' << "Chi_test" << ',';
    }

    std::cout << "Bowker_test"   << ',' << "Stuart_test"      << ',' << "Internal_Symmetry"
       << ',' << "Proposed_test" << ',' << "Alignment_Length" << endl;

    for (int i=0; i<size; i++) {
        for (int j=i+1; j<size; j++) {
            Sequence seq1 = sequences[i];
            Sequence seq2 = sequences[j];

            std::cout << "("<< seq1.id << "; " << seq2.id << ")" << ",";

            H = div_mat(seq1.seq, seq2.seq);
            n = H.sum(); //alignment length

            for (int i=0; i<4; i++) {
                seq1.nucl_freqs_al(i) = H.row(i).sum();
                seq2.nucl_freqs_al(i) = H.col(i).sum();
            }

            if (n == 0) {
                std::cout << "invalid matrix" << endl;
                continue;
            }

            //Saturation tests
            if (sat_tests) {
                double cas1 = sat_test_cas1(H, n, alignment.global_freqs);
                double cas2 = sat_test_cas2(seq1, seq2, H, n);
                double chi = chi_test_draft(seq1, seq2, H, n);


                std::cout << cas1 << ',';
                std::cout << cas2 << ',';
                std::cout << chi  << ',';
            }

            //Bowker Test
            m = get_m(H);
            B = get_B(H);
            stat_b = bowker_stat(m,B);

            std::cout << stat_b << ',';

            // Stuart test for marginal symmetry

            V = get_V(B);
            D = get_D(m);
            stat_s = stuart_stat(D,V);

            std::cout << stat_s << ',';

            // Internal symmetry
            int_sym = stat_b - stat_s;

            std::cout << int_sym << ',';


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

            std::cout << stat_q << ',' << n << std::endl;


        }


    }
}
