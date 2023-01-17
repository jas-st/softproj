#include<iostream>
#include <fstream>
#include <string.h>
#include "file_read.h"
#include <map>

#include "../lib/div_mat.h"
#include "../lib/stats.h"


using namespace std;


int main(int argc, char **argv) {
    string file_name;

    //read in flag argument -F for the file name with sequences
    int i = 1;
    while(i < argc) {
        if (strcmp( argv[i], "-F") == 0){
            file_name = argv[i+1];
        }
        i++;
    }

    //vector containing all of the sequences
    vector<string> sequences = seqs_read(file_name);
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
    int n = sequences[0].length();

    // print column names
    cout << "Bowker_test" << "       " << "Stuart_test" << "       " << "Proposed_test" << endl;

    for (int i=0; i<size; i++) {
        for (int j=i+1; j<size; j++) {

            H = div_mat(sequences[i],sequences[j]);

            //Bowker Test
            m = get_m(H);
            B = get_B(H);
            cout << bowker_stat(m,B) << "         ";

            // Stuart test for marginal symmetry
            V = get_V(B);
            D = get_D(m);

            cout << stuart_stat(D,V) << "         ";

            // Proposed test for quasy
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

                cout << stat/sqrt(var) << endl;
        }

    }


}
