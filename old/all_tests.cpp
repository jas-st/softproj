#include <iostream>
#include <string>
#include <map>

#include "div_mat.h"
#include "stats.h"
#include "file_read.h"

using namespace std;

int main(int argc, char **argv) {
  string file_name;
  //float alpha;

  //read in flag argument -F for the file name with sequences
  int i = 1;
  while(i < argc) {
      if (strcmp(argv[i], "-F") == 0){
          file_name = argv[i+1];
      } //else if (strcmp(argv[i], "-a")==0){
          //alpha=atof(argv[i+1]);
      //}
      i++;
  }

  //cout << alpha << endl;

  //vector containing all of the sequences
  vector<string> sequences = seqs_read(file_name);
  Matrix4d H;  //diversity matrix

  //std::map<float,float*> crit_val_bow;
  //std::map<float,float*> crit_val_stuart;
  //std::map<float,float*> crit_qs;
  Eigen::Matrix<double, 4, 4> P_hat; 
  Eigen::Matrix<double, 4, 4> V;
  Eigen::Matrix<double, 4, 1> d;
  Eigen::Matrix<double, 6, 6> B;
  Eigen::Matrix<double, 6, 1> m;
  Eigen::Matrix<double, 3, 3> VS;
  Eigen::Matrix<double, 3, 1> D;
  double bow_s;
  double stu_s;
  double qs_s;
  double stat;
  int index;

  int size = sequences.size();  //number of sequences
  int n = sequences[0].length();
  double var;

  //std::map<float, std::array<float, 2> > crit_val_bow;

  cout <<"Pair"<<"\t"<<"Bowker ts" << '\t' <<"Bowker nr"<< "\t" << "Stuart ts " << '\t' <<"Stuart nr"<< "\t" <<"IS ts" << "\t" <<"IS nr"<<"\t"<< "QS ts" <<"\t" <<"QS nr"<< "\t" << endl;
  for (int i = 0; i < size; i++) {
      for (int j = i+1; j < size; j++) {
        cout << "(" << i+1 << "," << j+1 << ")" << "\t ";

        //div matrix H
        H = div_mat(sequences[i], sequences[j]);

        // Bowker test for symmetry
        m = get_m(H);
        B = get_B(H);
        bow_s=bowker_stat(m,B);
        cout << bow_s << "\t";
        
        if (1.237 < bow_s) {
          if (bow_s < 14.449) {
            cout << 1 << "\t";
          } else {
            cout << 0 << "\t";
          }
        } else { 
          cout << 0 << "\t";
        }

        // Stuart test for margianl symmetry
        VS = get_V(B, index); 
        D = get_D(m, index);
        stu_s=stuart_stat(D,VS);
        cout << stu_s << "\t";
        
        if (0.216 < stu_s) {
          if (stu_s < 9.348) {
            cout << 1 << "\t";
          } else {
            cout << 0 << "\t";
          }
        } else { 
          cout << 0 << "\t";
        }

        //  internal symmetry
        cout << bow_s-stu_s << "\t";
        
        if (0.216 < bow_s-stu_s) {
          if (bow_s-stu_s < 9.348) {
            cout << 1 << "\t";
          } else {
            cout << 0 << "\t";
          }
        } else { 
          cout << 0 << "\t";
        }

        // Proposed test
        P_hat = H/(n);

        d = get_d(H);

        var = 4; //sum of variances

        stat = 0;

        for(int i=0; i<4; i++){
          for(int j=i; j<4; j++){
            if (i == j){
              V(i,j) = get_var(P_hat, n, i+1);
            }else{
              V(i,j) = get_covar(P_hat, n, i+1, j+1);
              V(j,i) = V(i,j);
            }
          }
        }
    
        for(int i=0; i<4; i++){
          stat += d(i,0)/sqrt(V(i,i));

          for(int j = 0; j<4; j++){
            if(i != j){
              var += 1/sqrt(V(i,i)) * 1/sqrt(V(j,j)) * V(i,j);
            }
          }
        }
        qs_s=stat/sqrt(var);
        cout << qs_s << "\t";
        
        if (-1.96 < qs_s) {
          if (qs_s < 1.96) {
            cout << 1 << endl;
          } else {
            cout << 0 << endl;
          }
        } else { 
          cout << 0 << endl;
        }
      }
  }
}
