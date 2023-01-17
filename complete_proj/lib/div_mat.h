#include <iostream>
#include <eigen-3.4.0/Eigen/Dense>
#include <eigen-3.4.0/Eigen/Eigenvalues>
using namespace Eigen;
using namespace std;

Matrix4d div_mat(string seq1, string seq2)
{
	Matrix4d out;
	string nuc1;
	string nuc2;
	map <string , int> dict;

	dict["A"] = 0;
	dict["C"] = 1;
	dict["G"] = 2;
	dict["T"] = 3;

	out << Eigen::Matrix<double, 4, 4>::Zero();

	for (int i = 0; i < seq1.size(); i++){
		if (seq1[i] != '-' && seq2[i] != '-') {
			nuc1 = seq1[i];
			nuc2 = seq2[i];
			out(dict[nuc1], dict[nuc2]) += 1;
		}

	}
	return(out);
}
