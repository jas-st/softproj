#include <iostream>
#include <eigen-3.4.0/Eigen/Dense>
#include <eigen-3.4.0/Eigen/Eigenvalues>
#include <eigen-3.4.0/Eigen/Core>

using namespace Eigen;
using namespace std;

//------------BOWKER STAT-------------------//
Eigen::Matrix<double, 6, 1> get_m(Matrix4d H)
{
	Eigen::Matrix<double, 6, 1> out;
	int k = 0;

	for (int i = 0; i < 4; i++){
		for (int j = i + 1; j < 4; j++){
			out(k) = abs(H(i,j)) - abs(H(j,i));
			k +=1;
		}
	}

	return(out);
}

Eigen::Matrix<double, 6, 6> get_B(Matrix4d H)
{
	Eigen::Matrix<double, 6, 6> out;
	int k = 0;
	out = Eigen::Matrix<double, 6, 6>::Zero();

	for (int i = 0; i < 4; i++){
		for (int j = 0; j < 4; j++){
			if (i < j){
				out(k,k) = H(i,j) + H(j,i);
				k += 1;
			}
		}
	}


	return(out);
}

double bowker_stat(Eigen::Matrix<double, 6, 1> m, Eigen::Matrix<double, 6, 6>  B)
{
	double stat;


	stat = m.transpose() * B.inverse() * m;



	return(stat);
}

//------------STUART STAT-------------------//

 Eigen::Matrix<double, 3, 3> get_V(Eigen::Matrix<double, 6, 6> B, int index = 1)
{
	Eigen::Matrix<double, 4, 6> C;
	Eigen::Matrix<double, 3, 6> C_index;
	Eigen::Matrix<double, 3, 3> out;
	int j = 0;

	C << 1,1,1,0,0,0,-1, 0, 0, 1, 1, 0,0, -1, 0, -1, 0, 1, 0, 0, -1, 0, -1, -1;

	for(int i = 0; i <4; i++){
		if(!(i==index)){
			C_index(j,all) = C(i,all);
			j += 1;
		}
	}

	out = C_index*B*C_index.transpose();

	return(out);

}

Eigen::Matrix<double, 3,1> get_D(Eigen::Matrix<double, 6, 1> m, int index = 1)
{
	Eigen::Matrix<double, 3, 1> out;
	Eigen::Matrix<double, 4, 1> D;
	Eigen::Matrix<double, 4, 6> C;
	int j;


	C << 1,1,1,0,0,0,-1, 0, 0, 1, 1, 0,0, -1, 0, -1, 0, 1, 0, 0, -1, 0, -1, -1;

	D = C*m;

	j = 0;
	for(int i = 0; i <4; i++){
		if(!(i==index)){
			out(j,all) = D(i,all);
			j += 1;
		}
	}

	return(out);
}

double stuart_stat(Eigen::Matrix<double, 3, 1> D, Eigen::Matrix<double, 3, 3>  V)
{
	double out;

	out = D.transpose()*V.inverse()*D;

	return(out);
}

//------------INTERNAL SYM STAT-------------------//

double intsym_stat(double bowk, double stu) {
	return (bowk - stu);
}

//------------QUASI-SYMMETRY STAT-------------------//

Eigen::Matrix<double, 4, 1> get_d(Eigen::Matrix<double, 4, 4> P)
{
	Eigen::Matrix<double, 4, 1> out;
	double d1, d2, d3, d4;

	d1 = P(0,1)* P(1,2)* P(2,0)- P(0,2)* P(2,1)* P(1,0);
	d2 = P(0,1)* P(1,3)* P(3,0)- P(0,3)* P(3,1)* P(1,0);
	d3 = P(0,2)* P(2,3)* P(3,0)- P(0,3)* P(3,2)* P(2,0);
	d4 = P(1,2)* P(2,3)* P(3,1)- P(1,3)* P(3,2)* P(2,1);

	out << d1, d2, d3, d4;

	return(out);

}

double get_expect_sqrt1(int N, double p1, double p2, double p3)
{
	double out;
	double tmp = log(p1*p2 + p1*p3 + p2*p3);
	double tmp2 = log(p1*p2*p3);
	out = log(1 + (N-3)*(p1+p2+p3) + exp(log(N-3)+log(N-4)+tmp) + exp(log(N-3)+log(N-4)+log(N-5)+tmp2));
	return(exp(log(N)+log(N-1)+log(N-2)+tmp2 + out));
}

Eigen::Matrix<double, 6, 1> get_p(Eigen::Matrix<double, 4, 4> P, int d1)
{
	Eigen::Matrix<double, 6, 1> out;

	if(d1 == 1){
		out << P(0,1), P(1,2), P(2,0), P(0,2), P(2,1), P(1,0);
	}else if(d1 == 2){
		out << P(0,1), P(1,3), P(3,0), P(0,3), P(3,1), P(1,0);
	}else if(d1 == 3){
		out << P(0,2), P(2,3), P(3,0), P(0,3), P(3,2), P(2,0);
	}else if(d1 == 4){
		out << P(1,2), P(2,3), P(3,1), P(1,3), P(3,2), P(2,1);
	}
	return(out);

}

double get_var(Eigen::Matrix<double, 4, 4> P, int N, int d1)
{
	Eigen::Matrix<double, 6, 1> tmp_p;
	double p1, p2, p3, p4, p5, p6;
	double out;
	double tmp, tmp1, tmp2, tmp3;
	tmp_p = get_p(P, d1);
	p1 = tmp_p(0,0);
	p2 = tmp_p(1,0);
	p3 = tmp_p(2,0);
	p4 = tmp_p(3,0);
	p5 = tmp_p(4,0);
	p6 = tmp_p(5,0);

	tmp = exp(log(N)+log(N-1)+log(N-2)+log(p1)+log(p2)+log(p3))-exp(log(N)+log(N-1)+log(N-2)+log(p4)+log(p5)+log(p6)); // E[a - b]
	tmp1 = get_expect_sqrt1(N,p1, p2, p3);// E[a^2]
	tmp2 = exp(log(N)+log(N-1)+log(N-2)+log(N-3)+log(N-4)+log(N-5)+log(p1)+log(p2)+log(p3)+log(p4)+log(p5)+log(p6));// E[ab]
	tmp3 = get_expect_sqrt1(N,p4, p5, p6);// E[b^2]
	out = tmp1 -2*tmp2 + tmp3 - tmp*tmp;

	return(out);
}

double get_covar_help(Eigen::Matrix<double, 4, 4> P, int N, int d1, int d2)
{
	Eigen::Matrix<double, 6, 1> tmp_p;
	double p1, p2, p3, p4, p5, p6;
	double b1, b2, b3, b4, b5, b6;
	double tmp1, tmp2, tmp3, tmp4, out;

	tmp_p = get_p(P, d1);
	p1 = tmp_p(0,0);
	p2 = tmp_p(1,0);
	p3 = tmp_p(2,0);
	p4 = tmp_p(3,0);
	p5 = tmp_p(4,0);
	p6 = tmp_p(5,0);


	tmp_p = get_p(P, d2);
	b1 = tmp_p(0,0);
	b2 = tmp_p(1,0);
	b3 = tmp_p(2,0);
	b4 = tmp_p(3,0);
	b5 = tmp_p(4,0);
	b6 = tmp_p(5,0);


	if(d1 == 1){
		if( d2 == 2){
			tmp1 = exp(log(N)+log(N-1)+log(N-2)+log(N-3)+log(N-4)+log(p2)+log(p3)+log(b1)+log(b2)+log(b3))*(1+(N-5)*p1);
			tmp2 = exp(log(N)+log(N-1)+log(N-2)+log(N-3)+log(N-4)+log(N-5)+log(p1)+log(p2)+log(p3)+log(b4)+log(b5)+log(b6));
			tmp3 = exp(log(N)+log(N-1)+log(N-2)+log(N-3)+log(N-4)+log(N-5)+log(p4)+log(p5)+log(p6)+log(b1)+log(b2)+log(b3));
			tmp4 = exp(log(N)+log(N-1)+log(N-2)+log(N-3)+log(N-4)+log(p4)+log(p5)+log(b4)+log(b5)+log(b6))*(1+(N-5)*p6);
		}else if(d2 == 3){
			tmp1 = exp(log(N)+log(N-1)+log(N-2)+log(N-3)+log(N-4)+log(N-5)+log(p1)+log(p2)+log(p3)+log(b1)+log(b2)+log(b3));
			tmp2 = exp(log(N)+log(N-1)+log(N-2)+log(N-3)+log(N-4)+log(p1)+log(p2)+log(b4)+log(b5)+log(b6))*(1+(N-5)*p3);
			tmp3 = exp(log(N)+log(N-1)+log(N-2)+log(N-3)+log(N-4)+log(p5)+log(p6)+log(b1)+log(b2)+log(b3))*(1+(N-5)*p4);
			tmp4 = exp(log(N)+log(N-1)+log(N-2)+log(N-3)+log(N-4)+log(N-5)+log(p4)+log(p5)+log(p6)+log(b4)+log(b5)+log(b6));
		}else if(d2 == 4){
			tmp1 = exp(log(N)+log(N-1)+log(N-2)+log(N-3)+log(N-4)+log(p1)+log(p3)+log(b1)+log(b2)+log(b3))*(1+(N-5)*p2);
			tmp2 = exp(log(N)+log(N-1)+log(N-2)+log(N-3)+log(N-4)+log(N-5)+log(p1)+log(p2)+log(p3)+log(b4)+log(b5)+log(b6));
			tmp3 = exp(log(N)+log(N-1)+log(N-2)+log(N-3)+log(N-4)+log(N-5)+log(p4)+log(p5)+log(p6)+log(b1)+log(b2)+log(b3));
			tmp4 = exp(log(N)+log(N-1)+log(N-2)+log(N-3)+log(N-4)+log(p4)+log(p6)+log(b4)+log(b5)+log(b6))*(1+(N-5)*p5);
		}
	}else if(d1 == 2){
		if(d2 == 3){
			tmp1 = exp(log(N)+log(N-1)+log(N-2)+log(N-3)+log(N-4)+log(p1)+log(p2)+log(b1)+log(b2)+log(b3))*(1+(N-5)*p3);
			tmp2 = exp(log(N)+log(N-1)+log(N-2)+log(N-3)+log(N-4)+log(N-5)+log(p1)+log(p2)+log(p3)+log(b4)+log(b5)+log(b6));
			tmp3 = exp(log(N)+log(N-1)+log(N-2)+log(N-3)+log(N-4)+log(N-5)+log(p4)+log(p5)+log(p6)+log(b1)+log(b2)+log(b3));
			tmp4 = exp(log(N)+log(N-1)+log(N-2)+log(N-3)+log(N-4)+log(p5)+log(p6)+log(b4)+log(b5)+log(b6))*(1+(N-5)*p4);
		}else if(d2 == 4){
			tmp1 = exp(log(N)+log(N-1)+log(N-2)+log(N-3)+log(N-4)+log(N-5)+log(p1)+log(p2)+log(p3)+log(b1)+log(b2)+log(b3));
			tmp2 = exp(log(N)+log(N-1)+log(N-2)+log(N-3)+log(N-4)+log(p1)+log(p3)+log(b4)+log(b5)+log(b6))*(1+(N-5)*p2);
			tmp3 = exp(log(N)+log(N-1)+log(N-2)+log(N-3)+log(N-4)+log(p4)+log(p6)+log(b1)+log(b2)+log(b3))*(1+(N-5)*p5);
			tmp4 = exp(log(N)+log(N-1)+log(N-2)+log(N-3)+log(N-4)+log(N-5)+log(p4)+log(p5)+log(p6)+log(b4)+log(b5)+log(b6));
		}
	}else if(d1 == 3){
		if(d2 == 4){
			tmp1 = exp(log(N)+log(N-1)+log(N-2)+log(N-3)+log(N-4)+log(p1)+log(p3)+log(b1)+log(b2)+log(b3))*(1+(N-5)*p2);
			tmp2 = exp(log(N)+log(N-1)+log(N-2)+log(N-3)+log(N-4)+log(N-5)+log(p1)+log(p2)+log(p3)+log(b4)+log(b5)+log(b6));
			tmp3 = exp(log(N)+log(N-1)+log(N-2)+log(N-3)+log(N-4)+log(N-5)+log(p4)+log(p5)+log(p6)+log(b1)+log(b2)+log(b3));
			tmp4 = exp(log(N)+log(N-1)+log(N-2)+log(N-3)+log(N-4)+log(p4)+log(p6)+log(b4)+log(b5)+log(b6))*(1+(N-5)*p5);
		}

	}
	out = tmp1 - tmp2 - tmp3 + tmp4;
	return(out);
}


double get_covar(Eigen::Matrix<double, 4, 4> P, int N, int d1, int d2)
{
	Eigen::Matrix<double, 6, 1> tmp_p;
	double p1, p2, p3, p4, p5, p6;
	double b1, b2, b3, b4, b5, b6;
	double out;
	double tmp, tmp1, tmp2, tmp3;
	double exp_d1, exp_d2, exp3, exp4;
	tmp_p = get_p(P, d1);
	p1 = tmp_p(0,0);
	p2 = tmp_p(1,0);
	p3 = tmp_p(2,0);
	p4 = tmp_p(3,0);
	p5 = tmp_p(4,0);
	p6 = tmp_p(5,0);

	exp_d1 = exp(log(N)+log(N-1)+log(N-2)+log(p1)+log(p2)+log(p3)) - exp(log(N)+log(N-1)+log(N-2)+log(p4)+log(p5)+log(p6));

	tmp_p = get_p(P, d2);
	b1 = tmp_p(0,0);
	b2 = tmp_p(1,0);
	b3 = tmp_p(2,0);
	b4 = tmp_p(3,0);
	b5 = tmp_p(4,0);
	b6 = tmp_p(5,0);

	exp_d2 = exp(log(N)+log(N-1)+log(N-2)+log(b1)+log(b2)+log(b3)) - exp(log(N)+log(N-1)+log(N-2)+log(b4)+log(b5)+log(b6));

	//cout << "get covar test comming:" << endl << exp_d1 << endl << endl;


	out = get_covar_help(P, N, d1, d2) - exp_d1*exp_d2;

	return(out);
}

double quasisym_stat(Eigen::Matrix<double, 4, 4> Phat, Eigen::Matrix<double, 4, 1> d, double n) {
	Eigen::Matrix<double, 4, 4> V2;
    double stat;
    double var;

	var = 4; //sum of variances
    stat = 0;

	for(int i=0; i<4; i++){
		for(int j=i; j<4; j++){
			if (i == j){
				V2(i,j) = get_var(Phat, n, i+1);
			}else{
				V2(i,j) = get_covar(Phat, n, i+1, j+1);
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

	return (stat/sqrt(var));
}
