#include <iostream>
#include <eigen-3.4.0/Eigen/Dense>
#include <eigen-3.4.0/Eigen/Eigenvalues>
#include <eigen-3.4.0/Eigen/Core>
#include <random>
using namespace Eigen;
using namespace std;

double my_rand_generator(int n){
	std::mt19937 generator (123);
  	std::uniform_real_distribution<double> dis(0.0, 1.0);
  	return(dis(generator));
}

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

Eigen::Matrix<double, 4,4> get_Dk(Eigen::Matrix<double, 4, 4> H, int k = 1)
{
	Eigen::Matrix<double, 4,4>  out = Eigen::Matrix<double, 4, 4>::Zero();
	for (int i = 0; i < 4; i++){
		out(i,i) = H(k,i)/H(i,k);
	}
	return(out);
}



Eigen::Matrix<double, 4,4> get_fDN(Eigen::Matrix<double, 4,4>  D, Eigen::Matrix<double, 4,4>  H)
{
	Eigen::Matrix<double, 4,4> out;
	/*float f;
	double fDN = 0;
	double fN = 0;
*/	out = D*H;

	/*for(int i = 0; i < 4; i++){
		for (int j = 0; j <4; j++){
			if ( i != j){
				fDN += out(i,j);
				fN += H(i,j);
			}
		}

	}

	//f = fDN/fN;
	f = fN/fDN;
	//f = 1.0;

	return(f*out);
	*/

	Eigen::Matrix<double, 4,4> D_new;

	Eigen::Matrix<double, 4,1> row_sum;
	Eigen::Matrix<double, 4,1> col_sum;

	row_sum << 0,0,0,0;
	col_sum = row_sum;





	for(int i = 0; i <4; i++){
		for(int j = 0; j <4; j++){
			if (i != j){
				row_sum(i,0) += out(i,j);
				col_sum(j,0) += out(i,j);
			}
		}
	}

	for(int i= 0; i<4; i++){
		D_new(i,i) = row_sum(i,0)/col_sum(i,0);
	}

	return(out*D_new);





}


double bowker_stat(Eigen::Matrix<double, 6, 1> m, Eigen::Matrix<double, 6, 6>  B)
{
	double stat;


	stat = m.transpose() * B.inverse() * m;



	return(stat);
}

double stuart_stat(Eigen::Matrix<double, 3, 1> D, Eigen::Matrix<double, 3, 3>  V)
{
	double out;

	out = D.transpose()*V.inverse()*D;

	return(out);
}

//----------------------------------------------------------------------------------------
Eigen::Matrix<double, 4, 4> get_Phat(Eigen::Matrix<double, 4, 4> H)
{
	Eigen::Matrix<double, 4, 1> row_sum;
	Eigen::Matrix<double, 4, 1> ones;
	Eigen::Matrix<double, 4, 4> D;
	D << Eigen::Matrix<double, 4, 4>::Zero();
	//float total_sum;

	ones << 1,1,1,1;
	/*for(int i = 0; i<4; i++){
		H(i,i) = 0;
	}*/

	row_sum = H*ones;

	for(int i = 0; i<4; i++){
		D(i,i) = row_sum(i,0);
	//	total_sum += row_sum(i,0);
	}


	//D = 1/total_sum * D;

	return(D.inverse()*H);

}

Eigen::Matrix<double, 1, 4> get_EV1(Eigen::Matrix<double, 4, 4> P, float n)
{
	Eigen::Matrix<double, 4, 1> out;
	EigenSolver<Matrix4d> es(P.transpose());
	Eigen::Matrix<double, 4, 4>D = es.pseudoEigenvalueMatrix();
	Eigen::Matrix<double, 4, 4> V = es.pseudoEigenvectors();
	int i = 0;


	while( abs(D(i,i)- n) >  0.01){
		i +=1 ;
	}

	out = V(all,i).transpose();

	// normalize
	float sum = 0;

	for (int j = 0; j < 4; j++){
		sum += out(j,0);
	}

	return(1/sum * out);

}


Eigen::Matrix<double, 4, 4> get_Shat(Eigen::Matrix<double, 4, 4> P, Eigen::Matrix<double, 1, 4> v)
{
	Eigen::Matrix<double, 4, 4>D;
	D << Eigen::Matrix<double, 4, 4>::Zero();

	for(int i = 0; i< 4; i++){
		D(i,i) = v(0,i);
	}


	return(D*P);
}

float get_nondiag(Eigen::Matrix<double, 4, 4> H)
{
	float out = 0;
	for (int i = 0; i<4; i++){
		for (int j = 0; j<4; j++){
			if (i != j){
				out += H(i,j);
			}
		}
	}

	return(out);

}

float get_minrowcolquot(Eigen::Matrix<double, 4, 4> H)
{
	float out = 1;
	Eigen::Matrix<double, 4, 1> row_sum;
	Eigen::Matrix<double, 4, 1> col_sum;
	Eigen::Matrix<double, 4, 1> ones;
	float quot;

	ones << 1,1,1,1;

	row_sum = H*ones;

	col_sum = (ones.transpose()*H).transpose();

	for (int i=0; i<4; i++){
		quot = row_sum(i,0)/col_sum(i,0);
		cout << quot << "\t";
		if (quot < out){
			out = quot;
		}
	}

	return(quot);



}

double bowker_stat(Eigen::Matrix<double, 4,4 > H, int n)
{
	double stat = 0;

	for (int i = 0; i <4; i++){
		for (int j = i+1; j <4; j++){
			if ( i != n){
				if (H(i,j) != 0 | H(j,i) != 0) {
					stat += pow((H(i,j)-H(j,i)),2)/(H(i,j)+H(j,i));
				}
			}

		}
	}





	return(stat);
}

Eigen::Matrix<double, 4, 4> get_Hhat(Eigen::Matrix<double, 4, 4> P_hat)
{
	EigenSolver<Matrix4d> es(P_hat);
	Eigen::Matrix<double, 4, 4> D = es.pseudoEigenvalueMatrix();
	Eigen::Matrix<double, 4, 4> V = es.pseudoEigenvectors();
// take ln of D
	for(int i = 0; i < 4; i++){
		D(i,i) = log(D(i,i));
	}
	return (V*D*V.inverse());
}


double naiv_test(Eigen::Matrix<double, 4, 4> H)
{
	double stat1;
	double stat2;
	double stat3;
	double stat4;
	double out;

	stat1 = pow(H(0,1)*H(1,2)*H(2,0)-H(0,2)*H(2,1)*H(1,0), 2)/(H(0,1)*H(1,2)*H(2,0)+H(0,2)*H(2,1)*H(1,0));
	stat2 = pow(H(0,1)*H(1,3)*H(3,0)-H(0,3)*H(3,1)*H(1,0), 2)/(H(0,1)*H(1,3)*H(3,0)+H(0,3)*H(3,1)*H(1,0));
	stat3 = pow(H(0,2)*H(2,3)*H(3,0)-H(0,3)*H(3,2)*H(2,0), 2)/(H(0,2)*H(2,3)*H(3,0)+H(0,3)*H(3,2)*H(2,0));
	stat4 = pow(H(1,2)*H(2,3)*H(3,1)-H(1,3)*H(3,2)*H(2,1), 2)/(H(1,2)*H(2,3)*H(3,1)+H(1,3)*H(3,2)*H(2,1));

	out = stat1+stat2+stat3+stat4;// - max(max(stat1,stat2),max(stat3,stat4));

	return(out);


}

Eigen::Matrix<double, 4, 4> get_diag(Eigen::Matrix<double, 1, 4> v)
{

	Eigen::Matrix<double, 4, 4>D;
	D << Eigen::Matrix<double, 4, 4>::Zero();

	for(int i = 0; i < 4; i++){
		D(i,i) = v(0,i);
	}
	return (D);
}

Eigen::Matrix<double, 4, 4> symm_orth_test(Eigen::Matrix<double, 4, 4> S)
{
	EigenSolver<Matrix4d> es(S);
	Eigen::Matrix<double, 4, 4> V = es.pseudoEigenvectors();
	Eigen::Matrix<double, 4, 4> V_norm;
	Eigen::Matrix<double, 4, 4> out;
	float length;

	// normalize Eigenvectors
	for(int i = 0; i<4; i++){
		length = 0;
		for(int j= 0; j<4; j++){
			length += V(j,i)*V(j,i);
		}
		for(int j= 0; j<4; j++){
			 V_norm(j,i) = V(j,i)/sqrt(length);
		}
	}

	return(V_norm.transpose()*V_norm);

}

float get_upper_tria(Eigen::Matrix<double, 4, 4> V)
{
	float out = 0;
	for(int i = 0; i<4; i++){
		for(int j =i+1; j<4; j++){
			out += V(i,j)*V(i,j);
		}
	}
	return(out);
}



double get_expect(Eigen::Matrix<double, 4, 4> P, int N,Eigen::Matrix<double, 3, 2> indices)
{
	double out;
	double p1 = P(0,1);
	double p2 = P(1,2);
	double p3 = P(2,0);
	out = N*(N-1)*(N-2)*p1*p2*p3;
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

/*double get_var1(Eigen::Matrix<double, 4, 4> P, int N,Eigen::Matrix<double, 3, 2> indices)
{
	double out;
	double exp = get_expect(P, N, indices);
	double exp_sqrt = get_expect_sqrt(P, N, indices);
	out = exp_sqrt - exp*exp;
	return(out);
}*/

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


Eigen::Matrix<double, 4, 4> get_Q(std::mt19937& generator, std::uniform_real_distribution<double> dis)
{
	Eigen::Matrix<double, 4, 1> d;
	Eigen::Matrix<double, 4, 1> row_sum;
	Eigen::Matrix<double, 4, 4> psi;
	Eigen::Matrix<double, 4, 4> S;
	Eigen::Matrix<double, 4, 4> out;
	double sum = 0;

	for(int i = 0; i<4; i++){
		d(i,0) = dis(generator);
		sum += d(i,0);
	}

	d = d/sum;
	psi = get_diag(d);

	for(int i = 0; i<4; i++){
		for(int j = i+1; j<4; j++){
			S(i,j) = dis(generator);
			S(j,i) = S(i,j);
		}
	}


	out = S*psi;
	row_sum = out.rowwise().sum();
	sum = 0;

	for(int i = 0; i<4; i++){
		out(i,i) = -row_sum(i,0);
		sum -= out(i,i)*d(i,0);
	}

	return(out/sum);

}

Eigen::Matrix<double, 4, 4> get_Q_given_pi(std::mt19937& generator, std::uniform_real_distribution<double> dis, Eigen::Matrix<double, 4, 1> d)
{

	Eigen::Matrix<double, 4, 1> row_sum;
	Eigen::Matrix<double, 4, 4> psi;
	Eigen::Matrix<double, 4, 4> S;
	Eigen::Matrix<double, 4, 4> out;
	double sum = 0;

	psi = get_diag(d);

	for(int i = 0; i<4; i++){
		for(int j = i+1; j<4; j++){
			S(i,j) = dis(generator);
			S(j,i) = S(i,j);
		}
	}


	out = S*psi;
	row_sum = out.rowwise().sum();
	sum = 0;

	for(int i = 0; i<4; i++){
		out(i,i) = -row_sum(i,0);
		sum -= out(i,i)*d(i,0);
	}

	return(out/sum);

}

Eigen::Matrix<double, 4, 4> get_Q_given_S(std::mt19937& generator, std::uniform_real_distribution<double> dis, Eigen::Matrix<double, 4, 4> S)
{

	Eigen::Matrix<double, 4, 1> d;
	Eigen::Matrix<double, 4, 1> row_sum;
	Eigen::Matrix<double, 4, 4> psi;
	Eigen::Matrix<double, 4, 4> out;
	double sum = 0;

	for(int i = 0; i<4; i++){
		d(i,0) = dis(generator);
		sum += d(i,0);
	}

	d = d/sum;
	psi = get_diag(d);

	out = S*psi;
	row_sum = out.rowwise().sum();
	sum = 0;

	for(int i = 0; i<4; i++){
		out(i,i) = -row_sum(i,0);
		sum -= out(i,i)*d(i,0);
	}

	return(out/sum);

}



//___________________________________________________________________________________________________________________________________________________________________________
