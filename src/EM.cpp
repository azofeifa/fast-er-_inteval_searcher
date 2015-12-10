#include "EM.h"
#include <cmath>
#include <algorithm>
using namespace std;
void get_distances(map<string, vector<segment> > query, map<string, vector<double> > & distances  ){
	typedef map<string, vector<segment> >::iterator it_type;
	for (it_type c= query.begin(); c!= query.end(); c++){
		vector<segment> q 	= c->second ;
		for (int i = 0 ; i < q.size(); i++){
			string INF 	= q[i].info+"|" ;
			int ON 		= q[i].overlaps.size();
			for (int o = 0; o < q[i].overlaps.size(); o++){
				double a 	= (q[i].start + q[i].stop) / 2.;
				double b 	= (q[i].overlaps[o].start + q[i].overlaps[o].stop) / 2.; 
				double d 	= a-b;
				distances[q[i].overlaps[o].info].push_back(d);
			}
		}
	}		
}

void BIN(vector<double> raw_data, int bins, double ** & X){
	//bin up the data
	X 	= new double*[bins];
	double min_x=raw_data[0], max_x=raw_data[0];
	int j 	= 0;
	for (int i = 0 ; i < max(int(raw_data.size()),bins) ; i++){
		if (j < bins){
			X[j] 	= new double[2];
		}
		j++;
		if (i < raw_data.size() and raw_data[i] < min_x){
			min_x 	= raw_data[i];
		}
		if (i < raw_data.size() and raw_data[i] > max_x){
			max_x 	= raw_data[i];
		}
	}
	double delta 	= (max_x - min_x) / bins;
	for (int i = 0 ; i < bins; i++){
		X[i][0] 		= min_x + i*delta;
		X[i][1] 		= 0;
	}
	for (int i = 0; i < raw_data.size(); i++){
		int j 	= 0;
		while (j < bins and X[j][0] < raw_data[i]){
			j++;
		}
		if (j < bins){
			X[j][1]+=1;
		}
	}
}

double norm_pdf(double x, double mu, double si){
	double vl  	= (1.0 / (sqrt(2*M_PI)*si  ))*exp(-pow(x-mu,2)/(2*pow(si,2))  );
	return vl;
}
double uni_pdf(double x, double a,double b){
	double vl 	= 1.0 / (b-a);
	return vl;
}
double depletion_func(double x, double MU, double SI, double a, double b){
	return uni_pdf(x, a,b)	 / (uni_pdf(x,a,b) + norm_pdf(x, MU, SI));
}

double EM(double ** X, int N, double & w, double & ll, double si, double mu, 
	int func_type, double C, double a, double b ){
	
	int T 	= 1000, t=0;
	bool converged = false;
	double prev_ll 	= 0;
	if (func_type==3){
		w 	= 1;
		double u;
		for (int i = 0; i < N;i++ ){
			u = uni_pdf(X[i][0], a, b);
			ll+=log(u)*X[i][1];
		}
		return 0;
	}
	double SI 	= si;
	while (not converged and t < T ){
		double EXN= 0, EXU=0, EX=0;
		ll=0;

		for (int i = 0; i < N;i++ ){
			double  n, u;
			u = uni_pdf(X[i][0], a, b)*(1-w);
			if (func_type==0){
				n = norm_pdf(X[i][0], mu, si)*w;
			}else if(func_type==1) {
				n = norm_pdf(X[i][0], mu, SI)*w;		
			}else if (func_type==2){
				n = depletion_func(X[i][0], mu, si, a,b)*w*(1.0 / C);				
			}
			if (func_type==1){
				EX+=pow(X[i][0]-mu,2 )*(n / (n+u))*X[i][1];
			}
			EXN+=((n / (n+u)  )*X[i][1]);
			EXU+=((u / (n+u)  )*X[i][1]);
			ll+=(log( n + u   )*X[i][1]) ;
		}
		if (func_type==1){
			SI 	= sqrt(EX / (EXN));
		}
		w 	= EXN / (EXN + EXU);
		if (abs(ll - prev_ll) < pow(10,-3) and t > 3) {
			converged=true;
		}
		prev_ll 	= ll;
		t++;
	}
	return SI;
}


double normal_constant(double MU, double SI, double A, double B){
	double step_size 	= 100000;
	double delta 		= (B-A) / step_size;
	double S 			= 0;
	double a,b;
	for (int i =0; i+1 < step_size; i++){ //simpsons method
		a 	= A + delta*(i);
		b 	= A + delta*(i+1);
		S+=((b-a)/6.)*(depletion_func(a, MU, SI, A, B) + 4*depletion_func((a+b)/2., MU, SI, A, B) +depletion_func(b, MU, SI, A, B));
	}
	return S;


}


map<string, vector<double> > get_stats(map<string, vector<segment>> query, 
	double a, double b , map<string, vector<double> > & distances){
	map<string, vector<double> > stats; 	

	//get_distances_by_motif
	get_distances(query, distances);
	typedef map<string, vector<double> >::iterator it_type;

	int BINS=200;

	double mu = 0 , si 	= 300;
	double C = normal_constant(mu, si, a, b);

	for (it_type m = distances.begin(); m!=distances.end(); m++){
		if (m->second.size() > 0){
			vector<double> current_stats(8);
			double ** X = NULL;
			BIN(m->second, BINS,X);
			double w_norm 	= 0.1, w_depletion=0.1, w_noise=1.0, ll_norm = 0, ll_depletion=0, ll_noise=0;
			EM(X, BINS, w_norm, ll_norm , si,mu, 0, C, a,b);
			double w_si 	= w_norm, ll_si 	= 0;
			double SI 		= EM(X, BINS, w_si, ll_si, si,mu, 1, C, a,b);


			EM(X, BINS, w_depletion, ll_depletion , si,mu, 2, C, a,b);
			EM(X, BINS, w_noise, ll_noise , si,mu, 3, C, a,b);
			if (X!=NULL){
				for (int i = 0; i < BINS; i++){
					delete X[i];
				}
			}
			current_stats[0]=w_norm,current_stats[1]=ll_norm;
			
			current_stats[2]=w_si,current_stats[3]=SI;

			current_stats[4]=w_depletion,current_stats[5]=ll_depletion;
			current_stats[6]=0.0,current_stats[7]=ll_noise;

			stats[m->first] 	= current_stats;

		}

	}
	return stats; 
}

