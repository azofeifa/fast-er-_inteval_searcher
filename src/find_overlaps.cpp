#include "load.h"
#include "find_overlaps.h"
#include <iostream>
#include <fstream>
#include <stdio.h>
#include "EM.h"
using namespace std;

void write_out(map<string, vector<segment>> query  , string OUT, string job_name, int upad, int MIN ){
	//write out raw distances
	ofstream FHW;
	FHW.open(OUT+ job_name+ "_raw_distances.bed");
	typedef map<string, vector<segment> >::iterator it_type;
	for (it_type c= query.begin(); c!= query.end(); c++){
		vector<segment> q 	= c->second ;
		for (int i = 0 ; i < q.size(); i++){
			string INF 	= q[i].info+"|" ;
			int ON 		= q[i].overlaps.size();
			for (int o = 0; o < q[i].overlaps.size(); o++){
				if (o+1 < ON){
					INF+=to_string(q[i].overlaps[o].start)+"-"+to_string(q[i].overlaps[o].stop)+":"+(q[i].overlaps[o].info)+",";
				}else{
					INF+=to_string(q[i].overlaps[o].start)+"-"+to_string(q[i].overlaps[o].stop)+":"+(q[i].overlaps[o].info);

				}
			}
			FHW<<q[i].chrom+"\t"+to_string(q[i].start) + "\t" + to_string(q[i].stop)+ "\t" + INF<<endl;
		}
	}		
	//write out EM results for each motif
	ofstream FHW_stats;
	FHW_stats.open(OUT+ job_name+ "_stats.tsv");

	double a 	= -upad, b = upad;
	map<string, vector<double> > distances ;
	map<string, vector<double> > stats 		= get_stats(query,a,b,distances, MIN) ;
	typedef map<string, vector<double> >::iterator it_type_2;
	for (it_type_2 m = stats.begin(); m!=stats.end(); m++){
		if (m->second.size()==8){
			FHW_stats<<m->first+"\t" + to_string(m->second[0]) + ","+ to_string(m->second[1]) + "\t";
			FHW_stats<<to_string(m->second[2]) + ","+ to_string(m->second[3]) + "\t";
			FHW_stats<<to_string(m->second[4]) + ","+ to_string(m->second[5]) + "\t";
			FHW_stats<<to_string(m->second[6]) + ","+ to_string(m->second[7]) + "\n";
		}
	}
	FHW_stats<<"#begin distances\n";
	for (it_type_2 m = distances.begin(); m!=distances.end(); m++){
		string line = "";
		FHW_stats<<m->first<<"\t";
		for (int i = 0 ; i < m->second.size(); i++){
			if (i+1 < m->second.size() ){
				line+=to_string(m->second[i])+ ",";
			}else{
				line+=to_string(m->second[i]);
			}
		}
		FHW_stats<<line<<endl;
		
	}



}

void compute_pairwise(vector<map<string, node *>> DBS, 
		vector<string> FILE_NAMES,  string out_directory){
	int N 	= DBS.size();
	int counts[N][N];
	//initialize to zero
	for (int i = 0; i < N; i++){
		for (int j = 0; j < N; j++){
			counts[i][j]=0;
		}
	}
	//make pairwise counts
	typedef map<string, node * >::iterator it_type;
	typedef map<string, vector<segment> >::iterator it_type_2;
	double percent = 0;
	int ct 	 	= 0;
	for (int i = 0; i < N; i++){
		map<string, node * > A =  DBS[i];
		map<string, vector<segment> > current;
		int i_ct 	= 0;
		for (it_type a_chrom = A.begin(); a_chrom != A.end(); a_chrom++ ){
			vector<segment> CURRENT;
			a_chrom->second->retrieve_nodes(CURRENT);
			i_ct+=int(CURRENT.size());
			current[a_chrom->first]=CURRENT;

		}	
		counts[i][i] 	= i_ct;
		#pragma omp parallel for  
		for (int j = 0; j < N; j++){
			if (j != i){
				map<string, node * > B 	= DBS[j];
				int ct 					= 0;
				for (it_type_2 a_chrom = current.begin(); a_chrom != current.end(); a_chrom++ ){
					if (B.find(a_chrom->first)!=B.end()  ){
						vector<segment> FINDS;
						for (int s = 0; s < a_chrom->second.size(); s++){
							B[a_chrom->first]->searchInterval(a_chrom->second[s].start, a_chrom->second[s].stop, FINDS);
						}
						ct+=int(FINDS.size());
					}
				}
				counts[i][j]=ct;
			}
		}
	}
	//write out counts  and file names
	ofstream OUT;
	OUT.open(out_directory + "pairwise_count_matrix.csv");
	for (int i = 0; i < N; i++){
		OUT<<to_string(i) + "," + FILE_NAMES[i] + "\n";
	}
	OUT<<"#------------------------------------------"<<endl;
	for (int i = 0 ; i < N ; i++){
		string line = "";
		for (int j = 0; j < N; j++){
			if (j+1 < N){
				line+=to_string(counts[i][j])+",";
			}else{
				line+=to_string(counts[i][j])+ "\n";
			}
		}
		OUT<<line;
	}
	OUT.close();



}

void search_overlaps(map<string, vector<segment>> query, vector<map<string, node *>> DBS, 
	string out, string job_name, int upad, int MIN ){
	typedef map<string, vector<segment> >::iterator it_type;
	typedef map<string, node * >::iterator it_type_2;

	for (it_type c= query.begin(); c!= query.end(); c++){
		vector<segment> q 	= c->second ;
		int N 				= q.size();
		#pragma omp parallel for  
		for (int i = 0 ; i < q.size(); i++){
			vector<segment> FINDS;
			for (int t = 0 ; t < DBS.size(); t++){
				if (DBS[t].find(c->first)!= DBS[t].end() ){
					DBS[t][c->first]->searchInterval(q[i].start, q[i].stop,FINDS);
				}
			}
			q[i].overlaps 	= FINDS;
		}
		query[c->first]=q;
	}
	write_out(query, out, job_name, upad, MIN);
	for (int t = 0; t < DBS.size(); t++){
		for (it_type_2 c=DBS[t].begin(); c!= DBS[t].end(); c++){
			delete DBS[t][c->first];
		}
	}

}