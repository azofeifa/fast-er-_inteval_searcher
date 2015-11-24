#include "load.h"
#include "find_overlaps.h"
#include <iostream>
#include <fstream>
#include <stdio.h>

using namespace std;

void write_out(map<string, vector<segment>> query  , string OUT, string job_name ){
	ofstream FHW;
	FHW.open(OUT+ job_name+ "_overlaps.bed");
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
}

void compute_pairwise(vector<map<string, node *>> DBS, 
		vector<string> FILE_NAMES, ofstream& FHW, string out_directory){
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
		if (double(i) / N > (percent + 0.05)){
			ct+=5;
			FHW<<to_string(ct) + "%,";
			FHW.flush();
			percent+=0.05;
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
							//printf("%d-%d,%d\n", a_chrom->second[s].start, a_chrom->second[s].stop, FINDS.size()  );
						}
						ct+=int(FINDS.size());
					}
				}
				counts[i][j]=ct;
			}
		}
	}
	FHW<<"done :)\n";
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
	string out, string job_name, ofstream& FHW ){
	typedef map<string, vector<segment> >::iterator it_type;
	typedef map<string, node * >::iterator it_type_2;

	for (it_type c= query.begin(); c!= query.end(); c++){
		vector<segment> q 	= c->second ;
		FHW<<"ID-ing overlaps for chromosome " + c->first+"...\n";
		FHW.flush();
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
		FHW<<"done...\n";
		FHW.flush();
		query[c->first]=q;

	
	}
	write_out(query, out, job_name);
	for (int t = 0; t < DBS.size(); t++){
		for (it_type_2 c=DBS[t].begin(); c!= DBS[t].end(); c++){
			delete DBS[t][c->first];
		}
	}

}