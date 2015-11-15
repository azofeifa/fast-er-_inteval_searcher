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

void search_overlaps(map<string, vector<segment>> query, vector<map<string, node>> DBS, 
	string out, string job_name, ofstream& FHW ){
	typedef map<string, vector<segment> >::iterator it_type;

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
					DBS[t][c->first].searchInterval(q[i].start, q[i].stop,FINDS);
				}
			}
			q[i].overlaps 	= FINDS;
		}
		FHW<<"done...\n";
		FHW.flush();
		query[c->first]=q;

	
	}
	write_out(query, out, job_name);


}