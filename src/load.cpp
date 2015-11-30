#include "load.h"
#include <stdio.h>
#include <dirent.h>
#include <string>
#include <map>
#include <vector>
#include "split.h"
#include <iostream>
#include <fstream>
#include <cctype>
#include <algorithm>
#include <omp.h>  
using namespace std;

segment::segment(){}
segment::segment(string c, int st, int sp, string inf, int i){
	chrom=c,info=inf, start=st, stop=sp, ID=i;
}
bool isNumeric(const string& input) {
    return all_of(input.begin(), input.end(), ::isdigit);
}


void load_DB(string filename, map<string, vector<segment> >& DB, int & i, string label, int upad, int pad ){
	ifstream FH(filename);
	map<string, vector<segment> > G;
	bool EXIT 	= false;
		
	if (FH){
		vector<string> lineArray;
		string line, chrom, info;
		int start, stop;
		string prevChrom="";
		while(getline(FH, line) and not EXIT){

			if (line.substr(0,1)!="#"){
				lineArray=splitter2(line, "\t");
				if (lineArray.size() > 2){
					if (isNumeric(lineArray[1]) and isNumeric(lineArray[2])  ){
						chrom 	= lineArray[0];
						start 	= stoi(lineArray[1]);
						stop 	= stoi(lineArray[2]);
						if (upad){
							double x 	= (start + stop)/2.;
							start 	= int(x - upad), stop = int(x + upad);
						}else if(pad){
							start-=pad, stop+=pad;
						}


						info 	= "";
						if (lineArray.size()>3){
							info 	= lineArray[3] + ":"+label;
						}
						segment S(chrom,start, stop, info, i);
						G[chrom].push_back(S);
						i+=1;
					}else if(isNumeric(lineArray[2]) and isNumeric(lineArray[3])  ){
						chrom 	= lineArray[1];
						start 	= stoi(lineArray[2]);
						stop 	= stoi(lineArray[3]);
						info 	= lineArray[0] + ":"+ label;
						segment S(chrom,start, stop, info, i);
						i+=1;
						G[chrom].push_back(S);
					}
					else{
						EXIT=true;
					}

				}else{
					EXIT 	= true;
				}
			}
		}
		if (EXIT){
			printf("ignoring non-bed formatted file: %s\n", filename.c_str() );
		}
	}else{
		EXIT=true;
		printf("couldn't open %s\n", filename.c_str() );
	}
	if (not EXIT){
		typedef map<string, vector<segment> >::iterator it_type;
		for (it_type c = G.begin(); c!=G.end(); c++){
			for (int i = 0; i < c->second.size(); i++){
				DB[c->first].push_back(c->second[i]);
			}
		}
	}	
}
node::node(){};

node::node(vector<segment> segments ){
	center 	= (double(segments[0].start)  + double(segments[segments.size()-1].stop)) / 2.;
	vector<segment> Left;
	vector<segment> Right;
	left=NULL, right=NULL;
	for (int i = 0 ; i < segments.size(); i++){
		if (segments[i].stop < center){
			Left.push_back(segments[i]);
		}
		else if (segments[i].start > center){
			Right.push_back(segments[i]);
		}
		else{
			current.push_back(segments[i]);
		}
	}
	if (Left.size() > 0){
		left 	= new node(Left);
	}
	if (Right.size() > 0){
		right 	= new node(Right);
	}
}

node::~node(){
	if (left != NULL){
		delete left;
	}
	if (right != NULL){
		delete right;
	}
}

void node::searchInterval(int st, int sp, vector<segment> & FINDS){
	for (int i = 0 ; i < current.size(); i++){
		if (sp > current[i].start and  st < current[i].stop  ){
			FINDS.push_back(current[i]);
		}
	}	
	if (sp >= center and right != NULL ){
		right->searchInterval(st, sp, FINDS);
	}
	if (st <= center and left !=NULL){
		left->searchInterval(st, sp, FINDS);
	}	
}

void node::retrieve_nodes(vector<segment> & saves){
	for (int i = 0; i < current.size(); i++){
		saves.push_back(current[i]);
	}
	if (right!= NULL){
		right->retrieve_nodes(saves);
	}
	if (left != NULL){
		left->retrieve_nodes(saves);		
	}

}


//sort database
void bubble_sort(map<string, vector<segment> > G, 
	map<string ,vector<segment>>& by_start ){
	typedef map<string, vector<segment> >::iterator it_type;
	for (it_type c = G.begin(); c!=G.end(); c++){
		bool GOOD=false;
		vector<segment> A 	= c->second;
		vector<segment> B 	= c->second;
		while (not GOOD){
			GOOD=true;
			//by starting position
			for (int i =1; i < A.size(); i++){
				if (A[i].start < A[i-1].start){
					segment cp 				= A[i-1];
					A[i-1] 			= A[i];
					A[i] 			= cp;
					GOOD=false;
				}
			}			
		}
		by_start[c->first]=A;
		
	}
}



vector<map<string, node * > > load_input_directory(string path, vector<string>&  FILE_NAMES, ofstream& FHW, int upad, int pad){
	struct dirent *entry;
	DIR *dp;

	vector<map<string, node *>> TS;
	dp = opendir(path.c_str());
	if (dp == NULL) {
		perror("opendir: Path does not exist or could not be read.");
		return TS;
	}
	int i 	= 0;
	typedef map<string, vector<segment> >::iterator it_type;
	while ((entry = readdir(dp))){
		string current_file_name 	= entry->d_name;
		map<string, node *> T;
		map<string, vector<segment> > DB;
		load_DB(path+current_file_name, DB, i, current_file_name, 0, 0);
		vector<string> chromosomes;
		for (it_type c = DB.begin(); c!= DB.end(); c++){
			chromosomes.push_back(c->first);
		}
		int N 	= chromosomes.size();
		if ( N > 0 ){
			FILE_NAMES.push_back(current_file_name);
			FHW<<"loading: " + current_file_name +"...";
			FHW.flush();
			#pragma omp parallel for  
			for (int c= 0; c <  N; c++){
				T[chromosomes[c]]= new node(DB[chromosomes[c]] );
			}
			FHW<<"done\n";
			FHW.flush();
			TS.push_back(T);
		}

	}

	
	closedir(dp);
	return TS;
}




