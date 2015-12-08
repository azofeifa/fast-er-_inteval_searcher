#include "read_in_parameters.h"
#include "load.h"
#include <omp.h>
#include <string>
#include <map>
#include "find_overlaps.h"
#include <string>
#include <iostream>
#include <fstream>
using namespace std;
int main(int argc, char* argv[]){
	params * P 	= new params();
    P 			= readInParameters(argv);
    P->display();

	string input_directory 	= P->p["-i"];
	string query_file 		= P->p["-q"];
	string out_directory 	= P->p["-o"];
	string job_name 		= P->p["-N"];
	int upad 				= stoi(P->p["-upad"]);
	int pad 				= stoi(P->p["-pad"]);
	int pairwise 			= stoi(P->p["-pairwise"]);
	vector<string> FILE_NAMES;
	vector<map<string, node *>> DBS	= load_input_directory(input_directory,  FILE_NAMES, upad, pad);
	if (not pairwise){

		map<string, vector<segment>> query;
		int q=0;
		load_DB(query_file, query, q,"", upad, pad);

		search_overlaps( query,  DBS, out_directory, job_name, upad );
	}else{
		compute_pairwise( DBS, FILE_NAMES ,  out_directory);
	}

	
	return 1;
}