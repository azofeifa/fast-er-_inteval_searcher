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
	string log_out 			= P->p["-log_out"];
	int pairwise 			= stoi(P->p["-pairwise"]);
	if (log_out.empty()){
		log_out 	= out_directory;
	}
	ofstream FHW;
	FHW.open(log_out+ "tmp-"+job_name+".log" );
	vector<string> FILE_NAMES;
	vector<map<string, node *>> DBS	= load_input_directory(input_directory,  FILE_NAMES, FHW);
	if (not pairwise){

		map<string, vector<segment>> query;
		int q=0;
		load_DB(query_file, query, q,"");

		search_overlaps( query,  DBS, out_directory, job_name, FHW );
	}else{
		FHW<<"(main) computing pairwise overlaps"<<endl;
		compute_pairwise( DBS, FILE_NAMES, FHW,  out_directory);
	}

	string name 	= log_out+ "tmp-"+job_name+".log" ;
	if( remove( name.c_str()) != 0 )
		perror( "Error deleting file" );
	
	return 1;
}