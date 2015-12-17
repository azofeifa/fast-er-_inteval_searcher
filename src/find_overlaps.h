#ifndef fo_H
#define fo_H
#include "load.h"
using namespace std;
void search_overlaps(map<string, vector<segment>>  , 
	vector<map<string, node * >>, string, string , int,int );
void compute_pairwise(vector<map<string, node *>> ,  
	vector<string> ,  string );	
#endif