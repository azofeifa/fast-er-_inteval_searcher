#include <iostream>
#include <fstream>
#include <unistd.h>
#include <string>
#include <vector>
#include <algorithm>
#include "read_in_parameters.h"
#include <stdio.h>
#include <ctype.h>

#include <stdio.h>
#include <time.h>
#include "split.h"
using namespace std;

params::params(){
	p["-rebuild"] 	="0";
	p["-pairwise"] 	= "0";
	p["-i"] 	= "";
	p["-q"] 	= "";
	p["-o"] 	= "";
	p["-pad"] 	= "0";
	p["-upad"] 	= "0";
	p["-N"] 	= "FIS_RUN";
}
const std::string currentDateTime() {
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    // Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
    // for more information about date/time format
    strftime(buf, sizeof(buf), "%m/%d/%Y %X", &tstruct);

    return buf;
}

void params::display(){
	cout<<"----------------------------------------------------------------"<<endl;
	
	cout<<"                Fast(er) Interval Searcher                      "<<endl;
	cout<<"Date Time   : "<<currentDateTime()<<endl;
	cout<<"-N          : "<<p["-N"]<<endl;

	cout<<"-i          : "<<p["-i"]<<endl;
	cout<<"-q          : "<<p["-q"]<<endl;
	cout<<"-o          : "<<p["-o"]<<endl;


	cout<<"Bugs/Questions? joseph[dot]azofeifa[at]colorado[dot]edu"<<endl;
	cout<<"----------------------------------------------------------------"<<endl;
}


void fillInOptions(char* argv[],params * P){
	string F 		= "";
	char * COM 		= "-";
	bool begin 		= true;
	bool GO_FORIT 	= false;
	bool INSERT 	= false;
	while (*argv){
		if ((*argv)[0] == COM[0] and GO_FORIT ){
			F 		=(*argv);
			INSERT 	=true;
		}else if(INSERT){
			P->p[F] 	= *argv;
			INSERT 		= false;
		}
		
		argv++;
		GO_FORIT=true;
	}

}


params * readInParameters( char* argv[]){	
	string userModParameter = "";
	params 	* P = new params;
	fillInOptions(argv, P);
	return P;
}


