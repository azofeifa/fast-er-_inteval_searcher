#include "split.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
using namespace std;

vector<string> splitter(string ELE, string D){
	int j = 0;
	vector<string> results;
	const char *d =D.c_str();
	while (not ELE.empty() and j != ELE.size()){
		if (ELE[j] == *d){
			results.push_back(ELE.substr(0,j));
			ELE=ELE.substr(j+1,ELE.size());
			j=0;
		}
		j++;
	}
	results.push_back(ELE.substr(0,j));

	return results;
}

vector<string> splitter2(string line, string delim){
	vector<string> tokens;
	istringstream iss(line);
	string token;
	while(std::getline(iss, token, '\t' )){   // but we can specify a different one
		tokens.push_back(token);
	}

	return tokens;
}


string strip(string ELE, string D){
	const char *d 	= D.c_str();
	string result 	= "";
	for (int i = 0; i < ELE.size(); i++){
		if (ELE[i]==*d){
			break;
		}else{
			result+=ELE[i];
		}
	}
	return result;
}

string join(vector<string> toBeJoined, string delim){
	typedef vector<string>::iterator vs_it;
	string result="";
	for (vs_it i =toBeJoined.begin(); i != toBeJoined.end(); i++){
		result=result + delim + *i;
	}
	return result;
}	
