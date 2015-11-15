#ifndef read_in_parameters_H
#define read_in_parameters_H

#include <iostream>
#include <fstream>
#include <unistd.h>
#include <string>
#include <vector>
#include <map>
using namespace std;


class params{
public:
	map<string, string> p;
	params();
	void display();
	void help();
	
	bool EXIT;
	
};


params * readInParameters(char**);

const std::string currentDateTime();
#endif