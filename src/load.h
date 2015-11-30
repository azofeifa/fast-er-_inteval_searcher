#ifndef load_H
#define load_H
#include <string>
#include <vector>
#include <map>
using namespace std;
class segment{
public:
	int start, stop, ID;
	string chrom,info;
	segment();
	segment(string, int,int,string, int);
	vector<segment> overlaps;
};

class node{
public:
	double center;
	int start, stop;
	node * left;
	node * right;
	vector<segment> current;
	void retrieve_nodes(vector<segment>&);
	node();
	~node();
	node(vector<segment>);
	void searchInterval(int, int, vector<segment>&);
};


vector<map<string, node *> > load_input_directory(string, vector<string> &, ofstream&, int, int);
void load_DB(string, map<string, vector<segment>>& , int &, string , int, int);

#endif
