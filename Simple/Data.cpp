#include "Data.h"
#include <fstream>
#include <iostream>

using namespace std;

Data Data::instance;

Data::Data()
{

}

void Data::load(const char* filename)
{
	fstream fin(filename, ios::in);

	t.clear(); Y.clear();
	double temp1, temp2;
	while(fin>>temp1 && fin>>temp2)
	{
		t.push_back(temp1);
		Y.push_back(temp2);
	}
	fin.close();

}

