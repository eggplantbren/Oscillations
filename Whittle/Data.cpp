/*
* Copyright (c) 2009, 2010, 2011, 2012 Brendon J. Brewer.
*
* This file is part of DNest3.
*
* DNest3 is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* DNest3 is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with DNest3. If not, see <http://www.gnu.org/licenses/>.
*/

#include "Data.h"
#include <fstream>
#include <iostream>
#include <algorithm>
#include <cmath>

using namespace std;

Data Data::instance;

Data::Data()
:loaded(false)
{

}

void Data::load(const char* filename)
{
	f.clear();
	y.clear();

	fstream fin(filename, ios::in);
	if(!fin)
	{
		cerr<<"# ERROR: Cannot open file "<<filename<<"."<<endl;
		return;
	}

	// Skip comment lines at the top of the file
	while(fin.peek() == '#')
		fin.ignore(1000000, '\n');

	// Read in data
	double temp1, temp2;
	while(fin>>temp1 && fin>>temp2)
	{
		f.push_back(temp1);
		y.push_back(temp2);
	}
	cout<<"# Loaded "<<f.size()<<" points from file "<<filename<<"."<<endl;
	fin.close();

	computeSummaries();
	loaded = true;
}

void Data::computeSummaries()
{
	fRange = *max_element(f.begin(), f.end())
			- *min_element(f.begin(), f.end());
}

