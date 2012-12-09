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

#include "Whittle.h"
#include "RandomNumberGenerator.h"
#include "Utils.h"
#include "Data.h"
#include <cmath>
#include <iostream>

using namespace std;
using namespace DNest3;

const int Whittle::maxNumComponents = 200;

Whittle::Whittle()
{

}

void Whittle::fromPrior()
{
	if(!Data::get_instance().get_loaded())
		cerr<<"# Warning: No data loaded!"<<endl;
	mockData.assign(Data::get_instance().get_N(), 0.);

	// Set limits on muAmplitudes
	minLogMu = log(1E-3*Data::get_instance().get_y_mean());
	maxLogMu = minLogMu + log(1E6);
	rangeLogMu = maxLogMu - minLogMu;

	// Set limits on frequencies
	minFreq = Data::get_instance().get_x_min();
	maxFreq = Data::get_instance().get_x_max();
	rangeFreq = maxFreq - minFreq;

	muAmplitudes = exp(minLogMu + rangeLogMu*randomU());
	numComponents = randInt(maxNumComponents + 1);

	frequencies.clear();
	amplitudes.clear();
	widths.clear();

	double A, f, width;
	for(int i=0; i<numComponents; i++)
	{
		A = -muAmplitudes*log(randomU());
		f = minFreq + rangeFreq*randomU();
		width = 3.*randomU();

		addComponent(A, f, width);
		amplitudes.push_back(A);
		frequencies.push_back(f);
		widths.push_back(width);
	}

}

double Whittle::perturb1()
{
	// Make a proposal for the new number of components
	int diff = static_cast<int>
			(round(maxNumComponents*pow(10., 1.5 - 6.*randomU())*randn()));
	if(diff == 0)
		diff = (randomU() <= 0.5)?(-1):(1);
	int proposal = numComponents + diff;
	proposal = mod(proposal, maxNumComponents + 1);
	int actual_diff = proposal - numComponents;

	if(actual_diff > 0)
	{
		double A, f, width;
		for(int i=0; i<actual_diff; i++)
		{
			A = -muAmplitudes*log(randomU());
			f = minFreq + rangeFreq*randomU();
			width = 0.3*randomU();

			addComponent(A, f, width);
			amplitudes.push_back(A);
			frequencies.push_back(f);
			widths.push_back(width);
			numComponents++;
		}
	}
	else if(actual_diff < 0)
	{
		int which;
		for(int i=0; i<-actual_diff; i++)
		{
			which = randInt(numComponents);
			addComponent(-amplitudes[which],
					frequencies[which], widths[which]);
			amplitudes.erase(amplitudes.begin() + which);
			frequencies.erase(frequencies.begin() + which);
			widths.erase(widths.begin() + which);
			numComponents--;
		}
	}

	staleness++;
	return 0.;
}

double Whittle::perturb2()
{
	double chance = pow(10., 0.5 - 4.*randomU());
	double scale = pow(10., 1.5 - 6.*randomU());

	int which = randInt(3);
	double temp;

	if(which == 0)
	{
		for(int i=0; i<numComponents; i++)
		{
			if(randomU() <= chance)
			{
				if(chance < 1.)
					addComponent(-amplitudes[i], frequencies[i],
						widths[i]);
				temp = 1. - exp(-amplitudes[i]/muAmplitudes);
				temp += scale*randn();
				temp = mod(temp, 1.);
				amplitudes[i] = -muAmplitudes*log(1. - temp);
				if(chance < 1.)
					addComponent(amplitudes[i], frequencies[i],
						widths[i]);
			}
		}
	}
	else if(which == 1)
	{
		for(int i=0; i<numComponents; i++)
		{
			if(randomU() <= chance)
			{
				if(chance < 1.)
					addComponent(-amplitudes[i], frequencies[i],
						widths[i]);

				frequencies[i] += scale*rangeFreq*randn();
				frequencies[i] = mod(frequencies[i] - minFreq, rangeFreq)
							+ minFreq;

				if(chance < 1.)
					addComponent(amplitudes[i], frequencies[i],
						widths[i]);
			}
		}
	}
	else
	{
		for(int i=0; i<numComponents; i++)
		{
			if(randomU() <= chance)
			{
				if(chance < 1.)
					addComponent(-amplitudes[i], frequencies[i],
						widths[i]);
				widths[i] += 3.*scale*randn();
				widths[i] = mod(widths[i], 3.);
				if(chance < 1.)
					addComponent(amplitudes[i], frequencies[i],
						widths[i]);
			}
		}
	}

	if(chance < 1.)
		staleness++;
	else
		calculateMockData();
	return 0.;
}

double Whittle::perturb3()
{
	double logH = 0.;

	double proposal = muAmplitudes;
	proposal = log(proposal);
	proposal += rangeLogMu*pow(10., 1.5 - 6.*randomU())*randn();
	proposal = mod(proposal - minLogMu, rangeLogMu) + minLogMu;
	proposal = exp(proposal);

	int which = randInt(2);
	if(which == 0)
	{
		double ratio = proposal/muAmplitudes;
		for(int i=0; i<numComponents; i++)
			amplitudes[i] *= ratio;
		for(size_t i=0; i<mockData.size(); i++)
			mockData[i] *= ratio;
		staleness++;
	}
	else
	{
		double logMu1 = log(muAmplitudes);
		double logMu2 = log(proposal);
		for(int i=0; i<numComponents; i++)
		{
			if(amplitudes[i] < 0.)
				cerr<<"# ERROR: Negative amplitude."<<endl;
			logH -= -logMu1 - amplitudes[i]/muAmplitudes;
			logH += -logMu2 - amplitudes[i]/proposal;
		}
	}

	muAmplitudes = proposal;
	return logH;
}

double Whittle::perturb()
{
	double logH = 0.;

	int which = randInt(3);
	if(which == 0)
	{
		logH = perturb1();
	}
	else if(which == 1)
	{
		logH = perturb2();
	}
	else if(which == 2)
	{
		logH = perturb3();
	}

	if(staleness > 1000)
		calculateMockData();

	return logH;
}

double Whittle::logLikelihood() const
{
	double logL = 0.;

	for(size_t i=0; i<mockData.size(); i++)
		logL += -log(mockData[i]) - Data::get_instance().get_y(i)/mockData[i];

	return logL;
}

void Whittle::calculateMockData()
{
	// Zero the mock data
	// Actually, put in known background
	for(size_t i=0; i<mockData.size(); i++)
		mockData[i] = 0.2;

	// Add each frequency
	for(int i=0; i<numComponents; i++)
		addComponent(amplitudes[i], frequencies[i], widths[i]);

	staleness = 0;
}

void Whittle::addComponent(double amplitude, double frequency, double width)
{
	if(amplitude == 0.)
		return;

	for(size_t i=0; i<mockData.size(); i++)
	{
		mockData[i] += amplitude/(1. + pow(
				(Data::get_instance().get_x(i) - frequency)/width
				, 2));
	}
}

void Whittle::print(std::ostream& out) const
{
	out<<numComponents<<' '<<muAmplitudes<<' '<<staleness<<' ';

	// Print amplitudes, use zero padding
	for(int i=0; i<numComponents; i++)
		out<<amplitudes[i]<<' ';
	for(int i=numComponents; i<maxNumComponents; i++)
		out<<0<<' ';

	// Print frequencies, use zero padding
	for(int i=0; i<numComponents; i++)
		out<<frequencies[i]<<' ';
	for(int i=numComponents; i<maxNumComponents; i++)
		out<<0<<' ';

	// Print widths, use zero padding
	for(int i=0; i<numComponents; i++)
		out<<widths[i]<<' ';
	for(int i=numComponents; i<maxNumComponents; i++)
		out<<0<<' ';
}

string Whittle::description() const
{
	string result("numComponents, muAmplitudes, staleness, ");
	result += "amplitudes, frequencies, widths.";
	return result;
}


