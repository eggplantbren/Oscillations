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
#include <gsl/gsl_sf_gamma.h>
using namespace std;
using namespace DNest3;

const int Whittle::maxNumComponents = 20;

Whittle::Whittle()
{

}

void Whittle::fromPrior()
{
	if(!Data::get_instance().get_loaded())
		cerr<<"# Warning: No data loaded!"<<endl;
	mockData.assign(Data::get_instance().get_N(), 0.);

	// Set limits on muAmplitudes
	minLogMu = log(1E-3*Data::get_instance().get_ySig());
	maxLogMu = minLogMu + log(1E4);
	rangeLogMu = maxLogMu - minLogMu;

	// Set limits on frequencies
	// Longest possible period = 10*data range
	// Shortest possible period = 1E-4*longest period
	minLogFreq = log(1./(10.*Data::get_instance().get_tRange()));
	maxLogFreq = log(1E4) + minLogFreq;
	rangeLogFreq = maxLogFreq - minLogFreq;

	muAmplitudes = exp(minLogMu + rangeLogMu*randomU());
	numComponents = randInt(maxNumComponents + 1);

	frequencies.clear();
	amplitudes.clear();
	phases.clear();
	double A, f, phi;
	for(int i=0; i<numComponents; i++)
	{
		A = -muAmplitudes*log(randomU());
		f = exp(minLogFreq + rangeLogFreq*randomU());
		phi = 2*M_PI*randomU();

		addComponent(A, f, phi);
		amplitudes.push_back(A);
		frequencies.push_back(f);
		phases.push_back(phi);
	}

	sigmaBoost = exp(log(1.) + log(100.)*randomU());
	nu = exp(log(0.5) + log(200.)*randomU());
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
		double A, f, phi;
		for(int i=0; i<actual_diff; i++)
		{
			A = -muAmplitudes*log(randomU());
			f = exp(minLogFreq + rangeLogFreq*randomU());
			phi = 2*M_PI*randomU();

			addComponent(A, f, phi);
			amplitudes.push_back(A);
			frequencies.push_back(f);
			phases.push_back(phi);
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
					frequencies[which], phases[which]);
			amplitudes.erase(amplitudes.begin() + which);
			frequencies.erase(frequencies.begin() + which);
			phases.erase(phases.begin() + which);
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
						phases[i]);
				temp = 1. - exp(-amplitudes[i]/muAmplitudes);
				temp += scale*randn();
				temp = mod(temp, 1.);
				amplitudes[i] = -muAmplitudes*log(1. - temp);
				if(chance < 1.)
					addComponent(amplitudes[i], frequencies[i],
						phases[i]);
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
						phases[i]);
				temp = log(frequencies[i]);
				temp += rangeLogFreq*scale*randn();
				temp = mod(temp - minLogFreq, rangeLogFreq)
					+ minLogFreq;
				frequencies[i] = exp(temp);
				if(chance < 1.)
					addComponent(amplitudes[i], frequencies[i],
						phases[i]);
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
						phases[i]);
				phases[i] += 2*M_PI*scale*randn();
				phases[i] = mod(phases[i], 2*M_PI);
				if(chance < 1.)
					addComponent(amplitudes[i], frequencies[i],
						phases[i]);
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

double Whittle::perturb4()
{
	double logH = 0.;
	int which = randInt(2);
	if(which == 0)
	{
		sigmaBoost = log(sigmaBoost);
		sigmaBoost += log(100.)*pow(10., 1.5 - 6.*randomU())*randn();
		sigmaBoost = mod(sigmaBoost - log(1.), log(100.)) + log(1.);
		sigmaBoost = exp(sigmaBoost);
	}
	else
	{
		nu = log(nu);
		nu += log(200.)*pow(10., 1.5 - 6.*randomU())*randn();
		nu = mod(nu - log(0.5), log(200.)) + log(0.5);
		nu = exp(nu);
	}
	return logH;
}

double Whittle::perturb()
{
	double logH = 0.;

	int which = randInt(4);
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
	else
	{
		logH = perturb4();
	}

	if(staleness > 1000)
		calculateMockData();

	return logH;
}

double Whittle::logLikelihood() const
{
	double logL = mockData.size()*
		(gsl_sf_lngamma((nu + 1.)/2) - gsl_sf_lngamma(nu/2)
		-0.5*log(M_PI*nu));

	for(size_t i=0; i<mockData.size(); i++)
	{
		double y = Data::get_instance().get_y(i);
		double sig = sigmaBoost*Data::get_instance().get_sig(i);
		logL += -log(sig);
		logL += -(nu+1.)/2*
			log(1. + 1./nu*pow((y - mockData[i])/sig, 2));
	}
	return logL;
}

void Whittle::calculateMockData()
{
	// Zero the mock data
	for(size_t i=0; i<mockData.size(); i++)
		mockData[i] = 0.;

	// Add each frequency
	for(int i=0; i<numComponents; i++)
		addComponent(amplitudes[i], frequencies[i], phases[i]);

	staleness = 0;
}

void Whittle::addComponent(double amplitude, double frequency, double phase)
{
	if(amplitude == 0.)
		return;

	for(size_t i=0; i<mockData.size(); i++)
		mockData[i] += amplitude*sin(2*M_PI*frequency*
				Data::get_instance().get_t(i) + phase);
}

void Whittle::print(std::ostream& out) const
{
	out<<numComponents<<' '<<muAmplitudes<<' ';
	out<<sigmaBoost<<' '<<nu<<' '<<staleness<<' ';

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

	// Print phases, use zero padding
	for(int i=0; i<numComponents; i++)
		out<<phases[i]<<' ';
	for(int i=numComponents; i<maxNumComponents; i++)
		out<<0<<' ';
}

string Whittle::description() const
{
	string result("numComponents, muAmplitudes, sigmaBoost, nu, staleness, ");
	result += "amplitudes, frequencies, phases.";
	return result;
}


