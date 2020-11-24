/*
 * predictor.cpp
 *
 *  Created on: Feb 16, 2020
 *      Author: ypmen
 */

#include <iostream>
#include <fstream>
#include <string.h>
#include <math.h>

#include <boost/algorithm/string.hpp>

#include "predictor.h"
#include "cheby.h"

using namespace std;

Predictor::Predictor()
{
	dispersion_constant = 0.;
	mjd_low = 0.;
	mjd_up = 0.;
	freq_low = 0.;
	freq_up = 0.;
	ncoef_t = 0;
	ncoef_f = 0;
	coef = NULL;
	dcoef = NULL;
}

Predictor::~Predictor()
{
	if (coef != NULL)
	{
		delete [] coef;
		coef = NULL;
	}

	if (dcoef != NULL)
	{
		delete [] dcoef;
		dcoef = NULL;
	}
}

void Predictor::set_coef(long double *coefs, int m, int n)
{
	if (coef == NULL and dcoef == NULL)
	{
		coef = new long double [m*n];
		dcoef = new long double [m*n];
		ncoef_f = m;
		ncoef_t = n;
	}

	if (m == ncoef_f and n == ncoef_t)
	{
		memcpy(coef, coefs, sizeof(long double)*m*n);
		memcpy(dcoef, coef, sizeof(long double)*m*n);
	}
	else
	{
		delete [] coef;
		delete [] dcoef;
		coef = new long double [m*n];
		dcoef = new long double [m*n];
		ncoef_f = m;
		ncoef_t = n;
	}

	compute_dcoef();
}

void Predictor::compute_dcoef()
{
	dcheby(dcoef, ncoef_f, ncoef_t);
}

long double Predictor::get_phase(long double mjd, long double freq)
{
	if (mjd>=mjd_low and mjd<=mjd_up)
	{
		long double x = -1. + 2.*(mjd-mjd_low)/(mjd_up-mjd_low);
		long double y = -1. + 2.*(freq-freq_low)/(freq_up-freq_low);
		long double phi = cheby2d_eval(coef, ncoef_f, ncoef_t, x, y) + dispersion_constant/(freq*freq);
		return phi;
	}
	else
	{
		cerr<<"Error: out of date"<<endl;
		return 0;
	}
}

double Predictor::get_pfold(long double mjd, long double freq)
{
	if (mjd>=mjd_low and mjd<=mjd_up)
	{
		long double x = -1. + 2.*(mjd-mjd_low)/(mjd_up-mjd_low);
		long double y = -1. + 2.*(freq-freq_low)/(freq_up-freq_low);
		return 1./(2./((mjd_up-mjd_low)*86400.)*cheby2d_eval(dcoef, ncoef_f, ncoef_t, x, y));
	}
	else
	{
		cerr<<"Error: out of date"<<endl;
		return 0;
	}
}

long double Predictor::get_phase_inverse(long double phase, long double freq)
{
	long double mjd = (mjd_low+mjd_up)/2.;
	int k = 0;
	long double dt = 1.;
	while (abs(dt)>1e-9)
	{
		dt = -(get_phase(mjd, freq)-phase)*get_pfold(mjd, freq);
		mjd += dt/86400.;
		if (++k > 8) break;
	}
	return mjd;
}

Predictors::Predictors()
{
	npred = 0;
	pred = NULL;
}

Predictors::Predictors(const string fname)
{
	read_t2pred(fname);
}

Predictors::~Predictors()
{
	if (pred != NULL)
	{
		delete [] pred;
		pred = NULL;
	}
}

int Predictors::get_nearest(long double mjd)
{
	int best_offset = 1e6;
	int inearest = -1;
	for (long int i=0; i<npred; i++)
	{
		long double mjd_low = pred[i].mjd_low;
		long double mjd_up = pred[i].mjd_up;
		if (mjd>=mjd_low and mjd<=mjd_up)
		{
			long double offset = abs((mjd_low+mjd_up)/2.-mjd);
			if (offset < best_offset)
			{
				inearest = i;
				best_offset = offset;
			}
		}
	}
	if (inearest == -1)
	{
		cerr<<"Error: out of date"<<endl;
	}
	return inearest;
}

long double Predictors::get_phase(long double mjd, long double freq)
{
	int k = get_nearest(mjd);
	return pred[k].get_phase(mjd, freq);
}

double Predictors::get_fphase(long double mjd, long double freq)
{
	long double phi = get_phase(mjd, freq);
	return phi-floor(phi);
}

double Predictors::get_pfold(long double mjd, long double freq)
{
	int k = get_nearest(mjd);
	return pred[k].get_pfold(mjd, freq);
}

long double Predictors::get_phase_inverse(long double phase, long double freq)
{
	long int k = 0;
	for (k=0; k<npred; k++)
	{
		if (phase >= pred[k].get_phase(pred[k].mjd_low, freq) and phase <= pred[k].get_phase(pred[k].mjd_up, freq))
			break;
	}

	if (k != npred)
	{
		return pred[k].get_phase_inverse(phase, freq);
	}
	else
	{
		cerr<<"out of phase"<<endl;
		return 0;
	}
}

long double Predictors::get_center_frequency()
{
	return (pred[0].freq_low+pred->freq_up)*0.5;
}

void Predictors::read_t2pred(const string fname)
{
	filename = fname;

	ifstream fpred(fname);
	string line;
	vector<string> items;
	int ncoeff_time = 0;
	int ncoeff_freq = 0;
	long double *coeffs = NULL;
	long int k = 0;
	Predictor *p = NULL;
	while (getline(fpred, line))
	{
		boost::split(items, line, boost::is_any_of(" "));
		if (items[0].compare("ChebyModelSet") == 0)
		{
			npred = stoi(items[1]);
			if (pred != NULL) delete [] pred;
			pred = new Predictor [npred];
			p = pred;
		}
		else if (items[0].compare("ChebyModel") == 0 and items[1].compare("BEGIN") == 0)
		{
			k = 0;
		}
		else if (items[0].compare("PSRNAME") == 0)
		{
			p->psrname = items[1];
		}
		else if (items[0].compare("SITENAME") == 0)
		{
			p->sitename = items[1];
		}
		else if (items[0].compare("TIME_RANGE") == 0)
		{
			p->mjd_low = stold(items[1]);
			p->mjd_up = stold(items[2]);
		}
		else if (items[0].compare("FREQ_RANGE") == 0)
		{
			p->freq_low = stold(items[1]);
			p->freq_up = stold(items[2]);
		}
		else if (items[0].compare("DISPERSION_CONSTANT") == 0)
		{
			p->dispersion_constant = stold(items[1]);
		}
		else if (items[0].compare("NCOEFF_TIME") == 0)
		{
			ncoeff_time = stoi(items[1]);
		}
		else if (items[0].compare("NCOEFF_FREQ") == 0)
		{
			ncoeff_freq = stoi(items[1]);
			coeffs = new long double [ncoeff_freq*ncoeff_time];
		}
		else if (items[0].compare("COEFFS") == 0)
		{
			for (long int i=1; i<items.size(); i++)
			{
				coeffs[(i-1)*ncoeff_time+k] = stold(items[i]);
			}
			k++;
		}
		else if (items[0].compare("ChebyModel") == 0 and items[1].compare("END") == 0)
		{
			p->set_coef(coeffs, ncoeff_freq, ncoeff_time);
			delete [] coeffs;
			coeffs = NULL;
			p++;
		}
	}
	p = NULL;
}
