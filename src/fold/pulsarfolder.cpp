/**
 * @author Yunpeng Men
 * @email ypmen@mpifr-bonn.mpg.de
 * @create date 2020-11-21 17:28:46
 * @modify date 2020-11-21 17:28:46
 * @desc [description]
 */

#include "dedisperse.h"
#include "pulsarfolder.h"
#include "utils.h"

using namespace Baseband;

PulsarFolder::PulsarFolder()
{
    npol = 0;
    nbin = 0;
    nsblk = 0;
    fold_period = 0.;
    
    nchan = 0;
    nsubint = 0;
    tsubint = 0.;
    offs_sub = 0.;
    period = 0.;
    freqref = 0.;
}

PulsarFolder::PulsarFolder(const std::string &fname) : predictor(fname)
{
    npol = 0;
    nbin = 0;
    nsblk = 0;
    fold_period = 0.;
    
    nchan = 0;
    nsubint = 0;
    tsubint = 0.;
    offs_sub = 0.;
    period = 0.;
    freqref = 0.;
}

PulsarFolder::~PulsarFolder()
{}

void PulsarFolder::prepare(DataBuffer<float> &databuffer)
{
    nchan = databuffer.nchans/npol;

    frequencies.resize(nchan, 0.);

    double fmax=-1e6, fmin=1e6;
    for (long int j=0; j<nchan; j++)
    {
        frequencies[j] = databuffer.frequencies[j];
        fmax = frequencies[j]>fmax ? frequencies[j]:fmax;
        fmin = frequencies[j]<fmin ? frequencies[j]:fmin;
    }

	if (predictor.empty())
	{
		freqref = (fmax+fmin)/2.;
	}
	else
	{
		freqref = predictor.get_center_frequency();
	}

    profiles.resize(npol*nchan*nbin, 0.);
    profilesTPF.resize(nbin*nchan*npol, 0.);

    hits.resize(nbin, 0);
}

void PulsarFolder::run(DataBuffer<float> &databuffer)
{
    double phi = 0.;
    double pfold = 0.;

    fill(profilesTPF.begin(), profilesTPF.end(), 0.);
    fill(hits.begin(), hits.end(), 0);

    long int i=0;
    while (i<databuffer.nsamples)
    {
        /**
         * @brief calculate block start phi 
         * 
         */
        if (i%nsblk == 0)
        {
            MJD block_epoch = start_epoch + (i/nsblk*nsblk+0.5)*databuffer.tsamp;
            if (!predictor.empty())
            {
                phi = predictor.get_fphase(block_epoch.to_day(), freqref);
                pfold = predictor.get_pfold(block_epoch.to_day(), freqref);
            }
            else
            {
                pfold = fold_period;
                long double tstart = block_epoch.to_second();
                long int nperiod = tstart/pfold;
                phi = (tstart - nperiod*pfold)/pfold;
            }
        }

        phi -= floor(phi);
        int ibin = int(phi*nbin);
        hits[ibin]++;
        phi += databuffer.tsamp/pfold;

        /**
         * @brief folding data
         * 
         */
        for (long int k=0; k<npol; k++)
        {
            for (long int j=0; j<nchan; j++)
            {
                profilesTPF[ibin*npol*nchan+k*nchan+j] += databuffer.buffer[i*npol*nchan+k*nchan+j];
            }
        }

        i++;
    }

    transpose_pad(&profiles[0], &profilesTPF[0], nbin, npol*nchan);    

    /**
     * @brief set zero to mean
     * 
     * 
     */
    for (long int ipol=0; ipol<npol; ipol++)
    {
        for (long int ichan=0; ichan<nchan; ichan++)
        {
            double mean = 0;
            long int count = 0;
            for (long int ibin=0; ibin<nbin; ibin++)
            {
                profiles[ipol*nchan*nbin+ichan*nbin+ibin] = profilesTPF[ibin*npol*nchan+ipol*nchan+ichan];
                if (hits[ibin] != 0)
                {
                    profiles[ipol*nchan*nbin+ichan*nbin+ibin] /= hits[ibin];
                    mean += profiles[ipol*nchan*nbin+ichan*nbin+ibin];
                    count++;
                }
            }
            mean /= count;
            for (long int ibin=0; ibin<nbin; ibin++)
            {
                if (hits[ibin] == 0)
                {
                    profiles[ipol*nchan*nbin+ichan*nbin+ibin] = mean;
                }
            }
        }
    }

    /**
     * @brief calculate tsubint,offs_sub,period
     * 
     */
    tsubint = databuffer.nsamples*databuffer.tsamp;
    MJD end_epoch = start_epoch + tsubint; 
    if (predictor.empty())
    {
        period = fold_period;

        long double med_epoch = (start_epoch + end_epoch).to_day()/2.;
        long double phase = fmod(med_epoch*86400., fold_period)/fold_period;
        epoch.format(med_epoch-phase*fold_period/86400.);
    }
    else
    {
        period = 0.;

        long double phase = floor(predictor.get_phase((start_epoch + end_epoch).to_day()/2., freqref));
        epoch.format(predictor.get_phase_inverse(phase, freqref));
    }

    nsubint++;
}
