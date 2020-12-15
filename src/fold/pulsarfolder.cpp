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

#include <eigen3/Eigen/Dense>

using namespace Baseband;

PulsarFolder::PulsarFolder()
{
    npol = 0;
    nbin = 0;
    nsblk = 0;
    fold_period = 0.;
    lambda = 0.;
    
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
    lambda = 0.;
    
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
}

void PulsarFolder::run(DataBuffer<float> &databuffer)
{
    double phi = 0.;
    double pfold = 0.;

    std::vector<float> profilesTPF(nbin*nchan*npol, 0.);
    std::vector<int> hits(nbin, 0);

    long int i=0;
    while (i<databuffer.nsamples)
    {
        /**
         * @brief calculate block start phi 
         * 
         */
        if (i%nsblk == 0)
        {
            MJD block_epoch_start = start_epoch + (i/nsblk*nsblk+0.5)*databuffer.tsamp;
            MJD block_epoch_mid = start_epoch + (i/nsblk*nsblk+0.5*nsblk)*databuffer.tsamp;
            if (!predictor.empty())
            {
                phi = predictor.get_fphase(block_epoch_start.to_day(), freqref);
                pfold = predictor.get_pfold(block_epoch_mid.to_day(), freqref);
            }
            else
            {
                pfold = fold_period;
                long double tstart = block_epoch_start.to_second();
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

void PulsarFolder::runLSM(DataBuffer<float> &databuffer)
{
    /**
     * @brief calculate tsubint,offs_sub,period
     * 
     */

    tsubint = databuffer.nsamples*databuffer.tsamp;
    MJD end_epoch = start_epoch + tsubint;
    long double med_epoch = (start_epoch + end_epoch).to_day()/2.;

    if (predictor.empty())
    {
        period = fold_period;

        long double phase = fmod(med_epoch*86400., fold_period)/fold_period;
        epoch.format(med_epoch-phase*fold_period/86400.);
    }
    else
    {
        period = predictor.get_pfold(med_epoch, freqref);

        long double phase = floor(predictor.get_phase(med_epoch, freqref));
        epoch.format(predictor.get_phase_inverse(phase, freqref));
    }

    double phi = 0.;
    double pfold = 0.;

    std::vector<double> mxWTW(nbin*nbin, 0.);
    std::vector<double> vWTd_T(nbin*npol*nchan, 0.);

    long int i=0;
    while (i<databuffer.nsamples)
    {
        /**
         * @brief calculate block start phi 
         * 
         */
        if (i%nsblk == 0)
        {
            MJD block_epoch_start = start_epoch + (i/nsblk*nsblk+0.5)*databuffer.tsamp;
            MJD block_epoch_mid = start_epoch + (i/nsblk*nsblk+0.5*nsblk)*databuffer.tsamp;
            if (!predictor.empty())
            {
                phi = predictor.get_fphase(block_epoch_start.to_day(), freqref);
                pfold = predictor.get_pfold(block_epoch_mid.to_day(), freqref);
            }
            else
            {
                pfold = fold_period;
                long double tstart = block_epoch_start.to_second();
                long int nperiod = tstart/pfold;
                phi = (tstart - nperiod*pfold)/pfold;
            }
        }

        double low_phi = phi;
        double high_phi = phi + databuffer.tsamp/pfold;
        phi = high_phi;

        if (low_phi > high_phi)
        {
            double tmp = low_phi;
            low_phi = high_phi;
            high_phi = tmp;
        }

        long int low_phin = floor(low_phi*nbin);
        long int high_phin = floor(high_phi*nbin);
        long int nphi = high_phin-low_phin+1;
       
        assert(nphi<=nbin);

        if (nphi == 1)
        {
            double vWli0 = 1.;

            vWli0 = (high_phi-low_phi)*nbin;

            long int l=low_phin%nbin;
            l = l<0 ? l+nbin:l;
            {
                mxWTW[l*nbin+l] += vWli0*vWli0;
            }

            for (long int j=0; j<npol*nchan; j++)
            {       
                vWTd_T[l*npol*nchan+j] += vWli0*databuffer.buffer[i*npol*nchan+j];
            }
        }
        else if (nphi == 2)
        {
            double vWli0=1., vWli1=1.;

            vWli0 = 1.-(low_phi*nbin-floor(low_phi*nbin));
            vWli1 = high_phi*nbin-floor(high_phi*nbin);

            long int l=low_phin%nbin;
            l = l<0 ? l+nbin:l;
            long int m = high_phin%nbin;
            m = m<0 ? m+nbin:m;

            mxWTW[l*nbin+l] += vWli0*vWli0;
            mxWTW[l*nbin+m] += vWli0*vWli1;
            mxWTW[m*nbin+l] += vWli1*vWli0;
            mxWTW[m*nbin+m] += vWli1*vWli1;

            for (long int j=0; j<npol*nchan; j++)
            {
                vWTd_T[l*npol*nchan+j] += vWli0*databuffer.buffer[i*npol*nchan+j];
                vWTd_T[m*npol*nchan+j] += vWli1*databuffer.buffer[i*npol*nchan+j];
            }
        }
        else
        {
            vector<double> vWli(nphi, 1.);
            vector<int> binplan(nphi, 0);

            vWli[0] = 1.-(low_phi*nbin-floor(low_phi*nbin));
            vWli[nphi-1] = high_phi*nbin-floor(high_phi*nbin);

            for (long int l=0; l<nphi; l++)
            {
                binplan[l] = (low_phin+l)%nbin;
                binplan[l] = binplan[l]<0 ? binplan[l]+nbin:binplan[l];
            }

            for (long int l=0; l<nphi; l++)
            {
                for (long int m=0; m<nphi; m++)
                {
                    mxWTW[binplan[l]*nbin+binplan[m]] += vWli[l]*vWli[m];
                }
            }

            for (long int l=0; l<nphi; l++)
            {
                for (long int j=0; j<npol*nchan; j++)
                {
                    vWTd_T[binplan[l]*npol*nchan+j] += vWli[l]*databuffer.buffer[i*npol*nchan+j];
                }
            }
        }

        i++;
    }

    for (long int l=0; l<nbin; l++)
    {
        for (long int m=0; m<nbin; m++)
        {
            mxWTW[l*nbin+m] /= (databuffer.nsamples*databuffer.tsamp/period);
        }
        mxWTW[l*nbin+l] += lambda;
    }

    vector<double> vWTd(npol*nchan*nbin, 0.);
    transpose_pad<double>(&vWTd[0], &vWTd_T[0], nbin, npol*nchan);

    Eigen::Map<Eigen::MatrixXd> eigen_mxWTW(mxWTW.data(), nbin, nbin);
    Eigen::Map<Eigen::MatrixXd> eigen_vWTd(vWTd.data(), nbin, npol*nchan);
    eigen_vWTd = eigen_mxWTW.colPivHouseholderQr().solve(eigen_vWTd);

    for (long int ipol=0; ipol<npol; ipol++)
    {
#ifdef _OPENMP
#pragma omp parallel for num_threads(num_threads)
#endif
        for (long int ichan=0; ichan<nchan; ichan++)
        {
            for (long int ibin=0; ibin<nbin; ibin++)
            {
                profiles[ipol*nchan*nbin+ichan*nbin+ibin] = vWTd[ipol*nchan*nbin+ichan*nbin+ibin];
            }
        }
    }

    nsubint++;
}
