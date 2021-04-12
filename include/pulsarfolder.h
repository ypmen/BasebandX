/**
 * @author Yunpeng Men
 * @email ypmen@mpifr-bonn.mpg.de
 * @create date 2020-11-21 17:12:36
 * @modify date 2020-11-21 17:12:36
 * @desc [description]
 */

#ifndef PULSARFOLDER_H
#define PULSARFOLDER_H

#include <vector>
#include "databuffer.h"
#include "predictor.h"
#include "mjd.h"

namespace Baseband
{
    /**
     * @brief folding input data chunk
     * 
     */
    class PulsarFolder
    {
    public:
        PulsarFolder();
        PulsarFolder(const std::string &fname);
        ~PulsarFolder();
        void prepare(DataBuffer<float> &databuffer);
#ifdef __AVX2__
        void run(std::vector<float, boost::alignment::aligned_allocator<float, 32>>::iterator data);
        void runLSM(std::vector<float, boost::alignment::aligned_allocator<float, 32>>::iterator data);
#else
        void run(std::vector<float>::iterator data);
        void runLSM(std::vector<float>::iterator data);
#endif
        void flush();
        void flushLSM();
    public:
        int npol;
        int nbin;
        int nsblk;
        MJD start_epoch;
        Predictors predictor;
        long double fold_period;
        double lambda;
    public:
        int nchan;
        int nsubint;
        double tsubint;
        double offs_sub;
        double period;
        double freqref;
        MJD epoch;
        std::vector<double> frequencies;
        std::vector<float> profiles;
        int cnt;
    private:
        double tsamp;
        double phi;
        double pfold;
        vector<double> mxWTW;
        vector<double> vWTd_T;
        vector<int> hits;
        vector<float> profilesTPF;
    };
}
#endif /* PULSARFOLDER_H */
