/**
 * @author Yunpeng Men
 * @email ypmen@pku.edu.cn
 * @create date 2020-11-02 10:54:12
 * @modify date 2020-11-02 10:54:12
 * @desc [description]
 */

#ifndef ARCHWRITER_H
#define ARCHWRITER_H

#include <string>
#include <vector>

#include "psrfits.h"
#include "mjd.h"
#include "integration.h"
#include "pulsarfolder.h"

using namespace std;

namespace Baseband
{
    class ArchiveWriter
    {
    public:
        ArchiveWriter();
        ~ArchiveWriter();
        void close();
        void prepare(PulsarFolder &pulsarfolder);
        void run(PulsarFolder &pulsarfolder);
    public:
        string rootname;
        string template_file;
        Integration::Mode mode;
        int ibeam;
        string src_name;
        string telescop;
        string ra;
        string dec;
        MJD start_mjd;
        double dm;
    public:
        int npol;
        int nchan;
        int nbin;
        double tbin;
        vector<double> frequencies;
    public:
        int nsubint;
    private:
        Integration it;
        Psrfits fits;
    };
}

#endif /* ARCHWRITER_H */
