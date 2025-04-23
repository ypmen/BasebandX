/**
 * @author Yunpeng Men
 * @email ypmen@pku.edu.cn
 * @create date 2020-11-02 16:24:33
 * @modify date 2020-11-02 16:24:33
 * @desc [description]
 */

#include "assert.h"

#include "archwriter.h"

using namespace Baseband;

ArchiveWriter::ArchiveWriter()
{
    mode = Integration::FOLD;
    ibeam = 1;
    npol = 0;
    nchan = 0;
    nbin = 0;
    nsubint = 0;
    tbin = 0.;
    dm = 0.;
}

ArchiveWriter::~ArchiveWriter()
{
    if (fits.fptr != NULL)
    {
        fits.close();
    }
}

void ArchiveWriter::close()
{
    if (fits.fptr != NULL)
    {
        fits.close();
    }    
}

void ArchiveWriter::prepare(PulsarFolder &pulsarfolder)
{
    assert(mode == Integration::FOLD);

    npol = pulsarfolder.npol;
    nchan = pulsarfolder.nchan;
    nbin = pulsarfolder.nbin;
    frequencies = pulsarfolder.frequencies;

    fits.primary.start_mjd = start_mjd;
    strcpy(fits.primary.src_name, src_name.c_str());
    strcpy(fits.primary.telesop, telescop.c_str());
    strcpy(fits.primary.ra, ra.c_str());
    strcpy(fits.primary.dec, dec.c_str());

    fits.filename = "!" + rootname + ".ar";
    strcpy(fits.primary.obs_mode, "PSR");
    strcpy(fits.primary.ibeam, to_string(ibeam).c_str());

    fits.primary.chan_dm = dm;
	fits.primary.obsfreq = 0.5 * (frequencies.back() + frequencies.front());
	fits.primary.obsbw = (frequencies.back() - frequencies.front()) / (frequencies.size() - 1) * frequencies.size();
	fits.primary.obsnchan = frequencies.size();

    fits.subint.mode = Integration::FOLD;
    fits.subint.dtype = Integration::SHORT;
    fits.subint.npol = npol;
    fits.subint.nchan = nchan;
    fits.subint.nbin = nbin;
    fits.subint.tbin = tbin;
    fits.subint.dm = dm;

    fits.parse_template(template_file);
    fits.primary.unload(fits.fptr);

    fits.t2predict.load(pulsarfolder.predictor.filename);
	fits.t2predict.unload(fits.fptr);

    fits.subint.unload_header(fits.fptr);

    it.mode = Integration::FOLD;
    it.dtype = Integration::SHORT;
}

void ArchiveWriter::run(PulsarFolder &pulsarfolder)
{
    it.load_data(&pulsarfolder.profiles[0], npol, nchan, nbin);
    it.load_frequencies(&frequencies[0], nchan);
    it.folding_period = pulsarfolder.fold_period;
    it.tsubint = pulsarfolder.tsubint;
    it.offs_sub = (pulsarfolder.epoch-start_mjd).to_second();
    fits.subint.unload_integration(fits.fptr, it);

    nsubint++;
}
