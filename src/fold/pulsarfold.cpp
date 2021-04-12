/**
 * @author Yunpeng Men
 * @email ypmen@mpifr-bonn.mpg.de
 * @create date 2020-11-22 12:53:59
 * @modify date 2020-11-22 12:53:59
 * @desc [description]
 */

#include <iostream>
#include <string.h>
#include <vector>
#include <boost/program_options.hpp>

#include "dedisperse.h"
#include "pulsarfolder.h"
#include "archwriter.h"

using namespace std;
using namespace boost::program_options;

unsigned int num_threads;

int main(int argc, const char *argv[])
{
	int verbose = 0;
    bool lsm = false;

    options_description desc{"Options"};
    desc.add_options()
            ("help,h", "Help")
            ("verbose,v", "Print debug information")
            ("threads,t", value<unsigned int>()->default_value(1), "Number of threads")
            ("template", value<string>()->default_value(""), "PSRFITS template file")
            ("t2predictor,P", value<vector<string>>()->multitoken()->composing(), "T2predictor file")
            ("parfile", value<vector<string>>()->multitoken()->composing(), "Tempo2 parfile")
            ("fold_period,c", value<vector<long double>>()->multitoken()->composing(), "Folding period")
            ("nbin,b", value<int>(), "Number of bins per period")
            ("tsubint,L", value<double>()->default_value(1), "Time length of integration")
            ("single,s", "Single pulse folding")
            ("nsub,n", value<int>()->default_value(0), "Output number of integration per file")
            ("dm,d", value<vector<double>>()->multitoken()->composing(), "Dispersion measure")
            ("telescope,k", value<string>(), "Telescope name")
            ("cont", "Input files are contiguous")
            ("trlsm", "Using Tikhonov-regularized least square method folding algorithm")
            ("lambda", value<double>()->default_value(1), "The regularization factor")
            ("input,f", value<vector<string>>()->multitoken()->composing(), "Input files");

    positional_options_description pos_desc;
    pos_desc.add("input", -1);
    command_line_parser parser{argc, argv};
    parser.options(desc).allow_unregistered().style(command_line_style::default_style | command_line_style::allow_short);
    parser.options(desc).positional(pos_desc).allow_unregistered();
    parsed_options parsed_options = parser.run();

    variables_map vm;
    store(parsed_options, vm);
    notify(vm);

    if (vm.count("help"))
    {
        std::cout << desc << '\n';
        return 0;
    }
    if (vm.count("verbose"))
    {
        verbose = 1;
    }
    if (vm.count("trlsm"))
    {
        lsm = true;
    }
    if (vm.count("input") == 0)
    {
        cerr<<"Error: no input file"<<endl;
        return -1;
    }
    if (vm.count("template") == 0)
    {
        cerr<<"Error: no template file"<<endl;
        return -1;
    }
    if (vm.count("parfile")==0 and vm.count("t2predictor") == 0 and vm.count("fold_period") == 0)
    {
        cerr<<"Error: no parfile or t2pred or folding period"<<endl;
        return -1;
    }

    num_threads = vm["threads"].as<unsigned int>();
    int nsub = vm["nsub"].as<int>();
    bool contiguous = vm.count("cont");
    vector<string> fnames = vm["input"].as<vector<string>>();

    int npulsar = vm.count("t2predictor")>vm.count("fold_period") ? vm.count("t2predictor"):vm.count("fold_period");
    npulsar = npulsar>vm.count("parfile") ? npulsar:vm.count("parfile");
    vector<double> dms(npulsar, 0.);
    if (vm.count("dm") == 1)
    {
        fill(dms.begin(), dms.end(), vm["dm"].as<vector<double>>()[0]);
    }
    else if (vm.count("dm") > 1)
    {
        dms = vm["dm"].as<vector<double>>();
    }

    /**
     * @brief read data info
     * 
     */
    long int npsf = fnames.size();
	Psrfits *psf = new Psrfits [npsf];
	for (long int i=0; i<npsf; i++)
	{
		psf[i].filename = fnames[i];
	}

	vector<MJD> tstarts;
	vector<MJD> tends;
	long int ntotal = 0;
	for (long int i=0; i<npsf; i++)
	{
		psf[i].open();
		psf[i].primary.load(psf[i].fptr);
		psf[i].load_mode();
		psf[i].subint.load_header(psf[i].fptr);
		ntotal += psf[i].subint.nsamples;
		tstarts.push_back(psf[i].primary.start_mjd);
		tends.push_back(psf[i].primary.start_mjd+psf[i].subint.nsamples*psf[i].subint.tbin);
		psf[i].close();
	}
	vector<size_t> idx = argsort(tstarts);
	for (long int i=0; i<npsf-1; i++)
	{
		if (abs((tends[idx[i]]-tstarts[idx[i+1]]).to_second())>0.5*psf[idx[i]].subint.tbin)
		{
			if (contiguous)
			{
				cerr<<"Warning: time not contiguous"<<endl;
			}
			else
			{
				cerr<<"Error: time not contiguous"<<endl;
				exit(-1);
			}
		}
	}

    psf[0].open();
	psf[0].primary.load(psf[0].fptr);
	psf[0].load_mode();
	psf[0].subint.load_header(psf[0].fptr);

	if (psf[0].mode != Integration::SEARCH)
	{
		cerr<<"Error: mode is not SEARCH"<<endl;
		exit(-1);
	}

    std::string src_name;
    std::string s_telescope;
    std::string s_ra;
    std::string s_dec;
    int ibeam = 1;
    /**
     * @brief read ibeam from data file
     * 
     */
    if (strcmp(psf[0].primary.ibeam, "") != 0)
			ibeam = stoi(psf[0].primary.ibeam);
    /**
     * @brief read src_name from data file
     * 
     */
	if (strcmp(psf[0].primary.src_name, "") != 0)
		src_name = psf[0].primary.src_name;

    if (vm.count("telescope"))
    {
        s_telescope = vm["telescope"].as<string>();
    }
    else
    {
        if (strcmp(psf[0].primary.telesop, "") != 0)
            s_telescope = psf[0].primary.telesop;
    }

    /**
     * @brief read ra and dec from data file
     * 
     */
	if (strcmp(psf[0].primary.ra, "") != 0)
    {
        s_ra = psf[0].primary.ra;
    }
    if (strcmp(psf[0].primary.dec, "") != 0)
    {
        s_dec = psf[0].primary.dec;
    }

	Integration it;
	psf[0].subint.load_integration(psf[0].fptr, 0, it);

	long int nchans = it.nchan;
    double tsamp = psf[0].subint.tbin;
    int nifs = it.npol;
    int nsblk = it.nsblk;

    MJD start_mjd = psf[idx[0]].primary.start_mjd;

    psf[0].close();

    /**
     * @brief initialize pipeline
     * 
     * 
     */
    vector<Baseband::PulsarFolder> folders(npulsar);
    vector<Baseband::ArchiveWriter> writers(npulsar);

    if (vm.count("t2predictor"))
    {
        std::vector<std::string> predfiles = vm["t2predictor"].as<std::vector<std::string>>();
    
        for (long int k=0; k<npulsar; k++)
        {
            folders[k].predictor.read_t2pred(predfiles[k]);
            folders[k].fold_period = 0.;
            writers[k].template_file = vm["template"].as<std::string>();
        }
    }
    else if (vm.count("parfile"))
    {
        std::vector<std::string> parfiles = vm["parfile"].as<std::vector<std::string>>();
        
        /**
         * @brief read DM from parfile
         * 
         */
        for (long int k=0; k<npulsar; k++)
        {
            string filename = parfiles[k];
            string line;
            ifstream ddplan(filename);
            int id = 0;
            while (getline(ddplan, line))
            {
                vector<string> items;
                boost::split(items, line, boost::is_any_of("\t "), boost::token_compress_on);
                if (items[0] == "DM")
                {
                    dms[k] = stod(items[1]);
                    break;
                }
            }
        }

        double mjd_start = (tstarts[idx[0]].to_day())-0.01;
        double mjd_end = (tends[idx.back()].to_day())+1;
        double fmin = it.frequencies[0];
        double fmax = it.frequencies[0];
        for (long int j=0; j<nchans; j++)
        {
            fmin = fmin<it.frequencies[j] ? fmin:it.frequencies[j];
            fmax = fmax>it.frequencies[j] ? fmax:it.frequencies[j];
        }
        double chanwidth = abs(it.frequencies[1]-it.frequencies[0]);
        fmin -= 0.5*chanwidth;
        fmax += 0.5*chanwidth;
        for (long int k=0; k<npulsar; k++)
        {
            string cmd = "tempo2 -npsr 1 -f " + parfiles[k] + " -pred ";
            cmd += "'";
            if (s_telescope.empty())
                cmd += "fake";
            else
                cmd += s_telescope;
            cmd += " ";
            cmd += to_string(mjd_start);
            cmd += " ";
            cmd += to_string(mjd_end);
            cmd += " ";
            cmd += to_string(fmin);
            cmd += " ";
            cmd += to_string(fmax);
            cmd += " ";
            cmd += "12 2";
            cmd += " ";
            cmd += "3600";
            cmd += "'";
            cmd += " > stdout.txt 2> stderr.txt";
            system(cmd.c_str());

            system("rm pred.tim");
            system(("mv t2pred.dat "+to_string(k)+".t2pred").c_str());

            folders[k].predictor.read_t2pred(to_string(k)+".t2pred");
            folders[k].fold_period = 0.;

            writers[k].template_file = vm["template"].as<std::string>();
        }
    }
    else
    {
        std::vector<long double> fold_periods = vm["fold_period"].as<std::vector<long double>>();

        for (long int k=0; k<npulsar; k++)
        {
            folders[k].fold_period = fold_periods[k];

            writers[k].template_file = vm["template"].as<std::string>();
        }
    }

    long int ndump = ceil(1./tsamp);
    ndump = std::min(ndump, ntotal);

    std::vector<int> blocksize(npulsar, ceil(vm["tsubint"].as<double>()/tsamp));

	DataBuffer<float> databuf(ndump, nifs*nchans);
	databuf.tsamp = tsamp;
	memcpy(&databuf.frequencies[0], it.frequencies, sizeof(double)*nchans);

    MJD mjd_tmp = start_mjd;
    char arname[4096];
	const char format[] = "%Y-%m-%d-%H:%M:%S";
	mjd_tmp.format();
	mjd_tmp.to_date(arname, 4096, format);

    std::vector<long double> phase(npulsar, 0.);

    for (long int k=0; k<npulsar; k++)
    {
        if (vm.count("nbin"))
        {
            folders[k].nbin = vm["nbin"].as<int>();
        }
        else
        {
            if (vm.count("t2predictor") or vm.count("parfile"))
            {
                double pfold = folders[k].predictor.get_pfold((tstarts[idx.front()]+tstarts[idx.back()]).dividedby2().to_day(), 0.5*(databuf.frequencies.front()+databuf.frequencies.back()));
                folders[k].nbin = std::pow(2, std::floor(std::log2(pfold/tsamp)));
            }
            else
            {
                folders[k].nbin = std::pow(2, std::floor(std::log2(folders[k].fold_period/tsamp)));
            }
        }

        folders[k].npol = nifs;
        folders[k].nsblk = nsblk;
        folders[k].lambda = vm["lambda"].as<double>();

        if (!contiguous)
            folders[k].start_epoch = psf[idx[0]].primary.start_mjd;
        else
            folders[k].start_epoch = psf[idx[0]].primary.start_mjd;

        folders[k].prepare(databuf);

        if (vm.count("single"))
        {
            if (vm.count("t2predictor") or vm.count("parfile"))
            {
                MJD block_epoch_start = folders[k].start_epoch + 0.5*tsamp;
                phase[k] = folders[k].predictor.get_phase(block_epoch_start.to_day(), folders[k].freqref);
                phase[k] += 1.;
                MJD block_epoch_end(folders[k].predictor.get_phase_inverse(phase[k], folders[k].freqref));
                
                blocksize[k] = (block_epoch_end-block_epoch_start).to_second()/tsamp;
            }
            else
            {
                blocksize[k] = folders[k].fold_period/tsamp;
            }
        }

        if (vm.count("t2predictor") or vm.count("parfile"))
        {
            long double phase = floor(folders[k].predictor.get_phase(mjd_tmp.to_day(), folders[k].freqref));
            mjd_tmp.format(folders[k].predictor.get_phase_inverse(phase, folders[k].freqref));
        }

        writers[k].start_mjd = mjd_tmp;
        writers[k].rootname = arname;

        writers[k].template_file = vm["template"].as<std::string>();
        writers[k].ibeam = ibeam;
        writers[k].src_name = src_name;
        writers[k].telescop = s_telescope;
        writers[k].ra = s_ra;
        writers[k].dec = s_dec;
        writers[k].dm = dms[k];
        
        writers[k].prepare(folders[k]);
    }

	long int ntot = 0;
	long int count = 0;
    long int bcnt1 = 0;
    long int intcnt = 0;

	for (long int idxn=0; idxn<npsf; idxn++)
	{
		long int n = idx[idxn];
        long int ns_psfn = 0;

		psf[n].open();
		psf[n].primary.load(psf[n].fptr);
		psf[n].load_mode();
		psf[n].subint.load_header(psf[n].fptr);

		for (long int s=0; s<psf[n].subint.nsubint; s++)
		{
			if (verbose)
			{
				cerr<<"\r\rfinish "<<setprecision(2)<<fixed<<tsamp*count<<" seconds ";
				cerr<<"("<<100.*count/ntotal<<"%)";
			}

			psf[n].subint.load_integration_data(psf[n].fptr, s, it);

			for (long int i=0; i<it.nsblk; i++)
			{
				count++;

                if (it.dtype == Integration::UINT8)
                {
                    for (long int k=0; k<nifs; k++)
                    {
                        for (long int j=0; j<nchans; j++)
                        {
                            databuf.buffer[bcnt1*nifs*nchans+k*nchans+j] = ((unsigned char *)(it.data))[i*nifs*nchans+k*nchans+j];
                        }
                    }
                }
                else if (it.dtype == Integration::FLOAT)
                {
                    for (long int k=0; k<nifs; k++)
                    {
                        for (long int j=0; j<nchans; j++)
                        {
                            databuf.buffer[bcnt1*nifs*nchans+k*nchans+j] = ((float *)(it.data))[i*nifs*nchans+k*nchans+j];
                        }
                    }
                }
                else
                {
                    cerr<<"Error: data type is not support"<<endl;
                }

                bcnt1++;
				ntot++;
                ns_psfn++;

				if (ntot%ndump == 0)
				{
                    for (long int ipsr=0; ipsr<npulsar; ipsr++)
                    {
                        int nleft = ndump;
                        while (nleft)
                        {
                            if (!lsm)
                                folders[ipsr].run(databuf.buffer.begin()+(databuf.nsamples-nleft)*databuf.nchans);
                            else
                                folders[ipsr].runLSM(databuf.buffer.begin()+(databuf.nsamples-nleft)*databuf.nchans);
                            
                            if (folders[ipsr].cnt == blocksize[ipsr])
                            {
                                if (!lsm)
                                    folders[ipsr].flush();
                                else
                                    folders[ipsr].flushLSM();

                                writers[ipsr].run(folders[ipsr]);

                                if (!contiguous)
                                    folders[ipsr].start_epoch = psf[n].primary.start_mjd+(ns_psfn-nleft+1)*psf[n].subint.tbin;
                                else
                                    folders[ipsr].start_epoch = psf[idx[0]].primary.start_mjd+(ntot-nleft+1)*psf[n].subint.tbin;

                                if (vm.count("single"))
                                {
                                    if (vm.count("t2predictor") or vm.count("parfile"))
                                    {
                                        MJD block_epoch_start = folders[ipsr].start_epoch + 0.5*tsamp;
                                        phase[ipsr] += 1.;
                                        MJD block_epoch_end(folders[ipsr].predictor.get_phase_inverse(phase[ipsr], folders[ipsr].freqref));
                                        
                                        blocksize[ipsr] = (block_epoch_end-block_epoch_start).to_second()/tsamp;
                                    }
                                    else
                                    {
                                        blocksize[ipsr] = folders[ipsr].fold_period/tsamp;
                                    }
                                }

                                if (++intcnt == nsub)
                                {
                                    mjd_tmp = psf[n].primary.start_mjd + (ns_psfn-nleft+1)*psf[n].subint.tbin;
                                    mjd_tmp.format();
                                    mjd_tmp.to_date(arname, 4096, format);
                                    for (long int ipsr=0; ipsr<npulsar; ipsr++)
                                    {
                                        writers[ipsr].close();
                                        
                                        if (vm.count("t2predictor") or vm.count("parfile"))
                                        {
                                            long double phase = floor(folders[ipsr].predictor.get_phase(mjd_tmp.to_day(), folders[ipsr].freqref));
                                            mjd_tmp.format(folders[ipsr].predictor.get_phase_inverse(phase, folders[ipsr].freqref));
                                        }

                                        writers[ipsr].start_mjd = mjd_tmp;
                                        writers[ipsr].rootname = arname;
                                        
                                        writers[ipsr].prepare(folders[ipsr]);
                                    }

                                    intcnt = 0;
                                }
                            }

                            nleft--;
                        }
                    }

                    bcnt1 = 0;
				}

                if (ns_psfn == psf[n].subint.nsamples)
				{
					goto next;
				}
			}
		}
        next:
		psf[n].close();
	}

    delete [] psf;

    for (long int ipsr=0; ipsr<npulsar; ipsr++)
    {
        writers[ipsr].close();
        if (vm.count("parfile"))
            system(("rm "+to_string(ipsr)+".t2pred").c_str());
    }

	if (verbose)
	{
		cerr<<"\r\rfinish "<<setprecision(2)<<fixed<<tsamp*count<<" seconds ";
		cerr<<"("<<100.*count/ntotal<<"%)"<<endl;
	}

	return 0;
}
