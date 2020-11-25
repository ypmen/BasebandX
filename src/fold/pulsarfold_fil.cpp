/**
 * @author Yunpeng Men
 * @email ypmen@mpifr-bonn.mpg.de
 * @create date 2020-11-25 10:52:55
 * @modify date 2020-11-25 10:52:55
 * @desc [description]
 */

#define FAST 1

#include <iostream>
#include <string.h>
#include <vector>
#include <boost/program_options.hpp>

#include "dedisperse.h"
#include "pulsarfolder.h"
#include "archwriter.h"
#include "filterbank.h"

using namespace std;
using namespace boost::program_options;

unsigned int num_threads;

int main(int argc, const char *argv[])
{
	int verbose = 0;

    options_description desc{"Options"};
    desc.add_options()
            ("help,h", "Help")
            ("verbose,v", "Print debug information")
            ("threads,t", value<unsigned int>()->default_value(1), "Number of threads")
            ("template", value<string>()->default_value(""), "PSRFITS template file")
            ("t2predictor,P", value<vector<string>>()->multitoken()->composing(), "T2predictor file")
            ("parfile", value<vector<string>>()->multitoken()->composing(), "Tempo2 parfile")
            ("fold_period,c", value<vector<double>>()->multitoken()->composing(), "Folding period")
            ("nbin,b", value<int>()->default_value(512), "Number of bins per period")
            ("tsubint,L", value<double>()->default_value(1), "Time length of integration")
            ("nsub,n", value<int>()->default_value(0), "Output number of integration per file")
            ("dm,d", value<vector<double>>()->multitoken()->composing(), "Dispersion measure")
            ("telescope,k", value<string>(), "Telescope name")
            ("cont", "Input files are contiguous")
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
	long int nfil = fnames.size();
	Filterbank *fil = new Filterbank [nfil];
	for (long int i=0; i<nfil; i++)
	{
		fil[i].filename = fnames[i];
	}

	vector<MJD> tstarts;
	vector<MJD> tends;
	long int ntotal = 0;
	for (long int i=0; i<nfil; i++)
	{
        fil[i].read_header();
        ntotal += fil[i].nsamples;
        MJD tstart(fil[i].tstart);
        tstarts.push_back(tstart);
        tends.push_back(tstart+fil[i].nsamples*fil[i].tsamp);
	}
	vector<size_t> idx = argsort(tstarts);
	for (long int i=0; i<nfil-1; i++)
	{
		if (abs((tends[idx[i]]-tstarts[idx[i+1]]).to_second())>0.5*fil[idx[i]].tsamp)
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

    std::string src_name;
    std::string s_telescope;
    std::string s_ra;
    std::string s_dec;
    int ibeam = 1;
    /**
     * @brief read ibeam from data file
     * 
     */
    if (fil[0].ibeam != 0)
        ibeam = fil[0].ibeam;
    /**
     * @brief read src_name from data file
     * 
     */
    if (strcmp(fil[0].source_name, "") != 0)
        src_name = fil[0].source_name;

    if (vm.count("telescope"))
    {
        s_telescope = vm["telescope"].as<string>();
    }
    else
    {
        get_telescope_name(fil[0].telescope_id, s_telescope);
    }

    /**
     * @brief read ra and dec from data file
     * 
     */
    double src_raj=0., src_dej=0.;
    if (fil[0].src_raj != 0.)
    {
        src_raj = fil[0].src_raj;
    }
    if (fil[0].src_dej != 0.)
    {
        src_dej = fil[0].src_dej;
    }

    get_s_radec(src_raj, src_dej, s_ra, s_dec);

	long int nchans = fil[0].nchans;
    double tsamp = fil[0].tsamp;
    int nifs = fil[0].nifs;
    int nsblk = 1024;

    MJD start_mjd = tstarts[idx[0]];

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

            folders[k].nbin = vm["nbin"].as<int>();
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
        double fmin = fil[0].frequency_table[0];
        double fmax = fil[0].frequency_table[0];
        for (long int j=0; j<nchans; j++)
        {
            fmin = fmin<fil[0].frequency_table[j] ? fmin:fil[0].frequency_table[j];
            fmax = fmax>fil[0].frequency_table[j] ? fmax:fil[0].frequency_table[j];
        }
        double chanwidth = abs(fil[0].frequency_table[1]-fil[0].frequency_table[0]);
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

            folders[k].nbin = vm["nbin"].as<int>();
            writers[k].template_file = vm["template"].as<std::string>();
        }
    }
    else
    {
        std::vector<double> fold_periods = vm["fold_period"].as<std::vector<double>>();

        for (long int k=0; k<npulsar; k++)
        {
            folders[k].fold_period = fold_periods[k];

            folders[k].nbin = vm["nbin"].as<int>();
            writers[k].template_file = vm["template"].as<std::string>();
        }
    }

    long int ndump = (int)(vm["tsubint"].as<double>()/tsamp);

	DataBuffer<float> databuf(ndump, nifs*nchans);
	databuf.tsamp = tsamp;
	memcpy(&databuf.frequencies[0], fil[0].frequency_table, sizeof(double)*nchans);

    MJD mjd_tmp = start_mjd;
    char arname[4096];
	const char format[] = "%Y-%m-%d-%H:%M:%S";
	mjd_tmp.format();
	mjd_tmp.to_date(arname, 4096, format);

    for (long int k=0; k<npulsar; k++)
    {
        folders[k].nbin = vm["nbin"].as<int>();
        folders[k].npol = nifs;
        folders[k].nsblk = nsblk;

        folders[k].prepare(databuf);

        if (vm.count("predictor") or vm.count("parfile"))
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
	for (long int idxn=0; idxn<nfil; idxn++)
	{
		long int n = idx[idxn];
        long int nseg = ceil(1.*fil[0].nsamples/nsblk);
        long int ns_filn = 0;

		for (long int s=0; s<nseg; s++)
		{
			if (verbose)
			{
				cerr<<"\r\rfinish "<<setprecision(2)<<fixed<<tsamp*count<<" seconds ";
				cerr<<"("<<100.*count/ntotal<<"%)";
			}

			fil[n].read_data(nsblk);
#ifdef FAST
			unsigned char *pcur = (unsigned char *)(fil[n].data);
#endif
			for (long int i=0; i<nsblk; i++)
			{
				count++;

				for (long int k=0; k<nifs; k++)
				{
					for (long int j=0; j<nchans; j++)
					{
						databuf.buffer[bcnt1*nifs*nchans+k*nchans+j] =  pcur[k*nchans+j];
					}
				}

                bcnt1++;
				ntot++;
                ns_filn++;

				if (ntot%ndump == 0)
				{
                    for (long int ipsr=0; ipsr<npulsar; ipsr++)
                    {
                        if (!contiguous)
                            folders[ipsr].start_epoch = (long double)fil[n].tstart+(long double)(ns_filn-ndump)*fil[n].tsamp/86400.;
                        else
                            folders[ipsr].start_epoch = (long double)fil[idx[0]].tstart+(long double)(ntot-ndump)*fil[n].tsamp/86400.;
                        
                        folders[ipsr].run(databuf);

                        writers[ipsr].run(folders[ipsr]);
                    }

                    bcnt1 = 0;

                    if (++intcnt == nsub)
                    {
                        mjd_tmp = (long double)fil[n].tstart + (long double)ns_filn*tsamp/86400.;
                        mjd_tmp.format();
                        mjd_tmp.to_date(arname, 4096, format);
                        for (long int ipsr=0; ipsr<npulsar; ipsr++)
                        {
                            writers[ipsr].close();
                            
                            if (vm.count("predictor") or vm.count("parfile"))
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

                if (ns_filn == fil[n].nsamples)
                {
                    goto next;
                }

				pcur += nifs*nchans;
			}
		}
        next:
		fil[n].close();
	}

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
