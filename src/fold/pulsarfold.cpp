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
#include "psrfitsreader.h"
#include "filterbankreader.h"

#include "logging.h"

using namespace std;
using namespace boost::program_options;

unsigned int num_threads;

int main(int argc, const char *argv[])
{
	init_logging();

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
			("ibeam,i", value<int>(), "Beam id")
            ("cont", "Input files are contiguous")
			("wts", "Apply DAT_WTS")
			("scloffs", "Apply DAT_SCL and DAT_OFFS")
			("zero_off", "Apply ZERO_OFF")
            ("trlsm", "Using Tikhonov-regularized least square method folding algorithm")
            ("lambda", value<double>()->default_value(1), "The regularization factor")
			("filterbank", "Input filterbank format data")
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
        BOOST_LOG_TRIVIAL(error)<<"no input file";
        return -1;
    }
    if (vm.count("template") == 0)
    {
        BOOST_LOG_TRIVIAL(error)<<"no template file";
        return -1;
    }
    if (vm.count("parfile")==0 and vm.count("t2predictor") == 0 and vm.count("fold_period") == 0)
    {
        BOOST_LOG_TRIVIAL(error)<<"no parfile or t2pred or folding period";
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

	bool apply_wts = false;
	bool apply_scloffs = false;
	bool apply_zero_off = false;

	if (vm.count("wts"))
		apply_wts = true;
	if (vm.count("scloffs"))
		apply_scloffs = true;
	if (vm.count("zero_off"))
		apply_zero_off = true;

	PSRDataReader * reader;

	if (vm.count("filterbank"))
		reader= new FilterbankReader;
	else
		reader= new PsrfitsReader;

	if (vm.count("ibeam"))
		reader->beam = std::to_string(vm["ibeam"].as<int>());
	if (vm.count("telescope"))
		reader->telescope = vm["telescope"].as<std::string>();
	reader->fnames = fnames;
	reader->sumif = false;
	reader->contiguous = contiguous;
	reader->verbose = verbose;
	reader->apply_scloffs = apply_scloffs;
	reader->apply_wts = apply_wts;
	reader->apply_zero_off = apply_zero_off;
	reader->check();
	reader->read_header();

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

        double mjd_start = (reader->mjd_starts[reader->idmap.front()].to_day())-0.01;
        double mjd_end = (reader->mjd_ends[reader->idmap.back()].to_day())+1;
        double fmin = 0.;
        double fmax = 0.;
		reader->get_fmin_fmax(fmin, fmax);
       
        for (long int k=0; k<npulsar; k++)
        {
            string cmd = "tempo2 -npsr 1 -f " + parfiles[k] + " -pred ";
            cmd += "'";
            if (reader->telescope.empty())
                cmd += "fake";
            else
                cmd += reader->telescope;
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

    long int ndump = ceil(1./reader->tsamp);
    ndump = std::min(ndump, (long int)reader->nsamples);

    std::vector<int> blocksize(npulsar, ceil(vm["tsubint"].as<double>()/reader->tsamp));

	DataBuffer<float> databuf(ndump, reader->nifs*reader->nchans);
	databuf.tsamp = reader->tsamp;
	databuf.frequencies = reader->frequencies;

    MJD mjd_tmp = reader->start_mjd;
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
                double pfold = folders[k].predictor.get_pfold((reader->mjd_starts[reader->idmap.front()]+reader->mjd_ends[reader->idmap.back()]).dividedby2().to_day(), 0.5*(databuf.frequencies.front()+databuf.frequencies.back()));
                folders[k].nbin = std::pow(2, std::floor(std::log2(pfold/reader->tsamp)));
            }
            else
            {
                folders[k].nbin = std::pow(2, std::floor(std::log2(folders[k].fold_period/reader->tsamp)));
            }
        }

        folders[k].npol = reader->nifs;
        folders[k].nsblk = reader->nsblk;
        folders[k].lambda = vm["lambda"].as<double>();
        folders[k].start_epoch = reader->start_mjd;
        folders[k].prepare(databuf);

        if (vm.count("single"))
        {
            if (vm.count("t2predictor") or vm.count("parfile"))
            {
                MJD block_epoch_start = folders[k].start_epoch + 0.5*reader->tsamp;
                phase[k] = folders[k].predictor.get_phase(block_epoch_start.to_day(), folders[k].freqref);
                phase[k] += 1.;
                MJD block_epoch_end(folders[k].predictor.get_phase_inverse(phase[k], folders[k].freqref));
                
                blocksize[k] = (block_epoch_end-block_epoch_start).to_second()/reader->tsamp;
            }
            else
            {
                blocksize[k] = folders[k].fold_period/reader->tsamp;
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
        writers[k].ibeam = std::stoi(reader->beam);
        writers[k].src_name = reader->source_name;
        writers[k].telescop = reader->telescope;
        writers[k].ra = reader->ra;
        writers[k].dec = reader->dec;
        writers[k].dm = dms[k];
        
        writers[k].prepare(folders[k]);
    }

	long int intcnt = 0;

	while (!reader->is_end)
	{
		if (reader->read_data(databuf, ndump) != ndump) break;

		size_t ns_psfn = reader->get_count_curfile();
		size_t ntot = reader->get_count();
		MJD start_mjd_curfile = reader->get_start_mjd_curfile();
		double tsamp_curfile = reader->get_tsamp_curfile();

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
                                    folders[ipsr].start_epoch = start_mjd_curfile + (ns_psfn - nleft + 1) * tsamp_curfile;
                                else
                                    folders[ipsr].start_epoch = reader->start_mjd + (ntot - nleft + 1) * tsamp_curfile;
					if (vm.count("single"))
					{
						if (vm.count("t2predictor") or vm.count("parfile"))
						{
							MJD block_epoch_start = folders[ipsr].start_epoch + 0.5*reader->tsamp;
							phase[ipsr] += 1.;
							MJD block_epoch_end(folders[ipsr].predictor.get_phase_inverse(phase[ipsr], folders[ipsr].freqref));
							
							blocksize[ipsr] = (block_epoch_end-block_epoch_start).to_second()/reader->tsamp;
						}
						else
						{
							blocksize[ipsr] = folders[ipsr].fold_period/reader->tsamp;
						}
					}

					if (++intcnt == nsub)
					{
						mjd_tmp = start_mjd_curfile + (ns_psfn - nleft + 1) * tsamp_curfile;
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
	}

    for (long int ipsr=0; ipsr<npulsar; ipsr++)
    {
        writers[ipsr].close();
        if (vm.count("parfile"))
            system(("rm "+to_string(ipsr)+".t2pred").c_str());
    }

	if (verbose)
	{
		cerr<<"\r\rfinish "<<setprecision(2)<<fixed<<reader->tsamp*reader->get_count()<<" seconds ";
		cerr<<"("<<100.*reader->get_count()/reader->nsamples<<"%)"<<endl;
	}

	return 0;
}
