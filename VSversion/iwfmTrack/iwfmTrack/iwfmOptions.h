#pragma once

#include <boost/program_options.hpp>

namespace po = boost::program_options;

template <typename T>
bool get_option(std::string optionName, po::variables_map &vm, T &out) {
	if(vm.count(optionName)) {
		out = vm[optionName].as<T>();
		return true;
	}
	else {
		std::cout << optionName << " option is missing" << std::endl;
		return false;
	}
}

struct iwfm_track_options {
	std::string config_file;
	std::string XYfile;
	std::string MSHfile;
	std::string STRATfile;
	std::string NRMLfile;
	std::string FACEZfile;
	std::string HFLOWfile;
	std::string VFLOWfile;
	std::string BCfile;
	std::string BCZfile;
	double radius;
	int Nminp;
	int Nsearch;
	int dir;
	std::string WELLSfile;
	std::string PARTICLEfile;
	std::string OUTfile;
	int Npart;
	bool bEndpnts;
	double threshold;
	double wellRadius;
	int maxSteps;
	double stepSize;
	double maxStepSize;
	double minStepSize;
	double tolStepSize;
	int pntPrecision;
	int velPrecision;
	int nThreads;
};


bool readInputParameters(int argc, char *argv[], iwfm_track_options& iwtropt) {
	// User command line options
	po::options_description commandLineOptions("User options");
	commandLineOptions.add_options()
		("version,v", "print version info")
		("help,h", "Get a list of options in the configuration file")
		("config,c", po::value<std::string >(), "Set configuration file")
		;

	po::variables_map vm_cmd;
	po::store(po::parse_command_line(argc, argv, commandLineOptions), vm_cmd);

	if (vm_cmd.size() == 0) {
		std::cout << " To run IWFMTrack specify the configuration file as" << std::endl;
		std::cout << "-c configfilename" << std::endl << std::endl;;
		std::cout << "Other command line options are:" << std::endl;
		std::cout << commandLineOptions << std::endl;
		return false;
	}

	if (vm_cmd.count("version")) {
		std::cout << "|------------------|" << std::endl;
		std::cout << "|    IWFM Trace    |" << std::endl;
		std::cout << "| Version : 0.1.00 |" << std::endl;
		std::cout << "|    by  giorgk    |" << std::endl;
		std::cout << "|------------------|" << std::endl;
		return false;
	}

	// User configuration options
	po::options_description config_opt("Config options");
	config_opt.add_options()
		("XY", "XY coordinate file")
		("MSH", "2D Mesh element file")
		("STRAT", "Stratificationfile")
		("NORMALS", "Horizontal velocity normals")
		("FACEZ", "Element Face elevations")
		("HFLOW", "Horizontal flow values")
		("VFLOW", "Vertical flow values")
		("BC", "Element barycenters")
		("BCZ", "Element barycenters elevation")
		("OUTFILE", "Output file")
		("WELLS", po::value<std::string>()->default_value(""), "Well information")
		("PARTICLES", po::value<std::string>()->default_value(""), "Initial location of particles")
		("RADIUS", po::value<double>()->default_value(3000.0), "Search radius during interpolation")
		("NSEARCH", po::value<int>()->default_value(20), "Number of search closest points ")
		("NMINPNTS", po::value<int>()->default_value(4), "Minimum number of points")
		("DIR", po::value<int>()->default_value(-1), "Direction of particle tracking")
		("NPART", po::value<int>()->default_value(100), "Number of particles per wells")
		("WELLRAD", po::value<double>()->default_value(10), "Distance from well to release particles")
		("ENDPNTS", po::value<bool>()->default_value(false), "Print the end points only")
		("THRES", po::value<double>()->default_value(0.01), "Threshold for inverse distance weight")
		("MAXSTEPS", po::value<int>()->default_value(1000), "Number of particles per wells")
		("STEPSIZE", po::value<double>()->default_value(5), "Initial Step size during particle tracking")
		("MAXSTEPSIZE", po::value<double>()->default_value(300), "Maximum step size during particle tracking")
		("MINSTEPSIZE", po::value<double>()->default_value(1), "Minimum step size during particle tracking")
		("TOLSTEPSIZE", po::value<double>()->default_value(1), "Tolerance of particle tracking")
		("PNTPRECISION", po::value<int>()->default_value(3), "Number of decimals in printing for points and age")
		("VELPRECISION", po::value<int>()->default_value(5), "Number of decimals in printing for velocity")
		("NTHREADS", po::value<int>()->default_value(0), "Number of threads")
		;
	

	if (vm_cmd.count("help")) {
		std::cout << "Copy paste the following list of options" << std::endl;
		std::cout << "into a configuration file without the dashes --" << std::endl;
		std::cout << "The options without default values are mandatory" << std::endl;
		std::cout << "All options are case sensitive" << std::endl;
		std::cout << "------------------------------" << std::endl;
		std::cout << config_opt << std::endl;
		return false;
	}

	po::variables_map vm_cfg;

	bool tf = true;
	if (vm_cmd.count("config")) {
		std::cout << "Reading options from " << vm_cmd["config"].as<std::string>().c_str() << std::endl;
		po::store(po::parse_config_file<char>(vm_cmd["config"].as<std::string>().c_str(), config_opt), vm_cfg);

		// Mandatory options
		tf = get_option<std::string>("XY", vm_cfg, iwtropt.XYfile);
		if (!tf) return false;

		tf = get_option<std::string>("MSH", vm_cfg, iwtropt.MSHfile);
		if (!tf) return false;

		tf = get_option<std::string>("STRAT", vm_cfg, iwtropt.STRATfile);
		if (!tf) return false;

		tf = get_option<std::string>("NORMALS", vm_cfg, iwtropt.NRMLfile);
		if (!tf) return false;

		tf = get_option<std::string>("FACEZ", vm_cfg, iwtropt.FACEZfile);
		if (!tf) return false;

		tf = get_option<std::string>("HFLOW", vm_cfg, iwtropt.HFLOWfile);
		if (!tf) return false;

		tf = get_option<std::string>("VFLOW", vm_cfg, iwtropt.VFLOWfile);
		if (!tf) return false;

		tf = get_option<std::string>("BC", vm_cfg, iwtropt.BCfile);
		if (!tf) return false;

		tf = get_option<std::string>("BCZ", vm_cfg, iwtropt.BCZfile);
		if (!tf) return false;
		
		tf = get_option<std::string>("OUTFILE", vm_cfg, iwtropt.OUTfile);
		if (!tf) return false;
		
		// Options with default values
		iwtropt.WELLSfile = vm_cfg["WELLS"].as<std::string>();
		iwtropt.PARTICLEfile = vm_cfg["PARTICLES"].as<std::string>();
		iwtropt.radius = vm_cfg["RADIUS"].as<double>();
		iwtropt.Nminp = vm_cfg["NMINPNTS"].as<int>();
		iwtropt.Nsearch = vm_cfg["NSEARCH"].as<int>();
		iwtropt.dir = vm_cfg["DIR"].as<int>();
		iwtropt.Npart = vm_cfg["NPART"].as<int>();
		iwtropt.bEndpnts = vm_cfg["ENDPNTS"].as<bool>();
		iwtropt.threshold = vm_cfg["THRES"].as<double>();
		iwtropt.wellRadius = vm_cfg["WELLRAD"].as<double>();
		iwtropt.maxSteps = vm_cfg["MAXSTEPS"].as<int>();
		iwtropt.stepSize = vm_cfg["STEPSIZE"].as<double>();
		iwtropt.maxStepSize = vm_cfg["MAXSTEPSIZE"].as<double>();
		iwtropt.minStepSize = vm_cfg["MINSTEPSIZE"].as<double>();
		iwtropt.tolStepSize = vm_cfg["TOLSTEPSIZE"].as<double>(); 
		iwtropt.pntPrecision = vm_cfg["PNTPRECISION"].as<int>();
		iwtropt.velPrecision = vm_cfg["VELPRECISION"].as<int>();
		iwtropt.nThreads = vm_cfg["NTHREADS"].as<int>();
	}
	else {
		std::cout << "IWFM Trace requires a configuration file with the following options:" << std::endl;
		std::cout << config_opt << std::endl;
		return false;
	}
	return tf;
}
