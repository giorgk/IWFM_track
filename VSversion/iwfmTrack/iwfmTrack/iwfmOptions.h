#pragma once

#include <boost/program_options.hpp>

namespace po = boost::program_options;

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
	po::options_description user_opt("User options");
	user_opt.add_options()
		("version,v", "print version string")
		("help,h", "produce help message")
		("config,c", po::value<std::string >(), "Main configuration file")
		;

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, user_opt), vm);

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
		("NTHREADS", po::value<int>()->default_value(6), "Number of threads")
		;
	po::variables_map vm1;

	if (vm.count("help")) {
		std::cout << user_opt << std::endl;
		return false;
	}

	if (vm.count("version")) {
		std::cout << "|------------------|" << std::endl;
		std::cout << "|    IWFM Trace    |" << std::endl;
		std::cout << "| Version : 0.0.04 |" << std::endl;
		std::cout << "|    by  giorgk    |" << std::endl;
		std::cout << "|------------------|" << std::endl;
		return false;
	}

	bool tf = true;
	if (vm.count("config")) {
		po::store(po::parse_config_file<char>(vm["config"].as<std::string>().c_str(), config_opt), vm1);
		if (vm1.count("XY")) {
			iwtropt.XYfile = vm1["XY"].as<std::string>();
		}
		else {
			std::cout << "XY option is missing" << std::endl;
			tf = false;
		}

		if (vm1.count("MSH")) {
			iwtropt.MSHfile = vm1["MSH"].as<std::string>();
		}
		else {
			std::cout << "MSH option is missing" << std::endl;
			tf = false;
		}

		if (vm1.count("STRAT")) {
			iwtropt.STRATfile = vm1["STRAT"].as<std::string>();
		}
		else {
			std::cout << "STRAT option is missing" << std::endl;
			tf = false;
		}

		if (vm1.count("NORMALS")) {
			iwtropt.NRMLfile = vm1["NORMALS"].as<std::string>();
		}
		else {
			std::cout << "NORMALS option is missing" << std::endl;
			tf = false;
		}

		if (vm1.count("FACEZ")) {
			iwtropt.FACEZfile = vm1["FACEZ"].as<std::string>();
		}
		else {
			std::cout << "FACEZ option is missing" << std::endl;
			tf = false;
		}

		if (vm1.count("HFLOW")) {
			iwtropt.HFLOWfile = vm1["HFLOW"].as<std::string>();
		}
		else {
			std::cout << "HFLOW option is missing" << std::endl;
			tf = false;
		}

		if (vm1.count("VFLOW")) {
			iwtropt.VFLOWfile = vm1["VFLOW"].as<std::string>();
		}
		else {
			std::cout << "VFLOW option is missing" << std::endl;
			tf = false;
		}

		if (vm1.count("BC")) {
			iwtropt.BCfile = vm1["BC"].as<std::string>();
		}
		else {
			std::cout << "BC option is missing" << std::endl;
			tf = false;
		}

		if (vm1.count("BCZ")) {
			iwtropt.BCZfile= vm1["BCZ"].as<std::string>();
		}
		else {
			std::cout << "BCZ option is missing" << std::endl;
			tf = false;
		}

		if (vm1.count("OUTFILE")) {
			iwtropt.OUTfile = vm1["OUTFILE"].as<std::string>();
		}
		else {
			std::cout << "OUTFILE option is missing" << std::endl;
			tf = false;
		}

		// Options with default values
		iwtropt.WELLSfile = vm1["WELLS"].as<std::string>();
		iwtropt.PARTICLEfile = vm1["PARTICLES"].as<std::string>();
		iwtropt.radius = vm1["RADIUS"].as<double>();
		iwtropt.Nminp = vm1["NMINPNTS"].as<int>();
		iwtropt.Nsearch = vm1["NSEARCH"].as<int>();
		iwtropt.dir = vm1["DIR"].as<int>();
		iwtropt.Npart = vm1["NPART"].as<int>();
		iwtropt.bEndpnts = vm1["ENDPNTS"].as<bool>();
		iwtropt.threshold = vm1["THRES"].as<double>();
		iwtropt.wellRadius = vm1["WELLRAD"].as<double>();
		iwtropt.maxSteps = vm1["MAXSTEPS"].as<int>();
		iwtropt.stepSize = vm1["STEPSIZE"].as<double>();
		iwtropt.maxStepSize = vm1["MAXSTEPSIZE"].as<double>();
		iwtropt.minStepSize = vm1["MINSTEPSIZE"].as<double>();
		iwtropt.tolStepSize = vm1["TOLSTEPSIZE"].as<double>(); 
		iwtropt.pntPrecision = vm1["PNTPRECISION"].as<int>();
		iwtropt.velPrecision = vm1["VELPRECISION"].as<int>();
		iwtropt.nThreads = vm1["NTHREADS"].as<int>();
		if (!tf) {
			std::cout << "IWFM Trace requires a configuration file with the following options:" << std::endl;
			std::cout << config_opt << std::endl;
			return false;
		}
	}
	else {
		std::cout << "IWFM Trace requires a configuration file with the following options:" << std::endl;
		std::cout << config_opt << std::endl;
		return false;
	}
	return tf;
}
