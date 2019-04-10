#ifndef IWFM_TRACK_OPTIONS_H
#define IWFM_TRACK_OPTIONS_H
#include <iostream>
#include <string>

#include <boost/program_options.hpp>

namespace po = boost::program_options;

/*!
 * \brief The iwfm_track_options struct holds all the parameters that this program needs
 */
struct iwfm_track_options{
    //! The main input configuration file
    std::string config_file;
    //! Te name of the Groundwater Zone budget file in hdf5 format.
    //! The program will read the mesh structure and the face flows from this file
    std::string GW_ZB_file;
    //! The IWFM Node file. Will read the node coordinates
    std::string Node_file;
    //! The IWFM stratigraphy file. Will read the layer elevations
    std::string Strat_file;
};

bool readInputParameters(int argc, char *argv[], iwfm_track_options& iwtropt){
    // Define the program options interface
    po::options_description iwfm_opt("All options");
    iwfm_opt.add_options()
            ("version,v", "print version string")
            ("help,h", "produce help message")
            ("gw-zbud,g", po::value<std::string >(), "GW_ZBudget.hdf output file")
            ("node,n", po::value<std::string >(), "IWFM node input file")
            ("strat,s", po::value<std::string >(), "C2Vsim stratigraphy file")
            ("config,c", po::value<std::string >(), "Main configuration file")
    ;

    // parse the input file
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, iwfm_opt), vm);

    if (vm.count("help")){
        std::cout << iwfm_opt << std::endl;
        return false;
    }

    if (vm.count("version")){
        std::cout << "|------------------|" << std::endl;
        std::cout << "|    IWFM Track    |" << std::endl;
        std::cout << "| Version : 0.0.01 |" << std::endl;
        std::cout << "|    by  giorgk    |" << std::endl;
        std::cout << "|__________________|" << std::endl;
        return 1;
    }

    // if there is a configuration file parse all options from it
    if (vm.count("config")){
        po::store(po::parse_config_file<char>(vm["config"].as< std::string >().c_str(), iwfm_opt),vm);
    }
    else{
        std::cout << "IWFM track requires a configuration file with the following options:" << std::endl;
        std::cout << iwfm_opt << std::endl;
        return false;
    }

    // Check that all required options are defined
    if (vm.count("node")){
        iwtropt.Node_file = vm["node"].as< std::string >();
    }
    else{
        std::cout << "node option is missing" << std::endl;
    }

    if (vm.count("strat")){
        iwtropt.Strat_file = vm["strat"].as< std::string >();
    }
    else{
        std::cout << "strat option is missing" << std::endl;
    }

    if (vm.count("gw-zbud")){
        iwtropt.GW_ZB_file = vm["gw-zbud"].as< std::string >();
    }
    else{
        std::cout << "gw-zbud option is missing" << std::endl;
    }

    return true;
}

#endif // IWFM_TRACK_OPTIONS_H
