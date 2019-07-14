#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>

#include <boost/program_options.hpp>

#include "BoucWen.hpp"

int main (int argc, char **argv)
{
     std::string inFileName, outFileName;
     std::vector<double> bwParameters;

     double DispTdT;
     BoucWen::NewtonRhapson::BoucWen_BaseClass<double> *BWen;
     
     // Arguments
     namespace po = boost::program_options;
     po::options_description desc("Allowed options");
     desc.add_options()
	  ("help,h", "produce help message")
	  ("input-file,i", po::value<std::string>(&inFileName), "input file")
	  ("output-file,o", po::value<std::string>(&outFileName), "output file")
	  ("params,p", po::value<std::vector<double>>(&bwParameters)->multitoken(), "BoucWen parameters")
	  ;

     po::variables_map vm;
     po::store(po::parse_command_line(argc, argv, desc), vm);
     po::notify(vm);    

     // Initialise BoucWen
     if (bwParameters.size() == 6) {
	  BWen = new BoucWen::NewtonRhapson::BoucWen<double>(bwParameters.data());
     } else if ( bwParameters.size() == 12) {
	  BWen = new BoucWen::NewtonRhapson::BoucWen_Deg<double>(bwParameters.data());
     } else if (bwParameters.size() == 18) {
	  BWen = new BoucWen::NewtonRhapson::BoucWen_BaberNoori<double>(bwParameters.data());
     } else {
	  std::cerr << "Wrong number of parameters for the BoucWewn model. Exiting." << std::endl;
	  std::exit(EXIT_FAILURE);
     }

     BWen->List();
	  
     // Handle files
     std::ifstream inFile(inFileName);
     std::ofstream outFile(outFileName);

     if (!inFile.is_open() || !outFile.is_open()){
	  if (!inFile.is_open()){
	       std::cerr << "Could not open: " << inFileName << ". Exiting." << std::endl;
	       std::exit(EXIT_FAILURE);
	  } else {
	       std::cerr << "Could not open: " << outFileName << ". Exiting." << std::endl;
	       std::exit(EXIT_FAILURE);
	  }
     }

     // Read until the end of file
     inFile >> DispTdT;
     while (!inFile.eof()){
	  outFile << DispTdT*5.0 << "\t" << BWen->Solver(DispTdT*5.0) << std::endl;
	  inFile >> DispTdT;
     }

     delete BWen;
     
     return 0;
}