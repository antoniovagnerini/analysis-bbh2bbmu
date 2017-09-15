#ifndef Analysis_Core_macro_config_h_
#define Analysis_Core_macro_config_h_

#include "boost/program_options.hpp"
#include <string>
#include <iostream>
#include <fstream>
#include <vector>

int macro_config(int argc, char * argv[]);

std::string cfg;

bool isMC;
bool signalregion;
std::string inputList;
std::string outputRoot;
std::string json;
int njets;
float btagwploose;
float btagwpmedium;
float btagwptight;
std::string btagwp;
std::string nonbtagwp;
std::string hltPath;

// int bWP = 1;
// float btagcut[3] = {0.46,0.84,0.92};
// // Cuts                                         // <<<<===== CMSDAS
// float ptmin[3]   = { 100.0, 100.0, 40.0 };
// float etamax[3]  = {   2.2,   2.2 , 2.2 };
// float btagmin[3] = { btagcut[bWP], btagcut[bWP], btagcut[bWP]};
// float nonbtag    = 0.46;
// float dRmin      = 1.;
// float detamax    = 1.55;

int macro_config(int argc, char * argv[])
{
   try
   {
      namespace po = boost::program_options;
      po::options_description desc("Options");
      desc.add_options()
         ("help,h","Show help messages")
         ("config,c",po::value<std::string>(&cfg),"Configuration file name");
      
      po::options_description config("Configuration");
      config.add_options()
         ("ntupleslist",po::value <std::string> (&inputList)->default_value("rootFileList.txt"),"File with list of ntuples")
         ("output",po::value <std::string> (&outputRoot)->default_value("histograms.root"),"Output root file")
         ("json",po::value <std::string> (&json)->default_value("json.txt"),"JSON file for data")
         ("nJets",po::value <int> (&njets)->default_value(3),"Minimum number of jets")
         ("isMC",po::value <bool> (&isMC)->default_value(true),"Flag for MC dataset")
         ("signalRegion",po::value <bool> (&signalregion)->default_value(true),"Flag for signal region")
         ("hltPath",po::value <std::string> (&hltPath),"HLT path name")
      
         ("btagWPLoose",po::value <float> (&btagwploose)->default_value(0.46),"BTag working point LOOSE")
         ("btagWPMedium",po::value <float> (&btagwpmedium)->default_value(0.84),"BTag working point MEDIUM")
         ("btagWPTight",po::value <float> (&btagwptight)->default_value(0.92),"BTag working point TIGHT");
      
      
      
      po::variables_map vm; 
      try
      {
         po::store(po::parse_command_line(argc, argv, desc), vm); // can throw
         // --help option
         
         if ( vm.count("help") )
         {
            std::cout << "SimpleBjetsAnalysis macro" << std::endl
                      << desc << std::endl
                      << config << std::endl;
            return 0;
         }
         po::notify(vm);
         
         std::ifstream cfg_s(cfg.c_str());
         po::store(po::parse_config_file(cfg_s, config), vm); // can throw
         if ( vm.count("config") )
         {

         }
         po::notify(vm);
      }
      catch(po::error& e)
      { 
         std::cerr << "ERROR: " << e.what() << std::endl << std::endl; 
         std::cerr << desc << std::endl; 
         return 1; 
      } 
      
   }
   catch(std::exception& e) 
   { 
      std::cerr << "ERROR: " << e.what() << std::endl << std::endl; 
      return 1; 
   } 

   
   return 0;
}

   
#endif
