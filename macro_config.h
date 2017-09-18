#ifndef Analysis_Core_macro_config_h_
#define Analysis_Core_macro_config_h_

#include "boost/program_options.hpp"
#include <string>
#include <iostream>
#include <fstream>
#include <vector>

int macro_config(int argc, char * argv[]);

std::string cfg_;

int nevtmax_;
bool isMC_;
bool signalregion_;
std::string inputlist_;
std::string outputRoot_;
std::string json_;
int njetsmin_;
int nbjetsmin_;
std::vector<float> jetsptmin_;
std::vector<float> jetsetamax_;
std::vector<float> jetsbtagmin_;

float drmin_;
float detamax_;
float btagwploose_;
float btagwpmedium_;
float btagwptight_;
float btagwp_;
float nonbtagwp_;


std::string hltPath_;
std::vector<std::string> triggerObjects_;

int macro_config(int argc, char * argv[])
{
   try
   {
      namespace po = boost::program_options;
      po::options_description desc("Options");
      desc.add_options()
         ("help,h","Show help messages")
         ("config,c",po::value<std::string>(&cfg_),"Configuration file name");
      
      po::options_description config("Configuration");
      config.add_options()
         ("ntuplesList",po::value <std::string> (&inputlist_)->default_value("rootFileList.txt"),"File with list of ntuples")
         ("nEventsMax",po::value <int> (&nevtmax_)->default_value(-1), "Maximum number of events")
         ("output",po::value <std::string> (&outputRoot_)->default_value("histograms.root"),"Output root file")
         ("json",po::value <std::string> (&json_)->default_value("json.txt"),"JSON file for data")
//      
         ("nJetsMin",po::value <int> (&njetsmin_)->default_value(3),"Minimum number of jets")
         ("nBJetsMin",po::value <int> (&nbjetsmin_)->default_value(0),"Minimum number of btgaged jets")
         ("jetsPtMin", po::value<std::vector<float> >(&jetsptmin_)->multitoken(),"Mimium pt of the jets")
         ("jetsEtaMax", po::value<std::vector<float> >(&jetsetamax_)->multitoken(),"Maximum |eta| of the jets")
         ("jetsBtagMin", po::value<std::vector<float> >(&jetsbtagmin_)->multitoken(),"Minimum btag of the jets; if < 0 -> reverse btag")
//      
         ("dRMin",po::value <float> (&drmin_)->default_value(1.0),"Minimum delta R between jets")
         ("dEtaMax",po::value <float> (&detamax_)->default_value(1.55),"Maximum delta eta between jets")
//      
         ("isMC",po::value <bool> (&isMC_)->default_value(true),"Flag for MC dataset")
         ("signalRegion",po::value <bool> (&signalregion_)->default_value(true),"Flag for signal region")
         ("hltPath",po::value <std::string> (&hltPath_),"HLT path name")
         ("triggerObjects", po::value<std::vector<std::string> >(&triggerObjects_)->multitoken(),"Trigger objects")
//      
         ("btagWPLoose",po::value <float> (&btagwploose_)->default_value(0.46),"BTag working point LOOSE")
         ("btagWPMedium",po::value <float> (&btagwpmedium_)->default_value(0.84),"BTag working point MEDIUM")
         ("btagWPTight",po::value <float> (&btagwptight_)->default_value(0.92),"BTag working point TIGHT")
//         
         ("btagWP",po::value <float> (&btagwp_)->default_value(0.8484),"Btag working point")
         ("nonbtagWP",po::value <float> (&nonbtagwp_)->default_value(0.46),"non-Btag working point");
      
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
            return 1;
         }
         po::notify(vm);
         
         std::ifstream cfg_s(cfg_.c_str());
         po::store(po::parse_config_file(cfg_s, config), vm); // can throw
         if ( vm.count("config") )
         {

         }
         po::notify(vm);
         if ( (int)jetsptmin_.size() != njetsmin_ )
         {
            std::cout << "Config Error *** Jet minimum pt were not defined or the definition does not match the minimum number of jets" <<std::endl;
            return -1;
         }
         if ( (int)jetsetamax_.size() != njetsmin_ )
         {
            std::cout << "Config Error *** Jet maximum |eta| were not defined or the definition does not match the minimum number of jets" <<std::endl;
            return -1;
         }
         if ( (int)jetsbtagmin_.size() != nbjetsmin_ && !jetsbtagmin_.empty() )
         {
            std::cout << "Config Error *** Jet minimum btag has been defined and does not match the minimum number of btagged jets" <<std::endl;
            return -1;
         }

      }
      catch(po::error& e)
      { 
         std::cerr << "ERROR: " << e.what() << std::endl << std::endl; 
         std::cerr << desc << std::endl; 
         return -1; 
      } 
      
   }
   catch(std::exception& e) 
   { 
      std::cerr << "ERROR: " << e.what() << std::endl << std::endl; 
      return -1; 
   } 
   
   return 0;
}

   
#endif
