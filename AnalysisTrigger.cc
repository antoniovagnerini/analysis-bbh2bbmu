#include <string>
#include <iostream>
#include <vector>

#include "TFile.h" 
#include "TFileCollection.h"
#include "TChain.h"
#include "TH1.h" 

#include "Analysis/Core/interface/Analysis.h"

using namespace std;
using namespace analysis;
using namespace analysis::tools;


// =============================================================================================   
int main(int argc, char * argv[])
{
   TH1::SetDefaultSumw2();  // proper treatment of errors when scaling histograms
   
   // Input files list
   std::string inputList = "rootfile.txt";
   Analysis analysis(inputList);
   
   // Physics Objects Collections
   analysis.triggerResults("MssmHbb/Events/TriggerResults");
   std::string hltPath = "HLT_Mu8_v";
   std::vector<std::string> l1seeds;
   l1seeds.push_back("L1_SingleMu3");
   l1seeds.push_back("L1_SingleMu5");
   l1seeds.push_back("L1_SingleMu7");
   
//    // Trigger objects
//    std::vector<std::string> jetTriggerObjects;
//    jetTriggerObjects.push_back("hltL1sDoubleJetC100");
//    jetTriggerObjects.push_back("hltDoubleJetsC100");
//    jetTriggerObjects.push_back("hltSingleBTagCSV0p84");
//    jetTriggerObjects.push_back("hltJetC350");
//    
//    std::string trgobj_path = "MssmHbb/Events/selectedPatTrigger/";
//    for ( auto & obj : jetTriggerObjects )
//       analysis.addTree<TriggerObject>(obj,trgobj_path+obj);
   
   // Analysis of events
   std::cout << "This analysis has " << analysis.size() << " events" << std::endl;
   int nevts = analysis.size();
   nevts = 1;
   for ( int i = 0 ; i < nevts ; ++i )
   {
      analysis.event(i);
      
      std::cout << "++++++    ENTRY  " << i;
      std::cout << std::endl;
      
      // Trigger results: fired and prescales
      // hltPath1
      int trg_fired = analysis.triggerResult(hltPath);
      std::map<std::string,int> l1ps = analysis.triggerPrescale(l1seeds);
      int hltps = analysis.triggerPrescale(hltPath);
      std::string s_accept = " fired ";
      if ( ! trg_fired ) s_accept = " did not fire ";
      std::cout << "The path " << hltPath << s_accept << "and has HLT PS = " << hltps << std::endl; 
      std::cout << "and its L1 seeds have the following prescales: " << endl;
      for ( auto & l1 : l1seeds )
         std::cout << "   " << l1 << ": " << l1ps[l1] << std::endl;
      
   }
   
//    
}

