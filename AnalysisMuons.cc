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
using namespace analysis::core;


// =============================================================================================   
int main(int argc, char * argv[])
{
   double ptmin = 0.;
   
   bool isMC = false;
   TH1::SetDefaultSumw2();  // proper treatment of errors when scaling histograms
   
   // Input files list
   std::string inputList = "rootFileList.txt";
   Analysis analysis(inputList);
   
   // Physics Objects Collections
   analysis.addTree<Muon> ("Muons","MssmHbb/Events/slimmedMuons");
   
//    // Trigger results
//    analysis.triggerResults("MssmHbb/Events/TriggerResults");
//    std::string hltPath = "HLT_DoubleJetsC100_DoubleBTagCSV_p014_DoublePFJetsC100MaxDeta1p6_v";
//    
//    // Trigger objects
//    std::vector<std::string> triggerObjects;
//    triggerObjects.push_back("hltL1sDoubleJetC100");
//    triggerObjects.push_back("hltDoubleJetsC100");
//    triggerObjects.push_back("hltBTagCaloCSVp014DoubleWithMatching");
//    triggerObjects.push_back("hltDoublePFJetsC100");
//    triggerObjects.push_back("hltDoublePFJetsC100MaxDeta1p6");
//    
//    std::string trgobj_path = "MssmHbb/Events/selectedPatTrigger/";
//    for ( auto & obj : triggerObjects )
//       analysis.addTree<TriggerObject>(obj,trgobj_path+obj);
//    
//    // Certified lumis
//    if ( ! isMC ) analysis.processJsonFile("json.txt");
   

   TFile hout("histograms_muons.root","recreate");
   
   std::map<std::string, TH1F*> h1;
   h1["n"]        = new TH1F("n" , "" , 30, 0, 30);
   h1["pt"]       = new TH1F("pt" , "" , 100, 0, 770);
   h1["eta"]      = new TH1F("eta" , "" , 100, -5.30, 3.30);
   h1["phi"]      = new TH1F("phi" , "" , 100, -3.7, 3.7);
   h1["q"]        = new TH1F("q", "", 4, -2, 2);
   h1["isPF"]       = new TH1F("isPF","",2,0,2);
   h1["isGlobal"]   = new TH1F("isGlobal","",2,0,2);
   h1["isTracker"]  = new TH1F("isTracker","",2,0,2);
   h1["isLoose"]    = new TH1F("isLoose","",2,0,2);
   h1["isMedium"]   = new TH1F("isMedium","",2,0,2);
   
   
   // Analysis of events
   std::cout << "This analysis has " << analysis.size() << " events" << std::endl;
   int nevents =  analysis.size();

//    int nFired = 0;
   for ( int i = 0 ; i < nevents ; ++i )
   {
      if ( i > 0 && i%100000 == 0 ) std::cout << i << " events processed..." << std::endl;
      
      analysis.event(i);
      
//       if ( !isMC )
//       {
//          if ( !analysis.selectJson() ) continue;
//       }

//       // Trigger results
//       // hltPath
//       int trgFired = analysis.triggerResult(hltPath);
//       if ( ! trgFired ) continue;
//       ++nFired;
            
//       std::cout << "++++++    ENTRY  " << i;
//       std::cout << std::endl;
      
      // Muons
      auto muons = analysis.collection<Muon>("Muons");
      int nmuons = 0;
      for ( int m = 0 ; m < muons->size() ; ++m )
      {
         Muon muon = muons->at(m);
         h1["pt"]         -> Fill(muon.pt());
         h1["eta"]        -> Fill(muon.eta());
         h1["phi"]        -> Fill(muon.phi());
         h1["q"]          -> Fill(muon.q());
         h1["isPF"]       -> Fill(muon.isPFMuon()     );
         h1["isGlobal"]   -> Fill(muon.isGlobalMuon() );
         h1["isTracker"]  -> Fill(muon.isTrackerMuon());
         h1["isLoose"]    -> Fill(muon.isLooseMuon()  );
         h1["isMedium"]   -> Fill(muon.isMediumMuon() );
         ++nmuons;
         
      }
      h1["n"] -> Fill(nmuons);
      
   }
   
   for (auto & ih1 : h1)
   {
      ih1.second -> Write();
   }
   
//    
}

