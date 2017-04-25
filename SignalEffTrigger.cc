#include <string>
#include <iostream>
#include <vector>

#include "TFile.h" 
#include "TFileCollection.h"
#include "TChain.h"
#include "TH1.h" 
#include "TGraphAsymmErrors.h" 

#include "Analysis/Core/interface/Analysis.h"

using namespace std;
using namespace analysis;
using namespace analysis::tools;


// =============================================================================================   
int main(int argc, char * argv[])
{
   TH1::SetDefaultSumw2();  // proper treatment of errors when scaling histograms
   
   // Input files list
   std::string inputList = "rootFileList.txt";
   Analysis analysis(inputList,"MssmHbb/Events/EventInfo");
   
   // Physics Objects Collections
   analysis.triggerResults("MssmHbb/Events/TriggerResults");
   std::string hltPath = "HLT_DoubleJetsC100_DoubleBTagCSV_p014_DoublePFJetsC100MaxDeta1p6_v";
   
   // BTag
//    analysis.addTree<JetTag> ("JetsTags","MssmHbbTrigger/Events/hltCombinedSecondaryVertexBJetTagsCalo");
   
   
   // Trigger objects
   std::vector<std::string> jetTriggerObjects;
   jetTriggerObjects.push_back("hltL1sDoubleJetC100");
   jetTriggerObjects.push_back("hltDoubleJetsC100");
   jetTriggerObjects.push_back("hltSingleBTagCSV0p84");
   jetTriggerObjects.push_back("hltJetC350");
   
//    std::string trgobj_path = "MssmHbbTrigger/Events/selectedPatTrigger/";
//    for ( auto & obj : jetTriggerObjects )
//       analysis.addTree<TriggerObject>(obj,trgobj_path+obj);
   
   TH1F * h_total = new TH1F("h_total","",70,0,70);
   TH1F * h_selected = new TH1F("h_selected","",70,0,70);
   
   // Analysis of events
   std::cout << "This analysis has " << analysis.size() << " events" << std::endl;
   for ( int i = 0 ; i < analysis.size() ; ++i )
   {
      analysis.event(i);
      
//      std::cout << "++++++    ENTRY  " << i;
//      std::cout << std::endl;
      
      int nPU = analysis.nPileup();
      
      h_total -> Fill(nPU);
      
      // Trigger results: fired and prescales
      // hltPath
      if ( ! analysis.triggerResult(hltPath) ) continue;
      
//       std::shared_ptr< Collection<JetTag> > jetstags   = analysis.collection<JetTag>("JetsTags");
//       std::vector<JetTag>  jetstagswp092;
//       for ( int jt = 0; jt < jetstags->size() ; ++jt )
//       {
//          JetTag jet = jetstags->at(jt);
//          if ( jet.pt() >= 80 )
//          {
//             if ( jet.btag() >= 0.92 ) jetstagswp092.push_back(jet);
//          }
//       }
      
//       if ( jetstagswp092.size() < 2 ) continue;
      
      h_selected -> Fill(nPU);
   }
   
   TGraphAsymmErrors * g_eff = new TGraphAsymmErrors(h_selected,h_total,"cl=0.683 b(1,1) mode");
   
   TFile * out = new TFile("signalEffTrigger.root","recreate");
   h_total -> Write();
   h_selected -> Write();
   g_eff -> Write();
   
   out->Close();
   
//    
}

