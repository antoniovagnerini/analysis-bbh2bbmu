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

void createHistograms();
void writeHistograms();

bool genPartPTordering(GenParticle p1, GenParticle p2) {return (p1.pt()>p2.pt());}
bool genJetPTordering(GenJet p1, GenJet p2) {return (p1.pt()>p2.pt());}

bool isSignal_;
std::map<std::string, TH1F*> h1_;

// =============================================================================================   
int main(int argc, char * argv[])
{
   TH1::SetDefaultSumw2();  // proper treatment of errors when scaling histograms
   
   isSignal_ = true;
   
   createHistograms();
   
   // Input files list
   std::string inputList = "rootFileList.txt";
   Analysis analysis(inputList);
   
   // Physics Objects Collections
   analysis.addTree<GenParticle> ("GenParticles","MssmHbb/Events/prunedGenParticles");
   analysis.addTree<GenJet> ("GenJets","MssmHbb/Events/slimmedGenJets");
   analysis.addTree<Jet> ("Jets","MssmHbb/Events/slimmedJetsReapplyJEC");

   
   // Analysis of events
   std::cout << "This analysis has " << analysis.size() << " events" << std::endl;
   int nevents = analysis.size();
   if ( nevents > 20000 ) nevents = 20000;
   for ( int i = 0 ; i < nevents ; ++i )
   {
      
      if ( i > 0 && i%10000 == 0 ) std::cout << i << " events processed..." << std::endl;
     
      analysis.event(i);
      
      // Generated Particles
      // The higgs (pdgId == 36) is also in such collection.
      // If higgs status=22, then hard-process intermediate (mass preserved)
      // If higgs status=62, then  particle with primordial kT included
      // The higgs daughters have initially status = 23, but they will undergo several steps 
      // until status=71, i.e. before hadronisation, which would be the best to use.
      // Thereofre, one can then try to follow the cascade down to the last b quarks
      // using delta_R matching. Maybe maybe some delta_pt/pt... but that's not done here.
      
      
      auto genParticles = analysis.collection<GenParticle>("GenParticles");
      GenParticle bquark;
      GenParticle bbarquark;
      std::vector<GenParticle> muons;
      
      for ( int j = 0 ; j < genParticles->size() ; ++j )
      {
         GenParticle gp = genParticles->at(j);
         if ( abs(gp.pdgId()) == 13 && gp.status() == 1 ) 
         {
            muons.push_back(gp);
         }
 
         if ( isSignal_ )
         {
            if ( gp.pdgId() == 36 && gp.status() > 60 )
            {
               h1_["h_hig_pt"]  -> Fill(gp.pt() );
               h1_["h_hig_eta"] -> Fill(gp.eta());
               h1_["h_hig_phi"] -> Fill(gp.phi());
            }
            if ( fabs(gp.pdgId()) == 5 && gp.higgsDaughter() ) 
            {
               if ( gp.pdgId() > 0 )  bquark = gp;
               if ( gp.pdgId() < 0 )  bbarquark = gp;
            }
         }
      }
      
      if ( isSignal_ )
      { 
      
         h1_["h_hig_d1_pt"]   -> Fill(bquark.pt() );
         h1_["h_hig_d1_eta"]  -> Fill(bquark.eta());
         h1_["h_hig_d1_phi"]  -> Fill(bquark.phi());
         
         h1_["h_hig_d2_pt"]   -> Fill(bbarquark.pt() );
         h1_["h_hig_d2_eta"]  -> Fill(bbarquark.eta());
         h1_["h_hig_d2_phi"]  -> Fill(bbarquark.phi());
         
         h1_["h_hig_dEta12"]  -> Fill(fabs(bquark.eta()-bbarquark.eta()));
         h1_["h_hig_dR12"]    -> Fill(bquark.p4().DeltaR(bbarquark.p4()));
         
      }
      
      
      // muons
      std::sort(muons.begin(),muons.end(),genPartPTordering);
      
      if ( muons.size() > 0 )
      {
         h1_["h_muon1_pt"]   -> Fill(muons[0].pt() );
         h1_["h_muon1_eta"]  -> Fill(muons[0].eta());
         h1_["h_muon1_phi"]  -> Fill(muons[0].phi());
      }
      
      
      
// ====================================================================      
      
      // Generated jets
      // Let's match with the Higgs daughters
      auto genJets = analysis.collection<GenJet>("GenJets");
      if ( genJets->size() > 1 )
      {
         GenJet j1 = genJets->at(0);
         GenJet j2 = genJets->at(1);
         
         TLorentzVector j12 = j1.p4() + j2.p4();
         double m12 = j12.M();
         
         h1_["h_jettrue_eta1"]  -> Fill(j1.eta());
         h1_["h_jettrue_eta2"]  -> Fill(j2.eta());  
         h1_["h_jettrue_pt1"]   -> Fill(j1.pt());
         h1_["h_jettrue_pt2"]   -> Fill(j2.pt());  
         h1_["h_jettrue_m12"]   -> Fill(m12);  
         h1_["h_jettrue_dEta12"]-> Fill(fabs(j1.eta()-j2.eta()));
         
      }
      
      
//       
// ====================================================================      
      
      
//       nb = 0;
//       nbbar = 0;
//       // reconstructed jets
//       // will now try to match the gen jets that are daughter of the higgs to the reconstructed jets
//       auto jets = analysis.collection<Jet>("Jets");
//       Jet bjet;
//       Jet bbarjet;
//       for ( int j = 0 ; j < jets->size() ; ++j )
//       {
//          Jet jet = jets->at(j);
//          if ( bgenjet.deltaR(jet) < 0.3 && nb == 0 )
//          {
//             bjet = jet ;
//             ++nb;
//          }
//          if ( bbargenjet.deltaR(jet) < 0.3 && nbbar == 0 )
//          {
//             bbarjet = jet ;
//             ++nbbar;
//          }
//          if ( nb > 0 && nbbar > 0 ) break;
//          
//       }
//       std::cout << "            I am a reco jet matched to a higgs daughter b quark! " << " : ";
//       std::cout << "pT  = "     << bjet.pt()      << ", ";
//       std::cout << "eta = "     << bjet.eta()     << ", ";
//       std::cout << "phi = "     << bjet.phi()     << ", " << std::endl;
//       
//       std::cout << "            I am a reco jet matched to a higgs daughter bbar quark! " << " : ";
//       std::cout << "pT  = "     << bbarjet.pt()      << ", ";
//       std::cout << "eta = "     << bbarjet.eta()     << ", ";
//       std::cout << "phi = "     << bbarjet.phi()     << ", " << std::endl;


//     std::cout << "=====" << std::endl;


  }
  
  writeHistograms();
   
//    
}

void createHistograms()
{
   // Higgs
   if ( isSignal_ )
   {
      h1_["h_hig_pt"]  = new TH1F("h_hig_pt","",100,0,200);
      h1_["h_hig_eta"] = new TH1F("h_hig_eta","",100,-10,10);
      h1_["h_hig_phi"] = new TH1F("h_hig_phi","",100,-3.15,3.15);
      
      h1_["h_hig_d1_pt"]  = new TH1F("h_hig_d1_pt","",100,0,200);
      h1_["h_hig_d1_eta"] = new TH1F("h_hig_d1_eta","",100,-10,10);
      h1_["h_hig_d1_phi"] = new TH1F("h_hig_d1_phi","",100,-3.15,3.15);
      
      h1_["h_hig_d2_pt"]  = new TH1F("h_hig_d2_pt","",100,0,200);
      h1_["h_hig_d2_eta"] = new TH1F("h_hig_d2_eta","",100,-10,10);
      h1_["h_hig_d2_phi"] = new TH1F("h_hig_d2_phi","",100,-3.15,3.15);
      
      h1_["h_hig_dEta12"]  = new TH1F("h_hig_dEta12","",100,0,10);
      h1_["h_hig_dR12"]    = new TH1F("h_hig_dR12","",100,0,5);
      
      
   }
   
   // Generated muon
   h1_["h_muon1_pt"]  = new TH1F("h_muon1_pt","",100,0,200);
   h1_["h_muon1_eta"] = new TH1F("h_muon1_eta","",100,-10,10);
   h1_["h_muon1_phi"] = new TH1F("h_muon1_phi","",100,-3.15,3.15);
   
   // gen jets
   h1_["h_jettrue_pt1"]   = new TH1F("h_jettrue_pt1","",200,0,200);
   h1_["h_jettrue_pt2"]   = new TH1F("h_jettrue_pt2","",200,0,200);
   
   h1_["h_jettrue_eta1"]   = new TH1F("h_jettrue_eta1","",100,-10,10);
   h1_["h_jettrue_eta2"]   = new TH1F("h_jettrue_eta2","",100,-10,10);
   h1_["h_jettrue_m12"]    = new TH1F("h_jettrue_m12","",200,0,2000);
   h1_["h_jettrue_dEta12"] = new TH1F("h_jettrue_dEta12","",100,0,10);
}

void writeHistograms()
{
   TFile * f_out = new TFile("analysis_generator.root","recreate");
   
   if ( isSignal_ )
   {
      f_out -> mkdir("higgs","higgs true info");
   }
   f_out -> mkdir("muons","muons true info");
   
   for (auto & ih1 : h1_)
   {
      f_out -> cd();
      if ( isSignal_ && TString(ih1.first).BeginsWith("h_hig") )  f_out -> cd("higgs");
      if ( TString(ih1.first).BeginsWith("h_muon") )              f_out -> cd("muons");
      
      ih1.second -> Write();
   }
   
   f_out -> Close();
}



