#include "boost/program_options.hpp"
#include "boost/algorithm/string.hpp"
#include <string>
#include <iostream>
#include <fstream>
#include <vector>

#include "TFile.h" 
#include "TFileCollection.h"
#include "TChain.h"
#include "TH1.h" 

#include "Analysis/Core/interface/Analysis.h"
#include "Analysis/Core/bin/macro_config.h"

using namespace std;
using namespace analysis;
using namespace analysis::core;

// =============================================================================================   
int main(int argc, char * argv[])
{
     
   if ( macro_config(argc, argv) != 0 ) return -1;
   
   TH1::SetDefaultSumw2();  // proper treatment of errors when scaling histograms
   TH2::SetDefaultSumw2();

   // Cuts
   std::string btagalgo = btagalgo_;
   float btagmin[3] = { btagwp_, btagwp_, btagwp_};
   
   // Input files list
   Analysis analysis(inputlist_);
   
   analysis.addTree<Jet> ("Jets","MssmHbb/Events/slimmedJetsPuppi");
   analysis.addTree<Muon>("Muons","MssmHbb/Events/slimmedMuons");
        
   for ( auto & obj : triggerObjectsJets_ )
       analysis.addTree<TriggerObject> (obj,Form("MssmHbb/Events/slimmedPatTrigger/%s",obj.c_str()));
   for ( auto & obj : triggerObjectsMuons_ )
       analysis.addTree<TriggerObject> (obj,Form("MssmHbb/Events/slimmedPatTrigger/%s",obj.c_str()));
   
   //L1 splitting of jet
   // std::string l1jetob = "hltL1Mu12er2p3Jet40er2p3dRMax0p4DoubleJet40er2p3dEtaMax1p6";
   // std::string path = "MssmHbb/Events/slimmedPatTrigger/";
   // analysis.addTree<TriggerObject>( l1jetob.c_str(),std::string(path + l1jetob).c_str());
    
   analysis.triggerResults("MssmHbb/Events/TriggerResults");

   
   if( !isMC_ ) analysis.processJsonFile(json_);
   
   std::string sr_s = "SR";
   if ( ! signalregion_ ) sr_s = "CR";
   boost::algorithm::replace_last(outputRoot_, ".root", "_"+sr_s+".root"); 
   
   TFile hout(outputRoot_.c_str(),"recreate");
   
   std::map<std::string, TH1F*> h1;
   std::map<std::string, TH2F*> h2;

   // Jets --------------------------------------------------------------------------------------

   h1["n"]            = new TH1F("n" ,         "" , 30, 0, 30);      // # jets 
   h1["n_btagsel"]    = new TH1F("n_btagsel" , "" , 30, 0, 30);      //         after btagsel
   h1["n_mujsel"]     = new TH1F("n_mujsel"  , "" , 30, 0, 30);      //         after mujsel
   h1["n_ptmin20"]    = new TH1F("n_ptmin20" , "" , 30, 0, 30);                                                 // # jet pT>20    
   h1["n_ptmin20_btagsel"] = new TH1F("n_ptmin20_btagsel" , "" , 30, 0, 30);                                    //             after btag sel 
   h1["n_ptmin20_mujsel"]  = new TH1F("n_ptmin20_mujsel"  , "" , 30, 0, 30);                                    //             after mujsel
   h2["n_ptmin20_phivseta"]         = new TH2F("n_ptmin20_phivseta","", 25, -2.5, 2.5, 20, -4.0, 4.0 );         // # jet pT>20 in bins of eta and phi 
   h2["n_ptmin20_btagsel_phivseta"] = new TH2F("n_ptmin20_btagsel_phivseta","", 25, -2.5, 2.5, 20, -4.0, 4.0 ); //                   & after btag sel

   // Muons -----------------------------------------------------------------------------------

   h1["n_muons"]           = new TH1F("n_muons"  , "" , 30, 0, 30);                           // # muons passing ID
   h2["n_muons_phivseta"]  = new TH2F("n_muons_phivseta","", 25, -2.5, 2.5, 20, -4.0, 4.0 );  //         passing kinematic

   //after dR selection
   h1["pt_mu"]         = new TH1F("pt_mu"    , "" , 100, 0, 300);  
   h1["eta_mu"]        = new TH1F("eta_mu   ", "" , 100, -5, 5);  
   h2["eta_phi_mu"]    = new TH2F("eta_phi_mu","", 25, -2.5, 2.5, 20, -4.0, 4.0 ); 
   h1["dR_muj"]        = new TH1F("dR_muj" , "" , 100, 0, 0.5);   
   h1["dR_muj0"]       = new TH1F("dR_muj0", "" , 100, 0, 0.5);
   h1["dR_muj1"]       = new TH1F("dR_muj1", "" , 100, 0, 0.5);

   //after offline  matching to
   h1["pt_mumatched"]  = new TH1F("pt_mumatched" , "" , 100, 0, 300); 
   h1["eta_mumatched"] = new TH1F("eta_mumatched", "" , 100, -5, 5);  
   h1["n_muinjet"]     = new TH1F("n_muinjet" , "" , 30, 0, 30);      // # muons matched in jet selected
   h1["n_mujets"]      = new TH1F("n_mujets" , "" , 5, 0.5, 5.5);     // # real mu-jets per event 

 
  // LEADING Jets --------------------------------------------------------------------
 
  for ( int i = 0 ; i < njetsmin_ ; ++i )
   {
      //no sel
      h1[Form("pt_%i",i)]         = new TH1F(Form("pt_%i",i) , "" , 100, 0, 1000);
      h1[Form("pt_%i_eta0p9",i)]      = new TH1F(Form("pt_%i_eta0p9",i) , "" , 100, 0, 1000);
      h1[Form("pt_%i_0p9eta1p5",i)]     = new TH1F(Form("pt_%i_0p9eta1p5",i) , "" , 100, 0, 1000);
      h1[Form("pt_%i_1p5eta2p0",i)]     = new TH1F(Form("pt_%i_1p5eta2p0",i) , "" , 100, 0, 1000);
      h1[Form("pt_%i_2p0eta",i)]      = new TH1F(Form("pt_%i_2p0eta",i) , "" , 100, 0, 1000);
      h1[Form("eta_%i",i)]        = new TH1F(Form("eta_%i",i) , "" , 100, -5, 5);
      h1[Form("phi_%i",i)]        = new TH1F(Form("phi_%i",i) , "" , 100, -4, 4);
      h1[Form("csv_%i",i)]        = new TH1F(Form("csv_%i",i) , "" , 200, 0, 1);
      h1[Form("deepcsv_%i",i)]    = new TH1F(Form("deepcsv_%i",i) , "" , 200, 0, 1);
      h2[Form("eta_phi_%i",i)]       = new TH2F(Form("eta_phi_%i",i) , "" , 25, -2.5, 2.5, 20, -4.0, 4.0 );
      
      // after btagsel sel 
      h1[Form("pt_%i_btagsel",i)]     = new TH1F(Form("pt_%i_btagsel",i) , "" , 100, 0, 1000);
      h1[Form("eta_%i_btagsel",i)]    = new TH1F(Form("eta_%i_btagsel",i) , "" , 100, -5, 5);
      h1[Form("phi_%i_btagsel",i)]    = new TH1F(Form("phi_%i_btagsel",i) , "" , 100, -4, 4);
      h1[Form("csv_%i_btagsel",i)]    = new TH1F(Form("csv_%i_btagsel",i) , "" , 200, 0, 1);  
      h1[Form("deepcsv_%i_btagsel",i)]= new TH1F(Form("deepcsv_%i_btagsel",i) , "" , 200, 0, 1);
      h2[Form("eta_phi_%i_btagsel",i)] = new TH2F(Form("eta_phi_%i_btagsel",i) , "" , 25, -2.5, 2.5, 20, -4.0, 4.0 );
      
      // after mujet sel
      h1[Form("pt_%i_mujsel",i)]     = new TH1F(Form("pt_%i_mujsel",i) , "" , 100, 0, 1000);
      h1[Form("eta_%i_mujsel",i)]    = new TH1F(Form("eta_%i_mujsel",i) , "" , 100, -5, 5);
      h1[Form("phi_%i_mujsel",i)]    = new TH1F(Form("phi_%i_mujsel",i) , "" , 100, -4, 4);
      h1[Form("csv_%i_mujsel",i)]    = new TH1F(Form("csv_%i_mujsel",i) , "" , 200, 0, 1);
      h1[Form("deepcsv_%i_mujsel",i)]= new TH1F(Form("deepcsv_%i_mujsel",i) , "" , 200, 0, 1);
 
   }

   // Dijet masss -----------------------------------------------------------------------------------

   h1["m12"]     = new TH1F("m12"     , "" , 50, 0, 1000);
   h1["m12_btagsel"] = new TH1F("m12_btagsel" , "" , 50, 0, 1000);
   h1["m12_mujsel"] = new TH1F("m12_mujsel" , "" , 50, 0, 1000);
   h1["m12_mujsel_loweta"] = new TH1F("m12_mujsel_loweta" , "" , 50, 0, 1000);
   h1["m12_mujsel_mediumeta"] = new TH1F("m12_mujsel_mediumeta" , "" , 50, 0, 1000);
   h1["m12_mujsel_higheta"] = new TH1F("m12_mujsel_higheta" , "" , 50, 0, 1000);
   h1["m12_mujsel_mixedeta"] = new TH1F("m12_mujsel_mixedeta" , "" , 50, 0, 1000);

   double mbb;
   double weight;
   TTree *tree = new TTree("MssmHbb_13TeV","");
   tree->Branch("mbb",&mbb,"mbb/D");
   tree->Branch("weight",&weight,"weight/D");
   
   // Analysis of events
   std::cout << "This analysis has " << analysis.size() << " events" << std::endl;
   
   // Cut flow
   // 0: triggered events
   // 1: 3+ idloose jets
   // 2: kinematics
   // 3: matched to online
   // 4: delta R
   // 5: delta eta
   // 6: btag (bbnb)
      
   int nsel[10] = { };  
   int nmatch_jet[10] = { };
   int nmatch_mu[10] = { };
   int nmatch_mujet = 0;   
   
   if ( nevtmax_ < 0 ) nevtmax_ = analysis.size();
   for ( int i = 0 ; i < nevtmax_ ; ++i )
   { 
      int njets = 0;
      int njets_btagsel = 0;
      int njets_mujsel = 0;
      int nmuons =0;
      bool goodEvent = true;
      bool muonJetEvent = false;   
       
      if ( i > 0 && i%100000==0 ) std::cout << i << "  events processed! " << std::endl;
     
      // if (analysis.run() > 304209 ) continue;
      analysis.event(i);
                  if (! isMC_ )
      {
         if (!analysis.selectJson() ) continue; // To use only goodJSonFiles
      }
	    
      // #0 TRIGGERED EVTs     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      int triggerFired = analysis.triggerResult(hltPath_);
      if ( !triggerFired ) continue;
      
      ++nsel[0];
      // =============== JET SELECTION =======================================================================
      // #1 ID LOOSE JET  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++      
      // match offline to online
      // analysis.match<Jet,TriggerObject>("Jets",triggerObjects_,0.5);
      
      // Jets - std::shared_ptr< Collection<Jet> >
      auto slimmedJets = analysis.collection<Jet>("Jets");
      auto slimmedMuons = analysis.collection<Muon>("Muons");
      // auto l1objs = analysis.collection<TriggerObject>(l1jetob);
   
      std::vector<Jet *> selectedJets;
      for ( int j = 0 ; j < slimmedJets->size() ; ++j )
      {
         if ( slimmedJets->at(j).idLoose() ) selectedJets.push_back(&slimmedJets->at(j));
      }
      if ( (int)selectedJets.size() < njetsmin_ ) continue;
      
      ++nsel[1];
      // #2 KINEMATICS SEL  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      // Kinematic selection - 3 leading jets
      for ( int j = 0; j < njetsmin_; ++j )
      {
         Jet * jet = selectedJets[j];
         if ( jet->pt() < jetsptmin_[j] || fabs(jet->eta()) > jetsetamax_[j] )
         {
            goodEvent = false;
            break;
         }
      }
      
      if ( ! goodEvent ) continue;
      
      ++nsel[2];
      // #3 DeltaR SEL +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      for ( int j1 = 0; j1 < njetsmin_-1; ++j1 )
      {
         const Jet & jet1 = *selectedJets[j1];
         for ( int j2 = j1+1; j2 < njetsmin_; ++j2 )
         {
            const Jet & jet2 = *selectedJets[j2];
            if ( jet1.deltaR(jet2) < drmin_ ) goodEvent = false;
         }
      }
      
      if ( ! goodEvent ) continue;
      
      ++nsel[3];
      // #4 Deta cut ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      
      if ( fabs(selectedJets[0]->eta() - selectedJets[1]->eta()) > detamax_ ) continue;
      
      ++nsel[4];

      // #5 BTAG SEL  ++++++++++++++++++++++++++++++++++++++++++

      // Fill histograms of JET SELECTION
      for ( int j = 0 ; j < (int)selectedJets.size() ; ++j )
      {
         if ( selectedJets[j]->pt() < 20. ) continue;
	 h2["n_ptmin20_phivseta"] -> Fill ( selectedJets[j]->eta(), selectedJets[j]->phi() );
         ++njets;
      }
      
      h1["n"] -> Fill(selectedJets.size());
      h1["n_ptmin20"] -> Fill(njets); // all jets but kinematic cuts only on the j012

      for ( int j = 0; j < njetsmin_; ++j )
      {
         Jet * jet = selectedJets[j];
         h1[Form("pt_%i",j)]   -> Fill(jet->pt());
	 if (fabs(jet->eta()) < 0.9)                                     h1[Form("pt_%i_eta0p9",j)]   -> Fill(jet->pt());
         else if ((fabs(jet->eta()) > 0.9) && (fabs(jet->eta()) < 1.5))  h1[Form("pt_%i_0p9eta1p5",j)]-> Fill(jet->pt());
	 else if ((fabs(jet->eta()) > 1.5) && (fabs(jet->eta()) < 2.0))  h1[Form("pt_%i_1p5eta2p0",j)]->Fill(jet->pt());
	 else if (fabs(jet->eta()) > 2.0)                                h1[Form("pt_%i_2p0eta",j)]   -> Fill(jet->pt());
         h1[Form("eta_%i",j)]  -> Fill(jet->eta());
         h1[Form("phi_%i",j)]  -> Fill(jet->phi());
         h1[Form("csv_%i",j)] -> Fill(jet->btag());
	 h1[Form("deepcsv_%i",j)] -> Fill(jet->btag("btag_deepb")+jet->btag("btag_deepbb"));
	 h2[Form("eta_phi_%i",j)] -> Fill(jet->eta(),jet->phi());
         
	 //btag selection
      	 float btagj; 
      	 if      ( btagalgo == "csv" )      btagj = jet->btag();
         else if ( btagalgo == "deepcsv" )  btagj = jet->btag("btag_deepb")+jet->btag("btag_deepbb");
      	 else return -1;
     
         if ( j < 2 && btagj < btagmin[j] )     goodEvent = false;
         if ( ! signalregion_ )
         {
            if ( j == 2 && btagj > nonbtagwp_ ) goodEvent = false; 
         }
         else
         {
            if ( j == 2 && btagj < btagmin[j] ) goodEvent = false; 
         }   
      }

      h1["m12"] -> Fill((selectedJets[0]->p4() + selectedJets[1]->p4()).M());

      if ( ! goodEvent ) continue;
     
      //+++++++++++++++++++++++++++++++++++++++++++
      // Fill histograms of B-JET SELECTION but no MATCHING
            
      for ( int j = 0 ; j < (int)selectedJets.size() ; ++j )
      {
         if ( selectedJets[j]->pt() < 20. ) continue;
	 h2["n_ptmin20_btagsel_phivseta"] -> Fill ( selectedJets[j]->eta(), selectedJets[j]->phi() );
         ++njets_btagsel;
      }
      h1["n_btagsel"] -> Fill(selectedJets.size());   // all jets in evt but leading bbnb tag
      h1["n_ptmin20_btagsel"] -> Fill(njets_btagsel); // all jets pt >20 but leading bbnb tag
      for ( int j = 0; j < njetsmin_; ++j )
      {
         Jet * jet = selectedJets[j];
         h1[Form("pt_%i_btagsel",j)]   -> Fill(jet->pt());
         h1[Form("eta_%i_btagsel",j)]  -> Fill(jet->eta());
         h1[Form("phi_%i_btagsel",j)]  -> Fill(jet->phi());
         h1[Form("csv_%i_btagsel",j)] -> Fill(jet->btag());
	 h1[Form("deepcsv_%i_btagsel",j)] -> Fill(jet->btag("btag_deepb")+jet->btag("btag_deepbb"));
	 h2[Form("eta_phi_%i_btagsel",j)] -> Fill(jet->eta(),jet->phi());

      }
      mbb = (selectedJets[0]->p4() + selectedJets[1]->p4()).M();
      if ( !signalregion_ )
      { 
         h1["m12_btagsel"] -> Fill(mbb);
	 //    weight = 1;
         //tree -> Fill();
      }
           
      ++nsel[5];

      //========== MUON SELECTION ============================================================================  
      
      //Muon ID sel +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      std::vector<Muon *> selectedMuons;
      for ( int m = 0 ; m < slimmedMuons->size() ; ++m )
      {	 

	//       if ( i > 0 && i%100==0 ) std::cout << i << "\t" << "slimmed" << "\t" <<slimmedMuons->at(m).pt() << std::endl;
       if ( slimmedMuons->at(m).isLooseMuon() )
	 { 
	 selectedMuons.push_back(&slimmedMuons->at(m));	    
	 ++nmuons;
	 }
      }
      if ( (int)selectedMuons.size() < 1 ) continue;
      h1["n_muons"] -> Fill(nmuons);

      //Muon-jet sel +++++++++++++++++++++++++++++++++++++++++++++++++++++
      std::vector<Muon *> selectedMuonsinJet;
            for (  size_t m = 0; m < selectedMuons.size(); ++m )      
       {
         Muon * muon = selectedMuons[m];
	 if ( muon->pt() < muonsptmin_[m] || fabs(muon->eta()) > muonsetamax_[m] ) continue;
	 
	 //	 if ( i > 0 && i%100==0 ) std::cout << i <<  "\t" << " kinem" <<  "\t" << muon->pt() << std::endl;
	 h2["n_muons_phivseta"]-> Fill(muon->eta(), muon->phi()) ;

	 float dR_muj0 = selectedJets[0]->deltaR(*muon) ;
	 float dR_muj1 = selectedJets[1]->deltaR(*muon) ;
	 
	 if ( dR_muj0 < drmax_  || dR_muj1 < drmax_ ) //at least 1 muon jet
            {
	      muonJetEvent  = true;

              h1["pt_mu"]  -> Fill( muon->pt());
              h1["eta_mu"] -> Fill( muon->eta());
              h2["eta_phi_mu"] -> Fill( muon->eta(), muon->phi());
	      h1["dR_muj"] -> Fill( (dR_muj0 < dR_muj1) ? dR_muj0 : dR_muj1 );
	      if       (dR_muj0 < drmax_)  h1["dR_muj0"] -> Fill( dR_muj0 );
	      else if  (dR_muj1 < drmax_)  h1["dR_muj1"] -> Fill( dR_muj1 );

      //	      if ( i > 0 && i%100==0 ) std::cout << i <<"\t" << " dR cut" << "\t" << muon->pt() << std::endl;
	      selectedMuonsinJet.push_back(muon);
	      break;
	      } 
       }     
      
      if ( ! muonJetEvent ) continue;

      ++nsel[6];       
      // #6 MATCHING OFFLINE to ONLINE  +++++++++++++++++++++++++++++++++++++++++++++++++++++
        
      //Jet L1 splitting 
	// for ( int L1 = 0 ; L1 < l1objs->size() ; ++L1 )
	//   {
        //     TriggerObject l1obj = l1objs->at(L1);
        //     if ( l1obj.type() == -99 ) analysis.match<Jet,TriggerObject>("Jets",l1jetob,0.5);    //l1jet seed is -99
	//   }

      analysis.match<Jet,TriggerObject>("Jets",triggerObjectsJets_,0.5);
      analysis.match<Muon,TriggerObject>("Muons",triggerObjectsMuons_,0.5);
      
      //DIJET loop -- NO MATCH if any of the 2 jets is not matching 
      //DR jet-mu  -- # mu-jets
      bool matched_jet[10] = {true,true,true,true,true,true,true,true,true,true}; //ntriggerobjects=10
      bool matched_mujet[2] = {false,false};  //njets =2 
      for ( int j = 0; j < 2; ++j )
      {
         Jet * jet = selectedJets[j];

	 //L1 JET
	 //	 if ( ! jet->matched(l1jetob) )                    matched_jet[0]    = false;
	 // if ( ! jet->matched(triggerObjectsJets_[io]) ) matched_jet[io+1] = false;

	 //DIJET loop	 
         for ( size_t io = 0; io < triggerObjectsJets_.size()-1 ; ++io ) //exclude hltBSoftMuonDiJet40Mu12L3FilterByDR
         {  
	    if ( ! jet->matched(triggerObjectsJets_[io]) ) matched_jet[io] = false;
         }
	 //DR jet-mu
	 if (jet->matched("hltBSoftMuonDiJet40Mu12L3FilterByDR")) matched_mujet[j] = true; 
      }

      //MUON loop -- MATCH if any of the muon matches  
      bool matched_mu[10] = {false,false,false,false,false,false,false,false,false,false};  //ntriggerobjects=10
      for ( size_t m = 0; m < selectedMuonsinJet.size(); ++m )
      {
         Muon * muon = selectedMuonsinJet[m];
	 	 
         for ( size_t io = 0; io < triggerObjectsMuons_.size() ; ++io )
         {       
	   if ( muon->matched(triggerObjectsMuons_[io]) ) matched_mu[io] = true;
         }
      }
   

      //---------------------------------------------------------------------------------------------------------------   
      //TRIGGER OBJECT JET & MUON loop 

      for ( size_t io = 0; io < triggerObjectsJets_.size() ; ++io ) // L1_jet included (+1) but jetbyDR excluded from list (-1)
      {
         if ( matched_jet[io] ) ++nmatch_jet[io];
         goodEvent = ( goodEvent && matched_jet[io] ); //set evt good if both jet matched & initially Goodevent
      }
      
      for ( size_t io = 0; io < triggerObjectsMuons_.size() ; ++io )  
      {
         if ( matched_mu[io] ) ++nmatch_mu[io];
         goodEvent = ( goodEvent && matched_mu[io] ); //set evt good if at least one goodMu matches & IF MATCHED DIJET
      }

      // if (! goodEvent) continue;

      // Mu-jets matching  
      int nmujets = 0;
      for ( int j = 0; j < 2; ++j ) 
	{
	  if ( matched_mujet[j] ) ++nmujets;
	  if ( matched_mujet[j] && nmujets==1 ) ++nmatch_mujet; //event with at least one mu-jet 
  	}

      if (!goodEvent || !(matched_mujet[0] || matched_mujet[1]) ) continue;


      //NEED TO ADD DOUBLE-MUJET CASE FOR EFFICIENCY CALC.
      ++nsel[7]; 

      
     //=============================================================================================== 
       // Fill histograms of passed MUJETS  selection
       for ( int j = 0 ; j < (int)selectedJets.size() ; ++j )
       {
          if ( selectedJets[j]->pt() < 20. ) continue;
          ++njets_mujsel;
       }

      h1["n_mujets"]  -> Fill(nmujets);                   // number real mujets per evt
      h1["n_muinjet"] -> Fill(selectedMuonsinJet.size()); // number of muons in dR cone
      h1["n_mujsel"]  -> Fill(selectedJets.size());      // all jets but leading 3 passing bbmunb selection
      h1["n_ptmin20_mujsel"] -> Fill(njets_mujsel);      //all jets pt > 20 but butleading 3 passing bbmunb selection 
  
      for ( int j = 0; j < njetsmin_; ++j )
      {
         Jet * jet = selectedJets[j];
         h1[Form("pt_%i_mujsel",j)]   -> Fill(jet->pt());
         h1[Form("eta_%i_mujsel",j)]  -> Fill(jet->eta());
         h1[Form("phi_%i_mujsel",j)]  -> Fill(jet->phi());
         h1[Form("csv_%i_mujsel",j)] -> Fill(jet->btag());
	 h1[Form("deepcsv_%i_mujsel",j)] -> Fill(jet->btag("btag_deepb")+jet->btag("btag_deepbb"));
      }
      mbb = (selectedJets[0]->p4() + selectedJets[1]->p4()).M();
      if ( !signalregion_ )
      { 
         h1["m12_mujsel"] -> Fill(mbb);
         weight = 1;
         tree -> Fill();
      }
      // both jets in low, medium or high eta region
      float eta_j0 = fabs(selectedJets[0]->eta());
      float eta_j1 = fabs(selectedJets[1]->eta());

      if      ( eta_j0 < 0.9  && eta_j1 < 0.9 )                                      h1["m12_mujsel_loweta"]   -> Fill(mbb);
      else if ( eta_j0 > 0.9  && eta_j0 < 1.5 && eta_j1 > 0.9 && eta_j1 < 1.5 )      h1["m12_mujsel_mediumeta"]-> Fill(mbb);
      else if ( eta_j0 > 1.5  && eta_j1 > 1.5 )                                      h1["m12_mujsel_higheta"]  -> Fill(mbb);
      else                                                                           h1["m12_mujsel_mixedeta"] -> Fill(mbb);

      for (  size_t m = 0; m < selectedMuonsinJet.size(); ++m )      
       {
         Muon * muon = selectedMuonsinJet[m];
	 //	 if ( i > 0 && i%100==0 ) std::cout << i  << "\t" << " matched" << "\t" << muon->pt() << std::endl;
	 //	 if ( i > 0 && i%100==0 ) std::cout << "______________________________________________________" << std::endl;
         h1["pt_mumatched"]  -> Fill( muon->pt() );
         h1["eta_mumatched"] -> Fill( muon->eta() );
         //h1["dR_muj"] -> Fill( (dR_muj0 < dR_muj1) ? dR_muj0 : dR_muj1 );
       } 

         
   } //end event loop
  
           
   for (auto & ih1 : h1)
   {
      ih1.second -> Write();
   }
   
   hout.Write();
   hout.Close();
   
// PRINT OUTS  ============================================================================================ 
   
   // Cut flow
   // 0: triggered events
   // 1: 3+ idloose jets
   // 2: matched to online
   // 3: kinematics
   // 4: delta R
   // 5: delta eta
   // 6: btag (bbnb)
   
   double fracAbs[10];
   double fracRel[10];
   std::string cuts[10];
   cuts[0] = "Triggered";
   cuts[1] = "Triple idloose-jet";
   cuts[2] = "Triple jet kinematics";
   cuts[3] = "Delta R(i;j)";
   cuts[4] = "Delta eta(j1;j2)";
   cuts[5] = "btagged (bbnb)";
   if ( signalregion_ ) cuts[5] = "btagged (bbb)";
   cuts[6] = "Muon jet sel";
   cuts[7] = "Match online tobj ";

   printf ("%-23s  %10s  %10s  %10s \n", std::string("Cut flow").c_str(), std::string("# events").c_str(), std::string("absolute").c_str(), std::string("relative").c_str() ); 
   for ( int i = 0; i < 8; ++i )
   {
      fracAbs[i] = double(nsel[i])/nsel[0];
      if ( i>0 )
         fracRel[i] = double(nsel[i])/nsel[i-1];
      else
         fracRel[i] = fracAbs[i];
      printf ("%-23s  %10d  %10.3f  %10.3f \n", cuts[i].c_str(), nsel[i], fracAbs[i], fracRel[i] ); 
   }

   // Trigger objects counts   
  std::cout << std::endl;
  printf ("%-40s  %10s \n", std::string("Trigger object dijet").c_str(), std::string("# events").c_str() );
  // printf ("%-40s  %10d \n", l1jetob.c_str(), nmatch_jet[0] );
  //   printf ("%-40s  %10d \n", triggerObjectsJets_[io].c_str(), nmatch_jet[io+1] ); 
  for ( size_t io = 0; io < triggerObjectsJets_.size()-1 ; ++io )
  {
    printf ("%-40s  %10d \n", triggerObjectsJets_[io].c_str(), nmatch_jet[io] ); 
  } 
  printf ("%-40s  %10d \n", "hltBSoftMuonDiJet40Mu12L3FilterByDR", nmatch_mujet ); 
  
    std::cout << std::endl;
   printf ("%-40s  %10s \n", std::string("Trigger object mu").c_str(), std::string("# events").c_str() ); 
  for ( size_t io = 0; io < triggerObjectsMuons_.size() ; ++io )
  {
     printf ("%-40s  %10d \n", triggerObjectsMuons_[io].c_str(),  nmatch_mu[io] ); 
  }     
   
   
   
      
   
} //end main

