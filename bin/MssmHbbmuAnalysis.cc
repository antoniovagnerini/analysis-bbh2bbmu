#include "boost/program_options.hpp"
#include "boost/algorithm/string.hpp"
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <ctime>

#include "TFile.h" 
#include "TFileCollection.h"
#include "TChain.h"
#include "TH1.h" 
#include "TRandom3.h"
#include "TMath.h"

#include "Analysis/Tools/interface/Analysis.h"
#include "Analysis/Tools/bin/macro_config.h"

using namespace std;
using namespace analysis;
using namespace analysis::tools;

// Functions declarations
void ScaleHistograms(const std::string & );
void CreateHistograms();
void WriteHistograms(const std::string & );

float btagMin(const std::string & btagwp);
void correctJetpt ( Jet & , const float & );
void copyJetpt    ( Jet & , Jet & );
double onlSFbtag ( double   , double   ) ;
double muonptWeight ( double  , double );
double jetTriggerScaleFactor ( double , double , double );

//float puWeight(float I, const string & var);
void applyPrescale ( const int& run, const double& random, float& prescaleEra, const int nsubsamples, int& window );

std::map<std::string, TH1F*> h1;
std::map<std::string, TH2F*> h2;

// auxiliar variables
int nGenTotal_;
float Efficiency_;
float lumi_  = 36500.;
int nsubsamples = 7;
std::vector<std::pair<int,int>> eraRanges = {{299337,302029},{302030,303434},{303435,304826},{304911,306462}}; //C,D,E,F                                              


// Systematic variations
std::vector<std::string> systematics = { "PU", "SFbtag", "onlSFbtag", "JER", "JES", "muID", "mu_trigeff", "jet_trigeff" };
//std::vector<std::string> systematics = { "CMS_PU", "CMS_scale_j", "CMS_res_j", "CMS_eff_b", "CMS_eff_bonl" };
std::vector<std::string> vars = { "up", "down" };


template <class T> struct greaterByPt {  bool operator() (const T& p1, const T& p2) const {return p1->pt() > p2->pt();} };

bool debug = false;

// =============================================================================================   
int main(int argc, char * argv[])
{
     
   if ( macro_config(argc, argv) != 0 ) return -1;
   
   TH1::SetDefaultSumw2();  // proper treatment of errors when scaling histograms
   TH2::SetDefaultSumw2();

   // Cuts
   std::string btagalgo = btagalgo_;
   //   float btagmin[3] = { btagwp_, btagwp_, btagwp_};
   float btagmin = btagMin(btagwp_);
   //   float ptj4max = ptmax_[0];

   // Input files list
   Analysis analysis(inputlist_);
  
   auto bsf_reader = analysis.btagCalibration(btagalgo_, btagsf_, btagwp_);
   auto bsf_reader_medium = analysis.btagCalibration(btagalgo_, btagsf_, "medium");
   auto jerinfo = analysis.jetResolutionInfo(jerpt_,jersf_);
   std::shared_ptr <PileupWeight> puweights = analysis.pileupWeights("/nfs/dust/cms/user/vagneria/store/CMSSW_9_4_2/src/Analysis/Tools/data/PileupWeight_Run2017CDEF_Mix_2017.root");
   std::shared_ptr <MuonIdWeight> muIDweights = analysis.muonIDWeights("/nfs/dust/cms/user/vagneria/store/CMSSW_9_4_2/src/Analysis/Tools/data/MuonTightID_B_F_weights.root");
   std::shared_ptr <JetTriggerWeight> jetTriggerweights = analysis.jetTriggerWeights("/nfs/dust/cms/user/vagneria/store/CMSSW_9_4_2/src/Analysis/Tools/data/JetTrigeff_weights.root");


   analysis.addTree<Jet>("Jets", Form("MssmHbb/Events/%s",jetsCol_.c_str()) );
   analysis.addTree<Muon>("Muons","MssmHbb/Events/slimmedMuons");

   if ( isMC_)
    {
     analysis.addTree<GenParticle>("GenParticles", "MssmHbb/Events/prunedGenParticles");
     analysis.addTree<GenJet> ("GenJets","MssmHbb/Events/slimmedGenJets");
    }

      
   for ( auto & obj : triggerObjectsJets_ )
       analysis.addTree<TriggerObject> (obj,Form("MssmHbb/Events/slimmedPatTrigger/%s",obj.c_str()));
   for ( auto & obj : triggerObjectsMuons_ )
       analysis.addTree<TriggerObject> (obj,Form("MssmHbb/Events/slimmedPatTrigger/%s",obj.c_str()));

  
   analysis.triggerResults("MssmHbb/Events/TriggerResults");
      
   if( !isMC_ ) analysis.processJsonFile(json_);
   
   //   float prescaleEra = 9.22;
   float prescaleEra = 7.0; //6.98 SR1 6.83 SR2 6.72 SR3                                                                                                                                  
   std::string sr_s = "SR";
   if ( ! signalregion_ ) sr_s = "CR";
   boost::algorithm::replace_last(outputRoot_, ".root", "_"+sr_s+".root"); 
   
   CreateHistograms();

   if ( isMC_ )
   {
      // Gen filter
      FilterResults genFilter = analysis.generatorFilter("MssmHbb/Metadata/GeneratorFilter");
      nGenTotal_ = genFilter.total;
      Efficiency_ = genFilter.efficiency;
      h1["nGen"] -> SetBinContent(1,nGenTotal_);
      h1["Eff"] -> SetBinContent(1,Efficiency_);
   }
   
   TFile * outFile = new TFile(outputRoot_.c_str(),"RECREATE");

   double mbb;
   double weight = 1;
   double mbbw[nsubsamples];

   TTree *tree = new TTree("MssmHbb_13TeV","");
   tree->Branch("mbb",&mbb,"mbb/D");
   tree->Branch("nsubsamples",&nsubsamples,"nsubsamples/I");
   tree->Branch("weight",&weight,"weight/D");

   TTree *tree0 = new TTree("MssmHbb_13TeV_0","");   
   tree0->Branch("mbb",&mbb,"mbb/D");
   tree0->Branch("weight",&weight,"weight/D");

   TTree *tree1 = new TTree("MssmHbb_13TeV_1","");   
   tree1->Branch("mbb",&mbb,"mbb/D");
   tree1->Branch("weight",&weight,"weight/D");

   TTree *tree2 = new TTree("MssmHbb_13TeV_2","");   
   tree2->Branch("mbb",&mbb,"mbb/D");
   tree2->Branch("weight",&weight,"weight/D");

   TTree *tree3 = new TTree("MssmHbb_13TeV_3","");   
   tree3->Branch("mbb",&mbb,"mbb/D");
   tree3->Branch("weight",&weight,"weight/D");

   TTree *tree4 = new TTree("MssmHbb_13TeV_4","");   
   tree4->Branch("mbb",&mbb,"mbb/D");
   tree4->Branch("weight",&weight,"weight/D");

   TTree *tree5 = new TTree("MssmHbb_13TeV_5","");   
   tree5->Branch("mbb",&mbb,"mbb/D");
   tree5->Branch("weight",&weight,"weight/D");

   TTree *tree6 = new TTree("MssmHbb_13TeV_6","");   
   tree6->Branch("mbb",&mbb,"mbb/D");
   tree6->Branch("weight",&weight,"weight/D");

   std::map<std::string, TTree *> m12_vars;

   for ( auto & sys : systematics)
     {
       for ( auto & var : vars )
	 {
	   m12_vars[("m12_"+sys+"_"+var).c_str()] = new TTree( ("MssmHbb_13TeV_"+ sys +"_" + var).c_str() ,"" );
	   m12_vars[("m12_"+sys+"_"+var).c_str()] ->Branch("mbb",&mbb,"mbb/D");
	   m12_vars[("m12_"+sys+"_"+var).c_str()] ->Branch("weight",&weight,"weight/D");
	   //m12_vars.push_back(m12_vars[("m12_"+sys+"_"+var).c_str()]);
	 }
     }

   // Analysis of events
   std::cout << "This analysis has " << analysis.size() << " events" << std::endl;
     
   int nsel[10] = { };  
   int nmatch_jet[10] = { };
   int nmatch_mu[10] = { };
   int nmatch_mujet = 0;   

   int nentries_m12 = 0;

   //   TRandom3 *R = new TRandom3( 243193272 );
   int seed = analysis.seed("seed.txt");
   TRandom3 *R = new TRandom3( seed );

 
   if ( nevtmax_ < 0 ) nevtmax_ = analysis.size();
   for ( int i = 0 ; i < nevtmax_ ; ++i )
   { 
      int njets = 0;
      int njets_btagsel = 0;
      int njets_mujsel = 0;
      int nmuons =0;
      int muj_i = 0; //muon-jet index
      int jet_i = 1; // jet index

      bool goodEvent = true;
      bool muonJetEvent = false;   

      int window = 0;

      double offbtag_weight       = 1;
      double offbtag_weight_up    = 1;
      double offbtag_weight_down  = 1;

      //double SFb_weight       = 1;
      double SFb_weight_up    = 1;
      double SFb_weight_down  = 1;
      //double SFc_weight       = 1;
      double SFc_weight_up    = 1;
      double SFc_weight_down  = 1;
      //double SFudsg_weight       = 1;
      double SFudsg_weight_up    = 1;
      double SFudsg_weight_down  = 1;

      double onlbtag_weight     = 1; //for double bjet matching
      double onlbtag_weight_up   = 1; 
      double onlbtag_weight_down = 1; 

      double PU_weight           = 1;
      double PU_weight_up        = 1;
      double PU_weight_down      = 1;

      double muID_weight           = 1;
      double muID_weight_up        = 1;
      double muID_weight_down      = 1;

      double mu_trigeff_weight           = 1;
      double mu_trigeff_weight_up        = 1;
      double mu_trigeff_weight_down      = 1;

      double jet_trigeff_weight           = 1;
      double jet_trigeff_weight_up        = 1;
      double jet_trigeff_weight_down      = 1;
      
      float genweight  = 1.;
      float evt_weight = 1.;

      if (debug) cout << i << "inside" << endl;
      if ( i > 0 && i%100000==0 ) std::cout << i << "  events processed! " << std::endl;
        
      analysis.event(i);
      if (! isMC_ )
      {
         if (!analysis.selectJson() ) continue; // To use only goodJSonFiles
       }

      // #0 TRIGGERED EVTs     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      
      if ( ! (isMC_ && signalregion_) ) // for data fire trigger first
     {
      int triggerFired = analysis.triggerResult(hltPath_);
      if ( !triggerFired ) continue;
     }
            
      if (debug)   cout <<  window << " " << prescaleEra << endl;
      ++nsel[0];

      // =============== GEN SELECTION =======================================================================

      if ( isMC_ ) 
	{
	  //Gen weight
	  genweight = analysis.genWeight()/fabs(analysis.genWeight());
	  evt_weight *= genweight;	  

	  h1["nentries"]->Fill((genweight+1.)/2.);
          h1["nentries_tot"]->Fill(1,evt_weight);

	  PU_weight         = puweights->weight(analysis.nTruePileup(),0);
	  PU_weight_up      = puweights->weight(analysis.nTruePileup(),2);
	  PU_weight_down    = puweights->weight(analysis.nTruePileup(),-2);

	  evt_weight *= PU_weight;	  
	}

      //  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++      
  
      std::vector<GenParticle *> higgsDaughters;
      std::vector<GenJet *>      selectedGenJets;  
      GenParticle * Higgs = new GenParticle;

      if (isMC_ )
     {
      auto genParticles = analysis.collection<GenParticle>("GenParticles");
      auto genJets = analysis.collection<GenJet>("GenJets");
                                                                                                                                                     \
      //Select good partons
      std::vector<GenParticle *> selectedGenParticles;
      std::vector<GenParticle *> bSpectators;
     
      for ( int j = 0 ; j < genParticles->size() ; ++j )
	{
	  GenParticle * gp = & genParticles->at(j) ;
	  if ( gp->pdgId() == 36 || gp->pdgId() == 25 ) Higgs = gp;              //Higgs particle 
	  //	  if ( gp->pt() < 20 || fabs(gp->eta()) > 2.3 )         continue;        //kinematic selection
	  if ( abs(gp->pdgId()) > 5 && abs(gp->pdgId()) != 21 ) continue;        // include only quark and gluon	   
	  selectedGenParticles.push_back(gp);	  
	  if ( abs(gp->pdgId()) == 5 && gp->status() == 23 ) 
	    bSpectators.push_back(gp);
	}

      if (debug) cout  << "selectedGenParticles "<< selectedGenParticles.size() << endl;

      //The last two entries are ALWAYS the Higgs daughters
      higgsDaughters.push_back(bSpectators.back());
      bSpectators.pop_back();
      higgsDaughters.push_back(bSpectators.back());
      bSpectators.pop_back();

      //Selected good genjets
      for ( int j = 0 ; j < genJets->size() ; ++j )
	{
	  GenJet * gj = & genJets->at(j) ;
	  if ( gj->pt() < 30 || fabs(gj->eta()) > 2.3 )         continue;        //kinematic selection
          selectedGenJets.push_back(gj);
	}      

      if (debug) cout  << "selectedJets  "<< selectedGenJets.size() << endl;

      //Before cascade
      //      std::vector<GenParticle *> bSpectators;
      // int n_bSpectators =0;

      /*
      for ( size_t j = 0 ; j < selectedGenParticles.size() ; ++j )
	{
           GenParticle * gp = selectedGenParticles[j];	 
	   const GenParticle & higgs = *Higgs;

	   //	   if ( gp->status() == 23 && abs(gp->pdgId()) == 5 ) //bquarks before ISR (bassociate & Higgs dau)
	   if ( abs(gp->pdgId()) == 5 ) //bquarks before ISR (bassociate & Higgs dau)                                                     
	   {
	       if ( gp->higgsDaughter() || n_bSpectators == 2 ) 
		 {
		   higgsDaughters.push_back(gp); // in LO 
		   if (debug) cout << "higgsDaughters  pT" << gp->pt() << "deltaR(H,d) =" << gp->deltaR(higgs) <<  endl;
		   
		 }
	       //	       
	       
	       else if ( n_bSpectators < 2 )
		 {
		   bSpectators.push_back(gp);		
		   if (debug) cout << "bSpectator " << n_bSpectators << " pT" << gp->pt() << "deltaR(H,a) =" << gp->deltaR(higgs) << endl;
		   ++n_bSpectators;
		   }
	     }
 	}
*/
      if (debug)
	{
	  cout  << "bSpectators "   << bSpectators.size() << endl;
	  cout  << "higgsDaughters "<< higgsDaughters.size() << endl;
	}

      std::sort(higgsDaughters.begin(),higgsDaughters.end(),greaterByPt<GenParticle *>());
      std::sort(   bSpectators.begin(),   bSpectators.end(),greaterByPt<GenParticle *>());

      //      GenParticle * bb = new GenParticle;
      if (higgsDaughters.size() == 2 &&  bSpectators.size()>0 )
	{
	  const GenParticle & h = *Higgs;
	  //	  bb->p4( higgsDaughters[0]->p4() + higgsDaughters[1]->p4());
	  //	  cout << "higgs cand pT= " << Higgs->pt() << " form di-bgen pT= "  << bb->pt() << endl; 
	  if (debug) cout << "higgs pt " << Higgs->pt() << endl;
	  h1["genpt_Higgs"]->Fill(Higgs->pt(), evt_weight);
	  h1["geneta_Higgs"]->Fill(Higgs->eta(), evt_weight);

	  h1["genpt_Higgsdau0"]-> Fill(higgsDaughters[0]->pt(), evt_weight);
          h1["genpt_Higgsdau1"]-> Fill(higgsDaughters[1]->pt(), evt_weight);
          h1["geneta_Higgsdau0"]-> Fill(higgsDaughters[0]->eta(), evt_weight);
          h1["geneta_Higgsdau1"]-> Fill(higgsDaughters[1]->eta(), evt_weight);

	  h1["genpt_bSpectator"] -> Fill( bSpectators[0]->pt(), evt_weight);
          h1["geneta_bSpectator"] -> Fill( bSpectators[0]->eta(), evt_weight);
	  /*
	  h1["deltaR_Ha1"]->Fill( bSpectators[0]->deltaR(h), evt_weight);
	  h1["deltaR_Hd1"]->Fill( higgsDaughters[0]->deltaR(h), evt_weight);
	  h1["deltaR_Hd2"]->Fill( higgsDaughters[1]->deltaR(h), evt_weight);
	  */
	}

      //      cout << "================================="<< endl;

      //      std::sort(higgsDaughters.begin(),higgsDaughters.end(),ptSort<GenParticle>);
      // std::sort(bSpectators.begin(),bSpectators.end(),ptSort<GenParticle>);
      //      std::sort(higgsDaughters.begin(),higgsDaughters.end(),greaterByPt<GenParticle *>());

      //Clearing vectors 
      selectedGenParticles.clear();
      bSpectators.clear();

     }                                                                                                               
      if (debug) cout << "particle jet selection" << endl;

      // =============== JET SELECTION =======================================================================

      // #1 ID LOOSE JET  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++             
     
      auto slimmedJets = analysis.collection<Jet>("Jets");
      auto slimmedMuons = analysis.collection<Muon>("Muons");
        
      if ( isMC_ )
	{
	  auto genJets = analysis.collection<GenJet>("GenJets");
	  slimmedJets->addGenJets(genJets);  // if not defined -> jerMatch = false, smearing will be done using the stochastic method only
	 }
      
      std::vector<Jet *> selectedJets;

      for ( int j = 0 ; j < slimmedJets->size() ; ++j )
      {

	Jet * jet = &slimmedJets->at(j) ;
	if ( !jet->pileupJetIdFullId("loose")) continue;  //PU ID
	if ( !jet->idTight() ) continue;                  //Jet ID

	if ( isMC_ ) 
	  { 
	    //JER
	    jet->jerInfo(*jerinfo,0.2); // this also performs matching to the added gen jets above, with delta R < 0.2 which is default and can be omitted
	 
	    /*	    if ( j < 3  && i%1000==0 )
	      {

		//	    std::cout << "diff get pT = " << gj->pt()-jet->pt() << endl;
	  
        	//GenJet * gj = & genJets->at(j) ;
	    //	    std::cout << "diff get pT = " << gj->pt()-jet->pt() << endl;
	    ///            std::cout << "    Jet #" << j << ": ";
		
	    std::cout << "    Jet #" << j << ": ";
	    std::cout << "diff get pT = " << gj->pt()-jet->pt() << endl;
	    std::cout << "pT  = "     << jet->pt()      << ", ";
	    std::cout << "eta = "     << jet->eta()     << ", ";
	    std::cout << "resolution = " << jet->jerPtResolution() << ", ";
	    std::cout << "JER SF  = "    << jet->jerSF()  << ", match = " << jet->jerMatch() << " jer corr = " << jet->jerCorrection() << " + ";
	    std::cout << jet->jerCorrection("up") << " - " << jet->jerCorrection("down") << std::endl; 
	    }*/
	    
	    if      ( jervar_ == "up"  ) correctJetpt( *jet, jet->jerCorrection("up") ); 
	    else if ( jervar_ == "down") correctJetpt( *jet, jet->jerCorrection("down") ); 
	    else    correctJetpt( *jet, jet->jerCorrection() );  //for variations jerCorrection("up","down")
	    
    	    //JES
	    if ( jecsigma_ != 0 ) correctJetpt( *jet, 1 + jecsigma_*jet->jecUncert()  );
	  }
  

	//Bjet energy regression Correction 
      	 float btagj = 0.; 
         if ( btagalgo == "deepcsv" )  btagj = jet->btag("btag_deepb")+jet->btag("btag_deepbb");
	 else if ( btagalgo == "deepflavour" )  btagj = jet->btag("btag_dfb")+jet->btag("btag_dfbb")+jet->btag("btag_dflepb");
	 if ( btagj > btagmin ) correctJetpt( *jet, jet->bRegCorr() );

	selectedJets.push_back(jet);

      }

      if ( (int)selectedJets.size() < njetsmin_ ) continue;
      
      /*      if ( isMC_ ) 
	{ 
	    //JER
	    //CHECK before kinemtic cuts JER & JES

	    //JER 
	    std::vector<std::string> vars = { "up", "down" }; // last entry to restore JER central value

	    Jet * j0 = selectedJets[0];
	    Jet * j1 = selectedJets[1];
	    Jet * j0_beforeVars = new Jet(j0->pt(), j0->eta(), j0->phi(), j0->e()); 
	    Jet * j1_beforeVars = new Jet(j1->pt(), j1->eta(), j1->phi(), j1->e()); 

	    //to cross-check only
	    Jet * j0_JES_up   = new Jet(j0->pt(), j0->eta(), j0->phi(), j0->e()); 
	    Jet * j0_JES_down = new Jet(j0->pt(), j0->eta(), j0->phi(), j0->e()); 
	    Jet * j1_JES_up   = new Jet(j1->pt(), j1->eta(), j1->phi(), j1->e()); 
	    Jet * j1_JES_down = new Jet(j1->pt(), j1->eta(), j1->phi(), j1->e()); 
	    
	    Jet * j0_JER_up   = new Jet(j0->pt(), j0->eta(), j0->phi(), j0->e()); 
	    Jet * j0_JER_down = new Jet(j0->pt(), j0->eta(), j0->phi(), j0->e()); 
	    Jet * j1_JER_up   = new Jet(j1->pt(), j1->eta(), j1->phi(), j1->e()); 
	    Jet * j1_JER_down = new Jet(j1->pt(), j1->eta(), j1->phi(), j1->e()); 

	  
	    for ( auto & var : vars )
	      {
		correctJetpt( *j0 , j0->jerCorrection(var.c_str())/ j0->jerCorrection() );
		correctJetpt( *j1 , j1->jerCorrection(var.c_str())/ j1->jerCorrection() );

		h1[Form("m12_mujsel_JER_%s_befkin",var.c_str())] -> Fill( (j0->p4() + j1->p4()).M() , evt_weight);
		h1[Form("pt_0_JER_%s_befkin",var.c_str())] -> Fill( j0->pt() , evt_weight);
		h1[Form("pt_1_JER_%s_befkin",var.c_str())] -> Fill( j1->pt() , evt_weight);
	
		//restore pT before variation
		if ( var =="up" )
		  { 
		    copyJetpt( *j0_JER_up, *j0);
		    copyJetpt( *j1_JER_up, *j1);
		  }

		if ( var =="down" ) 
		  {
		    copyJetpt( *j0_JER_down, *j0);
		    copyJetpt( *j1_JER_down, *j1);
		  }

		copyJetpt( *j0, *j0_beforeVars);
		copyJetpt( *j1, *j1_beforeVars);
	
	      }
	    	
	    //JES
	    std::vector<std::pair<string,int>> sigmas = { {"up",2}, {"down",-2} };

	    for ( auto & var : sigmas)
	      {
		correctJetpt( *j0, 1 + var.second*j0->jecUncert()  );
		correctJetpt( *j1, 1 + var.second*j1->jecUncert()  );
	      
	      h1[Form("m12_mujsel_JES_%s_befkin",var.first.c_str())] -> Fill( (j0->p4() + j1->p4()).M() , evt_weight); 
	      h1[Form("pt_0_JES_%s_befkin",var.first.c_str())] -> Fill( j0->pt() , evt_weight);
	      h1[Form("pt_1_JES_%s_befkin",var.first.c_str())] -> Fill( j1->pt() , evt_weight);
	
	      //restore pT before variation
	      if ( var.first =="up" )
		{ 
		  copyJetpt( *j0_JES_up, *j0);
		  copyJetpt( *j1_JES_up, *j1);
		}

	      if ( var.first =="down" ) 
		{
		  copyJetpt( *j0_JES_down, *j0);
		  copyJetpt( *j1_JES_down, *j1);
		  }

	      copyJetpt( *j0, *j0_beforeVars);
	      copyJetpt( *j1, *j1_beforeVars);     

	      }

	    h1["j0_JES_diff_befkin"]->Fill( j0_JES_up->pt()-j0_JES_down->pt() );
	    h1["j0_JER_diff_befkin"]->Fill( j0_JER_up->pt()-j0_JER_down->pt() );
	    h1["j1_JES_diff_befkin"]->Fill( j1_JES_up->pt()-j1_JES_down->pt() );
	    h1["j1_JER_diff_befkin"]->Fill( j1_JER_up->pt()-j1_JER_down->pt() );

	    h1["m12_JES_diff_befkin"]->Fill( (j0_JES_up->p4() + j1_JES_up->p4()).M() - (j0_JES_down->p4() + j1_JES_down->p4()).M()  );
	    h1["m12_JER_diff_befkin"]->Fill( (j0_JER_up->p4() + j1_JER_up->p4()).M() - (j0_JER_down->p4() + j1_JER_down->p4()).M()  );

	    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	    } //END variation */

      if (debug) cout << "correction jet selection" << endl;

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

      if (debug) cout << "kin jet selection" << endl;

      if (isMC_)  h1["nentries_kin"]->Fill(1,evt_weight);

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

      if (isMC_ ) h1["jet2_flavour"] ->Fill(selectedJets[2]->flavour(), evt_weight);

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
	 h1[Form("deepflavour_%i",j)] -> Fill(jet->btag("btag_dfb")+jet->btag("btag_dfbb")+jet->btag("btag_dflepb"));
	 h1[Form("deepcsv_%i",j)] -> Fill(jet->btag("btag_deepb")+jet->btag("btag_deepbb"));
	 h2[Form("eta_phi_%i",j)] -> Fill(jet->eta(),jet->phi());

	 //btag selection
      	 float btagj; 
 //      	 if      ( btagalgo == "csv" )      btagj = jet->btag();
         if ( btagalgo == "deepcsv" )  btagj = jet->btag("btag_deepb")+jet->btag("btag_deepbb");
	 else if ( btagalgo == "deepflavour" )  btagj = jet->btag("btag_dfb")+jet->btag("btag_dfbb")+jet->btag("btag_dflepb");

     	 else return -1;
	
	 if ( j < 2 && btagj < btagmin )     goodEvent = false;

	 /*	 if ( ! signalregion_ )
	   {
	     if ( j == 2 && btagj > nonbtagwp_ ) goodEvent = false; 
	   }
	 else
	   {
	     if ( j == 2 && btagj < btagmin ) goodEvent = false; 
	     } */  
	 
      }

      h1["m12"] -> Fill((selectedJets[0]->p4() + selectedJets[1]->p4()).M());
      
      if ( ! goodEvent ) continue;
     
      //cout << i << "jet btag sel" << endl;
      if (debug) cout << "b-jet selection" << endl;

     //+++++++++++++++++++++++++++++++++++++++++++
      // Fill histograms of B-JET SELECTION but no MATCHING
    
      //Apply offline Btag 
      if (isMC_)
	{

	  //	  h1["offbtag_sf_j0"] -> Fill(selectedJets[0]->btagSF (bsf_reader) );
	  //	  h1["offbtag_sf_j1"] -> Fill(selectedJets[1]->btagSF (bsf_reader) );
	  //	  h1["offbtag_sf_j0overj1"] -> Fill(selectedJets[0]->btagSF (bsf_reader)/selectedJets[1]->btagSF (bsf_reader) );

	  for (int j = 0 ; j < njetsmin_ ; ++j )
	    {
	      //    if (selectedJets[j]->flavour() == 5) offbtag_weight *= reader.eval(BTagEntry::FLAV_B,    selectedJets[j]->eta(), selectedJets[j]->pt(), 0);     	
	      offbtag_weight      *= selectedJets[j]->btagSF    (bsf_reader); //central  btagSFup method
	      offbtag_weight_up   *= selectedJets[j]->btagSFup  (bsf_reader);
	      offbtag_weight_down *= selectedJets[j]->btagSFdown(bsf_reader);

	      if (selectedJets[j]->flavour() == 5) 
		{
		  //		  SFb_weight      *= selectedJets[j]->btagSF    (bsf_reader); //central  btagSFup method
		  SFb_weight_up   *= selectedJets[j]->btagSFup  (bsf_reader);
		  SFb_weight_down *= selectedJets[j]->btagSFdown(bsf_reader);
		}

	      if (selectedJets[j]->flavour() == 4) 
		{
		  //		  SFc_weight      *= selectedJets[j]->btagSF    (bsf_reader); //central  btagSFup method
		  SFc_weight_up   *= selectedJets[j]->btagSFup  (bsf_reader);
		  SFc_weight_down *= selectedJets[j]->btagSFdown(bsf_reader);
		}

	      if (selectedJets[j]->flavour() != 4 && selectedJets[j]->flavour() != 5 ) 
		{
		  //		  SFudsg_weight      *= selectedJets[j]->btagSF    (bsf_reader); //central  btagSFup method
		  SFudsg_weight_up   *= selectedJets[j]->btagSFup  (bsf_reader);
		  SFudsg_weight_down *= selectedJets[j]->btagSFdown(bsf_reader);
		}

	    }
	  
	  evt_weight*= offbtag_weight;
	}


      //Fill histo
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
         h1[Form("deepflavour_%i_btagsel",j)] -> Fill(jet->btag("btag_dfb")+jet->btag("btag_dfbb")+jet->btag("btag_dflepb"));
	 h1[Form("deepcsv_%i_btagsel",j)] -> Fill(jet->btag("btag_deepb")+jet->btag("btag_deepbb"));
	 h2[Form("eta_phi_%i_btagsel",j)] -> Fill(jet->eta(),jet->phi());

      }

      float mbb_btagsel = (selectedJets[0]->p4() + selectedJets[1]->p4()).M();
      if ( !signalregion_ || isMC_ )
      { 
         h1["m12_btagsel"] -> Fill(mbb_btagsel);
         // weight = 1;
        // tree -> Fill();
      }
           
      ++nsel[5];

      int nTrueDaughters = 0;     
     if (isMC_ )
    {
      //      int nTrueDaughters = 0;
      for ( int j = 0; j < 2; ++j ) //dijet loop
	{
	  Jet & bjet = *selectedJets[j];
	  //	  cout << "passing_selection" << endl;

	  for ( size_t d =0; d < higgsDaughters.size(); d++ )
	    {
	      if ( bjet.deltaR(*higgsDaughters[d]) < 0.3 ) // 0.3 
		{ 
		  if (debug)
		    {
		      cout << "bjet " << j << " deltaR =" <<  bjet.deltaR(*higgsDaughters[d]) << endl;
		      cout << "bjet pt,eta,phi " << bjet.pt() << " " <<bjet.eta() << " " <<bjet.phi() << "Higgs daughter pt,eta,phi" << higgsDaughters[d]->pt()  << " " << higgsDaughters[d]->eta()  << " " << higgsDaughters[d]->phi()  << endl;
		    }
		  ++nTrueDaughters;
		  break;
		}
	    }
	}

      }

     if (debug) cout << "regression correction" << endl;
    
     
     //FSR correction - adding soft jets
       for ( size_t s = 3; s < selectedJets.size() ; ++s )  //soft jet loop - from 4th jet
            {
	      Jet & softjet = *selectedJets[s];
	      
	      //Cut on DeltaR softjet - b jet
	      float dRminsoft_bj = std::min({softjet.deltaR(*selectedJets[0]),softjet.deltaR(*selectedJets[1]),softjet.deltaR(*selectedJets[2])});
	      
	      if ( softjet.pt() < 20.0 ) continue;
	      if ( dRminsoft_bj > 0.8 )  continue;
	      //	      if ( softjet.qgLikelihood() > 0.5 ) continue;
	      //	      if (softjet.btag("btag_deepb")+softjet.btag("btag_deepbb") > 0.1522 )  continue;
	      for ( int j = 0; j < 3; ++j ) //dijet loop
		{
		  Jet & bjet = *selectedJets[j];
		  
		  if ( dRminsoft_bj != softjet.deltaR(bjet) ) continue;
	  
		  bjet.p4( bjet.p4() + softjet.p4() );
		  
		  //Remove bjet and soft jet and add corrected bjet
		  selectedJets.erase(selectedJets.begin()+j);
		  selectedJets.insert(selectedJets.begin()+j, &bjet );

		  h1["pt_softjet"] -> Fill(softjet.pt());                                                                                                               
                  h1["dR_jsoftjet"] -> Fill(softjet.deltaR(bjet));                                                                                                                        h1["rank_softjet"] -> Fill(s);                               
		}
    
		}
     
      if (debug) cout << "fsr correction" << endl;

      //    cout << "soft corrections " << isJetCorrected[0]  << " " << isJetCorrected[1] << " " << isJetCorrected[2] << endl;
      /*bool rerank_j1toj0 = false;
      bool rerank_j2 = false;
      
      if ( isJetCorrected[1] && selectedJets[1]->pt() > selectedJets[0]->pt()  ) rerank_j1toj0 = true;
      if ( isJetCorrected[2] && selectedJets[2]->pt() > selectedJets[1]->pt()  ) rerank_j2 = true; 
      
      if ( rerank_j1toj0 ) cout << i << " j0 pT " << selectedJets[0]->pt() <<" j1 pT " << selectedJets[1]->pt() << endl;
      if ( rerank_j2 )     cout << i << " j0 pT " << selectedJets[0]->pt() <<" j1 pT " << selectedJets[1]->pt() << " j2 pT " << selectedJets[2]->pt() << endl; 
      
      std::sort(selectedJets.begin(),selectedJets.end(),greaterByPt<Jet *>());
      
      if ( rerank_j1toj0 ) cout << i << " j0 pT " << selectedJets[0]->pt() << " j1 pT " << selectedJets[1]->pt() << endl; 
      if ( rerank_j2 )     cout << i << " j0 pT " << selectedJets[0]->pt() << " j1 pT " << selectedJets[1]->pt() << " j2 pT " << selectedJets[2]->pt() << endl;
      */

      // cout <<  "after correction m12 =" <<  (selectedJets[0]->p4() + selectedJets[1]->p4()).M() << endl; ;


      //========== MUON SELECTION ============================================================================  
      
      //Muon ID sel +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      std::vector<Muon *> selectedMuons;
      const char muID = muonsid_.at(0) ;
      for ( int m = 0 ; m < slimmedMuons->size() ; ++m )
      {	 

	//       if ( i > 0 && i%100==0 ) std::cout << i << "\t" << "slimmed" << "\t" <<slimmedMuons->at(m).pt() << std::endl;
	if ( (muID =='T' && slimmedMuons->at(m).isTightMuon()) || (muID =='M' && slimmedMuons->at(m).isMediumMuon())  )
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

	      // Muon jet vs non muon jet
	      h1["dR_muj"] -> Fill( (dR_muj0 < dR_muj1) ? dR_muj0 : dR_muj1 );
	      if       (dR_muj0 < drmax_)  h1["dR_muj0"] -> Fill( dR_muj0 );
	      else if  (dR_muj1 < drmax_) { 
		                            muj_i = 1;
	                               	    jet_i = 0;
	                       		    h1["dR_muj1"] -> Fill( dR_muj1 );
	                                  }

	      h1["pt_muj"] -> Fill( selectedJets[muj_i]->pt() );
	      h1["pt_j"] -> Fill( selectedJets[jet_i]->pt() );
	      h1["eta_muj"] -> Fill( selectedJets[muj_i]->eta() );
              h1["eta_j"] -> Fill( selectedJets[jet_i]->eta() );

      //	      if ( i > 0 && i%100==0 ) std::cout << i <<"\t" << " dR cut" << "\t" << muon->pt() << std::endl;
	      selectedMuonsinJet.push_back(muon);
	      break;
	      } 
       }     

      //Asymmetric pt cut
      if ( selectedJets[jet_i]->pt() < 50.0 ) continue;
      if ( selectedJets[muj_i]->pt() < 60.0 ) continue;
      
      if ( ! muonJetEvent ) continue;

      if (debug) cout << "muon selection" << endl;

      ++nsel[6];       
      // #6 MATCHING OFFLINE to ONLINE  +++++++++++++++++++++++++++++++++++++++++++++++++++++
   
      if ( isMC_ && signalregion_ ) // for MC apply trigger selection at the end                                                                                                       
	{
	  int triggerFired = analysis.triggerResult(hltPath_);
	  if ( !triggerFired ) continue;
	  ++nsel[7];
	}

    
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
	    //if ( io ==2 && i%50000  ) std::cout << "event # " << i << " j(" <<  j << ") matched " << matched_jet[io] <<  std::endl;
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
   
      //Matched BJETS and rank -- no muon requirement 

      std::vector<int> matched_jet_i;
      for ( size_t j = 0; j < selectedJets.size(); ++j )
	{
	  Jet * jet = selectedJets[j];
          bool matchedSelJet = true;
	  for ( size_t io = 0; io < triggerObjectsJets_.size()-1 ; ++io ) //exclude hltBSoftMuonDiJet40Mu12L3FilterByDR                                                                                      
	    {
	      if ( ! jet->matched(triggerObjectsJets_[io]) ) 
		{
		  matchedSelJet = false;
	          break;
		}
	    }

	  if ( matchedSelJet ) matched_jet_i.push_back(j);
	}

      if (debug) cout << "tob matching " << endl;

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

      //Jets matched to trigger obj
      h1["n_jetsmatched_mujsel"] -> Fill(matched_jet_i.size());     // # real mu-jets per event                                                                                  
      if (!goodEvent || !(matched_mujet[0] || matched_mujet[1]) ) continue;
    
      h1["n_jetsmatched_mujmatched"] -> Fill(matched_jet_i.size());     // # real mu-jets per event                                                                                                           

      //Trigger matching
      if ( ! (isMC_ && signalregion_) ) ++nsel[7]; 
      else ++nsel[8];


      Jet * thirdjet = selectedJets[2];
      
      //3rg btag jet selection
      float btag3j = 0; 
      if ( btagalgo == "deepcsv" )  btag3j = thirdjet->btag("btag_deepb")+thirdjet->btag("btag_deepbb");
      else if ( btagalgo == "deepflavour" )  btag3j = thirdjet->btag("btag_dfb")+thirdjet->btag("btag_dfbb")+thirdjet->btag("btag_dflepb");

      if ( ! signalregion_ )
	{
	  if ( btag3j > nonbtagwp_ ) goodEvent = false; 
	}

      else
	{
	  if ( btag3j < btagmin ) goodEvent = false; 
	}   

      if ( ! goodEvent ) continue;
   
      if ( ! (isMC_ && signalregion_) ) ++nsel[8];
      ++nsel[9]; 
         



      //MU ID weight for the muon in the jet
      //ONLINE BTAG WEIGHT
      //Mu trig eff weight

      if (isMC_) 
	{
	  //MuID
	  if ( selectedMuonsinJet[0]->pt() < 120. ) 
	    {
	      muID_weight         = muIDweights->weight(selectedMuonsinJet[0]->pt(), fabs(selectedMuonsinJet[0]->eta()), 0);
	      muID_weight_up      = muIDweights->weight(selectedMuonsinJet[0]->pt(), fabs(selectedMuonsinJet[0]->eta()), 2);
	      muID_weight_down    = muIDweights->weight(selectedMuonsinJet[0]->pt(), fabs(selectedMuonsinJet[0]->eta()),-2);
	      }
	    
	  // Muon pT trig eff 
	  mu_trigeff_weight         = muonptWeight(selectedMuonsinJet[0]->pt(), 0);
	  mu_trigeff_weight_up      = muonptWeight(selectedMuonsinJet[0]->pt(), 2);
	  mu_trigeff_weight_down    = muonptWeight(selectedMuonsinJet[0]->pt(),-2);

	  //onlbtag
	  onlbtag_weight      = onlSFbtag( selectedJets[0]->pt(),  0 )* onlSFbtag( selectedJets[1]->pt(),  0 ); //for double bjet matching
	  onlbtag_weight_up   = onlSFbtag( selectedJets[0]->pt(),  2 )* onlSFbtag( selectedJets[1]->pt(),  2 ); 
	  onlbtag_weight_down = onlSFbtag( selectedJets[0]->pt(), -2 )* onlSFbtag( selectedJets[1]->pt(), -2 ); 


	  //Jet trig eff for dijet matching

	  //BEFORE scale factor BINNED  
	  //jet_trigeff_weight  = jetTriggerweights->weight(selectedJets[0]->pt(), fabs(selectedJets[0]->eta()), 0);

	  jet_trigeff_weight         = jetTriggerScaleFactor(selectedJets[0]->pt(), fabs(selectedJets[0]->eta()), 0);
	  jet_trigeff_weight_up      = jetTriggerScaleFactor(selectedJets[0]->pt(), fabs(selectedJets[0]->eta()), 2);
	  jet_trigeff_weight_down    = jetTriggerScaleFactor(selectedJets[0]->pt(), fabs(selectedJets[0]->eta()),-2);

	  jet_trigeff_weight         *= jetTriggerScaleFactor(selectedJets[1]->pt(), fabs(selectedJets[1]->eta()), 0);
	  jet_trigeff_weight_up      *= jetTriggerScaleFactor(selectedJets[1]->pt(), fabs(selectedJets[1]->eta()), 2);
	  jet_trigeff_weight_down    *= jetTriggerScaleFactor(selectedJets[1]->pt(), fabs(selectedJets[1]->eta()),-2);

	  evt_weight*= onlbtag_weight*muID_weight*mu_trigeff_weight*jet_trigeff_weight;
	  //	   evt_weight*= onlbtag_weight*mu_trigeff_weight;

	  /*
	  float mbb_temp = (selectedJets[0]->p4() + selectedJets[1]->p4()).M();

	  evt_weight*= onlbtag_weight;
	  //	  cout << "----------------------------------------------" << endl;

	  //cout << "pT(jet0)= " << selectedJets[0]->pt() << " pt(jet1)= " << selectedJets[1]->pt() << "evt weight = " << evt_weight << endl;
	  //cout << "eta(jet0)= " << selectedJets[0]->eta() << " eta(jet1)= " << selectedJets[1]->eta() << "evt weight = " << evt_weight << endl;
	  //cout << "pT(muon)= " << selectedMuonsinJet[0]->pt() << " eta(muon)= " << selectedMuonsinJet[0]->eta() << endl;
	  
	  h1["pt_0_uncorr"]   -> Fill(selectedJets[0]->pt(), evt_weight);
	  h1["pt_1_uncorr"]   -> Fill(selectedJets[1]->pt(), evt_weight);
	  h1["m12_uncorr"]   -> Fill(mbb_temp, evt_weight);
	  cout << "before mass cut m12=" << mbb_temp << "weight = " << evt_weight << endl;

	  evt_weight*=mu_trigeff_weight;
	  h1["pt_0_mu_trigeff"]   -> Fill(selectedJets[0]->pt(), evt_weight);
	  h1["pt_1_mu_trigeff"]   -> Fill(selectedJets[1]->pt(), evt_weight);
	  h1["m12_mu_trigeff"]   -> Fill(mbb_temp, evt_weight);

	  //	  cout << "mueff w =" << mu_trigeff_weight << " " <<  evt_weight << endl;
	  evt_weight/=mu_trigeff_weight;
	  
	  evt_weight*=jet_trigeff_weight;
	  h1["pt_0_jet_trigeff"]   -> Fill(selectedJets[0]->pt(), evt_weight);
	  h1["pt_1_jet_trigeff"]   -> Fill(selectedJets[1]->pt(), evt_weight);
	  h1["m12_jet_trigeff"]   -> Fill(mbb_temp, evt_weight);
	  //          cout << "trigeff w =" << jet_trigeff_weight << " " <<  evt_weight << endl;
	  evt_weight/=jet_trigeff_weight;
	  
	  evt_weight*=muID_weight;
	  h1["pt_0_muID"]   -> Fill(selectedJets[0]->pt(), evt_weight);
	  h1["pt_1_muID"]   -> Fill(selectedJets[1]->pt(), evt_weight);
	  h1["m12_muID"]   -> Fill(mbb_temp, evt_weight);

	  //          cout << "muID w =" << muID_weight << " " <<  evt_weight << endl;
	  evt_weight/=muID_weight;
	  
	  */
	  //	  cout << selectedJets[0]->pt() << " " << selectedJets[1]->pt() << endl;

 
	}
   
     //=============================================================================================== 

      //Mass cut
      float mbb_mujsel = (selectedJets[0]->p4() + selectedJets[1]->p4()).M(); // set one mass window in TTree to fill

      // Mass cut
      if ( mbb_mujsel < 110.) continue;
      ++nentries_m12;

      //      if ( signalregion_ )
      //	{
	  if ( mbb_mujsel < 400 )                           h1["n_m12_sr1"] -> Fill (evt_weight);
	  if ( mbb_mujsel > 240  && mbb_mujsel < 575 )      h1["n_m12_sr2"] -> Fill (evt_weight);
	  if ( mbb_mujsel > 375  && mbb_mujsel < 800 )      h1["n_m12_sr3"] -> Fill (evt_weight);	  
//	}


//-----------------------------------------------------------------------------------------------------------------------------------------------
       // Fill histograms of passed MUJETS  selection
       for ( int j = 0 ; j < (int)selectedJets.size() ; ++j )
       {
          if ( selectedJets[j]->pt() < 20. ) continue;
          ++njets_mujsel;                             
       }
 

      h1["n_mujets"]  -> Fill(nmujets, evt_weight);                   // number real mujets per evt
      h1["n_muinjet"] -> Fill(selectedMuonsinJet.size(), evt_weight); // number of muons in dR cone
      h1["n_mujsel"]  -> Fill(selectedJets.size(), evt_weight);      // all jets but leading 3 passing bbmunb selection
      h1["n_ptmin20_mujsel"] -> Fill(njets_mujsel, evt_weight);      //all jets pt > 20 but butleading 3 passing bbmunb selection 
  
      for ( int j = 0; j < njetsmin_; ++j )
      {
         Jet * jet = selectedJets[j];
         h1[Form("pt_%i_mujsel",j)]   -> Fill(jet->pt(), evt_weight);
         h1[Form("eta_%i_mujsel",j)]  -> Fill(jet->eta(), evt_weight);
         h1[Form("phi_%i_mujsel",j)]  -> Fill(jet->phi(), evt_weight);
         h1[Form("deepflavour_%i_mujsel",j)] -> Fill(jet->btag("btag_dfb")+jet->btag("btag_dfbb")+jet->btag("btag_dflepb"), evt_weight);
	 h1[Form("deepcsv_%i_mujsel",j)] -> Fill(jet->btag("btag_deepb")+jet->btag("btag_deepbb"), evt_weight);
      }

      //PT muon categorization
      for (  size_t m = 0; m < selectedMuonsinJet.size(); ++m )      
       {
         Muon * muon = selectedMuonsinJet[m];
         h1["pt_mumatched"]  -> Fill( muon->pt() , evt_weight);
         h1["eta_mumatched"] -> Fill( muon->eta() , evt_weight);
       } 

      h1["pt_muj_mujsel"]  -> Fill( selectedJets[muj_i]->pt() , evt_weight);
      h1["pt_j_mujsel"]    -> Fill( selectedJets[jet_i]->pt() , evt_weight);
      h1["eta_muj_mujsel"] -> Fill( selectedJets[muj_i]->eta() , evt_weight);
      h1["eta_j_mujsel"]   -> Fill( selectedJets[jet_i]->eta() , evt_weight);
 

      //Prescale after mass cut
      if ( ! isMC_ && ! signalregion_ )  applyPrescale ( analysis.run(), R->Rndm(), prescaleEra, nsubsamples, window  );

      //      cout << "nentries = " << nentries_m12 << endl ;
      //cout << "window = " << window << endl ;

      //Initialise tree with different windows
      for ( int i=0; i< nsubsamples; i++)
	{
	  if (i != window) mbbw[i] = -1; // set other mass windows to -1 ( fixed nentries in Tree )
	  else mbbw[window] = mbb_mujsel;
	}

      //AFTER UNBLINDING
      h1["m12_mujsel"]  -> Fill(mbb_mujsel, evt_weight);  // unprescaled	 
      weight = evt_weight;
      mbb = mbb_mujsel; // tree with unprescaled  
      tree -> Fill();    


      if ( !signalregion_ || isMC_ )
      { 

	//	cout << "after mass cut m12=" << mbb_mujsel << "weight = " << evt_weight << endl;

	// 	 h1["m12_mujsel"]  -> Fill(mbb_mujsel, evt_weight);  // unprescaled	 
	// weight = evt_weight;
	// mbb = mbb_mujsel; // tree with unprescaled  
        // tree -> Fill();    
     
	 if     ( window ==0 ) tree0->Fill();
	 else if( window ==1 ) tree1->Fill();
         else if( window ==2 ) tree2->Fill();
         else if( window ==3 ) tree3->Fill();
         else if( window ==4 ) tree4->Fill();
         else if( window ==5 ) tree5->Fill();
         else if( window ==6 ) tree6->Fill();

	 if ( !isMC_ && ! signalregion_) h1[Form("m12_mujsel_%i",window)] -> Fill(mbbw[window], evt_weight);

      }

  
      if (isMC_)
	{
	  for ( int d =0 ; d <3; d++ )
	    if ( nTrueDaughters == d ) h1[Form("m12_mujsel_%id",d)] -> Fill(mbbw[window], evt_weight );
	  
	  // PU weight
	  double evt_weight_PU_up ;
	  double evt_weight_PU_down ;
	  
	  if ( PU_weight != 0 && evt_weight != 0 ) 
	    {
	      evt_weight_PU_up   = evt_weight / PU_weight * PU_weight_up   ; 
	      evt_weight_PU_down = evt_weight / PU_weight * PU_weight_down ;
	    }
	  else 
	    {
	      evt_weight_PU_up   =0;
	      evt_weight_PU_down =0;
	    }

	  //Offbtag
	  double evt_weight_offbtag_up   = evt_weight * offbtag_weight_up   / offbtag_weight; 
	  double evt_weight_offbtag_down = evt_weight * offbtag_weight_down / offbtag_weight;

	  //Onlbtag
	  double evt_weight_onlbtag_up   = evt_weight * onlbtag_weight_up   / onlbtag_weight; 
	  double evt_weight_onlbtag_down = evt_weight * onlbtag_weight_down / onlbtag_weight;

	  //muID weight
	  double evt_weight_muID_up   = evt_weight * muID_weight_up   / muID_weight; 
	  double evt_weight_muID_down = evt_weight * muID_weight_down / muID_weight;

	  //mupt trig eff weight
	  double evt_weight_mu_trigeff_up   = evt_weight * mu_trigeff_weight_up   / mu_trigeff_weight; 
	  double evt_weight_mu_trigeff_down = evt_weight * mu_trigeff_weight_down / mu_trigeff_weight;

          //jet trig eff weight
	  double evt_weight_jet_trigeff_up   = evt_weight * jet_trigeff_weight_up   / jet_trigeff_weight; 
	  double evt_weight_jet_trigeff_down = evt_weight * jet_trigeff_weight_down / jet_trigeff_weight;

	  //PU
	  h1["m12_mujsel_PU_up"]        -> Fill(mbb_mujsel, evt_weight_PU_up   );
	  h1["m12_mujsel_PU_down"]      -> Fill(mbb_mujsel, evt_weight_PU_down );

	  //Offline btag
	  h1["m12_mujsel_SFbtag_up"]    -> Fill(mbb_mujsel, evt_weight_offbtag_up  );
	  h1["m12_mujsel_SFbtag_down"]  -> Fill(mbb_mujsel, evt_weight_offbtag_down);

	  //Online btag
	  h1["m12_mujsel_onlSFbtag_up"]    -> Fill(mbb_mujsel, evt_weight_onlbtag_up  );
	  h1["m12_mujsel_onlSFbtag_down"]  -> Fill(mbb_mujsel, evt_weight_onlbtag_down);

	  //Muon ID
	  h1["m12_mujsel_muID_up"]        -> Fill(mbb_mujsel, evt_weight_muID_up   );
	  h1["m12_mujsel_muID_down"]      -> Fill(mbb_mujsel, evt_weight_muID_down );

	  //Muon trig eff
	  h1["m12_mujsel_mu_trigeff_up"]        -> Fill(mbb_mujsel, evt_weight_mu_trigeff_up   );
	  h1["m12_mujsel_mu_trigeff_down"]      -> Fill(mbb_mujsel, evt_weight_mu_trigeff_down );

	  //Jet trig eff
	  h1["m12_mujsel_jet_trigeff_up"]        -> Fill(mbb_mujsel, evt_weight_jet_trigeff_up   );
	  h1["m12_mujsel_jet_trigeff_down"]      -> Fill(mbb_mujsel, evt_weight_jet_trigeff_down );


	  //Offline btag CROSS CHECK by flavour
	  h1["m12_mujsel_SFb_up"]    -> Fill(mbb_mujsel, evt_weight *  SFb_weight_up   / offbtag_weight );
	  h1["m12_mujsel_SFb_down"]  -> Fill(mbb_mujsel, evt_weight *  SFb_weight_down / offbtag_weight );
	  h1["m12_mujsel_SFc_up"]    -> Fill(mbb_mujsel, evt_weight *  SFc_weight_up );
	  h1["m12_mujsel_SFc_down"]  -> Fill(mbb_mujsel, evt_weight *  SFc_weight_down );
	  h1["m12_mujsel_SFudsg_up"]    -> Fill(mbb_mujsel, evt_weight *  SFudsg_weight_up );
	  h1["m12_mujsel_SFudsg_down"]  -> Fill(mbb_mujsel, evt_weight *  SFudsg_weight_down );
	  

	  //Fill trees with evt_weight systematics: PU, offline and online btag
	  weight= evt_weight_PU_up;
	  m12_vars["m12_PU_up"]->Fill();
	  weight= evt_weight_PU_down;
	  m12_vars["m12_PU_down"]->Fill();

	  weight= evt_weight_offbtag_up;
	  m12_vars["m12_SFbtag_up"]->Fill();
	  weight= evt_weight_offbtag_down;
	  m12_vars["m12_SFbtag_down"]->Fill();

	  weight= evt_weight_onlbtag_up;
	  m12_vars["m12_onlSFbtag_up"]->Fill();
	  weight= evt_weight_onlbtag_down;
	  m12_vars["m12_onlSFbtag_down"]->Fill();

	  weight= evt_weight_muID_up;
	  m12_vars["m12_muID_up"]->Fill();
	  weight= evt_weight_muID_down;
	  m12_vars["m12_muID_down"]->Fill();

	  weight= evt_weight_mu_trigeff_up;
	  m12_vars["m12_mu_trigeff_up"]->Fill();
	  weight= evt_weight_mu_trigeff_down;
	  m12_vars["m12_mu_trigeff_down"]->Fill();

	  weight= evt_weight_jet_trigeff_up;
	  m12_vars["m12_jet_trigeff_up"]->Fill();
	  weight= evt_weight_jet_trigeff_down;
	  m12_vars["m12_jet_trigeff_down"]->Fill();


	  weight=evt_weight;
	  /*
	  cout << "XCHECK after all shit  m12=" << mbb_mujsel << "weight = " << evt_weight << endl;
	  cout << selectedJets[0]->pt() << "  "<<  selectedJets[1]->pt() << endl; 
	  cout << selectedJets[0]->jerMatch() << " " << selectedJets[1]->jerMatch() << endl;

	  cout << "----------------------------------------------" << endl;
	  */

	  //JER 
	  std::vector<std::string> vars = { "up", "down" }; // last entry to restore JER central value
	  Jet * j0 = selectedJets[0];
	  Jet * j1 = selectedJets[1];
	  Jet * j0_beforeVars = new Jet(j0->pt(), j0->eta(), j0->phi(), j0->e()); 
	  Jet * j1_beforeVars = new Jet(j1->pt(), j1->eta(), j1->phi(), j1->e()); 

	  //to cross-check only
	  /*	  Jet * j0_JES_up   = new Jet(j0->pt(), j0->eta(), j0->phi(), j0->e()); 
	  Jet * j0_JES_down = new Jet(j0->pt(), j0->eta(), j0->phi(), j0->e()); 
	  Jet * j1_JES_up   = new Jet(j1->pt(), j1->eta(), j1->phi(), j1->e()); 
	  Jet * j1_JES_down = new Jet(j1->pt(), j1->eta(), j1->phi(), j1->e()); 
	  */

	  Jet * j0_JER_up   = new Jet(j0->pt(), j0->eta(), j0->phi(), j0->e()); 
	  Jet * j0_JER_down = new Jet(j0->pt(), j0->eta(), j0->phi(), j0->e()); 
	  Jet * j1_JER_up   = new Jet(j1->pt(), j1->eta(), j1->phi(), j1->e()); 
	  Jet * j1_JER_down = new Jet(j1->pt(), j1->eta(), j1->phi(), j1->e()); 
	  
	  
	  for ( auto & var : vars )
	    {
	      correctJetpt( *j0 , j0->jerCorrection(var.c_str())/ j0->jerCorrection() );
	      correctJetpt( *j1 , j1->jerCorrection(var.c_str())/ j1->jerCorrection() );
	  
	      h1[Form("m12_mujsel_JER_%s",var.c_str())] -> Fill( (j0->p4() + j1->p4()).M() , evt_weight);
	      h1[Form("pt_0_JER_%s",var.c_str())] -> Fill( j0->pt() , evt_weight);
	      h1[Form("pt_1_JER_%s",var.c_str())] -> Fill( j1->pt() , evt_weight);
	      mbb = (j0->p4() + j1->p4()).M();
	      m12_vars[("m12_JER_"+var).c_str()]->Fill();

	      //restore pT before variation
	      	      if ( var =="up" )
		{ 
		  copyJetpt( *j0_JER_up, *j0);
		  copyJetpt( *j1_JER_up, *j1);
		}

	      if ( var =="down" ) 
		{
		  copyJetpt( *j0_JER_down, *j0);
		  copyJetpt( *j1_JER_down, *j1);
		  }

	      copyJetpt( *j0, *j0_beforeVars);
	      copyJetpt( *j1, *j1_beforeVars);
	    }
	  
	  //JES
	  std::vector<std::pair<string,int>> sigmas = { {"up",2}, {"down",-2} };

	  for ( auto & var : sigmas)
	    {
	      correctJetpt( *j0, 1 + var.second*j0->jecUncert()  );
	      correctJetpt( *j1, 1 + var.second*j1->jecUncert()  );
	      //	      cout << "JES " << var.first << " j0 = "<< j0->pt() << " j1 =  "<< j1->pt() << "m12 =  "<< (j0->p4() + j1->p4()).M() << endl;
	      
	      h1[Form("m12_mujsel_JES_%s",var.first.c_str())] -> Fill( (j0->p4() + j1->p4()).M() , evt_weight); 
	      h1[Form("pt_0_JES_%s",var.first.c_str())] -> Fill( j0->pt() , evt_weight);
	      h1[Form("pt_1_JES_%s",var.first.c_str())] -> Fill( j1->pt() , evt_weight);
	      mbb = (j0->p4() + j1->p4()).M();
	      m12_vars[("m12_JES_"+var.first).c_str()]->Fill();

	      //restore pT before variation
	      /*	      if ( var.first =="up" )
		{ 
		  copyJetpt( *j0_JES_up, *j0);
		  copyJetpt( *j1_JES_up, *j1);
		}

	      if ( var.first =="down" ) 
		{
		  copyJetpt( *j0_JES_down, *j0);
		  copyJetpt( *j1_JES_down, *j1);
		  }*/

	      copyJetpt( *j0, *j0_beforeVars);
	      copyJetpt( *j1, *j1_beforeVars);     

	   }
	  
	  //h1["j0_JES_diff"]->Fill( j0_JES_up->pt()-j0_JES_down->pt() );
	  h1["j0_JER_diff"]->Fill( j0_JER_up->pt()-j0_JER_down->pt() );
	  //h1["j1_JES_diff"]->Fill( j1_JES_up->pt()-j1_JES_down->pt() );
	  h1["j1_JER_diff"]->Fill( j1_JER_up->pt()-j1_JER_down->pt() );

	  //h1["m12_JES_diff"]->Fill( (j0_JES_up->p4() + j1_JES_up->p4()).M() - (j0_JES_down->p4() + j1_JES_down->p4()).M()  );
	  h1["m12_JER_diff"]->Fill( (j0_JER_up->p4() + j1_JER_up->p4()).M() - (j0_JER_down->p4() + j1_JER_down->p4()).M()  );
	  

	}


      /*      
      //4th jet veto optimization
     
      if ( !signalregion_ || isMC_ )
	{
	  if ( (int)selectedJets.size() > 3 && selectedJets[3]->pt() < 30. ) h1["m12_mujsel_j4veto30"] -> Fill(mbbw[window], evt_weight);
	  if ( (int)selectedJets.size() > 3 && selectedJets[3]->pt() < 35. ) h1["m12_mujsel_j4veto35"] -> Fill(mbbw[window], evt_weight);
	  if ( (int)selectedJets.size() > 3 && selectedJets[3]->pt() < 40. ) h1["m12_mujsel_j4veto40"] -> Fill(mbbw[window], evt_weight);
	}

      if ( !signalregion_ || isMC_ )
	{
          if ( selectedJets[2]->pt() < 60. ) h1["m12_mujsel_j3veto60"] -> Fill(mbbw[window], evt_weight);
	  if ( selectedJets[2]->pt() < 80. ) h1["m12_mujsel_j3veto80"] -> Fill(mbbw[window], evt_weight);
	  if ( selectedJets[2]->pt() < 100. ) h1["m12_mujsel_j3veto100"] -> Fill(mbbw[window], evt_weight);
	  if ( selectedJets[2]->pt() < 120. ) h1["m12_mujsel_j3veto120"] -> Fill(mbbw[window], evt_weight);
	}

        
     // Mujet optimization
      if ( !signalregion_ || isMC_ )
	{
	  
	  if ( window==1 )
	    {
	      //	if ( selectedJets[muj_i]->pt() > 45.0 ) h1["m12_mujsel_muj45"] -> Fill(mbbw[window], evt_weight);
	      if ( selectedJets[jet_i]->pt() > 47.0 && selectedJets[muj_i]->pt() > 55.0 ) h1["m12_mujsel_j47_muj55"] -> Fill(mbbw[window], evt_weight);
	      if ( selectedJets[jet_i]->pt() > 47.0 && selectedJets[muj_i]->pt() > 57.0 ) h1["m12_mujsel_j47_muj57"] -> Fill(mbbw[window], evt_weight);
	      if ( selectedJets[jet_i]->pt() > 47.0 && selectedJets[muj_i]->pt() > 60.0 ) h1["m12_mujsel_j47_muj60"] -> Fill(mbbw[window], evt_weight);
	      if ( selectedJets[jet_i]->pt() > 47.0 && selectedJets[muj_i]->pt() > 62.0 ) h1["m12_mujsel_j47_muj62"] -> Fill(mbbw[window], evt_weight);

	      if ( selectedJets[jet_i]->pt() > 50.0 && selectedJets[muj_i]->pt() > 55.0 ) h1["m12_mujsel_j50_muj55"] -> Fill(mbbw[window], evt_weight);
	      if ( selectedJets[jet_i]->pt() > 50.0 && selectedJets[muj_i]->pt() > 57.0 ) h1["m12_mujsel_j50_muj57"] -> Fill(mbbw[window], evt_weight);
	      if ( selectedJets[jet_i]->pt() > 50.0 && selectedJets[muj_i]->pt() > 60.0 ) h1["m12_mujsel_j50_muj60"] -> Fill(mbbw[window], evt_weight);
	      if ( selectedJets[jet_i]->pt() > 50.0 && selectedJets[muj_i]->pt() > 62.0 ) h1["m12_mujsel_j50_muj62"] -> Fill(mbbw[window], evt_weight);
	    }

	}

*/     
      Jet * jet12 = new Jet;
      jet12->p4( selectedJets[0]->p4() + selectedJets[1]->p4());
      h1["pt12_mujsel"]->Fill(jet12->pt(), evt_weight);
      //      h1["eta12_mujsel"]->Fill(jet12->eta(), evt_weight);
      
      if (isMC_ || ( ! isMC_ && ! signalregion_) ) h1["nentries_mujsel"]->Fill(1,evt_weight);

      /*      cout << "------------------------------------------------------------" << endl;
      cout << selectedJets[0]->pt() << " " << selectedJets[0]->eta() << " " << selectedJets[0]->phi() << endl ;
      cout << selectedJets[1]->pt() << " " << selectedJets[1]->eta() << " " << selectedJets[1]->phi() << endl ;
      cout << jet12->pt() << " " << jet12->eta() << " " << jet12->phi() << endl;
      */
      /*      
      if( jet12->pt() < mbbw[window]/2.) h1["m12_mujsel_pt12_m12_0p5"]-> Fill( mbbw[window] , evt_weight);
      if( jet12->pt() < mbbw[window]/3.) h1["m12_mujsel_pt12_m12_0p3"]-> Fill( mbbw[window] , evt_weight);
      if( selectedJets[2]->pt() < mbbw[window]/2.) h1["m12_mujsel_pt3_m12_0p5"]-> Fill( mbbw[window] , evt_weight);
      if( selectedJets[2]->pt() < mbbw[window]/3.) h1["m12_mujsel_pt3_m12_0p3"]-> Fill( mbbw[window] , evt_weight);
      */
                                                                                                        
      ScaleHistograms(inputlist_ );

      if (debug) cout << "all histos scaled and filled " << endl;

      
   //Clearing vectors 
   selectedGenJets.clear();
   higgsDaughters.clear();
   
   selectedJets.clear();
   selectedMuons.clear();
   selectedMuonsinJet.clear();
   matched_jet_i.clear();



   } //end event loop


   for ( auto & h : h1 )
     {
       h.second->Write();
     }
 
   tree->Write();
   tree0->Write();
   tree1->Write();
   tree2->Write();
   tree3->Write();
   tree4->Write();
   tree5->Write();
   tree6->Write();

   for ( auto & m12_var : m12_vars )
     {
       m12_var.second->Write();
     }

  
   outFile -> Close();
   delete outFile;
  
           
   
// PRINT OUTS  ============================================================================================ 

//   cout << nentries_m12 << endl;
   h1["n_m12_mujsel"]->SetBinContent(1,nentries_m12);
   
   double fracAbs[10];
   double fracRel[10];
   std::string cuts[10];
 
   cuts[0] = "Triggered with PS " + std::to_string(prescaleEra);
   if ( signalregion_ ) cuts[0] = "Triggered";
   if ( isMC_ && signalregion_ ) cuts[0] = "No Trigger";
   cuts[1] = "Triple idloose-jet";
   cuts[2] = "Triple jet kinematics";
   cuts[3] = "Delta R(i;j)";
   cuts[4] = "Delta eta(j1;j2)";
   cuts[5] = "btagged (bbnb)";
   if ( signalregion_ ) cuts[5] = "btagged (bbb)";
   cuts[6] = "Muon jet sel";
   cuts[7] = "Matched to online tobj ";
   cuts[8] = "4th jet veto ";
   if (isMC_ && signalregion_ )  
     {
       cuts[7] = "Triggered ";
       cuts[8] = "Matched to online tobj ";
       cuts[9] = "4th jet veto ";
     }

   printf ("%-23s  %10s  %10s  %10s \n", std::string("Cut flow").c_str(), std::string("# events").c_str(), std::string("absolute").c_str(), std::string("relative").c_str() ); 
   for ( int i = 0; i < 10; ++i ) // 9
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

//---------------------------------------------------------------------------------------------------


void CreateHistograms( )
{

  h1["nGen"]            = new TH1F("nGen","total number of generated events",1,0,1);
  h1["nentries"]        = new TH1F("nentries",       "total number of evts weighted"       ,2,-0.5,1.5);
  h1["nentries_tot"]    = new TH1F("nentries_tot",   "total number of evts total"          ,2,0,2);
  h1["nentries_kin"]    = new TH1F("nentries_kin",   "total number of evts after kin"      ,2,0,2);
  h1["nentries_mujsel"] = new TH1F("nentries_mujsel","total number of evts after all cuts" ,2,0,2);
  h1["n_m12_mujsel"] = new TH1F("n_m12_mujsel","# events in m12 range after mass cut",1,0,1);
  h1["Eff"]      = new TH1F("Eff","efficiency of filtering",1,0,1);

   // Jets --------------------------------------------------------------------------------------

   h1["n"]            = new TH1F("n" ,         "" , 30, 0, 30);      // # jets 
   h1["n_btagsel"]    = new TH1F("n_btagsel" , "" , 30, 0, 30);      //         after btagsel
   h1["n_mujsel"]     = new TH1F("n_mujsel"  , "" , 30, 0, 30);      //         after mujsel
   h1["n_ptmin20"]    = new TH1F("n_ptmin20" , "" , 30, 0, 30);                                                 // # jet pT>20    
   h1["n_ptmin20_btagsel"] = new TH1F("n_ptmin20_btagsel" , "" , 30, 0, 30);                                    //             after btag sel 
   h1["n_ptmin20_mujsel"]  = new TH1F("n_ptmin20_mujsel"  , "" , 30, 0, 30);                                    //             after mujsel
   h2["n_ptmin20_phivseta"]         = new TH2F("n_ptmin20_phivseta","", 25, -2.5, 2.5, 21, -3.15, 3.15 );         // # jet pT>20 in bins of eta and phi 
   h2["n_ptmin20_btagsel_phivseta"] = new TH2F("n_ptmin20_btagsel_phivseta","", 25, -2.5, 2.5, 21, -3.15, 3.15 ); //                   & after btag sel

   // Muons -----------------------------------------------------------------------------------

   h1["n_muons"]           = new TH1F("n_muons"  , "" , 30, 0, 30);                           // # muons passing ID
   h2["n_muons_phivseta"]  = new TH2F("n_muons_phivseta","", 25, -2.5, 2.5, 21, -3.15, 3.15 );  //         passing kinematic

   //after dR selection
   h1["pt_mu"]         = new TH1F("pt_mu"    , "" , 100, 0, 300);  
   h1["eta_mu"]        = new TH1F("eta_mu   ", "" , 100, -5, 5);  
   h2["eta_phi_mu"]    = new TH2F("eta_phi_mu","", 25, -2.5, 2.5, 21, -3.15, 3.15 ); 
   

   h1["dR_muj"]        = new TH1F("dR_muj" , "" , 100, 0, 0.5);   
   h1["dR_muj0"]       = new TH1F("dR_muj0", "" , 100, 0, 0.5);
   h1["dR_muj1"]       = new TH1F("dR_muj1", "" , 100, 0, 0.5);


   //after offline  matching to
   h1["pt_mumatched"]  = new TH1F("pt_mumatched" , "" , 100, 0, 300); 
   h1["eta_mumatched"] = new TH1F("eta_mumatched", "" , 100, -5, 5);  
   h1["n_muinjet"]     = new TH1F("n_muinjet" , "" , 30, 0, 30);      // # muons matched in jet selected
   h1["n_mujets"]      = new TH1F("n_mujets" , "" , 5, 0.5, 5.5);     // # real mu-jets per event 
   h1["n_jetsmatched_mujmatched"]  = new TH1F("n_jetsmatched_mujmatched" , "" , 5, 0.5, 5.5);     // # real mu-jets per event 
   h1["n_jetsmatched_mujsel"]      = new TH1F("n_jetsmatched_mujsel" , "" , 5, 0.5, 5.5);     // # real mu-jets per event


 
  // LEADING Jets --------------------------------------------------------------------
 
  for ( int i = 0 ; i < njetsmin_ ; ++i )
   {
      //no sel
      h1[Form("pt_%i",i)]         = new TH1F(Form("pt_%i",i) , "" , 400, 0, 1000);
      h1[Form("pt_%i_eta0p9",i)]      = new TH1F(Form("pt_%i_eta0p9",i) , "" , 400, 0, 1000);
      h1[Form("pt_%i_0p9eta1p5",i)]     = new TH1F(Form("pt_%i_0p9eta1p5",i) , "" , 400, 0, 1000);
      h1[Form("pt_%i_1p5eta2p0",i)]     = new TH1F(Form("pt_%i_1p5eta2p0",i) , "" , 400, 0, 1000);
      h1[Form("pt_%i_2p0eta",i)]      = new TH1F(Form("pt_%i_2p0eta",i) , "" , 400, 0, 1000);
      h1[Form("eta_%i",i)]        = new TH1F(Form("eta_%i",i) , "" , 100, -5, 5);
      h1[Form("phi_%i",i)]        = new TH1F(Form("phi_%i",i) , "" , 100, -4, 4);
      h1[Form("deepflavour_%i",i)]        = new TH1F(Form("deepflavour_%i",i) , "" , 400, 0, 1);
      h1[Form("deepcsv_%i",i)]    = new TH1F(Form("deepcsv_%i",i) , "" , 400, 0, 1);
      h2[Form("eta_phi_%i",i)]       = new TH2F(Form("eta_phi_%i",i) , "" , 25, -2.5, 2.5, 21, -3.15, 3.15 );
      
      // after btagsel sel 
      h1[Form("pt_%i_btagsel",i)]     = new TH1F(Form("pt_%i_btagsel",i) , "" , 400, 0, 1000);
      h1[Form("eta_%i_btagsel",i)]    = new TH1F(Form("eta_%i_btagsel",i) , "" , 100, -5, 5);
      h1[Form("phi_%i_btagsel",i)]    = new TH1F(Form("phi_%i_btagsel",i) , "" , 100, -4, 4);
      h1[Form("deepflavour_%i_btagsel",i)]    = new TH1F(Form("deepflavour_%i_btagsel",i) , "" , 400, 0, 1);  
      h1[Form("deepcsv_%i_btagsel",i)]= new TH1F(Form("deepcsv_%i_btagsel",i) , "" , 400, 0, 1);
      h2[Form("eta_phi_%i_btagsel",i)] = new TH2F(Form("eta_phi_%i_btagsel",i) , "" , 25, -2.5, 2.5, 21, -3.15, 3.15 );
      
      // after mujet sel
      h1[Form("pt_%i_mujsel",i)]     = new TH1F(Form("pt_%i_mujsel",i) , "" , 400, 0, 1000);
      h1[Form("eta_%i_mujsel",i)]    = new TH1F(Form("eta_%i_mujsel",i) , "" , 100, -5, 5);
      h1[Form("phi_%i_mujsel",i)]    = new TH1F(Form("phi_%i_mujsel",i) , "" , 100, -4, 4);
      h1[Form("deepflavour_%i_mujsel",i)]    = new TH1F(Form("deepflavour_%i_mujsel",i) , "" , 400, 0, 1);
      h1[Form("deepcsv_%i_mujsel",i)]= new TH1F(Form("deepcsv_%i_mujsel",i) , "" , 400, 0, 1);
 
   }

  h1["pt_muj"]        = new TH1F("pt_muj" , "" , 400, 0, 1000);
  h1["pt_j"]          = new TH1F("pt_j"   , "" , 400, 0, 1000);
  h1["eta_muj"]       = new TH1F("eta_muj" , "" , 100, -5, 5);
  h1["eta_j"]         = new TH1F("eta_j"   , "" , 100, -5, 5);
  h1["pt_muj_mujsel"]        = new TH1F("pt_muj_mujsel" , "" , 400, 0, 1000);
  h1["pt_j_mujsel"]          = new TH1F("pt_j_mujsel"   , "" , 400, 0, 1000);
  h1["eta_muj_mujsel"]       = new TH1F("eta_muj_mujsel" , "" , 100, -5, 5);
  h1["eta_j_mujsel"]         = new TH1F("eta_j_mujsel"   , "" , 100, -5, 5);
  /*
  h1["pt_j_mujsel_muj45"]        = new TH1F("pt_j_mujsel_muj45" , "" , 400, 0, 1000);
  h1["pt_j_mujsel_muj50"]        = new TH1F("pt_j_mujsel_muj50" , "" , 400, 0, 1000);
  */
  h1["pt_softjet"] = new TH1F("pt_softjet" ,   "" , 100, 0, 1000);
  h1["dR_jsoftjet"] = new TH1F("dR_jsoftjet" , "" , 30, 0.3, 1.0);
  h1["rank_softjet"] = new TH1F("rank_softjet" ,   "" , 100, 0, 10);

  /*  //m12 window optimization
  std::vector<int> m12_low;
  
  m12_low.push_back(80);
  m12_low.push_back(90);
  m12_low.push_back(100);

   for ( auto & cut : m12_low )
   {

     for ( int i = 0 ; i < njetsmin_ ; ++i )
       {
	 h1[Form("pt_%i_mujsel_m12low_%i",i,cut)]     = new TH1F(Form("pt_%i_mujsel_m12low_%i",i,cut) , "" , 400, 0, 1000);
       }

     h1[Form("pt_j_mujsel_m12low_%i",cut)]     = new TH1F(Form("pt_j_mujsel_m12low_%i",cut), "" , 400, 0, 1000);
     h1[Form("pt_muj_mujsel_m12low_%i",cut)]   = new TH1F(Form("pt_muj_mujsel_m12low_%i",cut), "" , 400, 0, 1000);
     h1[Form("pt_mumatched_m12low_%i",cut)]   = new TH1F(Form("pt_mumatched_m12low_%i",cut), "" , 100, 0, 300);   
   }
  */

   // Dijet mass -----------------------------------------------------------------------------------

   h1["n_m12_sr1"] = new TH1F("n_m12_sr1", "#entries in sr1"  ,2,0,2);
   h1["n_m12_sr2"] = new TH1F("n_m12_sr2", "#entries in sr2"  ,2,0,2);
   h1["n_m12_sr3"] = new TH1F("n_m12_sr3", "#entries in sr3"  ,2,0,2);

   h1["m12"]         = new TH1F("m12"     , "" , 455, 80, 1000);
   h1["m12_btagsel"] = new TH1F("m12_btagsel" , "" , 455, 80, 1000);
   h1["m12_mujsel"]  = new TH1F("m12_mujsel" , "" , 455, 80, 1000);

   //MC LO vs NLO 
   h1["pt12_mujsel"]         = new TH1F("pt12_mujsel"         , "" , 400, 0, 1200);
   h1["jet2_flavour"]        = new TH1F("jet2_flavour"        , "" , 10, 0, 10);

   h1["genpt_Higgs"]         = new TH1F("genpt_Higgs"         , "" , 400, 0, 1200);
   h1["geneta_Higgs"]        = new TH1F("geneta_Higgs"         , "" , 100, -5, 5);

   h1["genpt_Higgsdau0"]     = new TH1F("genpt_Higgsdau0"     , "" , 400, 0, 1200);
   h1["genpt_Higgsdau1"]     = new TH1F("genpt_Higgsdau1"     , "" , 400, 0, 1200);
   h1["geneta_Higgsdau0"]     = new TH1F("geneta_Higgsdau0"     , "" , 100, -5, 5);
   h1["geneta_Higgsdau1"]     = new TH1F("geneta_Higgsdau1"     , "" , 100, -5, 5);

   h1["genpt_bSpectator"]     = new TH1F("genpt_bSpectator"   , "" , 400, 0, 1200);
   h1["geneta_bSpectator"]    = new TH1F("geneta_bSpectator"  , "" , 100, -5, 5);

   h1["deltaR_Ha"]         = new TH1F("deltaR_Ha" , "" , 50, 0, 4.0);
   h1["deltaR_Hd1"]        = new TH1F("deltaR_Hd1", "" , 50, 0, 4.0);
   h1["deltaR_Hd2"]        = new TH1F("deltaR_Hd2", "" , 50, 0, 4.0);


   //   h1["offbtag_sf_j0"]           = new TH1F("offbtag_sf_j0" , "" , 50, 0, 0);
   //   h1["offbtag_sf_j1"]           = new TH1F("offbtag_sf_j1" , "" , 50, 0, 0);
   h1["j0_JES_diff"]     = new TH1F("j0_JES_diff" , "" , 150, -5, 30);
   h1["j0_JER_diff"]     = new TH1F("j0_JER_diff" , "" , 150, -5, 30);
   h1["j1_JES_diff"]     = new TH1F("j1_JES_diff" , "" , 150, -5, 30);
   h1["j1_JER_diff"]     = new TH1F("j1_JER_diff" , "" , 150, -5, 30);
   h1["m12_JES_diff"]     = new TH1F("m12_JES_diff" , "" , 150, -5, 30);
   h1["m12_JER_diff"]     = new TH1F("m12_JER_diff" , "" , 150, -5, 30);

   h1["j0_JES_diff_befkin"]     = new TH1F("j0_JES_diff_befkin" , "" , 150, -5, 30);
   h1["j0_JER_diff_befkin"]     = new TH1F("j0_JER_diff_befkin" , "" , 150, -5, 30);
   h1["j1_JES_diff_befkin"]     = new TH1F("j1_JES_diff_befkin" , "" , 150, -5, 30);
   h1["j1_JER_diff_befkin"]     = new TH1F("j1_JER_diff_befkin" , "" , 150, -5, 30);
   h1["m12_JES_diff_befkin"]     = new TH1F("m12_JES_diff_befkin" , "" , 150, -5, 30);
   h1["m12_JER_diff_befkin"]     = new TH1F("m12_JER_diff_befkin" , "" , 150, -5, 30);


   h1["pt_0_uncorr"]      =  new TH1F("pt_0_uncorr"  , "" , 455, 80, 1000);
   h1["pt_0_mu_trigeff"]  =  new TH1F("pt_0_mu_trigeff"  , "" , 455, 80, 1000);
   h1["pt_0_jet_trigeff"] =  new TH1F("pt_0_jet_trigeff" , "" , 455, 80, 1000);
   h1["pt_0_muID"]        =  new TH1F("pt_0_muID"        , "" , 455, 80, 1000);

   h1["pt_1_uncorr"]      =  new TH1F("pt_1_uncorr"  , "" , 455, 80, 1000);
   h1["pt_1_mu_trigeff"]  =  new TH1F("pt_1_mu_trigeff"  , "" , 455, 80, 1000);
   h1["pt_1_jet_trigeff"] =  new TH1F("pt_1_jet_trigeff" , "" , 455, 80, 1000);
   h1["pt_1_muID"]        =  new TH1F("pt_1_muID"        , "" , 455, 80, 1000);

   h1["m12_uncorr"]      =  new TH1F("m12_uncorr"  , "" , 455, 80, 1000);
   h1["m12_mu_trigeff"]  =  new TH1F("m12_mu_trigeff"  , "" , 455, 80, 1000);
   h1["m12_jet_trigeff"] =  new TH1F("m12_jet_trigeff" , "" , 455, 80, 1000);
   h1["m12_muID"]        =  new TH1F("m12_muID"        , "" , 455, 80, 1000);


/*
   h1["m12_mujsel_j4veto30"] = new TH1F("m12_mujsel_j4veto30" , "" , 455, 80, 1000);
   h1["m12_mujsel_j4veto35"] = new TH1F("m12_mujsel_j4veto35" , "" , 455, 80, 1000);
   h1["m12_mujsel_j4veto40"] = new TH1F("m12_mujsel_j4veto40" , "" , 455, 80, 1000);

   h1["m12_mujsel_j3veto60"] = new TH1F("m12_mujsel_j3veto60" , "" , 455, 80, 1000);
   h1["m12_mujsel_j3veto80"] = new TH1F("m12_mujsel_j3veto80" , "" , 455, 80, 1000);
   h1["m12_mujsel_j3veto100"] = new TH1F("m12_mujsel_j3veto100" , "" , 455, 80, 1000);
   h1["m12_mujsel_j3veto120"] = new TH1F("m12_mujsel_j3veto120" , "" , 455, 80, 1000);

   //scaled m12 cut
   h1["m12_mujsel_pt12_m12_0p5"] = new TH1F("m12_mujsel_pt12_m12_0p5" , "" , 455, 80, 1000);
   h1["m12_mujsel_pt12_m12_0p3"] = new TH1F("m12_mujsel_pt12_m12_0p3" , "" , 455, 80, 1000);
   h1["m12_mujsel_pt3_m12_0p5"] = new TH1F("m12_mujsel_pt3_m12_0p5" , "" , 455, 80, 1000);
   h1["m12_mujsel_pt3_m12_0p3"] = new TH1F("m12_mujsel_pt3_m12_0p3" , "" , 455, 80, 1000);

   h1["m12_mujsel_j47_muj55"] = new TH1F("m12_mujsel_j47_muj55" , "" , 455, 80, 1000);
   h1["m12_mujsel_j47_muj57"] = new TH1F("m12_mujsel_j47_muj57" , "" , 455, 80, 1000);
   h1["m12_mujsel_j47_muj60"] = new TH1F("m12_mujsel_j47_muj60" , "" , 455, 80, 1000);
   h1["m12_mujsel_j47_muj62"] = new TH1F("m12_mujsel_j47_muj62" , "" , 455, 80, 1000);

   h1["m12_mujsel_j50_muj55"] = new TH1F("m12_mujsel_j50_muj55" , "" , 455, 80, 1000);
   h1["m12_mujsel_j50_muj57"] = new TH1F("m12_mujsel_j50_muj57" , "" , 455, 80, 1000);
   h1["m12_mujsel_j50_muj60"] = new TH1F("m12_mujsel_j50_muj60" , "" , 455, 80, 1000);
   h1["m12_mujsel_j50_muj62"] = new TH1F("m12_mujsel_j50_muj62" , "" , 455, 80, 1000);


   h1["m12_mujsel_loweta"]    = new TH1F("m12_mujsel_loweta" , "" , 455, 80, 1000);
   h1["m12_mujsel_mediumeta"] = new TH1F("m12_mujsel_mediumeta" , "" , 455, 80, 1000);
   h1["m12_mujsel_higheta"]   = new TH1F("m12_mujsel_higheta" , "" , 455, 80, 1000);
   h1["m12_mujsel_mixedeta"]  = new TH1F("m12_mujsel_mixedeta" , "" , 455, 80, 1000);
   */

  for ( int i = 0 ; i < nsubsamples ; ++i )
     h1[Form("m12_mujsel_%i",i)]  = new TH1F(Form("m12_mujsel_%i",i) , "" , 455, 80, 1000);

   for ( int d =0 ; d <3; d++ )
     h1[Form("m12_mujsel_%id",d)] = new TH1F(Form("m12_mujsel_%id",d), "" , 455, 80, 1000);

       for (auto & var : vars)
	 {

	 for ( auto & sys : systematics)
	   	 {
		   h1[Form("m12_mujsel_%s_%s",sys.c_str(),var.c_str())] = new TH1F(Form("m12_mujsel_%s_%s",sys.c_str(),var.c_str()), "", 455, 80, 1000);
		   h1[Form("m12_mujsel_%s_%s_befkin",sys.c_str(),var.c_str())] = new TH1F(Form("m12_mujsel_%s_%s_befkin",sys.c_str(),var.c_str()), "", 455, 80, 1000);

		 }

	 h1[Form("pt_0_JER_%s",var.c_str())] = new TH1F(Form("pt_0_JER_%s",var.c_str()), "", 400, 0, 1000);
	 h1[Form("pt_0_JES_%s",var.c_str())] = new TH1F(Form("pt_0_JES_%s",var.c_str()), "", 400, 0, 1000);
	 h1[Form("pt_1_JER_%s",var.c_str())] = new TH1F(Form("pt_1_JER_%s",var.c_str()), "", 400, 0, 1000);
	 h1[Form("pt_1_JES_%s",var.c_str())] = new TH1F(Form("pt_1_JES_%s",var.c_str()), "", 400, 0, 1000);

	 h1[Form("pt_0_JER_%s_befkin",var.c_str())] = new TH1F(Form("pt_0_JER_%s_befkin",var.c_str()), "", 400, 0, 1000);
	 h1[Form("pt_0_JES_%s_befkin",var.c_str())] = new TH1F(Form("pt_0_JES_%s_befkin",var.c_str()), "", 400, 0, 1000);
	 h1[Form("pt_1_JER_%s_befkin",var.c_str())] = new TH1F(Form("pt_1_JER_%s_befkin",var.c_str()), "", 400, 0, 1000);
	 h1[Form("pt_1_JES_%s_befkin",var.c_str())] = new TH1F(Form("pt_1_JES_%s_befkin",var.c_str()), "", 400, 0, 1000);


	 }

   std::vector<string> flavours = { "b","c","udsg"};

   for ( auto & flavour : flavours )
     {
       for (auto & var : vars)
	 {
	   h1[Form("m12_mujsel_SF%s_%s",flavour.c_str(),var.c_str())] = new TH1F(Form("m12_mujsel_SF%s_%s",flavour.c_str(),var.c_str()), "", 455, 80, 1000);
	 }
     }

}


void ScaleHistograms( const std::string & inputlist )
{
 float sf = 1;

   /*   if ( process == "mssmhbb" ) // scaling with instantaneous luminosity 
   {

      float genLumi = nGenTotal_/SignalCrossSection(); //cross section
      sf = lumi_/genLumi;
      }*/
   /*
   if ( process == "qcd" || process == "ttbar" ) // scaling with instantaneous luminosity 
   { 
     Analysis analysis(inputlist);
     analysis.crossSections("MssmHbb/Metadata/CrossSections");
     float genLumi = nGenTotal_/analysis.crossSection(); //cross section
     sf = lumi_/genLumi;
     
     std::cout << "Xsec_ "     << analysis.crossSection() << std::endl;  
     std::cout << "Gen Lumi =" << genLumi                 << std::endl;
     std::cout << "sf = "      <<  sf                     << std::endl;
     
     }*/
   for ( auto & h : h1 )
   {
      if ( h.first != "nGen" )
      {
         h.second->Scale(sf);
      }
   }

}

/*

void WriteHistograms(const std::string & filename)
{
   
   TFile * outFile = new TFile(filename.c_str(),"RECREATE");
   
   for ( auto & h : h1 )
   {
      h.second->Write();
   }
   outFile -> Close();
   delete outFile;
}
*/


float btagMin(const string & wp)
{
  float min = -1000;
  if ( wp == "loose" )  min = btagwploose_;
  if ( wp == "medium" ) min = btagwpmedium_;
  if ( wp == "tight" )  min = btagwptight_;
   
  return min;
   
}

void correctJetpt ( Jet& jet , const float& corr ) 
{
  TLorentzVector corrjet ;
  corrjet.SetPtEtaPhiM( jet.pt()*corr, jet.eta(), jet.phi(), (jet.p4()).M());
  jet.p4( corrjet );
  
} 

void copyJetpt ( Jet& jet , Jet& newjet ) 
{
  TLorentzVector appo ;
  appo.SetPtEtaPhiE( newjet.pt() , newjet.eta(), newjet.phi(), newjet.e());
  jet.p4( appo );
} 


void applyPrescale ( const int& run, const double& random, float& prescaleEra, const int nsubsamples, int& window  )  
{
  //Get correct Prescale for each Era
  //  for ( size_t i = 0 ; i < eraRanges.size() ; i++ )
  //    if (run >= eraRanges[i].first && run <= eraRanges[i].second ) prescaleEra = 7.0;
  
  //Split bkg in nsubsamples equally likely distributions
  for ( int w=0 ; w < nsubsamples; w++ )
      {
	if ( w == nsubsamples -1 && random > (w+1)/prescaleEra ) continue; // if last window and Rand outside maximal window range skip  
	if ( w/prescaleEra < random  &&  random < (w+1)/prescaleEra ) 
	  {		  
	    window = w;
	    break;
	  }
	
      }
	
}

double onlSFbtag ( double pt , double nsigma ) 
{
  // SF = p0 +p1 * pT linear fit
  const double p0    =  0.852425;
  const double errp0 =  0.00713844;
  const double p1    = -0.0000616358;
  const double errp1 =  0.000055093;

  return (p0+nsigma*errp0) + (p1+nsigma*errp1)*pt;
  
} 

double muonptWeight ( double pt , double nsigma ) 
{
  // Scale factor muon pT trigger
  const double sf_pt12 =  0.94453947;
  const double sf_pt13 =  0.93695156;
  const double sf_pt16 =  0.96314276;
  const double sf_pt22 =  0.96445782;

  const double err_pt12 =  0.13653615;
  const double err_pt13 =  0.06519553;
  const double err_pt16 =  0.03708156;
  const double err_pt22 =  0.07786825;

  double weight;
  double uncert;

  if      ( pt > 11.5 && pt < 12.5 ) { weight = sf_pt12; uncert = err_pt12;}   
  else if ( pt > 12.5 && pt < 13.5 ) { weight = sf_pt13; uncert = err_pt13;}   
  else if ( pt > 13.5 && pt < 18.5 ) { weight = sf_pt16; uncert = err_pt16;}   
  else if ( pt > 18.5 && pt < 30.  ) { weight = sf_pt22; uncert = err_pt22;}   
  else { weight = sf_pt22; uncert = err_pt22;} 
  //  else { weight = sf_pt22; uncert = 2*err_pt22;} 

  return weight* ( 1.0 + nsigma*uncert ) ;
  
} 


double jetTriggerScaleFactor ( double x , double eta , double nsigma  )
{

  double p0 = -999, p1 = -999, p2=-999;
  double c00 = -999, c01 = -999, c02 = -999, c11 = -999, c12 = -999, c22 =-999;

  //Parameters and covariance matrix // cij = ci *cj 
  if ( fabs(eta) < 1.0)
    {
      p0 = 9.98377e-01;  p1 = 1.73918e+01; p2 = 3.19126e+01;
      c00 = 8.558e-08 ; c01 = -0.0005863; c02 = 0.0004388; c11 = 24.15 ; c12 = -15.41 ; c22 = 10.3;

    }
  else if ( fabs(eta) > 1.0 &&  fabs(eta) < 1.4 )
    {
      p0 = 9.99650e-01;  p1 = 1.95201e+01; p2 =3.29795e+01 ;
      c00 = 1.385e-07; c01 = -0.0002924; c02 = 0.0002412; c11 = 17.67; c12 = -10.83  ; c22 = 7.02 ;
    }

  else if ( fabs(eta) > 1.4 )
    {
      p0 = 9.98207e-01;  p1 =  3.51120e+01; p2 = 2.44015e+01;
      c00 = 9.755e-08 ; c01 = -0.0001355; c02 = 0.0001203; c11 = 10.8 ; c12 = -7.032; c22 = 4.954 ;
    }

  //Scale factor from erf fit
  double sf = p0*TMath::Erf((x-p1)/p2); 

  double t = (x-p1)/p2;
  double erf = TMath::Erf(t);
  double derf = (2.0/1.77245385091)*std::exp(-t*t); //derivative error function

  double df0 = erf;
  double df1 = -p0*derf/p2;
  double df2 = t*df1;

  double uncert = std::sqrt(c00*df0*df0 + 2*c01*df0*df1 + 2*c02*df0*df2 + c11*df1*df1+ 2*c12*df1*df2 + c22*df2*df2);  

  return sf * ( 1.0 + nsigma *uncert) ;
}








