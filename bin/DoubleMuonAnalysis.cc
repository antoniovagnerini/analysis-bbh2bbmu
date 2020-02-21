
#include "boost/program_options.hpp"
#include "boost/algorithm/string.hpp"
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include "TFile.h" 
#include "TFileCollection.h"
#include "TChain.h"
#include "TH1.h" 

#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TAxis.h"

//#include "tdrstyle.C"

#include "Analysis/Tools/interface/Analysis.h"
#include "Analysis/Tools/bin/macro_config.h"

using namespace std;
using namespace analysis;
using namespace analysis::tools;

static const double mJPsi = 3.10; //pdg value Y(1s) resonance

// =============================================================================================   
int main(int argc, char * argv[])
{
  
   if ( macro_config(argc, argv) != 0 ) return -1;
   
   TH1::SetDefaultSumw2();  // proper treatment of errors when scaling histograms
  
   // Input files list
   Analysis analysis(inputlist_);
   
   analysis.addTree<Muon>("Muons","MssmHbb/Events/slimmedMuons");

   for ( auto & obj : triggerObjectsMuons_ )
       analysis.addTree<TriggerObject> (obj,Form("MssmHbb/Events/slimmedPatTrigger/%s",obj.c_str()));

   analysis.triggerResults("MssmHbb/Events/TriggerResults");

   if( !isMC_ ) analysis.processJsonFile(json_);

   TFile hout(outputRoot_.c_str(),"recreate");
 
   std::map<std::string, TH1F*> h1;
 
   // Muons -----------------------------------------------------------------------------------

   h1["n_muons"]           = new TH1F("n_muons"      , "" , 30, 0, 30);  // # muons passing ID
   h1["n_goodmuons"]       = new TH1F("n_goodmuons"  , "" , 30, 0, 30);  // # muons passing kinematic sel

   h1["pt_mutag"]     = new TH1F("pt_mutag    " , "" , 100, 0, 300); //tag
   h1["eta_mutag"]    = new TH1F("eta_mutag   " , "" , 60, -3, 3);
   //   h1["rank_mutag"]   = new TH1F("rank_mutag  " , "" , 10, 0, 10);

   h1["dr_muL1"]  = new TH1F("dr_muL1" , "" , 80, 0.0, 0.4);
   h1["dr_muL3"]  = new TH1F("dr_muL3" , "" , 80, 0.0, 0.4);

   h1["pt_muprobe"]   = new TH1F("pt_muprobe  " , "" , 100, 0, 300); //probe
   h1["eta_muprobe"]  = new TH1F("eta_muprobe " , "" , 60, -3, 3);

   h1["pt_muprobefl"]   = new TH1F("pt_muprobefl  " , "" , 100, 0, 300); 
   h1["eta_muprobefl"]  = new TH1F("eta_muprobefl " , "" , 60, -3, 3);

   //M12 TagnProbe
   h1["m12_mutag"]           = new TH1F("m12_mutag"       , "" , 400, 0, 15); //m12  
   h1["m12_mutagnprobe"]     = new TH1F("m12_mutagnprobe" , "" , 400, 0, 15);
   h1["m12_mutagnprobefl"]     = new TH1F("m12_mutagnprobefl" , "" , 400, 0, 15);

   std::vector<std::pair<float,float> > pt_bins = {{11.5,12.5},{12.5,13.5},{13.5,20.5},{20.5,30.},{30.,50.}};


   //   if ( hltPath_ == "HLT_Mu8_v" ) pt_bins = {{5,6},{6,7},{7,7.75},{7.75,8.25},{8.25,9},{9,11},{11,30}};
   // else if (hltPath_=="HLT_SingleJet30_Mu12_SinglePFJet40_v") pt_bins = {{8,11},{11,13},{13,15},{15,50}};                       

   //   if      ( hltPath_ == "HLT_Mu8_v" && muonsid_ =="TIGHT" )  pt_bins = {{11.5,12.5},{12.5,13.5},{13.5,18.5},{18.5,30.}};
   //   if ( hltPath_ == "HLT_Mu8_v" ) pt_bins = {{11.5,12.5},{12.5,13.5},{13.5,20.5},{20.5,30.},{30.,50.} };                     
   //    else if (hltPath_=="HLT_SingleJet30_Mu12_SinglePFJet40_v") {{11.5,12.5},{12.5,13.5},{13.5,20.5},{20.5,30.},{30.,50.} };                     
    



   for ( auto & pt_bin : pt_bins ) 
     {
     h1[Form("m12_mutag_%.2fto%.2f"      , pt_bin.first, pt_bin.second)] = new TH1F(Form("m12_mutag_%.2fto%.2f"      , pt_bin.first, pt_bin.second) ,"", 400, 0, 15);
     h1[Form("m12_mutagnprobe_%.2fto%.2f", pt_bin.first, pt_bin.second)] = new TH1F(Form("m12_mutagnprobe_%.2fto%.2f", pt_bin.first, pt_bin.second) ,"", 400, 0, 15);
     h1[Form("m12_mutagnprobefl_%.2fto%.2f", pt_bin.first, pt_bin.second)] = new TH1F(Form("m12_mutagnprobefl_%.2fto%.2f", pt_bin.first, pt_bin.second) ,"", 400, 0, 15);
     }
 
   double m12_mutag;
   double m12_mutagnprobe;

   //edge array
   float edges[ pt_bins.size() +1 ]; //Nedges = Nbins+1 
   for ( size_t i = 0; i< pt_bins.size(); i++ )
    {                                                                                                            
     edges[i] = pt_bins[i].first;  //store lower bin edge
     std::cout<< "bin = " << i << " lowedge =" << edges[i]<< std::endl;
    }

   edges[pt_bins.size()] = pt_bins[pt_bins.size()-1].second; //store pt max as upper edge 
   std::cout<< "bin = " << pt_bins.size() << " highedge =" << edges[ pt_bins.size() ] <<std::endl;

   h1["h_num"] = new TH1F("h_num","",pt_bins.size(),edges);
   h1["h_den"] = new TH1F("h_den","",pt_bins.size(),edges);

   // Analysis of events
   std::cout << "This analysis has " << analysis.size() << " events" << std::endl;
       
   int nsel[10] = { };  
   //int nmatch_mu[10] = { };
   
   if ( nevtmax_ < 0 ) nevtmax_ = analysis.size();
   for ( int i = 0 ; i < nevtmax_ ; ++i )
   { 
      int nmuons =0;
      if ( i > 0 && i%100000==0 ) std::cout << i << "  events processed! " << std::endl;
     
      analysis.event(i);
      if (! isMC_ )
	{
         if (!analysis.selectJson() ) continue; // To use only goodJSonFiles
	}
   
      // #0 TRIGGERED EVTs     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      int triggerFired = analysis.triggerResult(hltPath_);
      if ( !triggerFired ) continue;
      
      ++nsel[0];
 
      //========== MUON SELECTION ============================================================================  
        
      auto slimmedMuons = analysis.collection<Muon>("Muons");
      
      //Muon ID sel +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      std::vector<Muon *> selectedMuons;
      const char muID = muonsid_.at(0) ;
      for ( int m = 0 ; m < slimmedMuons->size() ; ++m )
      {	 
	if ( (muID =='M' && slimmedMuons->at(m).isMediumMuon()) || (muID =='T' && slimmedMuons->at(m).isTightMuon())  )
	 { 
	   selectedMuons.push_back(&slimmedMuons->at(m));	    
	   ++nmuons;
	 }
      }

      if ( (int)selectedMuons.size() < nmuonsmin_ ) continue;
      h1["n_muons"] -> Fill(nmuons);

       ++nsel[1];

      //Muon kinematic sel +++++++++++++++++++++++++++++++++++++++++++++++++++++

       std::vector<Muon *> goodMuons;
       for (  size_t m = 0; m < selectedMuons.size(); ++m )      
       {
         Muon * muon = selectedMuons[m];
	 if ( muon->pt() < muonsptmin_[m] || fabs(muon->eta()) > muonsetamax_[m] ||  muon->pt() > muonsptmax_[m] ) continue;
	 goodMuons.push_back(muon);
       }

       if ( (int)goodMuons.size() < nmuonsmin_ ) continue;
       h1["n_goodmuons"] -> Fill((int)goodMuons.size());

      ++nsel[2];
    
 
     // Matching OFFLINE muons to ONLINE tobj +++++++++++++++++++++++++++++++++++      
           
      auto objL1 = analysis.collection<TriggerObject>(triggerObjectsMuons_[0]);
      auto objL3 = analysis.collection<TriggerObject>(triggerObjectsMuons_[1]);
       
      std::vector<TriggerObject *> triggerObjectL1 ;
      std::vector<TriggerObject *> triggerObjectL3 ;

      for ( int i = 0 ; i < objL1->size() ; ++i )
	{
	  TriggerObject * tobL1 = & objL1->at(i);
	  if (tobL1->pt()  < 12.0 ) continue;
	  if (tobL1->eta() > 2.3 ) continue;
          triggerObjectL1.push_back(tobL1) ;                                                                                                                                               
	  //	  triggerObjectL1.push_back(&objL1->at(i)) ;
	}

      for ( int j = 0 ; j < objL3->size() ; ++j )
	{
	  TriggerObject * tobL3 = & objL3->at(j);
          if (tobL3->pt()  < 12.0 ) continue;
          if (tobL3->eta() > 2.3  ) continue;
	  triggerObjectL3.push_back(tobL3) ;
	}


      std::vector<int> matchedL1_i ;                                    
      std::vector<int> matchedL3_i ;  
   
  	//analysis.match<Muon,TriggerObject>("Muons",triggerObjectsMuons_,drmax_);
      std::vector<Muon *> matchedMuons ;
      std::vector<int> matchedmu_i ;
    
      //std::cout << "--------------------------------------------------------" << std::endl;
      // MUON-TOBJ loop 
      for ( size_t m = 0; m < goodMuons.size(); ++m )
      {
         Muon * muon = goodMuons[m];	 	 
	 for ( size_t i = 0; i < triggerObjectL1.size() ; ++i )
	   {
	     const TriggerObject & tobL1 = *triggerObjectL1[i] ;
	     	    
	       for ( size_t j = 0; j < triggerObjectL3.size() ; ++j )
	       {
		 const TriggerObject & tobL3 = *triggerObjectL3[j] ;
		
		 if (muon->deltaR( tobL1 ) < drmax_ && muon->deltaR( tobL3 ) < 0.005 )
		   {
		     if ( std::find(matchedL1_i.begin(), matchedL1_i.end(), i) != matchedL1_i.end() ) continue; //if either L1 or L3 tobs already matched skip
		     if ( std::find(matchedL3_i.begin(), matchedL3_i.end(), j) != matchedL3_i.end() ) continue;
		     
		     matchedL1_i.push_back(i);
		     matchedL3_i.push_back(j);
		     matchedmu_i.push_back(m); //different muons matched to different tob
		     matchedMuons.push_back(muon);

		     h1["dr_muL1"]  -> Fill(muon->deltaR( tobL1 ) );
		     h1["dr_muL3"]  -> Fill(muon->deltaR( tobL3 ) );
		     //std::cout << "#l1obj = " << i << " #l3obj = " << j << " matched mu = " << m << std::endl;
		   }
	       }//end L3 loop
 
	   }//end L1 loop
      }

      //Muon-pair selection ++++++++++++++++++++++++++++++++++++++++++++++++++++
      std::vector<std::pair<int, int> >    TnProbePairs ;
      std::vector<std::pair<int, int> >  psTnProbePairs ;

      for ( size_t m1 = 0; m1 < matchedMuons.size(); ++m1 )
	{
	 const Muon & muon1  = *matchedMuons[m1];
	  int  nmupairs  =  0 ;
          int   probe_i  = -1 ;
          int psprobe_i  = -1 ;

	  //Tag selection
	  if ( muon1.pt() < ptmin_[0] ) continue; //tag tighter pT  
	  
	  // Tag & Probe pair selection
	  for ( size_t m2 = 0; m2 < goodMuons.size(); ++m2 )
	    {
	      const Muon & muon2 = *goodMuons[m2];

              if ( muon1.deltaR(muon2) > drmin_ && muon1.q()*muon2.q()< 0  )  
		{
		  ++nmupairs;
		  if ( nmupairs<2 )  probe_i = m2; 		  
		  if ( std::find(matchedmu_i.begin(), matchedmu_i.end(), m2) != matchedmu_i.end() ) // probe selection
                                      {       psprobe_i =m2;
				                    break;             }	  		   
		}
	    }
	  
	  if (probe_i   > -1) {  
	    TnProbePairs.push_back(std::make_pair(m1,probe_i));   
	   //std::cout << "tag muon = " << m1 << " probe mu = "   << probe_i  << std::endl; 
                            }
	  if (psprobe_i > -1) { 
            psTnProbePairs.push_back(std::make_pair(m1,psprobe_i));
            //std::cout << "tag muon = " << m1 << " psprobe mu = " << psprobe_i << std::endl;
                          }
	}
	

      //---------------------------------------------------------------------------------------------------------------   
      //TAG + probe
   
      if ( (int)TnProbePairs.size() < 1 ) continue;

      int tag_i   = TnProbePairs[0].first;
      int probe_i = TnProbePairs[0].second;
      m12_mutag = (goodMuons[tag_i]->p4() + goodMuons[probe_i]->p4()).M() ;

       h1["m12_mutag"]  -> Fill( m12_mutag ); 
       h1["pt_mutag"]   -> Fill( goodMuons[tag_i]->pt() );
       h1["eta_mutag"]  -> Fill( goodMuons[tag_i]->eta() );
       
       double pt_probe = goodMuons[probe_i]->pt();

       //       if ( hltPath_ == "HLT_Mu8_v" )
       //	 {
       for ( auto & pt_bin : pt_bins )
	 {
	   if ( pt_probe > pt_bin.first && pt_probe < pt_bin.second ) h1[Form("m12_mutag_%.2fto%.2f", pt_bin.first, pt_bin.second)]-> Fill(m12_mutag); 
	 }

       //std::cout << "pT probe = "  << pt_probe << std::endl;
       ++nsel[3];
                                            
       //TAG + PASSING probe +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++                                                                                                                             
       if ( (int)psTnProbePairs.size() > 0 )
       {
	 int psprobe_i = psTnProbePairs[0].second;
	 m12_mutagnprobe = (goodMuons[tag_i]->p4() + goodMuons[psprobe_i]->p4()).M() ;
	 
	 h1["m12_mutagnprobe"]  -> Fill( m12_mutagnprobe );
	 h1["pt_muprobe"]       -> Fill( goodMuons[psprobe_i]->pt() );
	 h1["eta_muprobe"]      -> Fill( goodMuons[psprobe_i]->eta() );
	 
	 double pt_psprobe = goodMuons[psprobe_i]->pt();
       
	 for ( auto & pt_bin : pt_bins )
	   {
	     if ( pt_psprobe > pt_bin.first && pt_psprobe < pt_bin.second ) h1[Form("m12_mutagnprobe_%.2fto%.2f", pt_bin.first, pt_bin.second)]-> Fill(m12_mutagnprobe);
	   }
       }

       //TAG + FAILING probe +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    
       else 
	 {
	   int flprobe_i = TnProbePairs[0].second;
	   m12_mutagnprobe = (goodMuons[tag_i]->p4() + goodMuons[flprobe_i]->p4()).M() ;

	   h1["m12_mutagnprobefl"]  -> Fill( m12_mutagnprobe );
	   h1["pt_muprobefl"]       -> Fill( goodMuons[flprobe_i]->pt() );
	   h1["eta_muprobefl"]      -> Fill( goodMuons[flprobe_i]->eta() );

	   double pt_flprobe = goodMuons[flprobe_i]->pt();
	  
	   for ( auto & pt_bin : pt_bins ) {
	     if ( pt_flprobe > pt_bin.first && pt_flprobe < pt_bin.second ) h1[Form("m12_mutagnprobefl_%.2fto%.2f", pt_bin.first, pt_bin.second)]-> Fill(m12_mutagnprobe);
	     }	  
	 }

       //std::cout << "pT passingprobe = "  << pt_flprobe << std::endl;
       //std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ "  << pt_psprobe << std::endl;

       ++nsel[4];
                                
   } //end event loop

   //Fill Trigger Efficiency

   TAxis *xaxis = h1["m12_mutag"]->GetXaxis();                                                                                                                                                             
   int low  = xaxis->FindBin( mJPsi - 0.2 );                                                                                                                                                               
   int high = xaxis->FindBin( mJPsi + 0.2 );                                                                                                                                                               
   double num  = 0;                                                                                                                                                                                        
   double den  = 0;                                                                                                                                                                                        
   double errnum  = 0;                                                                                                                                                                                     
   double errden  = 0;                                                                                                                                                                     
                                                                                                                                                                                          
   for ( size_t i = 0 ; i< pt_bins.size(); i++ )                                                                                                                                                           
     {
       //Add fit gauss+pol1
       den = h1[Form("m12_mutag_%.2fto%.2f"      , pt_bins[i].first, pt_bins[i].second)]->IntegralAndError(low,high,errden);      
       num = h1[Form("m12_mutagnprobe_%.2fto%.2f", pt_bins[i].first, pt_bins[i].second)]->IntegralAndError(low,high,errnum);                                                                               
       h1["h_num"]->SetBinContent(i+1, num);
       h1["h_den"]->SetBinContent(i+1, den);
       h1["h_num"]->SetBinError(i+1, errnum);
       h1["h_den"]->SetBinError(i+1, errden); 
       std::cout << "#bin = " << i+1 << " passing = " << num << " +/- " << errnum << " total = " << den << " +/- " << errden <<std::endl;
       
     }                                                                                                                                                                                                     

   //   setTDRStyle();

   //Plot trig eff 
   TCanvas *c1 = new TCanvas("c1","c1");                                                                                                                                                                   
   TGraphAsymmErrors * g_eff = new TGraphAsymmErrors(h1["h_num"],h1["h_den"],"cl=0.683 b(1,1) mode");
   g_eff->SetMaximum( 1.0 );                              
   g_eff->SetTitle( hltPath_.c_str() );
   g_eff->GetXaxis()->SetTitle("pt_{probe} /GeV");                                                                                                                                                            
   g_eff->GetYaxis()->SetTitle("trigger efficiency");
   g_eff->Draw("AP"); 
   c1->SaveAs(("Trig_eff_"+hltPath_+".pdf").c_str());  

   //Histos
   for (auto & ih1 : h1)   {
      ih1.second -> Write();
      
   }
   
   hout.Write();
   hout.Close();
   
// PRINT OUTS  ============================================================================================ 
   
   // Cut flow
   // 0: triggered events
   // 1: muon id selection
   // 2: kinematics
   // 3: tag matched
   // 4: probe matched
   
   double fracAbs[10];
   double fracRel[10];
   std::string cuts[10];
   cuts[0] = "Triggered";
   cuts[1] = "Muon id sel";
   cuts[2] = "Muon kinematic sel";
   cuts[3] = "Tag muon matched";
   cuts[4] = "Probe muon matched";

   printf ("%-23s  %10s  %10s  %10s \n", std::string("Cut flow").c_str(), std::string("# events").c_str(), std::string("absolute").c_str(), std::string("relative").c_str() ); 
   for ( int i = 0; i < 5; ++i )
   {
      fracAbs[i] = double(nsel[i])/nsel[0];
      if ( i>0 )
         fracRel[i] = double(nsel[i])/nsel[i-1];
      else
         fracRel[i] = fracAbs[i];
      printf ("%-23s  %10d  %10.3f  %10.3f \n", cuts[i].c_str(), nsel[i], fracAbs[i], fracRel[i] ); 
   }
   
/*   
   // Trigger objects counts   
   std::cout << std::endl;
   printf ("%-40s  %10s \n", std::string("Trigger object mu").c_str(), std::string("# events").c_str() ); 
   for ( size_t io = 0; io < triggerObjectsMuons_.size() ; ++io )
    {
     printf ("%-40s  %10d \n", triggerObjectsMuons_[io].c_str(),  nmatch_mu[io] ); 
    }     
*/ 
  
   
} //end main

