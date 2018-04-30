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
#include "TF1.h"
#include "TFormula.h"

#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TAxis.h"

#include "tdrstyle.C"

#include "Fit/Fitter.h"
#include "Fit/BinData.h"
#include "Fit/Chi2FCN.h"
#include "TH1.h"
#include "TList.h"
#include "Math/WrappedMultiTF1.h"
#include "HFitInterface.h"
#include "TCanvas.h"
#include "TStyle.h"


using namespace std;

static const double mJPsi = 3.10; //pdg value Y(1s) resonance


// =============================================================================================   

   // definition of shared parameter
// Matched Probe
int iparM[5] = { 0,      // amplitude (gaus)
		 1,      // mean      (gaus) 
		 2,      // sigma     (gaus)
		 3,      // trig eff  common
		 4,      //offset1
	   };

// Unmatched Probe
int iparU[5] = { 0,      // amplitude (gaus)                                                       
		 1,      // mean      (gaus)
		 2,      // sigma     (gaus)                                                
		 3,      // trig eff  common
		 //		 5,      // offset2 
		 4,
               };


struct GlobalChi2 {
  GlobalChi2(  ROOT::Math::IMultiGenFunction & f1,
	       ROOT::Math::IMultiGenFunction & f2) :
    fChi2_1(&f1), fChi2_2(&f2) {}
  
  // parameter vectors
  
  double operator() (const double *par) const {
    double p1[5];
    for (int i = 0; i < 5; ++i) p1[i] = par[iparM[i] ];
    
    double p2[5];
    for (int i = 0; i < 5; ++i) p2[i] = par[iparU[i] ];
    
    return (*fChi2_1)(p1) + (*fChi2_2)(p2);
  }
  
  const  ROOT::Math::IMultiGenFunction * fChi2_1;
  const  ROOT::Math::IMultiGenFunction * fChi2_2;
};

// =============================================================================================   
void Efficiency( string hltPath_, string era_, double halfwidth_ )
{
     
   TH1::SetDefaultSumw2();  // proper treatment of errors when scaling histograms
   std::map<std::string, TH1F*> h1;

   std::vector<std::pair<float,float> > pt_bins;
   if ( hltPath_ == "HLT_Mu8_v" ) pt_bins = {{7.75,8.25},{8.25,9},{9,11},{11,15},{15,50}};
   else if (hltPath_=="HLT_SingleJet30_Mu12_SinglePFJet40_v") pt_bins = {{11,13},{13,15},{15,50}};


   std::vector<std::string> eras = { era_ };
   
   TFile * f[10]; 

   for ( size_t i = 0 ; i < eras.size(); i++ ) 
     {
       if      ( hltPath_ == "HLT_Mu8_v" )                        f[i]= new TFile( ("ROOTFILES/DOUBLEMUON/histograms_2017"+eras[i]+".root").c_str(),"OLD");
       else if (hltPath_=="HLT_SingleJet30_Mu12_SinglePFJet40_v") f[i]= new TFile( ("ROOTFILES/BTAGCSV/histograms_2017"+eras[i]+".root").c_str(),"OLD");

        
       h1["m12_mutag"]       = (TH1F*) f[i]->Get("m12_mutag"); 
       h1["m12_mutagnprobe"] = (TH1F*) f[i]->Get("m12_mutagnprobe");

       
       //Edges variable bin size
       float edges[ pt_bins.size() +1 ]; //Nedges = Nbins+1 
       for ( size_t i = 0; i< pt_bins.size(); i++ )
	 {                                                                                                            
	   edges[i] = pt_bins[i].first;  //store lower bin edge
	   std::cout<< "bin = " << i << " lowedge =" << edges[i]<< std::endl;
	 }
       edges[pt_bins.size()] = pt_bins[pt_bins.size()-1].second; //store pt max as upper edge 

       h1["h_trig"] = new TH1F("h_num","",pt_bins.size(),edges);
       

       for ( auto & pt_bin : pt_bins ) 
	 {
	   h1[Form("m12_mutag_%.2fto%.2f"        , pt_bin.first, pt_bin.second)] = (TH1F*) f[i] ->Get(Form("m12_mutag_%.2fto%.2f"        , pt_bin.first, pt_bin.second)) ;
	   h1[Form("m12_mutagnprobe_%.2fto%.2f"  , pt_bin.first, pt_bin.second)] = (TH1F*) f[i] ->Get(Form("m12_mutagnprobe_%.2fto%.2f"  , pt_bin.first, pt_bin.second)) ;	   
           h1[Form("m12_mutagnprobefl_%.2fto%.2f", pt_bin.first, pt_bin.second)] = (TH1F*) f[i] ->Get(Form("m12_mutagnprobefl_%.2fto%.2f", pt_bin.first, pt_bin.second)) ;
	 }

       //Fill Trigger Efficiency
       
       const double hw =  halfwidth_;             
       TAxis *xaxis = h1["m12_mutag"]->GetXaxis();                                                                                                                                                    
       int nbins_fit  = xaxis->FindBin( mJPsi + hw ) - xaxis->FindBin( mJPsi - hw ) ;   
                                                                                                                                            
       for ( size_t i = 0 ; i< pt_bins.size(); i++ )                                                                                                                                                       
    
	 {
	   cout<< " ==============================================="<<endl;
	   cout<< " pT bin " << (pt_bins[i].first+pt_bins[i].second)/2 << endl;

	   //histos 

	   //	   TH1D * hM = new TH1D("hM","histo B",100,0,100);
	   //	   TH1D * hU = new TH1D("hU","histo S+B",100, 0,100);
	   string hps_ptbin = Form("m12_mutagnprobe_%.2fto%.2f"   , pt_bins[i].first, pt_bins[i].second);
           string hfl_ptbin = Form("m12_mutagnprobefl_%.2fto%.2f" , pt_bins[i].first, pt_bins[i].second);

	   TF1 * gausM = new TF1("gausM","    [3]*(gaus(0) + [4])", mJPsi -hw , mJPsi +hw ); // *[4]
	   //  TF1 * gausU = new TF1("gausU","(1-[3])*gaus(0) + [5]", mJPsi -hw , mJPsi +hw );
	   TF1 * gausU = new TF1("gausU","(1-[3])*(gaus(0) + [4])", mJPsi -hw , mJPsi +hw );

	   gausM->SetParameters( 10., mJPsi , 0.5*hw , 0.5 , -1.); //-0.1
	   /*    gausM->SetParLimits(0, 0.01, 1000. );
           gausM->SetParLimits(1, mJPsi - hw , mJPsi + hw );
           gausM->SetParLimits(2, 0.2*hw, 2*hw );
           gausM->SetParLimits(3, 0.1, 1.0 ); 
           gausM->SetParLimits(4, -20., 20.);
	   */
           gausU->SetParameters( 10., mJPsi , 0.5*hw , 0.5 , -1.);
	   /*  gausU->SetParLimits(0, 0.01, 1000. );
           gausU->SetParLimits(1, mJPsi - hw , mJPsi + hw );
           gausU->SetParLimits(2, 0.2*hw, 2*hw );
           gausU->SetParLimits(3, 0.1, 1.0 );
           gausU->SetParLimits(5, -20., 20.);
	   */

	   //Fit
	   ROOT::Math::WrappedMultiTF1 wgausM (*gausM,1);
	   ROOT::Math::WrappedMultiTF1 wgausU (*gausU,1);

	   ROOT::Fit::DataOptions opt;
	   ROOT::Fit::DataRange rangeM;
	   // set the data range
	   rangeM.SetRange( mJPsi -hw , mJPsi +hw );
	   ROOT::Fit::BinData dataM(opt,rangeM);
	   ROOT::Fit::FillData(dataM, h1[hps_ptbin]);

	   ROOT::Fit::DataRange rangeU;
	   rangeU.SetRange( mJPsi -hw , mJPsi +hw );
	   ROOT::Fit::BinData dataU(opt,rangeU);
	   ROOT::Fit::FillData(dataU, h1[hfl_ptbin]);

	   ROOT::Fit::Chi2Function chi2_M (dataM, wgausM);
	   ROOT::Fit::Chi2Function chi2_U (dataU, wgausU);

	   GlobalChi2 globalChi2(chi2_M, chi2_U);

	   ROOT::Fit::Fitter fitter;
	   const int Npar = 5; //6
	   double par0[Npar] = { 10., mJPsi , 0.5*hw , 0.5 , -1. }; // add -0.1

	   // create before the parameter settings in order to fix or set range on them
	   fitter.Config().SetParamsSettings(5,par0); //6

	   //	   fitter.Config().ParSettings(1).Fix(); // fix m0
	   fitter.Config().ParSettings(1).SetLimits(mJPsi - hw , mJPsi + hw);
	   fitter.Config().ParSettings(2).SetLimits( 0.2*hw, 2*hw );
	   fitter.Config().ParSettings(3).SetLimits( 0.1, 1.0 );
	   fitter.Config().ParSettings(4).SetLimits( -20., 20. );
	   //	   fitter.Config().ParSettings(5).SetLimits( -40., 40. );

	   fitter.Config().MinimizerOptions().SetPrintLevel(0);
	   fitter.Config().SetMinimizer("Minuit2","Migrad");

	   // fit FCN function directly
	   // (specify optionally data size and flag to indicate that is a chi2 fit)

	   fitter.FitFCN(Npar,globalChi2,0,dataM.Size()+dataU.Size(),true); //6 before
	   ROOT::Fit::FitResult result = fitter.Result();
	   result.Print(std::cout);

	   //Trigger results no BKG subtraction
 
	   //cout << "trigeff= " <<  << " +/- " <<  << endl;

	   setTDRStyle();

	   TCanvas * c1 = new TCanvas("Simfit","Simultaneous fit of two histograms",
                             10,10,700,700);
	   c1->Divide(1,2);
	   c1->cd(1);
	   gStyle->SetOptFit(1111);

	   gausM->SetFitResult( result, iparM);
	   gausM->SetRange(rangeM().first, rangeM().second);
	   gausM->SetLineColor(kBlue);
	   h1[hps_ptbin]->GetListOfFunctions()->Add(gausM);
	   h1[hps_ptbin]->GetXaxis()->SetRangeUser( mJPsi - hw -0.1,  mJPsi + hw + 0.1);
	   h1[hps_ptbin]->Draw();
	   
	   c1->cd(2);
	   gausU->SetFitResult( result, iparU);
	   gausU->SetRange(rangeU().first, rangeU().second);
	   gausU->SetLineColor(kRed);
	   h1[hfl_ptbin]->GetListOfFunctions()->Add(gausU);
	   h1[hfl_ptbin]->GetXaxis()->SetRangeUser( mJPsi - hw -0.1,  mJPsi + hw + 0.1);
	   h1[hfl_ptbin]->Draw();

	   c1->SaveAs(("PLOTS/"+hps_ptbin+era_+".png").c_str());
	   
	   const std::vector<double> * pars = & result.Parameters();   
	   const std::vector<double> * errs = & result.Errors();

	   h1["h_trig"]->SetBinContent(i+1, pars->at(3) );
           h1["h_trig"]->SetBinError  (i+1, errs->at(3) );

       	 }                      
       
     }
    //Plot trig eff
   
   setTDRStyle();
      
   //Plot trig eff 
   TCanvas *c = new TCanvas("c","c");                                                                                                                                                                   
   h1["h_trig"]->SetMaximum( 1.0 );                              
   h1["h_trig"]->SetTitle( hltPath_.c_str() );
   h1["h_trig"]->GetXaxis()->SetTitle("pt_{probe} /GeV");                                                                                                                                                            
   h1["h_trig"]->GetYaxis()->SetTitle("trigger efficiency");
   h1["h_trig"]->Draw(); 
   c->SaveAs(("PLOTS/"+hltPath_+era_+".png").c_str());  
     
     
   
}

