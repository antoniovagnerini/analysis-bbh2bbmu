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
#include "Fit/PoissonLikelihoodFCN.h"
#include "TList.h"
#include "Math/WrappedMultiTF1.h"
#include "HFitInterface.h"
#include "TStyle.h"


using namespace std;

static const double mJPsi = 3.0969; //pdg value Y(1s) resonance

/* useful typedefs

typedef Chi2FCN< ROOT::Math::IMultiGenFunction > 	ROOT::Fit::Chi2Function
typedef PoissonLikelihoodFCN<ROOT::Math::IMultiGenFunction, ROOT::Math::IParamMultiFunction> PoissonLLFunction;
typedef PoissonLikelihoodFCN<ROOT::Math::IMultiGradFunction, ROOT::Math::IParamMultiFunction> PoissonLLGradFunction;

*/

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

//Global Chi2
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

//Global Likelihood
struct GlobalLikelihood {
  GlobalLikelihood(  ROOT::Math::IMultiGenFunction & f1,
		     ROOT::Math::IMultiGenFunction & f2) :
    fLikelihood_1(&f1), fLikelihood_2(&f2) {}
  
  // parameter vectors
  
  double operator() (const double *par) const {
    double p1[5];
    for (int i = 0; i < 5; ++i) p1[i] = par[iparM[i] ];
    
    double p2[5];
    for (int i = 0; i < 5; ++i) p2[i] = par[iparU[i] ];
    
    return (*fLikelihood_1)(p1) + (*fLikelihood_2)(p2);
  }
  
  const  ROOT::Math::IMultiGenFunction * fLikelihood_1;
  const  ROOT::Math::IMultiGenFunction * fLikelihood_2;
};

// =============================================================================================   
void Efficiency( string hltPath_, string era_, double halfwidth_, string estimator_ )
{
   setTDRStyle();

   TH1::SetDefaultSumw2();  // proper treatment of errors when scaling histograms
   std::map<std::string, TH1F*> h1;

   std::vector<std::pair<float,float> > pt_bins;
   if ( hltPath_ == "HLT_Mu8_v" ) pt_bins = {{7.75,8.25},{8.25,9},{9,11},{11,30}};
   else if (hltPath_=="HLT_SingleJet30_Mu12_SinglePFJet40_v") pt_bins = {{11,13},{13,15},{15,50}};


   std::vector<std::string> eras = { era_ };
   
   TFile * f[10]; 

   for ( size_t i = 0 ; i < eras.size(); i++ ) 
     {
       //       if      ( hltPath_ == "HLT_Mu8_v" )                        f[i]= new TFile( ("ROOTFILES/DOUBLEMUON/histograms_2017"+eras[i]+".root").c_str(),"OLD");
       if      ( hltPath_ == "HLT_Mu8_v" )                        f[i]= new TFile( ("ROOTFILES/DOUBLEMUON/RERECO/histograms_2017"+eras[i]+".root").c_str(),"OLD");    
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

	   string hps_ptbin = Form("m12_mutagnprobe_%.2fto%.2f"   , pt_bins[i].first, pt_bins[i].second);
           string hfl_ptbin = Form("m12_mutagnprobefl_%.2fto%.2f" , pt_bins[i].first, pt_bins[i].second);

	   TF1 * gausM = new TF1("gausM","    [3]*(gaus(0) + [4])", mJPsi -hw , mJPsi +hw ); // *[4]
	   //  TF1 * gausU = new TF1("gausU","(1-[3])*gaus(0) + [5]", mJPsi -hw , mJPsi +hw );
	   TF1 * gausU = new TF1("gausU","(1-[3])*(gaus(0) + [4])", mJPsi -hw , mJPsi +hw );

	   gausM->SetParameters( 10., mJPsi , 0.5*hw , 0.5 , -1.); //-0.1
           gausU->SetParameters( 10., mJPsi , 0.5*hw , 0.5 , -1.);

	   //Fit options
	   ROOT::Math::WrappedMultiTF1 wgausM (*gausM,1);
	   ROOT::Math::WrappedMultiTF1 wgausU (*gausU,1);

	   ROOT::Fit::DataOptions opt;
	   ROOT::Fit::DataRange rangeM;
	   rangeM.SetRange( mJPsi -hw , mJPsi +hw );
	   ROOT::Fit::BinData dataM(opt,rangeM);
	   ROOT::Fit::FillData(dataM, h1[hps_ptbin]);

	   ROOT::Fit::DataRange rangeU;
	   rangeU.SetRange( mJPsi -hw , mJPsi +hw );
	   ROOT::Fit::BinData dataU(opt,rangeU);
	   ROOT::Fit::FillData(dataU, h1[hfl_ptbin]);

	   //Chi2 method
	   ROOT::Fit::Chi2Function chi2_M (dataM, wgausM);
	   ROOT::Fit::Chi2Function chi2_U (dataU, wgausU);
	   GlobalChi2 globalChi2(chi2_M, chi2_U);
	  	   
	   //Likelihood method
	   ROOT::Fit::PoissonLLFunction Likelihood_M (dataM, wgausM);
	   ROOT::Fit::PoissonLLFunction Likelihood_U (dataU, wgausU);
	   GlobalLikelihood globalLikelihood(Likelihood_M, Likelihood_U) ;
	   	   
	   ROOT::Fit::Fitter fitter;
	   const int Npar = 5; //6
	   double par0[Npar] = { 10., mJPsi , 0.5*hw , 0.5 , 2. }; // add -0.1

	   // Parameter settings 
	   fitter.Config().SetParamsSettings(5,par0); //6

	   fitter.Config().ParSettings(1).Fix(); // fix m0
	   //fitter.Config().ParSettings(1).SetLimits(mJPsi - hw , mJPsi + hw);
	   fitter.Config().ParSettings(2).SetLimits( 0.2*hw, 2*hw );
	   fitter.Config().ParSettings(3).SetLimits( 0.1, 1.0 );
	   fitter.Config().ParSettings(4).SetLimits( -20., 20. );
	   //	   fitter.Config().ParSettings(5).SetLimits( -40., 40. );

	   fitter.Config().MinimizerOptions().SetPrintLevel(0);
	   fitter.Config().SetMinimizer("Minuit2","Migrad");

	   // Simultaneous Fit 
	   // (specify optionally data size and flag to indicate that is a chi2 fit)

	   if ( estimator_ == "Chi2" )      fitter.FitFCN(Npar,globalChi2      ,0,dataM.Size()+dataU.Size(),true);
	   else if ( estimator_ == "LL" )   fitter.FitFCN(Npar,globalLikelihood,0,dataM.Size()+dataU.Size(),false); 
	    
	   ROOT::Fit::FitResult result = fitter.Result();
	   result.Print(std::cout);

	   //PLOT SIM FIT
	   
	   TCanvas * c1 = new TCanvas("Simfit","Simultaneous fit of two histograms", 1400,1000);
	   c1->Divide(2,1);
	   c1->cd(1);
	   //	   gStyle->SetOptFit(1111);

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

	   c1->SaveAs(("PLOTS/"+hps_ptbin+era_+estimator_+".png").c_str());
	   
	   const std::vector<double> * pars = & result.Parameters();   
	   const std::vector<double> * errs = & result.Errors();

	   h1["h_trig"]->SetBinContent(i+1, pars->at(3) );
           h1["h_trig"]->SetBinError  (i+1, errs->at(3) );

       	 }                      
       
     }
    //Plot trig eff
 
      
   //Plot trig eff 
   TCanvas *c = new TCanvas("c","c");    
   h1["h_trig"]->GetXaxis()->SetRangeUser( 0.,30.);                                                                                                                          
   h1["h_trig"]->SetMaximum( 1.0 );                              
   h1["h_trig"]->SetTitle( hltPath_.c_str() );
   h1["h_trig"]->GetXaxis()->SetTitle("p_{Tprobe} [GeV]");                                                                                                                                                            
   h1["h_trig"]->GetYaxis()->SetTitle("trigger efficiency");
   h1["h_trig"]->Draw(); 
   c->SaveAs(("PLOTS/"+hltPath_+era_+estimator_+".png").c_str());  
          
   
}

