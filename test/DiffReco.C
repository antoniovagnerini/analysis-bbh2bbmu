#include "tdrstyle.C"
#include <iostream>
#include <cstring>
#include <string>

#include "TH1.h"
#include "TH2.h"
#include "TGaxis.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "Fit/Fitter.h"
#include "TStyle.h"
#include <iostream>

#include "HbbStylesNew.cc"

using namespace std;

int DiffReco( string era_, string var_, string varunit_ , string xtitle_, float xlow_ ,float xhigh_ , float ylow_ ,float yhigh_ , float yRlow_ ,float yRhigh_  )
{

  HbbStylesNew style;
  style.SetStyle();
  gStyle->SetOptStat(0);

  setTDRStyle();

  std::string var = var_;
  std::string varunit = varunit_;
  std::string xtitle = xtitle_;
  
  const float xlow  = xlow_;
  const float xhigh = xhigh_;
  const float ylow  = ylow_;
  const float yhigh = yhigh_;
  const float yRlow  = yRlow_;
  const float yRhigh = yRhigh_;

  float lumi = 36.5; // all 2017 C to F

  if ( era_ == "C" ) lumi = 9.58;
  if ( era_ == "D" ) lumi = 4.22;
  if ( era_ == "E" ) lumi = 9.26;
  if ( era_ == "F" ) lumi = 13.5;
  if ( era_ == "CtoE" ) lumi = 23.1;
 

  TFile * f1 = new TFile( ("ROOTFILES/PROMPTRECO/histograms_2017"+era_+".root").c_str(),"OLD");
  TFile * f2 = new TFile( ("ROOTFILES/RERECO/histograms_2017"+era_+".root").c_str(),"OLD");

  TH1F * hist_prompt = (TH1F*)f1->Get(var.c_str());
  TH1F * hist_rereco = (TH1F*)f2->Get(var.c_str());

  std::vector<double> lines1 = {1.0, 1.5, 2.0};
  std::vector<double> linesF = {1.0, 2.0, 3.0, 4.0};

  TCanvas* c1 = style.MakeCanvas("c1","",700,700);
  //c1 -> SetLogy();
  style.InitHist(hist_prompt, xtitle.c_str(),"Entries / 20 GeV ",kBlack,0);
  style.InitHist(hist_rereco, xtitle.c_str(),"Entries / 20 GeV ",kRed,0);

  //Ratio plot
  auto rp = new TRatioPlot(hist_rereco,hist_prompt);
  rp -> SetH1DrawOpt("E");
  rp -> SetH2DrawOpt("E");
  rp -> SetGridlines(lines1);
  rp -> Draw();
  
  rp -> GetLowerRefGraph() -> SetLineWidth(2);
  rp -> GetLowerRefGraph() -> SetMarkerSize(1);
  
  rp -> SetRightMargin(0.05);
  rp -> SetUpTopMargin(0.055); // 0.03
  rp -> SetLeftMargin(0.13);
  rp -> SetSeparationMargin(0.07);
  rp -> SetLowBottomMargin(0.4); 
  rp -> SetLowTopMargin(0.0);
  
  rp -> GetLowerRefYaxis() -> SetTitle("ReReco / Prompt");
  rp -> GetLowerRefYaxis() -> SetRangeUser(yRlow,yRhigh);
  rp -> GetUpperRefYaxis() -> SetRangeUser(ylow,yhigh);
  rp -> GetUpperRefXaxis() -> SetRangeUser(xlow,xhigh);

  TGaxis::SetMaxDigits(3);
    
  rp -> GetLowerRefYaxis() -> SetTitleSize(0.035);
  rp -> GetLowerRefYaxis() -> SetTitleOffset(1.7);
  //  rp -> GetUpperRefYaxis() -> SetTitleSize(0.05);  
  rp -> GetUpperRefYaxis() -> SetTitleOffset(1.7); //1.3
  rp -> GetUpperRefYaxis() -> SetTitleSize(0.035);
  rp -> GetLowerRefXaxis() -> SetTitleSize(0.05);
  rp -> GetLowerRefXaxis() -> SetTitleOffset(1.15);
  rp -> GetLowerRefXaxis() -> SetLabelSize(0.035);
  rp -> GetLowerRefXaxis() -> SetLabelOffset(0.02);
  rp -> GetLowerRefYaxis() -> SetLabelSize(0.035);
  rp -> GetLowerRefYaxis() -> SetLabelOffset(0.02);
  rp -> GetUpperRefYaxis() -> SetLabelSize(0.035);
  rp -> GetUpperRefYaxis() -> SetLabelOffset(0.02);

  rp -> GetLowYaxis() -> SetNdivisions(505);
  rp -> GetUpperPad() -> cd();
  
  TLegend* leg = new TLegend(0.58,0.63,0.98,0.93);
  style.SetLegendStyle(leg);
  leg -> AddEntry(hist_prompt,("Era "+era_+ ", Prompt").c_str(),"L");
  leg -> AddEntry(hist_rereco,("Era "+era_+ ", ReReco").c_str(),"L");
  leg -> Draw("SAME");
  
  CMSPrelim( Form("%.1f fb^{-1} (13 TeV)", lumi ) , 0.15, 0.78);

  c1 -> Update();
  c1 -> SaveAs(("PLOTS/CompReco_"+ var_+ era_ + ".png").c_str());

  //pad1 -> SetLogy();
  //c1 -> SaveAs(("PLOTS/CompReco"+ era_ + "Log.png").c_str());

  return 0;
}

