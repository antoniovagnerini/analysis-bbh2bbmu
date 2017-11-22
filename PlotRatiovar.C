#include "tdrstyle.C"
#include <iostream>
#include <cstring>
#include <string>

using namespace std;
//int PlotRatiovar( string var_  )
int PlotRatiovar( string var_, string varunit_ , float xlow_ ,float xhigh_ , float ylow_ ,float yhigh_ , float yRlow_ ,float yRhigh_  )
{
  setTDRStyle();
  float relative_size = 0.49;

  std::string var = var_;
  std::string varunit = varunit_;
  // const std::string varlegend; 
  
  const float xlow  = xlow_;
  const float xhigh = xhigh_;
  const float ylow  = ylow_;
  const float yhigh = yhigh_;
  const float yRlow  = yRlow_;
  const float yRhigh = yRhigh_;


  TFile * f1 = new TFile( "ROOTFILES/histograms_2017C-v.root","OLD");
  TFile * f2 = new TFile( "ROOTFILES/histograms_2017D-v1.root","OLD");
  TFile * f3 = new TFile( "ROOTFILES/histograms_2017E-oldinj.root","OLD");
  TFile * f3bis = new TFile( "ROOTFILES/histograms_2017E-newinj.root","OLD");
  TFile * f4 = new TFile( "ROOTFILES/histograms_2017F-v1.root","OLD");

  TH1F * var_0C = (TH1F*)f1->Get(var.c_str());
  TH1F * var_0D = (TH1F*)f2->Get(var.c_str());
  TH1F * var_0E = (TH1F*)f3->Get(var.c_str());
  TH1F * var_0Ebis = (TH1F*)f3bis->Get(var.c_str());
  TH1F * var_0F = (TH1F*)f4->Get(var.c_str());

  std::cout << var << " "  << var.length() <<std::endl; 
  float legx;  
  if (var.length() < 10) legx = 0.715; //0.655
  else legx = 0.715 + (var.length()-10)*0.02; //0.655
  
  //TLegend * leg1 = new TLegend(0.1926503,0.550778,legx,0.8816568);
  TLegend * leg1 = new TLegend(0.3819599,0.6039126,legx,0.9347905);
  // TLegend * leg2 = new TLegend(0.7561247,0.7255245,0.9443207,0.9460731);  
  TLegend * leg2 = new TLegend(0.7120,0.6313878,0.9777283,0.9433835);
  TCanvas * c1 = new TCanvas("c1","",900,1000); //600,700
  
  //PAD1
  TPad* pad1 = new TPad("pad_top","Dijet m_{12} bbmunb distribution",0,1 - relative_size,1,1);
  // pad1->SetTopMargin(gStyle->GetPadTopMargin() * 1. / (relative_size));
  //pad1->SetBottomMargin(0.05); //0.03
  pad1->Range(0.5621891,1.838102,1.008955,5.356941);
  pad1->SetLeftMargin(0.08463252);
  pad1->SetRightMargin(0.02004454);
  pad1->SetBottomMargin(0.04600896);
  pad1->Draw();
  pad1->cd();

  TH1F *ratio =  (TH1F*)var_0D->Clone("ratio");
  TH1F *ratio1 = (TH1F*)var_0E->Clone("ratio1");
  TH1F *ratio1bis = (TH1F*)var_0Ebis->Clone("ratio1bis");
  TH1F *ratio2 = (TH1F*)var_0F->Clone("ratio2");

  var_0C -> GetYaxis()->SetTitle("entries");
  var_0C -> SetLineColor(kBlue+1); 
  var_0C -> GetXaxis() -> SetRangeUser(xlow,xhigh); // x axis range
  var_0C -> GetYaxis() -> SetRangeUser(ylow,yhigh); // y axis range
  var_0C -> DrawCopy("");

  var_0D -> SetLineColor(kRed);
  var_0D -> Draw("same");
  var_0E -> SetLineColor(kGreen+2);
  var_0E -> Draw("same");
  var_0Ebis -> SetLineColor(kGreen-1);
  var_0Ebis -> Draw("same");
  var_0F -> SetLineColor(kOrange);
  var_0F -> Draw("same");
   
  leg1->AddEntry(var_0C,(var + " Era2017C").c_str(),"lp");
  leg1->AddEntry(var_0D,(var + " Era2017D").c_str(),"lp");
  leg1->AddEntry(var_0E,(var + " Era2017E oldIS").c_str(),"lp");
  leg1->AddEntry(var_0Ebis,(var + " Era2017E newIS").c_str(),"lp");
  leg1->AddEntry(var_0F,(var + " Era2017F").c_str(),"lp");
  leg1->SetBorderSize(1);
  //leg1->SetFontSize(18);

  leg1 -> Draw();
  c1->cd();
  
  //Ratio plot
  TPad *pad2 = new TPad("pad2","pad2",0,0.05,1,relative_size);
  // pad2-> SetTopMargin(0.11); //0.04
  //pad2->SetBottomMargin(gStyle->GetPadBottomMargin() * 0.8 / relative_size);
  pad2->SetLeftMargin(0.08351893);
  pad2->SetRightMargin(0.02004454);
  pad2->SetTopMargin(0.04854763);
  pad2->SetBottomMargin(0.214497);
  pad2->Draw();
  pad2->cd();

  //Ratio plot 1 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  //RATIO OPTIONS
  ratio->Sumw2();
  ratio->Divide(var_0C);
  ratio->SetMarkerStyle(20);
  //ratio->SetMarkerSize(0.5);
  ratio->SetLineColor(kRed);
  ratio->SetMarkerColor(kRed);

  //RATIO X and Y axis style
  ratio -> GetXaxis() -> SetRangeUser(xlow,xhigh);
  ratio -> GetYaxis() -> SetRangeUser(yRlow,yRhigh);

 // Y axis ratio plot settings
  ratio->GetXaxis()->SetTitle( (var + "/" + varunit).c_str() );   // !!!!!!!!!!!!!!!!!!!!!!!!
  ratio->GetYaxis()->SetTitle("Era 2017X / Era 2017C");
  ratio->GetYaxis()->SetNdivisions(505);
  ratio->GetYaxis()->SetTitleSize(21); //22
  ratio->GetYaxis()->SetTitleFont(43);
  ratio->GetYaxis()->SetTitleOffset(1.55);
  ratio->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  ratio->GetYaxis()->SetLabelSize(20);

  // X axis ratio plot settings
  ratio->GetXaxis()->SetTitleSize(21); //22
  ratio->GetXaxis()->SetTitleFont(43);
  ratio->GetXaxis()->SetTitleOffset(4.);
  ratio->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  ratio->GetXaxis()->SetLabelSize(20);

  ratio->DrawCopy("");

  //Ratio plot 1 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++=                                                                                                                     
  ratio1->Sumw2();
  ratio1->Divide(var_0C);
  ratio1->SetMarkerStyle(21);
  ratio1-> SetLineColor(kGreen+2);
  ratio1-> SetMarkerColor(kGreen+2);

  ratio1->GetYaxis()->SetNdivisions(505);
  ratio1->GetYaxis()->SetTitleSize(21);                                                                                                                                                           
  ratio1->GetYaxis()->SetTitleFont(43);
  ratio1->GetYaxis()->SetTitleOffset(1.55);
  ratio1->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)                                                                                                                
  ratio1->GetYaxis()->SetLabelSize(20);
  ratio1->GetXaxis()->SetTitleSize(21); //22                                                                                                                                                      
  ratio1->GetXaxis()->SetTitleFont(43);
  ratio1->GetXaxis()->SetTitleOffset(4.);
  ratio1->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)                                                                                                                         
  ratio1->GetXaxis()->SetLabelSize(20);
  ratio1->Draw("same");
  fixOverlay();

  //Ratio plot 1bis +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++=                                                                                                                     
  ratio1bis->Sumw2();
  ratio1bis->Divide(var_0C);
  ratio1bis->SetMarkerStyle(22);
  ratio1bis-> SetLineColor(kGreen-1);
  ratio1bis-> SetMarkerColor(kGreen-1);

  ratio1bis->GetYaxis()->SetNdivisions(505);
  ratio1bis->GetYaxis()->SetTitleSize(21);                                                                                                                                                           
  ratio1bis->GetYaxis()->SetTitleFont(43);
  ratio1bis->GetYaxis()->SetTitleOffset(1.55);
  ratio1bis->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)                                                                                                                
  ratio1bis->GetYaxis()->SetLabelSize(20);
  ratio1bis->GetXaxis()->SetTitleSize(21); //22                                                                                                                                                      
  ratio1bis->GetXaxis()->SetTitleFont(43);
  ratio1bis->GetXaxis()->SetTitleOffset(4.);
  ratio1bis->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)                                                                                                                         
  ratio1bis->GetXaxis()->SetLabelSize(20);
  ratio1bis->Draw("same");
  fixOverlay();

  //Ratio plot 2 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++=                                                                                                                     
  ratio2->Sumw2();
  ratio2->Divide(var_0C);
  ratio2->SetMarkerStyle(23);
  ratio2-> SetLineColor(kOrange);
  ratio2-> SetMarkerColor(kOrange);

  ratio2->GetYaxis()->SetNdivisions(505);
  ratio2->GetYaxis()->SetTitleSize(21);                                                                                                                                                           
  ratio2->GetYaxis()->SetTitleFont(43);
  ratio2->GetYaxis()->SetTitleOffset(1.55);
  ratio2->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)                                                                                                                
  ratio2->GetYaxis()->SetLabelSize(20);
  ratio2->GetXaxis()->SetTitleSize(21); //22                                                                                                                                                      
  ratio2->GetXaxis()->SetTitleFont(43);
  ratio2->GetXaxis()->SetTitleOffset(4.);
  ratio2->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)                                                                                                                         
  ratio2->GetXaxis()->SetLabelSize(20);
  ratio2->Draw("same");
  fixOverlay();

  //LEGEND
  leg2->AddEntry(ratio,  "EraD / EraC","pl");
  leg2->AddEntry(ratio1, "EraE oldIS/ EraC","pl");
  leg2->AddEntry(ratio1bis, "EraE newIS/ EraC","pl");
  leg2->AddEntry(ratio2, "EraF / EraC","pl");
  leg2->SetBorderSize(1);
  leg2->SetTextSize(0.05);
  leg2 -> Draw();
  c1->cd();
 
  c1 -> SaveAs(("PLOTS/"+ var + ".pdf").c_str());
  c1 -> SaveAs(("PLOTS/"+ var + ".png").c_str());

  pad1->cd();
  pad1 -> SetLogy();
  c1->cd();
  c1 -> SaveAs(("PLOTS/"+ var + "Log.pdf").c_str());
  c1 -> SaveAs(("PLOTS/"+ var + "Log.png").c_str());  
   
  return 0;
}

