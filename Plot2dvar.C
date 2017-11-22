#include "tdrstyle.C"

void Plot2dvar(string var_)
{

  // std::string var = "n_ptmin20_btagsel_phivseta";
  std::string var = var_;
  //  setTDRStyle();

  gStyle->SetOptStat(0);

  std::vector<std::string> eras;
  eras.push_back("C-v");
  eras.push_back("D-v1");
  eras.push_back("E-oldinj");
  eras.push_back("E-newinj");
  eras.push_back("F-v1");


  TFile * f[20]; // signal

  for ( size_t e = 0; e < eras.size(); e++  )
   {
     f[e] = new TFile(Form("ROOTFILES/histograms_2017%s.root",eras[e].c_str()),"OLD");
  
     TH2F * hvar = (TH2F*) f[e]->Get(var.c_str());  
     std::cout << var << " "  << var.length() <<std::endl; 

     TCanvas *c1 = new TCanvas("c1", "c1",2067,115,1310,869); //2067,113
     c1->SetFillColor(0);
     c1->SetBorderMode(0);
     c1->SetBorderSize(2);
     c1->SetTheta(89.88808);
     c1->SetPhi(-0.06576451);
     c1->SetFrameBorderMode(0);
   
     hvar -> GetXaxis()->SetTitle("Eta");
     hvar -> GetYaxis()->SetTitle("Phi");
     hvar -> GetYaxis()->SetTitleOffset(1.05);
     hvar -> GetXaxis()->SetRangeUser(-2.3,2.3);
     hvar -> GetYaxis()->SetRangeUser(-3.,3.);
     hvar ->Draw("colz");
   

     c1 -> SaveAs(("PLOTS/"+ var + eras[e] + ".png").c_str());

     f[e]-> Close();
   }
   
  

}
