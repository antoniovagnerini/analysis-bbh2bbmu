/*
 * Histograms.cpp
 *
 *  Created on: 11 jan. 2017 г.
 *      Author: vagnerini
 */

#include "Analysis/Semilep/interface/Histograms.h"

Histograms::Histograms() : size_(0) {
// TODO Auto-generated constructor stub
}

Histograms::~Histograms() {
  //	std::cout<<"I'm at ~Histograms"<<std::endl;
	// TODO Auto-generated destructor stub
	//Clean memory:
	histo_.clear();
	histo2D_.clear();
}

//Initialise histos
void Histograms::Make(const int &size /*, const bool & lowM*/) 
{
  size_ = size;
  TH1::SetDefaultSumw1();
  TH1::SetDefaultSumw2(); //add TH2
  
  //Declare All Basic Histograms
  histo_["TotalNumberOfGenEvents"]            = new TH1D("TotalNumberOfGenEvents",           "Total number of generated events",1,0,5.e+08);
  histo_["NumberOfGenEvents_afterMHat"]       = new TH1D("NumberOfGenEvents_afterMHat",      "Total number of generated events after mHat cut",1,0,5.e+08);
  histo_["NumberOfGenEvents_afterMHat_rewPU"] = new TH1D("NumberOfGenEvents_afterMHat_rewPU","Total number of generated events after mHat cut and PU reweighting",1,0,5.e+08);

}


//Gets
std::map<std::string, TH1* >& Histograms::getHisto() {
	return histo_;
}

std::map<std::string, TH2*>& Histograms::getHisto2D() {
	return histo2D_;
}
