/*
 * Histograms.h
 *
 *  Created on: 10 jan. 2017 .
 *      Author: vagnerini
 */

#ifndef ANALYSIS_SEMILEP_HISTOGRAMS_H_
#define ANALYSIS_SEMILEP_HISTOGRAMS_H_

#include <iostream>		// standard in/out
#include <memory> 		// for std::shared_ptr
#include <string>
#include <map>
#include <array>
#include "/usr/include/root/TMath.h"

#include "/usr/include/root/TH1.h"
#include "/usr/include/root/TH2.h"
#include "/usr/include/root/TEfficiency.h"

#include "stdlib.h"

class Histograms {
	public:
		Histograms();
		virtual ~Histograms();

		//Initialise
		void Make(const int &size = 100);

		//Gets
		std::map<std::string, TH1*>& getHisto();
		std::map<std::string, TH2*>& getHisto2D();

	protected:

		int size_;
		std::map<std::string,TH1* > histo_;
		std::map<std::string,TH2* > histo2D_;

};

#endif /* ANALYSIS_SEMILEP_HISTOGRAMS_H_ */
