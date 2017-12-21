#!/bin/csh -f

# PlotRatio var( varname, varunit ,xlow ,xhigh ,ylow ,yhigh ,yRlow ,yRhigh ) 
#deepCSV
root -b -q -l 'PlotRatiovar.C("deepcsv_0","",0.0,1.0,100.,30000.,0.25,1.0)'
root -b -q -l 'PlotRatiovar.C("deepcsv_1","",0.0,1.0,100.,30000.,0.25,1.0)'
root -b -q -l 'PlotRatiovar.C("deepcsv_2","",0.0,1.0,100.,30000.,0.25,1.0)'

root -b -q -l 'PlotRatiovar.C("deepcsv_0_btagsel","",0.6,1.0,100.,30000.,0.25,1.0)'
root -b -q -l 'PlotRatiovar.C("deepcsv_1_btagsel","",0.6,1.0,100.,30000.,0.25,1.0)'
root -b -q -l 'PlotRatiovar.C("deepcsv_2_btagsel","",0.0,0.255,100.,30000.,0.25,1.0)'

#pt
root -b -q -l 'PlotRatiovar.C("pt_0","GeV",0.,1000.,10.,190000.,0.25,1.0)'
root -b -q -l 'PlotRatiovar.C("pt_1","GeV",0.,1000.,10.,190000.,0.25,1.0)'
root -b -q -l 'PlotRatiovar.C("pt_2","GeV",0.,1000.,10.,190000.,0.25,1.0)'

#pt for different eta ranges
root -b -q -l 'PlotRatiovar.C("pt_0","GeV",0.,1000.,10.,190000.,0.25,1.0)'
root -b -q -l 'PlotRatiovar.C("pt_0_eta0p9","GeV",0.,1000.,10.,190000.,0.25,1.0)'
root -b -q -l 'PlotRatiovar.C("pt_0_0p9eta1p5","GeV",0.,1000.,10.,190000.,0.25,1.0)'
root -b -q -l 'PlotRatiovar.C("pt_0_1p5eta2p0","GeV",0.,1000.,10.,190000.,0.25,1.0)'
root -b -q -l 'PlotRatiovar.C("pt_0_2p0eta","GeV",0.,1000.,10.,190000.,0.25,1.0)'

root -b -q -l 'PlotRatiovar.C("pt_1","GeV",0.,1000.,10.,400000.,0.25,1.0)'
root -b -q -l 'PlotRatiovar.C("pt_1_eta0p9","GeV",0.,1000.,10.,400000.,0.25,1.0)'
root -b -q -l 'PlotRatiovar.C("pt_1_0p9eta1p5","GeV",0.,1000.,10.,40000.,0.25,1.0)'
root -b -q -l 'PlotRatiovar.C("pt_1_1p5eta2p0","GeV",0.,1000.,10.,40000.,0.25,1.0)'
root -b -q -l 'PlotRatiovar.C("pt_1_2p0eta","GeV",0.,1000.,10.,40000.,0.25,1.0)'

root -b -q -l 'PlotRatiovar.C("pt_2","GeV",0.,1000.,10.,400000.,0.25,1.0)'
root -b -q -l 'PlotRatiovar.C("pt_2_eta0p9","GeV",0.,1000.,10.,400000.,0.25,1.0)'
root -b -q -l 'PlotRatiovar.C("pt_2_0p9eta1p5","GeV",0.,1000.,10.,40000.,0.25,1.0)'
root -b -q -l 'PlotRatiovar.C("pt_2_1p5eta2p0","GeV",0.,1000.,10.,40000.,0.25,1.0)'
root -b -q -l 'PlotRatiovar.C("pt_2_2p0eta","GeV",0.,1000.,10.,40000.,0.25,1.0)'

#eta
root -b -q -l 'PlotRatiovar.C("eta_0","",-2.5,2.5,4000.,80000.,0.25,1.0)'
root -b -q -l 'PlotRatiovar.C("eta_1","",-2.5,2.5,4000.,80000.,0.25,1.0)'
root -b -q -l 'PlotRatiovar.C("eta_2","",-2.5,2.5,4000.,80000.,0.25,1.0)'

#phi
root -b -q -l 'PlotRatiovar.C("phi_0","",-3.3,3.3,500.,30000.,0.35,0.9)'
root -b -q -l 'PlotRatiovar.C("phi_1","",-3.3,3.3,500.,30000.,0.35,0.9)'
root -b -q -l 'PlotRatiovar.C("phi_2","",-3.3,3.3,500.,30000.,0.35,0.9)'

#m12
root -b -q -l 'PlotRatiovar.C("m12_mujsel",          "GeV",0.,1000.,10.,80000.,0.25,1.0)'
root -b -q -l 'PlotRatiovar.C("m12_mujsel_loweta",   "GeV",0.,1000.,10.,80000.,0.25,1.0)'
root -b -q -l 'PlotRatiovar.C("m12_mujsel_mediumeta","GeV",0.,1000.,10.,60000.,0.25,1.0)'
root -b -q -l 'PlotRatiovar.C("m12_mujsel_higheta",  "GeV",0.,1000.,10.,20000.,0.25,1.0)'
root -b -q -l 'PlotRatiovar.C("m12_mujsel_mixedeta", "GeV",0.,1000.,10.,80000.,0.25,1.0)'

#2D histos
root -b -q -l 'Plot2dvar.C("n_ptmin20_phivseta")'
root -b -q -l 'Plot2dvar.C("n_ptmin20_btagsel_phivseta")'
root -b -q -l 'Plot2dvar.C("eta_phi_mu")'
root -b -q -l 'Plot2dvar.C("n_muons_phivseta")'
