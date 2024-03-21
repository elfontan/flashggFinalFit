// ------------------------------------------------------- //
// Macro to plot the correlation matrix between parameters //
// ------------------------------------------------------- //

#include "RooRealVar.h"
#include "RooDataSet.h"
#include <RooDataHist.h>
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooAddPdf.h"
#include <RooBernstein.h>
#include "RooFitResult.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TText.h"
#include "TAttLine.h"
#include "TLegend.h"
#include "RooPlot.h"
#include "TFile.h"
#include "TStyle.h"
#include "TH2.h"
#include "TMatrixDSym.h"
#include <RooChi2Var.h>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <string>

using namespace std;
using namespace RooFit;

void plot_corr(){
  //double massindex = 0;
  //const vector<int> v_massindex = {0,   1,   2,   3,   4,   5,   6,   7,   8,   9,  10,   11,   12,  13 };
  const vector<double> v_mass     = {5., 10., 15., 20., 25., 30., 35., 40., 45., 50., 55.,  60.,  65., 70.}; 

  TFile* sig_file = NULL; 
  
  for (int i=0; i<14; ++i){
    //INPUT FILE WITH HISTOGRAMS TO FIT SIGNAL
    sig_file=TFile::Open(Form("new_simpleFits_largerRange_initialParamsFrom_paramModel/sig_ggh_M%d.root", int(v_mass.at(i))));
    cout << "Signal file: " << sig_file->GetName() << endl;
    
    // Get the histograms
    TH2D* h_corr=(TH2D*)sig_file->Get("correlation_matrix");
    
    //-------------------------
    // Plotting
    //-------------------------	
    TCanvas c_sig("c_sig", "c_sig", 1200, 1000);                
    gStyle->SetOptStat(00000);
    h_corr->SetTitle(Form("Signal ggH (m = %d GeV)", int(v_mass.at(i))));
    h_corr->Draw("colz");
  
    /*
      TLegend *leg = new TLegend(0.88,0.88,0.6,0.76);
      leg->SetLineColor(kWhite);
      leg->SetFillColor(0);
      leg->SetTextFont(42);
      leg->SetTextSize(0.03);
      leg->AddEntry(sig_frame->findObject("SignalModel"),"ggH Signal Model","l");
      leg->Draw();
    */
    
    c_sig.SaveAs(Form("new_simpleFits_largerRange_initialParamsFrom_paramModel/h_corr_ggH_mass%d.png", int(v_mass.at(i))));
  }
}
