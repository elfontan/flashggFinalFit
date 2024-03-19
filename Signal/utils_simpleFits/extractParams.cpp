// ---------------------------------------------------------- //
// Macro to plot extract lists of parameter from RooFitResult //
// ---------------------------------------------------------- //

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

void extractParams(){
  //double massindex = 0;
  //const vector<int> v_massindex = {0,   1,   2,   3,   4,   5,   6,   7,   8,   9,  10,   11,   12,  13 };
  //const vector<double> v_mass     = {5., 10., 15., 20., 25., 30., 35., 40., 45., 50., 55.,  60.,  65., 70.}; 
  const vector<double> v_mass     = { 10., 15., 20., 25., 30., 35., 40., 45., 50., 55.,  60.,  65., 70.}; 

  //INPUT FILE WITH HISTOGRAMS TO FIT SIGNAL
  TFile* sig_file = NULL; 
  RooFitResult* result = {};
  vector<float> v_frac1 = {};  
  vector<float> v_frac2 = {};  
  vector<float> v_mean1 = {};  
  vector<float> v_mean2 = {};  
  vector<float> v_mean3 = {};  
  vector<float> v_sig1 = {};  
  vector<float> v_sig2 = {};  
  vector<float> v_sig3 = {};  

  vector<float> v_frac1_err = {};  
  vector<float> v_frac2_err = {};  
  vector<float> v_mean1_err = {};  
  vector<float> v_mean2_err = {};  
  vector<float> v_mean3_err = {};  
  vector<float> v_sig1_err = {};  
  vector<float> v_sig2_err = {};  
  vector<float> v_sig3_err = {};  
  
  for (int i=0; i<14; ++i){
    sig_file=TFile::Open(Form("tests_freezingParams_signalModel/sig_ggh_M%d.root", int(v_mass.at(i))));  
    //sig_file=TFile::Open(Form("simpleFits_initialParamsFrom_paramModel/sig_ggh_M%d.root", int(v_mass.at(i))));  
    //sig_file=TFile::Open(Form("simpleFits_results/sig_ggh_M%d.root", int(v_mass.at(i))));  
    cout << "Signal file: " << sig_file->GetName() << endl;
    // Get the RooFitResult Object
    result = (RooFitResult*)sig_file->Get("fitresult_sig_model_data_sig");

    RooArgList list = result->floatParsFinal();
    cout << " Printing list \n";
    list.Print("s");

    RooRealVar* frac1  = (RooRealVar*)(&list[0]); 
    v_frac1.push_back(frac1->getVal());
    v_frac1_err.push_back(frac1->getError());
    RooRealVar* frac2  = (RooRealVar*)(&list[1]); 
    v_frac2.push_back(frac2->getVal());
    v_frac2_err.push_back(frac2->getError());
    RooRealVar* sig1  = (RooRealVar*)(&list[2]); 
    v_sig1.push_back(sig1->getVal());
    v_sig1_err.push_back(sig1->getError());
    RooRealVar* sig2  = (RooRealVar*)(&list[3]); 
    v_sig2.push_back(sig2->getVal());
    v_sig2_err.push_back(sig2->getError());
    RooRealVar* sig3  = (RooRealVar*)(&list[4]); 
    v_sig3.push_back(sig3->getVal());
    v_sig3_err.push_back(sig3->getError());
    /*
    RooRealVar* mean1  = (RooRealVar*)(&list[2]); 
    v_mean1.push_back(mean1->getVal());
    v_mean1_err.push_back(mean1->getError());
    RooRealVar* mean2  = (RooRealVar*)(&list[3]); 
    v_mean2.push_back(mean2->getVal());
    v_mean2_err.push_back(mean2->getError());
    RooRealVar* mean3  = (RooRealVar*)(&list[4]); 
    v_mean3.push_back(mean3->getVal());
    v_mean3_err.push_back(mean3->getError());
    RooRealVar* sig1  = (RooRealVar*)(&list[5]); 
    v_sig1.push_back(sig1->getVal());
    v_sig1_err.push_back(sig1->getError());
    RooRealVar* sig2  = (RooRealVar*)(&list[6]); 
    v_sig2.push_back(sig2->getVal());
    v_sig2_err.push_back(sig2->getError());
    RooRealVar* sig3  = (RooRealVar*)(&list[7]); 
    v_sig3.push_back(sig3->getVal());
    v_sig3_err.push_back(sig3->getError());
    */
    //cout << "print frac1\n" << frac1->getVal() << endl;
    //cout << "print err frac1\n" << frac1->getError() << endl;
  }

  cout << "frac1 = ";
  for (int p=0; p<14; ++p){
    cout << v_frac1.at(p) << ", ";
  }  
  cout << endl;

  cout << "frac2 = ";
  for (int p=0; p<14; ++p){
    cout << v_frac2.at(p) << ", ";
  }  
  cout << endl;
  /*
  cout << "mean1 = ";
  for (int p=0; p<14; ++p){
    cout << v_mean1.at(p) << ", ";
  }  
  cout << endl;

  cout << "mean2 = ";
  for (int p=0; p<14; ++p){
    cout << v_mean2.at(p) << ", ";
  }  
  cout << endl;

  cout << "mean3 = ";
  for (int p=0; p<14; ++p){
    cout << v_mean3.at(p) << ", ";
  }  
  cout << endl;
  */
  cout << "sigma1 = ";
  for (int p=0; p<14; ++p){
    cout << v_sig1.at(p) << ", ";
  }  
  cout << endl;

  cout << "sigma2 = ";
  for (int p=0; p<14; ++p){
    cout << v_sig2.at(p) << ", ";
  }  
  cout << endl;

  cout << "sigma3 = ";
  for (int p=0; p<14; ++p){
    cout << v_sig3.at(p) << ", ";
  }  
  cout << endl;

  cout << "frac1_err = ";
  for (int p=0; p<14; ++p){
    cout << v_frac1_err.at(p) << ", ";
  }  
  cout << endl;

  cout << "frac2_err = ";
  for (int p=0; p<14; ++p){
    cout << v_frac2_err.at(p) << ", ";
  }  
  cout << endl;
  /*
  cout << "mean1_err = ";
  for (int p=0; p<14; ++p){
    cout << v_mean1_err.at(p) << ", ";
  }  
  cout << endl;

  cout << "mean2_err = ";
  for (int p=0; p<14; ++p){
    cout << v_mean2_err.at(p) << ", ";
  }  
  cout << endl;

  cout << "mean3_err = ";
  for (int p=0; p<14; ++p){
    cout << v_mean3_err.at(p) << ", ";
  }  
  cout << endl;
  */
  cout << "sigma1_err = ";
  for (int p=0; p<14; ++p){
    cout << v_sig1_err.at(p) << ", ";
  }  
  cout << endl;

  cout << "sigma2_err = ";
  for (int p=0; p<14; ++p){
    cout << v_sig2_err.at(p) << ", ";
  }  
  cout << endl;

  cout << "sigma3_err = ";
  for (int p=0; p<14; ++p){
    cout << v_sig3_err.at(p) << ", ";
  }  
  cout << endl;

    //-------------------------
    // Plotting
    //-------------------------	
    //TCanvas c_sig("c_sig", "c_sig", 1200, 1000);                
    //gStyle->SetOptStat(00000);
    //h_corr->SetTitle(Form("Signal ggH (m = %d GeV)", int(v_mass.at(massindex))));
    //h_corr->Draw("colz");
    
    /*
      TLegend *leg = new TLegend(0.88,0.88,0.6,0.76);
      leg->SetLineColor(kWhite);
      leg->SetFillColor(0);
      leg->SetTextFont(42);
      leg->SetTextSize(0.03);
      leg->AddEntry(sig_frame->findObject("SignalModel"),"ggH Signal Model","l");
      leg->Draw();
    */
    
    //c_sig.SaveAs(Form("h_mass%d.png", int(v_mass.at(massindex))));
}
