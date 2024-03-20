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

void dcb_extractParams(){
  //double massindex = 0;
  //const vector<int> v_massindex = {0,   1,   2,   3,   4,   5,   6,   7,   8,   9,  10,   11,   12,  13 };
  //const vector<double> v_mass     = {5., 10., 15., 20., 25., 30., 35., 40., 45., 50., 55.,  60.,  65., 70.}; 
  const vector<double> v_mass     = { 10., 15., 20., 25., 30., 35., 40., 45., 50., 55.,  60.,  65., 70.}; 

  //INPUT FILE WITH HISTOGRAMS TO FIT SIGNAL
  TFile* sig_file = NULL; 
  RooFitResult* result = {};

  // Define parameters for Double Crystal Ball and Gaussian                                                                                                         
  vector<float> v_frac = {};  
  vector<float> v_mean = {};  
  vector<float> v_sigma_dcb = {};  
  vector<float> v_sigma_gaus = {};  
  vector<float> v_a1 = {};  
  vector<float> v_a2 = {};  
  vector<float> v_n1 = {};  
  vector<float> v_n2 = {};  

  vector<float> v_frac_err = {};  
  vector<float> v_mean_err = {};  
  vector<float> v_sigma_dcb_err = {};  
  vector<float> v_sigma_gaus_err = {};  
  vector<float> v_a1_err = {};  
  vector<float> v_a2_err = {};  
  vector<float> v_n1_err = {};  
  vector<float> v_n2_err = {};  
  
  for (int i=0; i<13; ++i){
    sig_file=TFile::Open(Form("tests_signalModel_cat1/sig_ggh_M%d.root", int(v_mass.at(i))));  
    cout << "Signal file: " << sig_file->GetName() << endl;
    // Get the RooFitResult Object
    result = (RooFitResult*)sig_file->Get("fitresult_sig_model_data_sig");

    RooArgList list = result->floatParsFinal();
    cout << " Printing list \n";
    list.Print("s");
 
    RooRealVar* a1  = (RooRealVar*)(&list[0]); 
    v_a1.push_back(a1->getVal());
    v_a1_err.push_back(a1->getError());
    RooRealVar* a2  = (RooRealVar*)(&list[1]); 
    v_a2.push_back(a2->getVal());
    v_a2_err.push_back(a2->getError());
    RooRealVar* frac  = (RooRealVar*)(&list[2]); 
    v_frac.push_back(frac->getVal());
    v_frac_err.push_back(frac->getError());
    RooRealVar* mean  = (RooRealVar*)(&list[3]); 
    v_mean.push_back(mean->getVal());
    v_mean_err.push_back(mean->getError());
    RooRealVar* n1  = (RooRealVar*)(&list[4]); 
    v_n1.push_back(n1->getVal());
    v_n1_err.push_back(n1->getError());
    RooRealVar* n2  = (RooRealVar*)(&list[5]); 
    v_n2.push_back(n2->getVal());
    v_n2_err.push_back(n2->getError());
    RooRealVar* sigma_dcb  = (RooRealVar*)(&list[6]); 
    v_sigma_dcb.push_back(sigma_dcb->getVal());
    v_sigma_dcb_err.push_back(sigma_dcb->getError());
    RooRealVar* sigma_gaus  = (RooRealVar*)(&list[7]); 
    v_sigma_gaus.push_back(sigma_gaus->getVal());
    v_sigma_gaus_err.push_back(sigma_gaus->getError());

  }

  cout << "a1 = ";
  for (int p=0; p<13; ++p){
    cout << v_a1.at(p) << ", ";
  }  
  cout << endl;
  cout << "a2 = ";
  for (int p=0; p<13; ++p){
    cout << v_a2.at(p) << ", ";
  }  
  cout << endl;
  cout << "frac = ";
  for (int p=0; p<13; ++p){
    cout << v_frac.at(p) << ", ";
  }  
  cout << endl;
  cout << "mean = ";
  for (int p=0; p<13; ++p){
    cout << v_mean.at(p) << ", ";
  }  
  cout << endl;
  cout << "n1 = ";
  for (int p=0; p<13; ++p){
    cout << v_n1.at(p) << ", ";
  }  
  cout << endl;
  cout << "n2 = ";
  for (int p=0; p<13; ++p){
    cout << v_n2.at(p) << ", ";
  }  
  cout << endl;
  cout << "sigma_dcb = ";
  for (int p=0; p<13; ++p){
    cout << v_sigma_dcb.at(p) << ", ";
  }  
  cout << endl;
  cout << "sigma_gaus = ";
  for (int p=0; p<13; ++p){
    cout << v_sigma_gaus.at(p) << ", ";
  }  
  cout << endl;

  cout << "a1_err = ";
  for (int p=0; p<13; ++p){
    cout << v_a1_err.at(p) << ", ";
  }  
  cout << endl;
  cout << "a2_err = ";
  for (int p=0; p<13; ++p){
    cout << v_a2_err.at(p) << ", ";
  }  
  cout << endl;
  cout << "frac_err = ";
  for (int p=0; p<13; ++p){
    cout << v_frac_err.at(p) << ", ";
  }  
  cout << endl;
  cout << "mean_err = ";
  for (int p=0; p<13; ++p){
    cout << v_mean_err.at(p) << ", ";
  }  
  cout << endl;
  cout << "n1_err = ";
  for (int p=0; p<13; ++p){
    cout << v_n1_err.at(p) << ", ";
  }  
  cout << endl;
  cout << "n2_err = ";
  for (int p=0; p<13; ++p){
    cout << v_n2_err.at(p) << ", ";
  }  
  cout << endl;
  cout << "sigma_dcb_err = ";
  for (int p=0; p<13; ++p){
    cout << v_sigma_dcb_err.at(p) << ", ";
  }  
  cout << endl;
  cout << "sigma_gaus_err = ";
  for (int p=0; p<13; ++p){
    cout << v_sigma_gaus_err.at(p) << ", ";
  }  
  cout << endl;



}
