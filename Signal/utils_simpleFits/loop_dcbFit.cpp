// -------------------------------------------------------------- //
// Macro for the signal modelling of the low mass diphoton search //
// -------------------------------------------------------------- //

#include "RooRealVar.h"
#include "RooDataSet.h"
#include <RooDataHist.h>
#include "RooGaussian.h"
#include "RooCBShape.h"
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
#include "../../../HiggsAnalysis/CombinedLimit/interface/RooDoubleCBFast.h"

using namespace std;
using namespace RooFit;

void loop_dcbFit(TString year="2018", TString cat="cat1"){
  //const vector<int> v_massindex = {0,   1,   2,   3,   4,   5,   6,   7,   8,   9,  10,   11,   12,  13 };
  //const vector<double> v_mass     = {5., 10., 15., 20., 25., 30., 35., 40., 45., 50., 55.,  60.,  65., 70.}; 
  //const vector<int> v_massindex= {0,   1,   2,   3,   4,   5,   6,   7,   8,   9,  10,   11,   12 };
  const vector<double> v_mass   = {10., 15., 20., 25., 30., 35., 40., 45., 50., 55.,  60.,  65., 70.}; 

  // Iterate over each mass index
  for (size_t massindex = 0; massindex < v_mass.size(); ++massindex) {  
    //INPUT FILE WITH HISTOGRAMS TO FIT SIGNAL
    TFile* sig_file = NULL; 
    //sig_file=TFile::Open(Form("/eos/user/e/elfontan/DiPhotonAnalysis/Oct2022_xsec1_massPoints/signal_histos_cat0/nBins140-rebinFor_5_10/histogram_ggH_M%d.root", int(v_mass.at(massindex))));    
    if (cat=="cat0"){
      std::cout << "--- CATEGORY 0" << std::endl;
      sig_file=TFile::Open(Form("/eos/user/e/elfontan/DiPhotonAnalysis/diphotonBDT/SignalModeling/signal_histos_cat0/histogram_ggH_M%d.root", int(v_mass.at(massindex))));
    }else if (cat=="cat1"){
      std::cout << "--- CATEGORY 1" << std::endl;
      sig_file=TFile::Open(Form("/eos/user/e/elfontan/DiPhotonAnalysis/diphotonBDT/SignalModeling/signal_histos_cat1/histogram_ggH_M%d.root", int(v_mass.at(massindex))));
    }else {
      // Handle the case when neither category is selected
      cout << "Please select a category to process." << endl;
    }
    
    cout << "Signal file: " << sig_file->GetName() << endl;
    
    // Get the histograms
    //TH1D* cat_sig=(TH1D*)sig_file->Get(Form("ggh_M%d_cat0", int(v_mass.at(massindex))));
    TH1D* cat_sig=(TH1D*)sig_file->Get(Form("ggh_M%d_cat1", int(v_mass.at(massindex))));
    double massLow  =  cat_sig->GetXaxis()->GetXmin();
    double massHigh =  cat_sig->GetXaxis()->GetXmax();
    double n_bins = cat_sig->GetSize()-2;
    double fit_min = cat_sig->GetMean()-3*cat_sig->GetRMS();
    double fit_max = cat_sig->GetMean()+3*cat_sig->GetRMS();
    int fit_nbins = int((fit_max - fit_min)/((massHigh - massLow)/n_bins));
    
    vector<float> v_mean_dcb; 
    vector<float> v_sigma_dcb;
    vector<float> v_a1_dcb;
    vector<float> v_a2_dcb;
    vector<float> v_n1_dcb;
    vector<float> v_n2_dcb;
    vector<float> v_dm_dcb;
    vector<float> v_sigma_gaus;
    vector<float> v_frac;
    vector<float> v_fea;
    
    // Compute mass point and define ROOFit variables
    RooRealVar CMS_Hgg_mass("CMS_Hgg_mass", "CMS_Hgg_mass", fit_min, fit_max);
    //EF RooRealVar CMS_Hgg_mass("CMS_Hgg_mass", "CMS_Hgg_mass", massLow, massHigh);
    
    
    // Define parameters for Double Crystal Ball and Gaussian
    RooRealVar mean_dcb("mean_dcb", "mean_dcb", v_mass[massindex], massLow, massHigh);
    RooRealVar sigma_dcb("sigma_dcb", "sigma_dcb", 0.4, 0.01, 4.);
    RooRealVar a1("a1", "a1", 0.5, 0.1, 5.);
    RooRealVar n1("n1", "n1", 10.0, 0.1, 20.);
    RooRealVar a2("a2", "a2", 0.3, 0.1, 5.);
    RooRealVar n2("n2", "n2", 10.0, 0.1, 25.);
    RooRealVar frac("frac", "frac", 0.2, 0., 1.);
    RooRealVar sigma_gaus("sigma_gaus", "sigma_gaus", 0.3, 0.05, 1.5);
    
    // -----------------------
    // Define the signal model
    // -----------------------
    RooDataHist data_sig("data_sig", "", RooArgList(CMS_Hgg_mass), cat_sig);                                                                    
    RooRealVar sig_norm("sig_norm", "",cat_sig->Integral());
    
    // Double Crystal Ball PDF
    RooDoubleCBFast doubleCB("doubleCB", "doubleCB", CMS_Hgg_mass, mean_dcb, sigma_dcb, a1, n1, a2, n2);
    // Gaussian
    RooGaussian gauss("gauss", "gauss", CMS_Hgg_mass, mean_dcb, sigma_gaus);
    
    // Define the signal model as a sum of Double Crystal Ball and Gaussian
    //RooAddPdf sig_model("sig_model", "sig_model", RooArgList(doubleCB, gauss), RooArgList(frac));
    RooAddPdf sig_model("sig_model", "sig_model", RooArgList(doubleCB, gauss), RooArgList(frac), true); 	
    

    //-------------------------                                                                                                                                   
    // Fitting                                                                                                                                                    
    //-------------------------                                                                                                                                   
    RooFitResult* result = sig_model.fitTo(data_sig, Range(fit_min, fit_max), Save());
    //RooFitResult* result = sig_model.fitTo(data_sig, Save());                                                                                                   
    cout << "#########################" << endl;
    cout << "##### PRINT RESULTS #####" << endl;
    cout << "#########################" << endl;
    result->Print("v");
    TH2* hcorr = result->correlationHist();

    /*TFile* sig_ggh_out = NULL;
    if (cat=="cat0"){
      sig_ggh_out = new TFile(Form("tests_signalModel_cat0/sig_ggh_M%d.root", int(v_mass.at(massindex))), "RECREATE");
    }
    else if (cat=="cat1"){
      sig_ggh_out = new TFile(Form("tests_signalModel_cat1/sig_ggh_M%d.root", int(v_mass.at(massindex))), "RECREATE");
    }
    */
    TFile* sig_ggh_out = new TFile(Form("tests_signalModel_cat1/sig_ggh_M%d.root", int(v_mass.at(massindex))), "RECREATE");
    //sig_ggh_out->cd();
    result->Write();
    hcorr->Write();
    sig_ggh_out->Close();

     // Each RooRealVar is set to be a constant
    mean_dcb.setConstant(kTRUE);
    sigma_dcb.setConstant(kTRUE);
    a1.setConstant(kTRUE);
    n1.setConstant(kTRUE);
    a2.setConstant(kTRUE);
    n2.setConstant(kTRUE);
    sigma_gaus.setConstant(kTRUE);
    frac.setConstant(kTRUE);
  

    //-------------------------
    // RooPlot
    //-------------------------	
    RooPlot *sig_frame = CMS_Hgg_mass.frame();                                                 
    cout << "-----------------------------" << endl;
    cout << "# TH1D Histo bins: " << n_bins << endl;
    cout << "# NBins data_sig = " << (data_sig.numEntries() -2) << endl;
    // Print number of floating parameters in the PDF                                                                          
    //cout << "# Number of floating parameters: " << result->floatParsFinal().getSize() << endl;                                                        
    
    sig_frame->SetTitle("");
    sig_frame->GetXaxis()->SetTitle("Mass [GeV]");
    sig_frame->GetYaxis()->SetTitle("Events/0.1");
    data_sig.plotOn(sig_frame);                                  
    if (cat=="cat0"){
      sig_model.plotOn(sig_frame, RooFit::Name("SignalModel"), LineColor(kAzure-4));
    }
    else if (cat=="cat1"){
      sig_model.plotOn(sig_frame, RooFit::Name("SignalModel"), LineColor(kAzure-1));
    }
    
    std::cout << "################" << std::endl;
    std::cout << "##### CHI2 #####" << std::endl;
    std::cout << "################" << std::endl;
    double chi2_over_ndof = sig_frame->chiSquare(0); //tell the chiSquare method how many free parameters you have         
    //double chi2_over_ndof = sig_frame->chiSquare(result->floatParsFinal().getSize()); //tell the chiSquare method how many free parameters you have         
    cout << "Chi2/ndof = " << chi2_over_ndof << endl;
    cout << "-----------------------------" << endl;
    
    
    //-------------------------
    // Plotting
    //-------------------------	
    TCanvas c_sig("c_sig", "c_sig", 1200, 1000);                
    sig_frame->Draw("");
    
    TLegend *leg = new TLegend(0.88,0.88,0.6,0.76);
    //EF leg->SetHeader("#chi^{2};/N_{dof} = ...");
    leg->SetHeader(Form("#chi^{2}/N_{dof} = %.3f", chi2_over_ndof));
    leg->SetLineColor(kWhite);
    leg->SetFillColor(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.03);
    leg->AddEntry(sig_frame->findObject("SignalModel"),"ggH Signal Model","l");
    leg->Draw();
    //TLegendEntry *header = (TLegendEntry*)leg->GetListOfPrimitives()->First();
    //header->SetTextSize(.03);
    
    if (cat=="cat0"){
      c_sig.SaveAs(Form("/eos/user/e/elfontan/www/LowMassDiPhoton/diphotonBDT/ParametricBDT/TensorFlow/SignalModeling/newNtuples_newTraining_NN0p907/dcb_standalone/sig_ggH_cat0_"+year+"_dcb_mass%d.png", int(v_mass.at(massindex))));
    }else if (cat=="cat1"){
      c_sig.SaveAs(Form("/eos/user/e/elfontan/www/LowMassDiPhoton/diphotonBDT/ParametricBDT/TensorFlow/SignalModeling/newNtuples_newTraining_NN0p907/dcb_standalone/sig_ggH_cat1_"+year+"_dcb_mass%d.png", int(v_mass.at(massindex))));
    }else {
      // Handle the case when neither category is selected
      cout << "Please select a category to process." << endl;
    }
    
    /*
    //-------------------------
  // Save into ROO workspace
  //-------------------------
  RooWorkspace dpworkspace("dpworkspace", "");
  dpworkspace.import(sig_model);
  dpworkspace.writeToFile(Form("output/dpWorkspace"+year+suff+"_%d_new_wgt.root",i));
    */

  }
}
