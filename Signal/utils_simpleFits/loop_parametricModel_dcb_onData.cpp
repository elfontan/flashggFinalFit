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

void loop_parametricModel_dcb_onData(TString year="2018", TString cat="cat0"){
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
    TH1D* cat_sig=(TH1D*)sig_file->Get(Form("ggh_M%d_%s", int(v_mass.at(massindex)), cat.Data()));
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
    
    // Compute mass point and define ROOFit variables
    RooRealVar CMS_Hgg_mass("CMS_Hgg_mass", "CMS_Hgg_mass", fit_min, fit_max);
    //EF RooRealVar CMS_Hgg_mass("CMS_Hgg_mass", "CMS_Hgg_mass", massLow, massHigh);
    
    // Linear interpolation from simple fit results with mean values fixed - with (first line) and without (second line) error bars considered) 	
    if (cat=="cat0"){
      /*
      v_mean_dcb = {};
      v_sigma_dcb = {};
      v_a1_dcb = {};
      v_a2_dcb = {};
      v_n1_dcb = {};
      v_n2_dcb = {};
      v_sigma_gaus = {};
      v_frac = {};
      v_dm_dcb = {};
      */
      /*v_mean_dcb = {10.000446301331976, 14.991267825476825, 19.982089350000024, 24.972910875454545, 29.96373239904642, 34.95455392450094, 39.94537544995546, 44.936196975409985, 49.927018500864506, 54.917840018868446, 59.90866154432297, 64.89948306977749, 69.89030459523201};
      v_sigma_dcb = {0.6068516373634338, 0.6211022138595581, 0.6353528499603271, 0.6496034860610962, 0.6638540625572205, 0.6781046986579895, 0.6923553347587585, 0.7066059112548828, 0.7208565473556519, 0.7351071834564209, 0.7493577599525452, 0.7636083960533142, 0.7778589725494385};
      v_a1_dcb = {1.3843529224395752, 1.3202941417694092, 1.2562353610992432, 1.1921764612197876, 1.1281176805496216, 1.064058780670166, 1.0, 0.9359411597251892, 0.8718823194503784, 0.8078235387802124, 0.7437646985054016, 0.6797058582305908, 0.61564701795578};
      v_a2_dcb = {1.4081695079803467, 1.4501131772994995, 1.492056965827942, 1.5340006351470947, 1.5759443044662476, 1.6178879737854004, 1.6598316431045532, 1.7017754316329956, 1.7437191009521484, 1.7856627702713013, 1.827606439590454, 1.8695502281188965, 1.9114938974380493};
      v_n1_dcb = {35.0125846862793, 35.0134162902832, 35.014251708984375, 35.01508331298828, 35.01591873168945, 35.01675033569336, 35.017581939697266, 35.01841735839844, 35.019248962402344, 35.020084381103516, 35.02091598510742, 35.021751403808594, 35.0225830078125};
      v_n2_dcb = {19.04859733581543, 19.54859733581543, 20.04859733581543, 20.54859733581543, 21.04859733581543, 21.54859733581543, 22.04859733581543, 22.54859733581543, 23.04859733581543, 23.54859733581543, 24.04859733581543, 24.54859733581543, 25.04859733581543};
      v_sigma_gaus = {0.1117207258939743, 0.14166198670864105, 0.1716032326221466, 0.20154449343681335, 0.2314857542514801, 0.26142701506614685, 0.2913682758808136, 0.32130953669548035, 0.3512507975101471, 0.38119202852249146, 0.4111332893371582, 0.44107455015182495, 0.4710158109664917};
      v_frac = {0.17271964251995087, 0.1927834302186966, 0.21284721791744232, 0.23291102051734924, 0.25297480821609497, 0.2730385959148407, 0.2931023836135864, 0.31316620111465454, 0.33322998881340027, 0.353293776512146, 0.3733575642108917, 0.39342135190963745, 0.41348516941070557};*/
      v_mean_dcb = {};
      v_sigma_dcb = {};
      v_a1_dcb = {};
      v_a2_dcb = {};
      v_n1_dcb = {};
      v_n2_dcb = {};
      v_sigma_gaus = {};
      v_frac = {};      
    } else if (cat=="cat1"){
      /*v_mean_dcb = {9.94032847136259, 14.936102859675884, 19.93187725543976, 24.92765164375305, 29.923426039516926, 34.91920042783022, 39.91497481614351, 44.91074921190739, 49.90652360022068, 54.902297995984554, 59.89807238429785, 64.89384678006172, 69.88962116837502};
      v_sigma_dcb = {1.0406465530395508, 1.018246054649353, 0.9958455562591553, 0.9734450578689575, 0.9510445594787598, 0.928644061088562, 0.9062435626983643, 0.8838430643081665, 0.8614425659179688, 0.8390420079231262, 0.8166415095329285, 0.7942410111427307, 0.771840512752533};
      v_a1_dcb = {1.4002058506011963, 1.3338029384613037, 1.2674001455307007, 1.200997233390808, 1.1345943212509155, 1.0681915283203125, 1.00178861618042, 0.9353857636451721, 0.8689829111099243, 0.8025800585746765, 0.7361771464347839, 0.6697742938995361, 0.6033714413642883};
      v_a2_dcb = {1.1081123352050781, 1.100339651107788, 1.0925670862197876, 1.0847944021224976, 1.0770217180252075, 1.069249153137207, 1.061476469039917, 1.053703784942627, 1.0459312200546265, 1.0381585359573364, 1.0303858518600464, 1.022613286972046, 1.0148406028747559};
      v_n1_dcb = {23.625934600830078, 23.637266159057617, 23.648595809936523, 23.659927368164062, 23.67125701904297, 23.682588577270508, 23.693920135498047, 23.705249786376953, 23.716581344604492, 23.7279109954834, 23.739242553710938, 23.750572204589844, 23.761903762817383};
      v_n2_dcb = {27.144695281982422, 27.644695281982422, 28.144695281982422, 28.644695281982422, 29.144695281982422, 29.644695281982422, 30.144695281982422, 30.644695281982422, 31.144695281982422, 31.644695281982422, 32.14469528198242, 32.64469528198242, 33.14469528198242};
      v_sigma_gaus = {0.1940184086561203, 0.24374337494373322, 0.29346832633018494, 0.34319329261779785, 0.39291825890541077, 0.4426432251930237, 0.4923681914806366, 0.5420931577682495, 0.5918181538581848, 0.6415430903434753, 0.6912680864334106, 0.7409930229187012, 0.7907180190086365};
      v_frac = {0.01507380697876215, 0.08722556382417679, 0.159377321600914, 0.23152907192707062, 0.3036808371543884, 0.37583258748054504, 0.44798433780670166, 0.5201361179351807, 0.5922878384590149, 0.6644396185874939, 0.7365913987159729, 0.8087431192398071, 0.8808948993682861};*/
      v_mean_dcb = {};
      v_sigma_dcb = {};
      v_a1_dcb = {};
      v_a2_dcb = {};
      v_n1_dcb = {};
      v_n2_dcb = {};
      v_sigma_gaus = {};
      v_frac = {};


    } else {
      // Handle the case when neither category is selected
      cout << "Please select a category to process." << endl;
    }
    
    // Define parameters for Double Crystal Ball and Gaussian
    RooRealVar mean_dcb("mean_dcb", "mean_dcb", v_mass[massindex], massLow, massHigh);
    RooRealVar sigma_dcb("sigma_dcb", "sigma_dcb", 0.4, 0.1, 2.);
    RooRealVar a1("a1", "a1", 2.0, 0.1, 10.);
    RooRealVar n1("n1", "n1", 1.0, 0.1, 20.);
    RooRealVar a2("a2", "a2", 1.0, 0.1, 10.);
    RooRealVar n2("n2", "n2", 1.0, 0.1, 50.);
    RooRealVar frac("frac", "frac", 0.2, 0.01, 0.9);
    RooRealVar sigma_gaus("sigma_gaus", "sigma_gaus", 0.2, 0.05, 1.5);
    
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
    
    std::cout << "PARAMETERS USED: " << std::endl;
    std::cout << "mean_dcb: " << v_mean_dcb.at(massindex) << std::endl;
    //std::cout << "mean0: " << v_mass.at(massindex) + v_dm0.at(massindex) << std::endl;
    //std::cout << "mean0: " << v_dm0.at(massindex) << std::endl;
    
  
    /* NOTE about the mean: best fit from flashggFinalFit given as dm0, dm1, dm2
       RooFormulaVar::mean_g0_GG2H_2018_UntaggedTag_0_13TeV[ actualVars=(MH,dm_g0_GG2H_2018_UntaggedTag_0_13TeV) formula="(@0+@1)" ] = 70.0027
       RooFormulaVar::mean_g1_GG2H_2018_UntaggedTag_0_13TeV[ actualVars=(MH,dm_g1_GG2H_2018_UntaggedTag_0_13TeV) formula="(@0+@1)" ] = 69.4859
       RooFormulaVar::mean_g2_GG2H_2018_UntaggedTag_0_13TeV[ actualVars=(MH,dm_g2_GG2H_2018_UntaggedTag_0_13TeV) formula="(@0+@1)" ] = 67.0678*/
    
    // Each RooRealVar is set to be a constant 
    mean_dcb.setConstant(kTRUE);
    sigma_dcb.setConstant(kTRUE);
    a1.setConstant(kTRUE);
    n1.setConstant(kTRUE);
    a2.setConstant(kTRUE);
    n2.setConstant(kTRUE);
    sigma_gaus.setConstant(kTRUE);
    frac.setConstant(kTRUE);
   //mean1.setVal(v_dm0.at(massindex));
    //mean2.setVal(v_dm1.at(massindex));
    //mean3.setVal(v_dm2.at(massindex));
  
    mean_dcb.setVal(v_mean_dcb.at(massindex));
    sigma_dcb.setVal(v_sigma_dcb.at(massindex));
    a1.setVal(v_a1_dcb.at(massindex));
    a2.setVal(v_a2_dcb.at(massindex));
    n1.setVal(v_n1_dcb.at(massindex));
    n2.setVal(v_n2_dcb.at(massindex));
    sigma_gaus.setVal(v_sigma_gaus.at(massindex));
    frac.setVal(v_frac.at(massindex));
    
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
      sig_model.plotOn(sig_frame, RooFit::Name("gauss"), RooFit::Components("gauss"), LineColor(kCyan+1), LineStyle(kDashed));
      sig_model.plotOn(sig_frame, RooFit::Name("doubleCB"), RooFit::Components("doubleCB"), LineColor(kBlue-9), LineStyle(kDashed));
      sig_model.plotOn(sig_frame, RooFit::Name("SignalModel"), LineColor(kBlue));
   }
    else if (cat=="cat1"){
      sig_model.plotOn(sig_frame, RooFit::Name("gauss"), RooFit::Components("gauss"), LineColor(kRed-4), LineStyle(kDashed));
      sig_model.plotOn(sig_frame, RooFit::Name("doubleCB"), RooFit::Components("doubleCB"), LineColor(kRed+1), LineStyle(kDashed));
      sig_model.plotOn(sig_frame, RooFit::Name("SignalModel"), LineColor(kOrange-3));
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
    

    /*TLegend *legend = new TLegend(0.17,0.66,0.5,0.76,NULL,"brNDC");                                                                                          
    legend->SetFillColor(0);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextSize(0.04);
    //legend->AddEntry(h201618,"Full reconstruction","l"); */

    TLatex *latex2 = new TLatex();
    latex2->SetNDC();
    latex2->SetTextSize(0.4*c_sig.GetTopMargin());
    latex2->SetTextFont(42);
    latex2->SetTextAlign(31);// # align right                                                                                                                       
    latex2->DrawLatex(0.9, 0.92,"13 TeV");                                                                                                    
    //latex2->DrawLatex(0.9, 0.92,"96.6 fb^{-1} (13 TeV)");                                                                

    latex2->SetTextSize(0.55*c_sig.GetTopMargin());
    latex2->SetTextFont(62);
    latex2->SetTextAlign(11);// # align right                                                                                                                       
    latex2->DrawLatex(0.11, 0.92, "CMS"); // Out of the canvas        
    latex2->SetTextFont(42);
    latex2->SetTextSize(0.4*c_sig.GetTopMargin());
    latex2->DrawLatex(0.2, 0.92, "Simulation Preliminary"); // Out of the canvas        
    //latex2->DrawLatex(0.18, 0.82, "CMS"); //NORM                                                                                                                  

    latex2->SetTextSize(0.5*c_sig.GetTopMargin());
    latex2->SetTextFont(52);
    latex2->SetTextAlign(11);

    TLegend *leg = new TLegend(0.88,0.88,0.6,0.66);
    //EF leg->SetHeader("#chi^{2};/N_{dof} = ...");
    leg->SetHeader(Form("ggH MC (#chi^{2}/N_{dof} = %.3f)", chi2_over_ndof));
    leg->SetLineColor(kWhite);
    leg->SetFillColor(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.03);
    leg->AddEntry(sig_frame->findObject("SignalModel"),"Total signal fit","l");
    leg->AddEntry(sig_frame->findObject("gauss"),"Gaussian","l");
    leg->AddEntry(sig_frame->findObject("doubleCB"),"Double Crystal Ball","l");
    leg->Draw();
    //TLegendEntry *header = (TLegendEntry*)leg->GetListOfPrimitives()->First();
    //header->SetTextSize(.03);
    
    if (cat=="cat0"){
      c_sig.SaveAs(Form("/eos/user/e/elfontan/www/LowMassDiPhoton/diphotonBDT/ParametricBDT/TensorFlow/SignalModeling/newNtuples_newTraining_NN0p907/TESTS/sig_ggH_cat0_"+year+"_dcb_mass%d.png", int(v_mass.at(massindex))));
      //c_sig.SaveAs(Form("/eos/user/e/elfontan/www/LowMassDiPhoton/diphotonBDT/ParametricBDT/TensorFlow/SignalModeling/newNtuples_newTraining_NN0p907/TESTS/mNom40_tunedParams/sig_ggH_cat0_"+year+"_dcb_mass%d.png", int(v_mass.at(massindex))));
    }else if (cat=="cat1"){
      c_sig.SaveAs(Form("/eos/user/e/elfontan/www/LowMassDiPhoton/diphotonBDT/ParametricBDT/TensorFlow/SignalModeling/newNtuples_newTraining_NN0p907/TESTS/sig_ggH_cat1_"+year+"_dcb_mass%d.png", int(v_mass.at(massindex))));
      //c_sig.SaveAs(Form("/eos/user/e/elfontan/www/LowMassDiPhoton/diphotonBDT/ParametricBDT/TensorFlow/SignalModeling/newNtuples_newTraining_NN0p907/TESTS/mNom40_tunedParams/sig_ggH_cat1_"+year+"_dcb_mass%d.png", int(v_mass.at(massindex))));
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
