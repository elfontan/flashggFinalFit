// -------------------------------------------------------------- //
// Macro for the signal modelling of the low mass diphoton search //
// -------------------------------------------------------------- //

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

void signalModel_freezingSomeParams(TString year="2018"){
  const int massindex = 0;
  //const vector<int> v_massindex = {0,   1,   2,   3,   4,   5,   6,   7,   8,   9,  10,   11,   12,  13 };
  const vector<double> v_mass     = {5., 10., 15., 20., 25., 30., 35., 40., 45., 50., 55.,  60.,  65., 70.}; 
  
	//INPUT FILE WITH HISTOGRAMS TO FIT SIGNAL
        TFile* sig_file = NULL; 
        //sig_file=TFile::Open(Form("/eos/user/e/elfontan/DiPhotonAnalysis/Oct2022_xsec1_massPoints/signal_histos_cat0/nBins140-rebinFor_5_10/histogram_ggH_M%d.root", int(v_mass.at(massindex))));    
        sig_file=TFile::Open(Form("/eos/user/e/elfontan/DiPhotonAnalysis/Oct2022_xsec1_massPoints/signal_histos_cat0/histogram_ggH_M%d.root", int(v_mass.at(massindex))));    
	cout << "Signal file: " << sig_file->GetName() << endl;

	// Get the histograms
	TH1D* cat_sig=(TH1D*)sig_file->Get(Form("ggh_M%d_cat0", int(v_mass.at(massindex))));
	double massLow  =  cat_sig->GetXaxis()->GetXmin();
	double massHigh =  cat_sig->GetXaxis()->GetXmax();
	double n_bins = cat_sig->GetSize()-2;
	double fit_min = cat_sig->GetMean()-3*cat_sig->GetRMS();
	double fit_max = cat_sig->GetMean()+3*cat_sig->GetRMS();
	int fit_nbins = int((fit_max - fit_min)/((massHigh - massLow)/n_bins)); 
	
	// Compute mass point and define ROOFit variables
	RooRealVar CMS_Hgg_mass("CMS_Hgg_mass", "CMS_Hgg_mass", fit_min, fit_max);
	//RooRealVar CMS_Hgg_mass("CMS_Hgg_mass", "CMS_Hgg_mass", massLow, massHigh);

        // FlashggFinalFit: linear interpolation from simple fit results                                                                   
        const vector<float> v_frac0 = {
          0.6709164285714286, 0.6564129560439561, 0.6419094835164835, 0.6274060109890109, 0.6129025384615384, 0.5983990659340659, 0.5838955934065932, 0.5693921208791207, 0.5548886483516482, 0.5403851758241757, 0.5258817032967031, 0.5113782307692305, 0.496874758241758, 0.48237128571428545
        };
        const vector<float> v_frac1 = {
          0.3971574411428576, 0.41470643503296745, 0.4322554289230773, 0.4498044228131871, 0.467353416703297, 0.4849024105934068, 0.5024514044835167, 0.520000398373264, 0.5375493922637363, 0.5550983861538461, 0.572647380043956, 0.5901963739340659, 0.6077453678241757, 0.6252943617142854
        };
        const vector<float> v_dm0 = {
          4.997394285714296, 9.993632527472538, 14.989870769230782, 19.986109010989026, 24.98234725274727, 29.978585494505513, 34.97482373626375, 39.97106197802199, 44.967300219780235, 49.963538461538484, 54.959776703296725, 59.95601494505497, 64.95225318681322, 69.94849142857146
        };
        const vector<float> v_dm1 = {
          4.453346285714301, 9.427265978021996, 14.401185670329689, 19.375105362637385, 24.349025054945077, 29.32294474725277, 34.296864439560466, 39.27078413186816, 44.24470382417585, 49.21862351648355, 54.19254320879124, 59.16646290109893, 64.14038259340663, 69.11430228571432
        };
        const vector<float> v_dm2 = {
          5.693821142857154, 10.612516351648363, 15.531211560439573, 20.449906769230783, 25.368601978021992, 30.287297186813202, 35.20599239560441, 40.124687604395625, 45.04338281318683, 49.962078021978044, 54.88077323076925, 59.799468439560464, 64.71816364835166, 69.63685885714287
        };
        const vector<float> v_sig0 = {
          0.11679265714285728, 0.15592185274725287, 0.1950510483516485, 0.2341802439560441, 0.2733094395604397, 0.3124386351648353, 0.35156783076923087, 0.39069702637362647, 0.42982622197802206, 0.46895541758241766, 0.5080846131868133, 0.547213808791209, 0.5863430043956045, 0.6254722000000001
        };
        const vector<float> v_sig1 = {
          1.5329935714285718, 1.5006269890109893, 1.4682604065934066, 1.4358938241758241, 1.4035272417582416, 1.3711606593406591, 1.3387940769230766, 1.3064274945054941, 1.2740609120879114, 1.241694329670329, 1.2093277472527464, 1.1769611648351639, 1.1445945824175814, 1.1122279999999987
        };
        const vector<float> v_sig2 = {
          1.7244708285714272, 1.6520574813186801, 1.5796441340659328, 1.5072307868131856, 1.4348174395604385, 1.3624040923076912, 1.289990745054944, 1.2175773978021969, 1.1451640505494496, 1.0727507032967023, 1.0003373560439552, 0.927924008791208, 0.8555106615384608, 0.7830973142857136
        };

	// Good set of initial parameters for simple fits using all mass points	
	//RooRealVar mean1("mean1", "mean1", v_mass[massindex], massLow, massHigh);
	//RooRealVar mean2("mean2", "mean2", v_mass[massindex], massLow, massHigh);
	//RooRealVar mean3("mean3", "mean3", v_mass[massindex], massLow, massHigh);
	RooRealVar sigma1("sigma1", "sigma1", 0.1, 0.1, 1.); 
	RooRealVar sigma2("sigma2", "sigma2", 0.4, 0.1, 4.); 
	RooRealVar sigma3("sigma3", "sigma3", 0.4, 0.1, 7.); 
	RooRealVar frac1("frac1", "frac1", 0.45, 0., 1.2);
	RooRealVar frac2("frac2", "frac2", 0.45, 0., 1.);	

	RooRealVar mean1("mean1", "mean1", v_dm0.at(massindex));
	RooRealVar mean2("mean2", "mean2", v_dm1.at(massindex));
	RooRealVar mean3("mean3", "mean3", v_dm2.at(massindex));
	//RooRealVar sigma1("sigma1", "sigma1", v_sig0.at(massindex)); 
	//RooRealVar sigma2("sigma2", "sigma2", v_sig1.at(massindex)); 
	//RooRealVar sigma3("sigma3", "sigma3", v_sig2.at(massindex)); 
	//RooRealVar frac1("frac1", "frac1", v_frac0.at(massindex));
	//RooRealVar frac2("frac2", "frac2", v_frac1.at(massindex));	

	std::cout << "mean0: " << v_dm0.at(massindex) << std::endl;
	std::cout << "mean1: " << v_dm1.at(massindex) << std::endl;
	std::cout << "mean2: " << v_dm2.at(massindex) << std::endl;
	std::cout << "sig0: " << v_sig0.at(massindex) << std::endl;
	std::cout << "sig1: " << v_sig1.at(massindex) << std::endl;
	std::cout << "sig2: " << v_sig2.at(massindex) << std::endl;
	std::cout << "frac0: " << v_frac0.at(massindex) << std::endl;
	std::cout << "frac1: " << v_frac1.at(massindex) << std::endl;

	/*
	RooRealVar mean1("mean1", "mean1", v_mass.at(massindex) + v_dm0.at(massindex), massLow, massHigh);
	RooRealVar mean2("mean2", "mean2", v_mass.at(massindex) + v_dm1.at(massindex), massLow, massHigh);
	RooRealVar mean3("mean3", "mean3", v_mass.at(massindex) + v_dm2.at(massindex), massLow, massHigh);
	RooRealVar sigma1("sigma1", "sigma1", v_sig0.at(massindex), 0.1, 3.); 
	RooRealVar sigma2("sigma2", "sigma2", v_sig1.at(massindex), 0.3, 3.); 
	RooRealVar sigma3("sigma3", "sigma3", v_sig2.at(massindex), 0.3, 10.); 
	RooRealVar frac1("frac1", "frac1", v_frac0.at(massindex), 0., 2.);
	RooRealVar frac2("frac2", "frac2", v_frac1.at(massindex), 0., 1.);
	*/

	// -----------------------
	// Define the signal model
	// -----------------------
	RooDataHist data_sig("data_sig", "", RooArgList(CMS_Hgg_mass), cat_sig);        
	RooRealVar sig_norm("sig_norm", "",cat_sig->Integral());

	//cout << "RooDataHist bins: " << data_sig << endl;

	//RooGaussian sig_model("sig_model", "sig_model", CMS_Hgg_mass, mean, sigma);
	RooGaussian gauss1("gauss1", "gauss1", CMS_Hgg_mass, mean1, sigma1);
	RooGaussian gauss2("gauss2", "gauss2", CMS_Hgg_mass, mean2, sigma2);
	RooGaussian gauss3("gauss3", "gauss3", CMS_Hgg_mass, mean3, sigma3);

	// Recursive fractions: Each coefficient is interpreted as the fraction of the left-hand component of the i-th recursive sum (as it is used in flashggFInalFit)
	RooAddPdf sig_model("sig_model", "sig_model", RooArgList(gauss1,gauss2,gauss3), RooArgList(frac1,frac2), true); 	
	//RooAddPdf sig_model("sig_model", "sig_model", RooArgList(gauss1,gauss2,gauss3), RooArgList(frac1,frac2));

	/* NOTE about the mean: best fit from flashggFinalFit given as dm0, dm1, dm2
	   RooFormulaVar::mean_g0_GG2H_2018_UntaggedTag_0_13TeV[ actualVars=(MH,dm_g0_GG2H_2018_UntaggedTag_0_13TeV) formula="(@0+@1)" ] = 70.0027
	   RooFormulaVar::mean_g1_GG2H_2018_UntaggedTag_0_13TeV[ actualVars=(MH,dm_g1_GG2H_2018_UntaggedTag_0_13TeV) formula="(@0+@1)" ] = 69.4859
	   RooFormulaVar::mean_g2_GG2H_2018_UntaggedTag_0_13TeV[ actualVars=(MH,dm_g2_GG2H_2018_UntaggedTag_0_13TeV) formula="(@0+@1)" ] = 67.0678*/

	//-------------------------
	// Fitting
	//-------------------------
	//sig_model.fitTo(data_sig);
	RooFitResult* result = sig_model.fitTo(data_sig, Range(fit_min, fit_max), Save());
	//RooFitResult* result = sig_model.fitTo(data_sig, Save());
	cout << "#########################" << endl;
	cout << "##### PRINT RESULTS #####" << endl;
	cout << "#########################" << endl;
	result->Print("v");
	TH2* hcorr = result->correlationHist();
	TFile* sig_ggh_out = new TFile(Form("tests_freezingParams_signalModel/sig_ggh_M%d.root", int(v_mass.at(massindex))), "RECREATE");
	result->Write();
	hcorr->Write();
	sig_ggh_out->Close();

	// Each RooRealVar is set to be a constant after the result of the fit
	mean1.setConstant(kTRUE);
	mean2.setConstant(kTRUE);
	mean3.setConstant(kTRUE);
	sigma1.setConstant(kTRUE);
	sigma2.setConstant(kTRUE);
	sigma3.setConstant(kTRUE);
	frac1.setConstant(kTRUE);
	frac2.setConstant(kTRUE);
	
	//-------------------------
	// RooPlot
	//-------------------------	
	RooPlot *sig_frame = CMS_Hgg_mass.frame();                                             

	std::cout << "################" << std::endl;
	std::cout << "##### CHI2 #####" << std::endl;
	std::cout << "################" << std::endl;
	//RooChi2Var chi2_sig_model("chi2_sig_model","chi2", sig_model, data_sig, RooFit::DataError(RooAbsData::Expected));
	//std::cout << chi2_sig_model.getVal() << std::endl;
	//double chi2_over_ndof = chi2_sig_model.getVal() / (fit_nbins - result->floatParsFinal().getSize());
	//std::cout << "Chi2/ndof = " << chi2_over_ndof << std::endl;
	
	sig_frame->SetTitle("");
	sig_frame->GetXaxis()->SetTitle("Mass [GeV]");
	sig_frame->GetYaxis()->SetTitle("Events/0.1");
	data_sig.plotOn(sig_frame);                                  
	sig_model.plotOn(sig_frame, RooFit::Name("SignalModel"), LineColor(kRed));

	cout << "-----------------------------" << endl;
	cout << "# TH1D Histo bins: " << n_bins << endl;
	cout << "# NBins data_sig = " << (data_sig.numEntries() -2) << endl;
	// Print number of floating parameters in the PDF
	//cout << "# Number of floating parameters: " << result->floatParsFinal().getSize() << endl;
	double chi2_over_ndof = sig_frame->chiSquare(result->floatParsFinal().getSize()); //tell the chiSquare method how many free parameters you have
	cout << "Chi2/ndof = " << chi2_over_ndof << endl;
	cout << "-----------------------------" << endl;

	//-------------------------
	// Plotting
	//-------------------------	
	TCanvas c_sig("c_sig", "c_sig", 1200, 1000);                
	sig_frame->Draw("");

	TLegend *leg = new TLegend(0.88,0.88,0.6,0.76);
	leg->SetHeader(Form("#chi^{2}/N_{dof} = %.4f", chi2_over_ndof));
	//leg->SetHeader(Form("#chi^{2} = %.2f", chi2_sig_model.getVal()));
	leg->SetLineColor(kWhite);
	leg->SetFillColor(0);
	leg->SetTextFont(42);
	leg->SetTextSize(0.03);
	leg->AddEntry(sig_frame->findObject("SignalModel"),"ggH Signal Model","l");
	leg->Draw();
	//TLegendEntry *header = (TLegendEntry*)leg->GetListOfPrimitives()->First();
	//header->SetTextSize(.03);

	// Add text to frame
	/*
	TLatex *latex = new TLatex(); // prepare text in LaTeX format
	latex->SetTextSize(0.035);
	latex->SetNDC();
	latex->DrawLatex(0.6, 0.9, Form("#frac{#chi^{2}}{N_{dof}} = %.2f", chi2_sig_model.getVal())); 
	sig_frame->addObject(latex);
	*/
	c_sig.SaveAs(Form("tests_freezingParams_signalModel/sig_ggH_cat0_"+year+"_trigauss_mass%d.png", int(v_mass.at(massindex))));

	/*
	//-------------------------
	// Save into ROO workspace
	//-------------------------
	RooWorkspace dpworkspace("dpworkspace", "");
	dpworkspace.import(sig_model);
	dpworkspace.writeToFile(Form("output/dpWorkspace"+year+suff+"_%d_new_wgt.root",i));
	*/
}

// ----------------------------------------
// Additional notes and tests!
// ----------------------------------------

	// A d d   b o x   w i t h   p d f   p a r a m e t e r s 
	// -----------------------------------------------------
	// Left edge of box starts at 55% of Xaxis)
	//sig_model.paramOn(sig_frame,Layout(0.85)) ;
	
	// A d d   b o x   w i t h   d a t a   s t a t i s t i c s
	// ------------------------------------------------------- 
	// X size of box is from 55% to 99% of Xaxis range, top of box is at 80% of Yaxis range)
	//data_sig.statOn(sig_frame,Layout(0.6,0.9,0.9)) ;
	

	/*
	  RooRealVar mean1("mean1", "mean1", 65, 60, 80);
	  RooRealVar sigma1("sigma1", "sigma1", 0.3, 0.001, 1);
	  RooRealVar mean2("mean2", "mean2", 70, 60, 80);
	  RooRealVar sigma2("sigma2", "sigma2", 0.3, 0.001, 1);
	  RooRealVar mean3("mean3", "mean3", 75, 60, 80);
	  RooRealVar sigma3("sigma3", "sigma3", 1, 0.001, 3);
	  RooRealVar frac1("frac1", "frac1", 0.45, 0, 1);
	  RooRealVar frac2("frac2", "frac2", 0.45, 0, 1);
	  
	  RooRealVar mean1("mean1", "mean1", 25, 23, 27); #25
	  RooRealVar sigma1("sigma1", "sigma1", 1., 0.01, 1.); #25
	  RooRealVar mean1("mean1", "mean1", 20., 18., 22.); #20
	  RooRealVar mean1("mean1", "mean1", 15., 13., 27.); #15
	  RooRealVar mean3("mean1", "mean3", 10., 9., 11.); #10
	  RooRealVar sigma1("sigma1", "sigma1", 0.2, 0.01, 1.);#10
	  RooRealVar mean1("mean1", "mean1", 5., 4.5, 5.5);
	*/
