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

void signalModel(TString year="2018"){
  const int massindex = 13;
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
	//RooRealVar nBins("nBins", "nBins", 10);

        const vector<float> v_frac0 = {
          0.9190168976783752, 0.9034169912338257, 0.8878171443939209, 0.8722172975540161, 0.8566174507141113, 0.8410176038742065, 0.825417697429657, 0.8098178505897522, 0.7942180037498474, 0.7786181569099426, 0.7630183100700378, 0.7474184036254883, 0.7318185567855835, 0.7162187099456787
        };
        const vector<float> v_frac1 = {
          0.4959886074066162, 0.5209783911705017, 0.5459681749343872, 0.5709580183029175, 0.595947802066803, 0.6209375858306885, 0.645927369594574, 0.6709171533584595, 0.695906937122345, 0.7208967804908752, 0.7458865642547607, 0.7708763480186462, 0.7958661317825317, 0.8208559155464172
        };
        const vector<float> v_dm0 = {
          -0.09266149252653122, -0.08413475006818771, -0.07560800760984421, -0.0670812651515007, -0.058554526418447495, -0.05002778396010399, -0.04150104150176048, -0.03297429904341698, -0.02444756031036377, -0.015920817852020264, -0.00739407679066062, 0.0011326646199449897, 0.009659405797719955, 0.018186146393418312
        };
        const vector<float> v_dm1 = {
          -1.4578336477279663, -1.4222248792648315, -1.3866159915924072, -1.351007103919983, -1.3153983354568481, -1.2797894477844238, -1.2441805601119995, -1.2085717916488647, -1.1729629039764404, -1.1373540163040161, -1.1017452478408813, -1.066136360168457, -1.0305275917053223, -0.994918704032898
        };
        const vector<float> v_dm2 = {
          -0.23807260394096375, -0.19001522660255432, -0.1419578641653061, -0.09390048682689667, -0.04584311693906784, 0.0022142515517771244, 0.05027162283658981, 0.09832899272441864, 0.14638635516166687, 0.1944437325000763, 0.24250109493732452, 0.29055845737457275, 0.3386158347129822, 0.3866732120513916
        };
        const vector<float> v_sig0 = {
          0.22452326118946075, 0.26481401920318604, 0.3051047623157501, 0.3453955054283142, 0.3856862783432007, 0.42597702145576477, 0.46626776456832886, 0.5065585374832153, 0.546849250793457, 0.5871400237083435, 0.62743079662323, 0.6677215099334717, 0.7080122828483582, 0.7483029961585999
        };
        const vector<float> v_sig1 = {
          0.9823096990585327, 0.9824902415275574, 0.982670783996582, 0.9828513264656067, 0.9830319285392761, 0.9832124710083008, 0.9833930134773254, 0.9835735559463501, 0.9837540984153748, 0.9839346408843994, 0.9841151833534241, 0.9842957258224487, 0.9844762682914734, 0.9846568703651428
        };
        const vector<float> v_sig2 = {
          8.659833908081055, 8.659462928771973, 8.65909194946289, 8.658720016479492, 8.65834903717041, 8.657978057861328, 8.65760612487793, 8.657235145568848, 8.656864166259766, 8.656493186950684, 8.656121253967285, 8.655750274658203, 8.655379295349121, 8.655007362365723
        };
	// Good set of initial parameters for simple fits using all mass points
	
	RooRealVar mean1("mean1", "mean1", v_mass[massindex], massLow, massHigh);
	RooRealVar mean2("mean2", "mean2", v_mass[massindex], massLow, massHigh);
	RooRealVar mean3("mean3", "mean3", v_mass[massindex], massLow, massHigh);
	RooRealVar sigma1("sigma1", "sigma1", 0.2, 0.1, 3.); 
	RooRealVar sigma2("sigma2", "sigma2", 0.1, 0.3, 8.); 
	RooRealVar sigma3("sigma3", "sigma3", 0.1, 0.3, 10.); 
	RooRealVar frac1("frac1", "frac1", 0.45, 0., 1.2);
	RooRealVar frac2("frac2", "frac2", 0.45, 0., 1.);
	
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
	TFile* sig_ggh_out = new TFile(Form("tests_signalModel/sig_ggh_M%d.root", int(v_mass.at(massindex))), "RECREATE");
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
	c_sig.SaveAs(Form("tests_signalModel/sig_ggH_cat0_"+year+"_trigauss_mass%d.png", int(v_mass.at(massindex))));

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
