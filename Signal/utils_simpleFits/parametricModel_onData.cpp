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

void parametricModel_onData(TString year="2018"){
  const int massindex = 5;
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
	//EF RooRealVar CMS_Hgg_mass("CMS_Hgg_mass", "CMS_Hgg_mass", massLow, massHigh);

	// Linear interpolation from simple fit results with mean values fixed - with (first line) and without (second line) error bars considered) 	
	const vector<float> v_frac0 = {
	  0.7196685671806335, 0.7083272337913513, 0.6969859004020691, 0.6856445670127869, 0.6743032336235046, 0.6629619002342224, 0.5661643147468567, 0.5628259181976318, 0.559487521648407, 0.5561491250991821, 0.5528107285499573, 0.5494723320007324, 0.5461339354515076, 0.5427955389022827
	};
	const vector<float> v_frac1 = {
	  0.7160645127296448, 0.6962210536003113, 0.676377534866333, 0.6565340757369995, 0.636690616607666, 0.6168471574783325, 0.5660428404808044, 0.5315342545509338, 0.49702566862106323, 0.462517112493515, 0.4280085265636444, 0.3934999406337738, 0.3589913845062256, 0.324482798576355
	};
	const vector<float> v_dm0 = {
	  -0.03291655704379082, -0.04268718138337135, -0.05245780944824219, -0.06222843378782272, -0.07199905812740326, -0.0817696824669838, 0.038349710404872894, 0.038380302488803864, 0.03841089457273483, 0.0384414866566658, 0.03847207501530647, 0.03850266709923744, 0.03853325918316841, 0.03856385126709938
	};
	const vector<float> v_dm1 = {
	  -0.3752392530441284, -0.3684430718421936, -0.3616468906402588, -0.35485073924064636, -0.34805455803871155, -0.34125837683677673, -0.6976693868637085, -0.6839905381202698, -0.670311689376831, -0.6566327810287476, -0.6429539322853088, -0.6292750835418701, -0.6155961751937866, -0.6019173264503479
	};
	const vector<float> v_dm2 = {
	  0.041547831147909164, 0.04119031876325607, 0.04083280265331268, 0.04047529026865959, 0.0401177741587162, 0.03976025804877281, -0.21418417990207672, -0.2641841769218445, -0.31418418884277344, -0.36418417096138, -0.41418418288230896, -0.4641841650009155, -0.5141841769218445, -0.5641841888427734
	};
	const vector<float> v_sig0 = {
	  0.09639949351549149, 0.16336040198802948, 0.23032130300998688, 0.29728221893310547, 0.36424311995506287, 0.43120402097702026, 0.4471903145313263, 0.47340548038482666, 0.499620646238327, 0.5258358120918274, 0.5520510077476501, 0.5782661437988281, 0.6044813394546509, 0.6306964755058289
	};
	const vector<float> v_sig1 = {
	  0.41178637742996216, 0.5703612565994263, 0.7289361357688904, 0.8875110149383545, 1.0460858345031738, 1.2046607732772827, 1.0499999523162842, 1.0, 0.949999988079071, 0.8999999761581421, 0.8500000238418579, 0.800000011920929, 0.75, 0.699999988079071
	};
	const vector<float> v_sig2 = {
	  0.46699094772338867, 0.4169909954071045, 0.3669910430908203, 0.31699109077453613, 0.26699113845825195, 0.21699120104312897, 0.319436252117157, 0.5129575133323669, 0.7064787745475769, 0.8999999761581421, 1.093521237373352, 1.287042498588562, 1.480563759803772, 1.674085021018982
	};

	RooRealVar mean1("mean1", "mean1", v_mass[massindex], massLow, massHigh);
	RooRealVar mean2("mean2", "mean2", v_mass[massindex], massLow, massHigh);
	RooRealVar mean3("mean3", "mean3", v_mass[massindex], massLow, massHigh);
	RooRealVar sigma1("sigma1", "sigma1", 0.1, 0.1, 1.); //Extracting values with mean fixed 
	RooRealVar sigma2("sigma2", "sigma2", 0.4, 0.1, 4.); //Extracting values with mean fixed
	RooRealVar sigma3("sigma3", "sigma3", 0.4, 0.1, 7.); //Extracting values with mean fixed
	//RooRealVar sigma1("sigma1", "sigma1", 0.2, 0.1, 1.);
	//RooRealVar sigma2("sigma2", "sigma2", 0.1, 0.3, 8.);
	//RooRealVar sigma3("sigma3", "sigma3", 0.1, 0.3, 10.);
	RooRealVar frac1("frac1", "frac1", 0.45, 0., 1.2);
	RooRealVar frac2("frac2", "frac2", 0.45, 0., 1.);

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

	
	std::cout << "PARAMETERS USED: " << std::endl;
	std::cout << "mean0: " << v_mass.at(massindex) + v_dm0.at(massindex) << std::endl;
	std::cout << "mean1: " << v_mass.at(massindex) + v_dm1.at(massindex) << std::endl;
	std::cout << "mean2: " << v_mass.at(massindex) + v_dm2.at(massindex) << std::endl;
	//std::cout << "mean0: " << v_dm0.at(massindex) << std::endl;
	//std::cout << "mean1: " << v_dm1.at(massindex) << std::endl;
	//std::cout << "mean2: " << v_dm2.at(massindex) << std::endl;
	std::cout << "sig0: " << v_sig0.at(massindex) << std::endl;
	std::cout << "sig1: " << v_sig1.at(massindex) << std::endl;
	std::cout << "sig2: " << v_sig2.at(massindex) << std::endl;
	std::cout << "frac0: " << v_frac0.at(massindex) << std::endl;
	std::cout << "frac1: " << v_frac1.at(massindex) << std::endl;	       
	
	
	/* NOTE about the mean: best fit from flashggFinalFit given as dm0, dm1, dm2
	   RooFormulaVar::mean_g0_GG2H_2018_UntaggedTag_0_13TeV[ actualVars=(MH,dm_g0_GG2H_2018_UntaggedTag_0_13TeV) formula="(@0+@1)" ] = 70.0027
	   RooFormulaVar::mean_g1_GG2H_2018_UntaggedTag_0_13TeV[ actualVars=(MH,dm_g1_GG2H_2018_UntaggedTag_0_13TeV) formula="(@0+@1)" ] = 69.4859
	   RooFormulaVar::mean_g2_GG2H_2018_UntaggedTag_0_13TeV[ actualVars=(MH,dm_g2_GG2H_2018_UntaggedTag_0_13TeV) formula="(@0+@1)" ] = 67.0678*/
	
	// Each RooRealVar is set to be a constant
	mean1.setConstant(kTRUE);
	mean2.setConstant(kTRUE);
	mean3.setConstant(kTRUE);
	sigma1.setConstant(kTRUE);
	sigma2.setConstant(kTRUE);
	sigma3.setConstant(kTRUE);
	frac1.setConstant(kTRUE);
	frac2.setConstant(kTRUE);
	//mean1.setVal(v_dm0.at(massindex));
	//mean2.setVal(v_dm1.at(massindex));
	//mean3.setVal(v_dm2.at(massindex));
	mean1.setVal(v_mass.at(massindex) + v_dm0.at(massindex));
	mean2.setVal(v_mass.at(massindex) + v_dm1.at(massindex));
	mean3.setVal(v_mass.at(massindex) + v_dm2.at(massindex));
	sigma1.setVal(v_sig0.at(massindex));
	sigma2.setVal(v_sig1.at(massindex));
	sigma3.setVal(v_sig2.at(massindex));
	frac1.setVal(v_frac0.at(massindex));
	frac2.setVal(v_frac1.at(massindex));
	
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
	sig_model.plotOn(sig_frame, RooFit::Name("SignalModel"), LineColor(kBlue));

	
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
	//EF leg->SetHeader("#chi^{2}/N_{dof} = ...");
	leg->SetHeader(Form("#chi^{2}/N_{dof} = %.3f", chi2_over_ndof));
	leg->SetLineColor(kWhite);
	leg->SetFillColor(0);
	leg->SetTextFont(42);
	leg->SetTextSize(0.03);
	leg->AddEntry(sig_frame->findObject("SignalModel"),"ggH Signal Model","l");
	leg->Draw();
	//TLegendEntry *header = (TLegendEntry*)leg->GetListOfPrimitives()->First();
	//header->SetTextSize(.03);

	c_sig.SaveAs(Form("/eos/user/e/elfontan/www/LowMassDiPhoton/Feb2023_test/modelling/twoRanges_flashggFinalFit/Mar30_sig_ggH_cat0_"+year+"_trigauss_mass%d.png", int(v_mass.at(massindex))));

	/*
	//-------------------------
	// Save into ROO workspace
	//-------------------------
	RooWorkspace dpworkspace("dpworkspace", "");
	dpworkspace.import(sig_model);
	dpworkspace.writeToFile(Form("output/dpWorkspace"+year+suff+"_%d_new_wgt.root",i));
	*/
}


	/*
	// FlashggFinalFit: linear interpolation from simple fit results 	 
	const vector<float> v_frac0 = {
	  0.6709164285714286, 0.6564129560439561, 0.6419094835164835, 0.6274060109890109, 0.6129025384615384, 0.5983990659340659, 0.5838955934065932, 0.5693921208791207, 0.5548886483516482, 0.5403851758241757, 0.5258817032967031, 0.5113782307692305, 0.496874758241758, 0.48237128571428545
	};
	const vector<float> v_frac1 = {
	  0.3971574411428576, 0.41470643503296745, 0.4322554289230773, 0.4498044228131871, 0.467353416703297, 0.4849024105934068, 0.5024514044835167, 0.5200003983736264, 0.5375493922637363, 0.5550983861538461, 0.572647380043956, 0.5901963739340659, 0.6077453678241757, 0.6252943617142854
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
	*/
