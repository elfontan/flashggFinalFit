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
	//EF RooRealVar CMS_Hgg_mass("CMS_Hgg_mass", "CMS_Hgg_mass", massLow, massHigh);

	// Linear interpolation from simple fit results with mean values fixed - with (first line) and without (second line) error bars considered) 	
	const vector<float> v_frac0 = {
	  0.8681726633568282, 0.832518282442386, 0.796863901527944, 0.7612095206135018, 0.7255551396990597, 0.6899007587846175, 0.6542463778701755, 0.6185919969557334, 0.5829376160412914, 0.5472832351268492, 0.511628854212407, 0.47597447329796494, 0.4403200923835228, 0.4046657114690807
	  //0.7796216285714281, 0.7541717956043952, 0.7287219626373622, 0.7032721296703293, 0.6778222967032963, 0.6523724637362635, 0.6269226307692305, 0.6014727978021975, 0.5760229648351646, 0.5505731318681316, 0.5251232989010988, 0.49967346593406575, 0.4742236329670328, 0.4487737999999999
	};
	const vector<float> v_frac1 = {
	  //0.5740844439816892, 0.5449233610721715, 0.5157622781626537, 0.486601195253136, 0.45744011234361825, 0.4282790294341005, 0.39911794652458277, 0.36995686361506497, 0.3407957807055472, 0.3116346977960295, 0.28247361488651174, 0.253312531976994, 0.22415144906747625, 0.19499036615795845
	  0.6311067999999996, 0.606760193406593, 0.5824135868131864, 0.5580669802197799, 0.5337203736263733, 0.5093737670329666, 0.48502716043956, 0.4606805538461534, 0.4363339472527468, 0.41198734065934023, 0.38764073406593363, 0.363294127472527, 0.3389475208791204, 0.3146009142857138
	};
	const vector<float> v_dm0 = {
	  //4.991263255157044, 9.986574217800262, 14.98188518044348, 19.9771961430867, 24.97250710572992, 29.967818068373138, 34.96312903101636, 39.958439993659574, 44.95375095630279, 49.949061918946015, 54.94437288158923, 59.93968384423245, 64.93499480687566, 69.93030576951888
	  4.997394285714296, 9.993632527472538, 14.989870769230782, 19.986109010989026, 24.98234725274727, 29.978585494505513, 34.97482373626375, 39.97106197802199, 44.967300219780235, 49.963538461538484, 54.959776703296725, 59.95601494505497, 64.95225318681322, 69.94849142857146
	};
	const vector<float> v_dm1 = {
	  //4.999999999995275, 9.95868046303297, 14.917360926070664, 19.87604138910836, 24.834721852146053, 29.793402315183748, 34.75208277822144, 39.71076324125914, 44.66944370429683, 49.628124167334526, 54.58680463037222, 59.545485093409916, 64.5041655564476, 69.4628460194853
	  4.453346285714301, 9.427265978021996, 14.401185670329689, 19.375105362637385, 24.349025054945077, 29.32294474725277, 34.296864439560466, 39.27078413186816, 44.24470382417585, 49.21862351648355, 54.19254320879124, 59.16646290109893, 64.14038259340663, 69.11430228571432
	};
	const vector<float> v_dm2 = {
	  //5.002094631622204, 10.006132152607107, 15.010169673592008, 20.014207194576915, 25.018244715561817, 30.02228223654672, 35.02631975753162, 40.030357278516526, 45.03439479950143, 50.03843232048633, 55.04246984147123, 60.04650736245613, 65.05054488344105, 70.05458240442594
	  5.693821142857154, 10.612516351648363, 15.531211560439573, 20.449906769230783, 25.368601978021992, 30.287297186813202, 35.20599239560441, 40.124687604395625, 45.04338281318683, 49.962078021978044, 54.88077323076925, 59.799468439560464, 64.71816364835166, 69.63685885714287
	};
	const vector<float> v_sig0 = {
	  0.14777673426893626, 0.18997813114682544, 0.2321795280247146, 0.27438092490260374, 0.3165823217804929, 0.35878371865838204, 0.40098511553627125, 0.4431865124141604, 0.48538790929204956, 0.5275893061699387, 0.5697907030478279, 0.611992099925717, 0.6541934968036062, 0.6963948936814954
	  //0.14750411428571442, 0.1896552505494507, 0.23180638681318694, 0.2739575230769232, 0.31610865934065946, 0.35825979560439575, 0.40041093186813204, 0.44256206813186827, 0.4847132043956045, 0.5268643406593407, 0.5690154769230771, 0.6111666131868133, 0.6533177494505495, 0.6954688857142859
	};
	const vector<float> v_sig1 = {
	  0.18679170685846533, 0.3006461641563175, 0.4145006214541697, 0.5283550787520219, 0.6422095360498741, 0.7560639933477262, 0.8699184506455785, 0.9837729079434306, 1.097627365241283, 1.211481822539135, 1.3253362798369872, 1.4391907371348394, 1.5530451944326915, 1.6668996517305439
	  //0.4258284857142858, 0.5048287076923077, 0.5838289296703296, 0.6628291516483517, 0.7418293736263736, 0.8208295956043956, 0.8998298175824175, 0.9788300395604395, 1.0578302615384614, 1.1368304835164833, 1.2158307054945054, 1.2948309274725274, 1.3738311494505493, 1.4528313714285712
	};
	const vector<float> v_sig2 = {
	  0.23889849976540423, 0.28269520700437284, 0.32649191424334145, 0.37028862148231007, 0.4140853287212787, 0.4578820359602473, 0.5016787431992159, 0.5454754504381845, 0.5892721576771531, 0.6330688649161218, 0.6768655721550904, 0.7206622793940589, 0.7644589866330276, 0.8082556938719963
	  //0.5389454571428572, 0.5700658813186814, 0.6011863054945056, 0.6323067296703297, 0.6634271538461539, 0.6945475780219781, 0.7256680021978023, 0.7567884263736264, 0.7879088505494507, 0.8190292747252748, 0.850149698901099, 0.8812701230769231, 0.9123905472527474, 0.9435109714285714
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
	//std::cout << "mean0: " << v_mass.at(massindex) + v_dm0.at(massindex) << std::endl;
	//std::cout << "mean1: " << v_mass.at(massindex) + v_dm1.at(massindex) << std::endl;
	//std::cout << "mean2: " << v_mass.at(massindex) + v_dm2.at(massindex) << std::endl;
	std::cout << "mean0: " << v_dm0.at(massindex) << std::endl;
	std::cout << "mean1: " << v_dm1.at(massindex) << std::endl;
	std::cout << "mean2: " << v_dm2.at(massindex) << std::endl;
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
	mean1.setVal(v_dm0.at(massindex));
	mean2.setVal(v_dm1.at(massindex));
	mean3.setVal(v_dm2.at(massindex));
	//mean1.setVal(v_mass.at(massindex) + v_dm0.at(massindex));
	//mean2.setVal(v_mass.at(massindex) + v_dm1.at(massindex));
	//mean3.setVal(v_mass.at(massindex) + v_dm2.at(massindex));
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

	c_sig.SaveAs(Form("tests_freezingParams_signalModel/closure/sig_ggH_cat0_"+year+"_trigauss_mass%d.png", int(v_mass.at(massindex))));

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
