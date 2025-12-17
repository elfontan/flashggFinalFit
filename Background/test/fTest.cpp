#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>

#include "boost/program_options.hpp"
#include "boost/lexical_cast.hpp"

#include "TFile.h"
#include "TMath.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "RooWorkspace.h"
#include "RooDataSet.h"
#include "RooHist.h"
#include "RooAbsData.h"
#include "RooGaussian.h"
#include "RooAbsPdf.h"
#include "RooAddPdf.h"
#include "RooArgSet.h"
#include "RooFitResult.h"
#include "RooMinimizer.h"
#include "RooMsgService.h"
#include "RooDataHist.h"
#include "RooExtendPdf.h"
#include "RooChi2Var.h"
#include "TRandom3.h"
#include "TLatex.h"
#include "TMacro.h"
#include "TH1F.h"
#include "TH1I.h"
#include "TArrow.h"
#include "TKey.h"

#include "RooCategory.h"
#include "HiggsAnalysis/CombinedLimit/interface/RooMultiPdf.h"
#include "HiggsAnalysis/CombinedLimit/interface/RooDoubleCBFast.h"
#include "RooProdPdf.h"
#include "RooGenericPdf.h"

#include "../interface/PdfModelBuilder.h"
#include <Math/PdfFuncMathCore.h>
#include <Math/ProbFunc.h>
#include <iomanip>
#include "boost/program_options.hpp"
#include "boost/algorithm/string/split.hpp"
#include "boost/algorithm/string/classification.hpp"
#include "boost/algorithm/string/predicate.hpp"

#include "../../tdrStyle/tdrstyle.C"
#include "../../tdrStyle/CMS_lumi.C"

using namespace std;
using namespace RooFit;
using namespace boost;

namespace po = program_options;

bool runFtestCheckWithToys=false;

// Mass range configuration
float mN = 2.75;
float sigma = 0.025;
int nsigma = 10;
float mN_low  = 80;
float mN_high = 800;

// not configured in main
int nBinsForFit  = 4*(mN_high-mN_low); // kept baseline values for Hgg 
int nBinsForPlot = 4*(mN_high-mN_low);  // ""

RooRealVar *intLumi_ = new RooRealVar("IntLumi","hacked int lumi", 1000.);

TRandom3 *RandomGen = new TRandom3();

RooAbsPdf* getPdf(PdfModelBuilder &pdfsModel, string type, int order, const char* ext=""){
/* Construct the pdf model using the PdfModelBuilder */
  
  if (type=="Bernstein") return pdfsModel.getBernstein(Form("%s_bern%d",ext,order),order); 
  else if (type=="Exponential") return pdfsModel.getExponentialSingle(Form("%s_exp%d",ext,order),order); 
  else if (type=="ExponentialSum") return pdfsModel.getExponential(Form("%s_expsum%d",ext,order),order); 
  else if (type=="PowerLaw") return pdfsModel.getPowerLawSingle(Form("%s_pow%d",ext,order),order); 
  else if (type=="PowerLawSingle") return pdfsModel.getPowerLawSingle(Form("%s_powsing%d",ext,order),order); 
  else if (type=="PowerLawSum") return pdfsModel.getPowerLaw(Form("%s_pow%d",ext,order),order); 
  else if (type=="PowerLawGeneric") return pdfsModel.getPowerLawGeneric(Form("%s_powgen%d",ext,order),order); 
  else if (type=="Laurent") return pdfsModel.getLaurentSeries(Form("%s_lau%d",ext,order),order); 
  else if (type=="Chebychev") return pdfsModel.getChebychev(Form("%s_cheb%d",ext,order),order); 
  else {
    cerr << "[ERROR] -- getPdf() -- type " << type << " not recognised." << endl;
    return NULL;
  }
}
#include "RooChi2Var.h"

struct Chi2Result { double chi2, chi2red; int ndof; double pval; };

// Restituisce una pdf estesa pronta per il χ² (possibilmente wrappando quella in ingresso)
// Se trova 'bkg_norm' e/o 'z_norm' li somma e li usa come yield totale; altrimenti scala sui sideband.
std::unique_ptr<RooAbsPdf> makeExtendedForGOF(RooAbsPdf* pdf, RooRealVar* mass,
                                              RooAbsData* data, const char* rangeList) {
  // Se è già extended (RooAddPdf con yields o RooExtendPdf), usa così com’è
  if (pdf->InheritsFrom("RooExtendPdf")) {
    return std::unique_ptr<RooAbsPdf>((RooAbsPdf*)pdf); // ATTENZIONE: non deletare il pdf originale!
  }
  // Se è RooAddPdf con yields, è comunque extended: non wrappare
  if (pdf->InheritsFrom("RooAddPdf")) {
    // heuristica: se ha una lista di coefficienti che NON somma a 1 è extended; ma non abbiamo accesso semplice qui
    // più semplice: prova a chiamare expectedEvents; se >0, trattala come extended.
    double nexp = pdf->expectedEvents(RooArgSet(*mass));
    if (nexp > 0) return std::unique_ptr<RooAbsPdf>((RooAbsPdf*)pdf);
  }

  // Altrimenti costruiamo una RooExtendPdf con una norma sensata
  // 1) prova a leggere una 'norm' dal modello (es. bkg_norm, z_norm)
  double norm_from_fit = 0.0;
  std::unique_ptr<RooArgSet> pars(pdf->getParameters((const RooArgSet*)0));
  TIterator* it = pars->createIterator();
  while (RooAbsArg* a = (RooAbsArg*)it->Next()) {
    TString n = a->GetName();
    if (n.Contains("bkg_norm") || n.Contains("z_norm") || n.Contains("nbkg")) {
      if (auto* rr = dynamic_cast<RooRealVar*>(a)) norm_from_fit += rr->getVal();
    }
  }
  delete it;

  RooRealVar* normVar = nullptr;
  if (norm_from_fit > 0) {
    normVar = new RooRealVar("__gof_norm","__gof_norm", norm_from_fit, 0., 1e15);
    normVar->setConstant(kTRUE);
  } else {
    // 2) fallback: scala in modo che attesi nei sideband = osservati nei sideband
    std::unique_ptr<RooAbsReal> frac_pdf(
      pdf->createIntegral(RooArgSet(*mass), RooFit::NormSet(RooArgSet(*mass)),
                          RooFit::Range(rangeList))
    );
    RooDataHist* dh = dynamic_cast<RooDataHist*>(data);
    std::unique_ptr<RooDataHist> dh_owner;
    if (!dh) { dh_owner.reset(new RooDataHist("__gof_dh","__gof_dh", RooArgSet(*mass), *data)); dh = dh_owner.get(); }
    double n_sb = dh->sumEntries(nullptr, rangeList);
    double f_sb = frac_pdf ? frac_pdf->getVal() : 0.;
    double norm_val = (f_sb>0) ? (n_sb / f_sb) : n_sb; // robust
    normVar = new RooRealVar("__gof_norm","__gof_norm", norm_val, 0., 1e15);
    normVar->setConstant(kTRUE);
  }
  return std::unique_ptr<RooAbsPdf>(new RooExtendPdf("__gof_ext","__gof_ext", *pdf, *normVar));
}

Chi2Result computeChi2ReducedExtended(RooRealVar* mass, RooAbsPdf* pdf_in, RooAbsData* data,
                                      const char* rangeList) {
    // Step 1: assicurati di avere RooDataHist
    RooDataHist* dh = dynamic_cast<RooDataHist*>(data);
    std::unique_ptr<RooDataHist> dh_owner;
    if (!dh) {
        dh_owner.reset(new RooDataHist("__chi2_dh","__chi2_dh", RooArgSet(*mass), *data));
        dh = dh_owner.get();
    }

    // Step 2: decidi quale pdf usare
    RooAbsPdf* pdf_ext = nullptr;
    std::unique_ptr<RooExtendPdf> owned_ext; // tiene in vita se la creiamo qui

    if (pdf_in->InheritsFrom("RooExtendPdf")) {
        pdf_ext = pdf_in; // già extended → usiamo direttamente
    } else {
        // creo la norm dal fit (se disponibile) o dai sideband
        double norm_val = pdf_in->expectedEvents(RooArgSet(*mass));
        if (norm_val <= 0) {
            // calcola scala dai sideband
            auto frac_pdf = std::unique_ptr<RooAbsReal>(
                pdf_in->createIntegral(RooArgSet(*mass), RooFit::NormSet(RooArgSet(*mass)),
                                       RooFit::Range(rangeList))
            );
            double n_sb = dh->sumEntries(nullptr, rangeList);
            double f_sb = frac_pdf ? frac_pdf->getVal() : 1.0;
            norm_val = (f_sb > 0) ? (n_sb / f_sb) : n_sb;
        }
        auto normVar = new RooRealVar("__gof_norm","__gof_norm", norm_val, 0., 1e15);
        normVar->setConstant(kTRUE);
        owned_ext.reset(new RooExtendPdf("__gof_ext","__gof_ext", *pdf_in, *normVar));
        pdf_ext = owned_ext.get();
    }

    // Step 3: calcola il chi² (Poisson)
    RooChi2Var chi2var("__chi2","__chi2", *pdf_ext, *dh,
                       RooFit::Range(rangeList),
                       RooFit::DataError(RooAbsData::Poisson));
    double chi2 = chi2var.getVal();

    // Step 4: conta bin usati
    auto tmp = std::unique_ptr<RooPlot>(mass->frame());
    dh->plotOn(tmp.get(), RooFit::CutRange(rangeList));
    auto* h = dynamic_cast<RooHist*>(tmp->getObject(int(tmp->numItems()-1)));
    int nbins_used = h ? h->GetN() : 0;

    // Step 5: conta parametri liberi della pdf_ext
    std::unique_ptr<RooArgSet> allPars(pdf_ext->getParameters(*dh));
    RooAbsCollection* floats = allPars->selectByAttrib("Constant", kFALSE);
    int nfloat = floats ? floats->getSize() : 0;
    delete floats;

    int ndof = nbins_used - nfloat;
    if (ndof < 1) ndof = 1;
    double chi2red = chi2 / ndof;
    double pval = TMath::Prob(chi2, ndof);

    return {chi2, chi2red, ndof, pval};
}


void runFit(RooAbsPdf *pdf, RooAbsData *data, double *NLL, int *stat_t, int MaxTries, double *chi2_out=nullptr, bool blindSignalRegion=false){
/* Basic fitting routine, fit is not extended */

  int ntries=0;
  RooArgSet *params_test = pdf->getParameters((const RooArgSet*)(0));
  
  data->Print("v");
  params_test->Print("v");
  int stat=1;
  double minnll=10e8;
  
  while (stat!=0){
    if (ntries>=MaxTries) break;
    RooFitResult *fitTest;
    
    // Choose fit range based on blinding flag
    if (blindSignalRegion) {
      std::cout << "[INFO] Fitting with signal region [100-140] GeV blinded" << std::endl;
      fitTest = pdf->fitTo(*data,RooFit::Save(1),RooFit::Minimizer("Minuit2","minimize"),RooFit::Strategy(0),RooFit::PrintLevel(-1),RooFit::Optimize(1),RooFit::SumW2Error(kFALSE),RooFit::Range("blind_low,blind_high"));
    } else {
      std::cout << "[INFO] Fitting with full mass range" << std::endl;
      fitTest = pdf->fitTo(*data,RooFit::Save(1),RooFit::Minimizer("Minuit2","minimize"),RooFit::Strategy(0),RooFit::PrintLevel(-1),RooFit::Optimize(1),RooFit::SumW2Error(kFALSE));
    }
    stat = fitTest->status();
    minnll = fitTest->minNll();
    if (stat!=0) params_test->assignValueOnly(fitTest->randomizePars());
    ntries++; 
  }
  *stat_t = stat;
  *NLL = minnll;

  
  // Calculate chi2 if requested
  if (chi2_out && stat == 0) {
     // Normalizzazione: conta solo eventi nei sidebands
    double Nvis = -1;
    RooArgSet *final_params = pdf->getParameters((const RooArgSet*)(0));
    RooRealVar* bkg_norm = nullptr;
    TIterator* iter = final_params->createIterator();
    RooAbsArg* arg;
    while ((arg = (RooAbsArg*)iter->Next())) {
      if (TString(arg->GetName()).Contains("bkg_norm")) {
        bkg_norm = (RooRealVar*)arg;
      }
    }
    delete iter;
    if (bkg_norm) {
      // Calcola l'integrale della PDF di background nel range di fit e nel range di plot
      //double fitIntegral = pdf->createIntegral(*mass, NormSet(*mass), Range("blind_low,blind_high"))->getVal();
      double normVal = bkg_norm->getVal();
      Nvis = normVal;// * (fitIntegral );
      std::cout << "[INFO] Using scaled PDF integral for Nvis: " << Nvis << " events (plot range)" << std::endl;
    } else {
      if (blindSignalRegion) {
        Nvis = data->sumEntries("1", "blind_low,blind_high");
        std::cout << "[INFO] Blinded fit: Nvis = " << Nvis << " events in sidebands" << std::endl;
      } else {
        Nvis = data->sumEntries();
        std::cout << "[INFO] Unblinded fit: Nvis = " << Nvis << " events total" << std::endl;
      }
    }
    // Create a temporary plot to calculate chi2 properly
    RooRealVar* mass = (RooRealVar*)data->get()->first();
    RooPlot* tempPlot = mass->frame();
    data->plotOn(tempPlot, Name("data"),Range("blind_low,blind_high"));
    
    // DEBUG: Check if PDF is RooAddPdf with explicit normalizations
    bool isAddPdf = pdf->InheritsFrom("RooAddPdf");
    cout << "[DEBUG] PDF type: " << pdf->ClassName() << ", isAddPdf: " << isAddPdf << endl;
    cout << "[DEBUG] Data entries: " << data->sumEntries() << endl;
    
    // Plot for chi2 calculation - let RooFit handle natural normalization
    pdf->plotOn(tempPlot, Name("pdf"),Range("blind_low,blind_high"));
    cout << "[DEBUG] Using natural normalization for chi2 calculation" << endl;
    
    int np = pdf->getParameters(*data)->getSize();
    *chi2_out = tempPlot->chiSquare("pdf", "data", np);
    cout << "[DEBUG] Chi2 calculation: np=" << np << ", chi2=" << *chi2_out << endl;
    delete tempPlot;

  }
  
  // DEBUG: Print final parameters after fit
  if (stat == 0) {
    std::cout << "=== FINAL PARAMETERS AFTER FIT ===" << std::endl;
    RooArgSet *final_params = pdf->getParameters((const RooArgSet*)(0));
    final_params->Print("v");
    
    // Print specific Z-related parameters if they exist
    RooRealVar* z_norm = nullptr;
    RooRealVar* bkg_norm = nullptr;
    
    TIterator* iter = final_params->createIterator();
    RooAbsArg* arg;
    while ((arg = (RooAbsArg*)iter->Next())) {
      if (TString(arg->GetName()).Contains("z_norm")) {
        z_norm = (RooRealVar*)arg;
      }
      if (TString(arg->GetName()).Contains("bkg_norm")) {
        bkg_norm = (RooRealVar*)arg;
      }
    }
    delete iter;
    
    if (z_norm) {
      std::cout << "Z normalization: " << z_norm->getVal() << " +/- " << z_norm->getError() << " events" << std::endl;
      
      if (bkg_norm) {
        std::cout << "Background normalization: " << bkg_norm->getVal() << " +/- " << bkg_norm->getError() << " events" << std::endl;
        double total = z_norm->getVal() + bkg_norm->getVal();
        std::cout << "Z fraction: " << (z_norm->getVal()/total)*100.0 << "%" << std::endl;
      }
    } else {
      // Fallback: look for coefficient approach
      RooRealVar* z_coeff = nullptr;
      TIterator* iter2 = final_params->createIterator();
      RooAbsArg* arg2;
      while ((arg2 = (RooAbsArg*)iter2->Next())) {
        if (TString(arg2->GetName()).Contains("z_coeff")) {
          z_coeff = (RooRealVar*)arg2;
          break;
        }
      }
      delete iter2;
      
      if (z_coeff) {
        std::cout << "Z coefficient: " << z_coeff->getVal() << " +/- " << z_coeff->getError() << std::endl;
      }
    }
    std::cout << "===================================" << std::endl;
  }
}

double getProbabilityFtest(double chi2, int ndof,RooAbsPdf *pdfNull, RooAbsPdf *pdfTest, RooRealVar *mass, RooAbsData *data, std::string name){
/* Get probability of the f-test. Currently toys are not used and only simple TMath::Prob is used. */
 
  double prob_asym = TMath::Prob(chi2,ndof);
  if (!runFtestCheckWithToys) return prob_asym;

  int ndata = data->sumEntries();
   
  RooFitResult *fitNullData;// fit the pdfs to the data and keep this fit Result (for randomizing)
  RooFitResult *fitTestData;
  
  fitNullData = pdfNull->fitTo(*data,RooFit::Save(1),RooFit::Strategy(0)
    ,RooFit::Minimizer("Minuit2","minimize"),RooFit::PrintLevel(-1),RooFit::Optimize(1));
  fitTestData = pdfTest->fitTo(*data,RooFit::Save(1),RooFit::Strategy(0)
    ,RooFit::Minimizer("Minuit2","minimize"),RooFit::PrintLevel(-1),RooFit::Optimize(1)); 

  // Ok we want to check the distribution in toys then 
  // Step 1, cache the parameters of each pdf so as not to upset anything 
  RooArgSet *params_null = pdfNull->getParameters((const RooArgSet*)(0));
  RooArgSet preParams_null;
  params_null->snapshot(preParams_null);
  RooArgSet *params_test = pdfTest->getParameters((const RooArgSet*)(0));
  RooArgSet preParams_test;
  params_test->snapshot(preParams_test);
 
  int ntoys = 100;  // Reduced for faster testing
  TCanvas *can = new TCanvas();
  can->SetLogy();
  TH1F toyhist(Form("toys_fTest_%s.pdf",pdfNull->GetName()),";Chi2;",60,-2,10);
  TH1I toyhistStatN(Form("Status_%s.pdf",pdfNull->GetName()),";FitStatus;",8,-4,4);
  TH1I toyhistStatT(Form("Status_%s.pdf",pdfTest->GetName()),";FitStatus;",8,-4,4);

  TGraph *gChi2 = new TGraph();
  gChi2->SetLineColor(kGreen+2);
  double w = toyhist.GetBinWidth(1);

  int ipoint=0;

  for (int b=0;b<0/*toyhist.GetNbinsX()*/;b++){
    double x = toyhist.GetBinCenter(b+1);
    if (x>0){
      gChi2->SetPoint(ipoint,x,(ROOT::Math::chisquared_pdf(x,ndof)));
      ipoint++;
    }
  }
  int npass =0; int nsuccesst =0;
  // mass->setBins(nBinsForFit); // Binning now taken from workspace
  for (int itoy = 0 ; itoy < ntoys ; itoy++){

    params_null->assignValueOnly(preParams_null);
    params_test->assignValueOnly(preParams_test);
    RooDataSet *binnedtoy = pdfNull->generate(RooArgSet(*mass),ndata,0,1);

    int stat_n=1;
    int stat_t=1;
    int ntries = 0;
    double nllNull,nllTest;
    // Iterate on the fit 
    int MaxTries = 2;
    while (stat_n!=0){
      if (ntries>=MaxTries) break;
      RooFitResult *fitNull;
      fitNull = pdfNull->fitTo(*binnedtoy,RooFit::Save(1),RooFit::Strategy(0)
                                              ,RooFit::Minimizer("Minuit2","minimize"),RooFit::Minos(0),RooFit::Hesse(0),RooFit::PrintLevel(-1),RooFit::Optimize(1));
      //,RooFit::Optimize(0));

      nllNull = fitNull->minNll();
      stat_n = fitNull->status();
      if (stat_n!=0) params_null->assignValueOnly(fitNullData->randomizePars());
      ntries++; 
    }
  
    ntries = 0;
    while (stat_t!=0){
      if (ntries>=MaxTries) break;
      RooFitResult *fitTest;
      fitTest = pdfTest->fitTo(*binnedtoy,RooFit::Save(1),RooFit::Strategy(0)
                                              ,RooFit::Minimizer("Minuit2","minimize"),RooFit::Minos(0),RooFit::Hesse(0),RooFit::PrintLevel(-1),RooFit::Optimize(1));
      nllTest = fitTest->minNll();
      stat_t = fitTest->status();
      if (stat_t!=0) params_test->assignValueOnly(fitTestData->randomizePars()); 
      ntries++; 
    }
       
    toyhistStatN.Fill(stat_n);
    toyhistStatT.Fill(stat_t);

    if (stat_t !=0 || stat_n !=0) continue;
    nsuccesst++;
    double chi2_t = 2*(nllNull-nllTest);
    if (chi2_t >= chi2) npass++;
        toyhist.Fill(chi2_t);

  } // end loop over toys

  double prob=0;
  if (nsuccesst!=0)  prob = (double)npass / nsuccesst;
  toyhist.Scale(1./(w*toyhist.Integral()));
  toyhist.Draw();
  TArrow lData(chi2,toyhist.GetMaximum(),chi2,0);
  lData.SetLineWidth(2);
  lData.Draw();
  gChi2->Draw("L");
  TLatex *lat = new TLatex();
  lat->SetNDC();
  lat->SetTextFont(42);
  lat->DrawLatex(0.1,0.91,Form("Prob (asymptotic) = %.4f (%.4f)",prob,prob_asym));
  std::cout << "debug get probability f test: " << name.c_str() << std::endl; 
  can->SaveAs(name.c_str());

  TCanvas *stas =new TCanvas();
  toyhistStatN.SetLineColor(2);
  toyhistStatT.SetLineColor(1); 
  TLegend *leg = new TLegend(0.2,0.6,0.4,0.87); leg->SetFillColor(0);
  leg->SetTextFont(42);
  leg->AddEntry(&toyhistStatN,"Null Hyp","L");
  leg->AddEntry(&toyhistStatT,"Test Hyp","L");
  toyhistStatN.Draw();
  toyhistStatT.Draw("same");
  leg->Draw();
  stas->SaveAs(Form("%s_fitstatus.pdf",name.c_str()));
  //reassign params
  params_null->assignValueOnly(preParams_null);
  params_test->assignValueOnly(preParams_test);

  delete can; delete stas;
  delete gChi2;
  delete leg;
  delete lat;

  // Still return the asymptotic prob (usually its close to the toys one)
  return prob_asym;

}

double getGoodnessOfFit(RooRealVar *mass, RooAbsPdf *mpdf, RooAbsData *data, std::string name){
/* Get goodness of fit, based on chi-square, using binned dataset and fitted pdf 
   use toys or chi-square distributions depending on avg number of events in bin */

  double prob;
  int ntoys = 50;   // Reduced for faster testing
  // Routine to calculate the goodness of fit. 
  name+="_gofTest.pdf";
  RooRealVar norm("norm","norm",data->sumEntries(),0,10E6);
  //norm.removeRange();

  RooExtendPdf *pdf = new RooExtendPdf("ext","ext",*mpdf,norm);

  // get The Chi2 value from the data
  RooPlot *plot_chi2 = mass->frame();
  data->plotOn(plot_chi2,Binning(nBinsForFit),Name("data"));
  pdf->plotOn(plot_chi2,Name("pdf")); // Let RooFit handle natural normalization

  int np = pdf->getParameters(*data)->getSize();

  double chi2 = plot_chi2->chiSquare("pdf","data",np);
  std::cout << "[INFO] Calculating GOF for pdf " << pdf->GetName() << ", using " <<np << " fitted parameters" <<std::endl;

  // The first thing is to check if the number of entries in any bin is < 5 
  // if so, we don't rely on asymptotic approximations
 
  if ((double)data->sumEntries()/nBinsForFit < 5 ){

    std::cout << "[INFO] Running toys for GOF test " << std::endl;
    // store pre-fit params 
    RooArgSet *params = pdf->getParameters(*data);
    RooArgSet preParams;
    params->snapshot(preParams);
    int ndata = data->sumEntries();
 
    int npass =0;
    std::vector<double> toy_chi2;
    for (int itoy = 0 ; itoy < ntoys ; itoy++){
      std::cout << "[INFO] " <<Form("\t.. %.1f %% complete\r",100*float(itoy)/ntoys) << std::flush;
      params->assignValueOnly(preParams);
      int nToyEvents = RandomGen->Poisson(ndata);
      RooDataSet *binnedtoy = pdf->generate(RooArgSet(*mass),nToyEvents,0,1);
      //RooDataHist *binnedtoy = pdf->generateBinned(RooArgSet(*mass),nToyEvents,0,1);
//      pdf->fitTo(*binnedtoy,RooFit::Minimizer("Minuit2","minimize"),RooFit::Minos(0),RooFit::Hesse(0),RooFit::PrintLevel(-1),RooFit::Strategy(0)); 

      RooPlot *plot_t = mass->frame();
      binnedtoy->plotOn(plot_t);
      pdf->fitTo(*binnedtoy,RooFit::Save(1),RooFit::Strategy(0)
                                          ,RooFit::Minimizer("Minuit2","minimize"),RooFit::Minos(0),RooFit::Hesse(0),RooFit::PrintLevel(-1),RooFit::Optimize(1));
      pdf->plotOn(plot_t); // Let RooFit handle natural normalization
      double chi2_t = plot_t->chiSquare(np);
      if( chi2_t>=chi2) npass++;
      toy_chi2.push_back(chi2_t*(nBinsForFit-np));
      delete plot_t;
    }
    std::cout << "[INFO] complete" << std::endl;
    prob = (double)npass / ntoys;

    TCanvas *can = new TCanvas();
    double medianChi2 = toy_chi2[(int)(((float)ntoys)/2)];
    double rms = TMath::Sqrt(medianChi2);

    TH1F toyhist(Form("gofTest_%s.pdf",pdf->GetName()),";Chi2;",50,medianChi2-5*rms,medianChi2+5*rms);
    for (std::vector<double>::iterator itx = toy_chi2.begin();itx!=toy_chi2.end();itx++){
      toyhist.Fill((*itx));
    }
    toyhist.Draw();

    TArrow lData(chi2*(nBinsForFit-np),toyhist.GetMaximum(),chi2*(nBinsForFit-np),0);
    lData.SetLineWidth(2);
    lData.Draw();
    can->SaveAs(name.c_str());

    // back to best fit   
    params->assignValueOnly(preParams);
  } else {
    prob = TMath::Prob(chi2*(nBinsForFit-np),nBinsForFit-np);
  }
  std::cout << "[INFO] GOF Chi2 in Observed =  " << chi2*(nBinsForFit-np) << std::endl;
  std::cout << "[INFO] GOF Chi2 in Observed =  " << chi2 << std::endl;
  std::cout << "[INFO] GOF p-value  =  " << prob << std::endl;
  delete pdf;
  return prob;

}

void plotComponents(RooRealVar *mass, RooAbsPdf *pdf, RooAbsPdf *zModel, RooAbsData *data, string name,vector<string> flashggCats_, int status, double *prob, double chi2FromFit=-1, bool fitWasBlinded=false){
/* Plot pdf vs data with separate Z component and ratio plot */
    
  // Calculate goodness of fit
  *prob = getGoodnessOfFit(mass,pdf,data,name);
  
  RooPlot *plot = mass->frame();
  
  // Blinding: escludi solo la SR [100,135] GeV
  mass->setRange("blind_low", mass->getMin(), 100);   // Sideband basso
  mass->setRange("blind_high", 135, mass->getMax()); // Sideband alto

  // Fit: usa solo sidebands
  mass->setRange("fit_low", mass->getMin(), 100);     // Low sideband: min-100 GeV
  mass->setRange("fit_high", 135, mass->getMax());   // High sideband: 135-max GeV

  // Plot data con blinding SR se richiesto
  if (fitWasBlinded) {
    data->plotOn(plot, MarkerStyle(20), MarkerSize(0.8), CutRange("blind_low,blind_high"));
  } else {
    data->plotOn(plot, MarkerStyle(20), MarkerSize(0.8));
  }

  // Normalizzazione: conta solo eventi nei sidebands
  double Nvis = -1;
  RooArgSet *final_params = pdf->getParameters((const RooArgSet*)(0));
  RooRealVar* bkg_norm = nullptr;
  TIterator* iter = final_params->createIterator();
  RooAbsArg* arg;
  while ((arg = (RooAbsArg*)iter->Next())) {
    if (TString(arg->GetName()).Contains("bkg_norm")) {
      bkg_norm = (RooRealVar*)arg;
    }
  }
  delete iter;
  if (bkg_norm) {
    // Calcola l'integrale della PDF di background nel range di fit e nel range di plot
    double fitIntegral = pdf->createIntegral(*mass, NormSet(*mass), Range("blind_low,blind_high"))->getVal();
    double normVal = bkg_norm->getVal();
    Nvis = normVal * (fitIntegral );
    std::cout << "[INFO] Using scaled PDF integral for Nvis: " << Nvis << " events (plot range)" << std::endl;
  } else {
    if (fitWasBlinded) {
      Nvis = data->sumEntries("1", "blind_low,blind_high");
      std::cout << "[INFO] Blinded fit: Nvis = " << Nvis << " events in sidebands" << std::endl;
    } else {
      Nvis = data->sumEntries();
      std::cout << "[INFO] Unblinded fit: Nvis = " << Nvis << " events total" << std::endl;
    }
  }

  // Plot PDF su tutto il range con normalizzazione naturale
  std::cout << "[INFO] Plotting function over full mass range with natural normalization" << std::endl;
  pdf->plotOn(plot,LineColor(kRed),LineWidth(2),Name("total"),Range(mass->getMin(),mass->getMax()),NormRange("blind_low,blind_high"),
                     Normalization(Nvis, RooAbsReal::NumEvent));

  // Componenti: stessa cosa
  if (pdf->InheritsFrom("RooAddPdf")) {
    RooAddPdf* addPdf = (RooAddPdf*)pdf;
    RooArgList comps = addPdf->pdfList();
    if (comps.getSize() >= 2) {
      RooAbsPdf* bkgPdf = (RooAbsPdf*)comps.at(0);
      if (bkgPdf) {
        // Usa sempre la normalizzazione del fit (bkg_norm->getVal()) e NormRange coerente
        pdf->plotOn(plot,Components(*bkgPdf),LineColor(kBlue),LineStyle(kDashed),LineWidth(2),Name("background"),
                   NormRange("blind_low,blind_high"),
                   Normalization(Nvis, RooAbsReal::NumEvent),Range(mass->getMin(),mass->getMax()));
      }
      
      // Plot Z component with proper normalization
      if (zModel) {
        // Anche la Z viene normalizzata con Nvis e NormRange coerente
        pdf->plotOn(plot,Components(*zModel),LineColor(kGreen+2),LineWidth(2),Name("Z"),
                   NormRange("blind_low,blind_high"),
                   Normalization(Nvis, RooAbsReal::NumEvent));
      }
    }
  }

  // Use chi2 from fit if provided
  double chi2;
  if (chi2FromFit > 0) {
    chi2 = chi2FromFit;
    cout << "[INFO] Using chi2 from fit: " << chi2 << endl;
  } else {
    int np = pdf->getParameters(*data)->getSize();
    chi2 = plot->chiSquare(np);
    cout << "[WARNING] Using chi2 from plot (may be unreliable): " << chi2 << endl;
  }

  TCanvas *canv = new TCanvas("canv","canv",800,800);
  canv->Divide(1,2);
  
  // Upper plot: data and fits
  canv->cd(1);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.02);
  gPad->SetPad(0.01,0.3,0.99,0.99);
  plot->GetXaxis()->SetTitleSize(0.04);
  plot->GetXaxis()->SetTitle("");
  plot->GetXaxis()->SetLabelSize(0);
  plot->GetYaxis()->SetTitleSize(0.05);
  plot->GetYaxis()->SetLabelSize(0.04);
  plot->GetYaxis()->SetTitleOffset(1.2);
  plot->GetYaxis()->SetTitle("Events");
  
  plot->SetTitle("");
  plot->SetMaximum(plot->GetMaximum()*1.3);
  plot->Draw();
  
  // Add legend
  TLegend *leg = new TLegend(0.65,0.65,0.89,0.89);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->AddEntry("Graph","Data","pe");
  leg->AddEntry("total","Total fit","l");
  leg->AddEntry("background","Background","l");
  leg->AddEntry("Z","Z resonance","l");
  leg->Draw();
  
  TLatex *lat = new TLatex();
  lat->SetNDC();
  lat->SetTextFont(42);
  lat->SetTextSize(0.034);
  lat->DrawLatex(0.15,0.94,Form("#chi^{2}/ndof = %.3f, Prob = %.2f",chi2,*prob));

  // Lower plot: ratio (Data-Background)/Background to highlight Z contribution
  canv->cd(2);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.15);
  gPad->SetPad(0.01,0.01,0.99,0.3);
  gPad->SetGridy();
  
  RooPlot *ratioPlot = mass->frame();
  
  // Create pull plot with conditional blinding
  if (fitWasBlinded) {
    data->plotOn(ratioPlot,MarkerStyle(20),MarkerSize(0.8),CutRange("blind_low,blind_high"));
    pdf->plotOn(ratioPlot,LineColor(kRed),LineWidth(2),NormRange("blind_low,blind_high"),
                     Normalization(Nvis, RooAbsReal::NumEvent)); // Let RooFit handle natural normalization
  } else {
    data->plotOn(ratioPlot,MarkerStyle(20),MarkerSize(0.8));
    pdf->plotOn(ratioPlot,LineColor(kRed),LineWidth(2)); // Let RooFit handle natural normalization
  }
  
  // Calculate pulls excluding signal region
  RooHist* hpull = ratioPlot->pullHist("", "", true); // Calculate pulls
  
  RooPlot *pullPlotFinal = mass->frame();
  pullPlotFinal->addPlotable(hpull,"P");
  
  pullPlotFinal->GetYaxis()->SetNdivisions(504);
  pullPlotFinal->GetYaxis()->SetLabelSize(0.12);
  pullPlotFinal->GetYaxis()->SetTitleSize(0.12);
  pullPlotFinal->GetYaxis()->SetTitleOffset(0.5);
  pullPlotFinal->GetYaxis()->SetTitle("Pull");
  pullPlotFinal->GetXaxis()->SetTitleSize(0.12);
  pullPlotFinal->GetXaxis()->SetLabelSize(0.12);
  pullPlotFinal->GetXaxis()->SetTitle("m_{ll} (GeV)");
  pullPlotFinal->SetTitle("");
  pullPlotFinal->SetMaximum(3);
  pullPlotFinal->SetMinimum(-3);
  pullPlotFinal->Draw();

  canv->SaveAs(Form("%s_components.pdf",name.c_str()));
  canv->SaveAs(Form("%s_components.png",name.c_str()));

  delete canv;
  delete lat;
  delete leg;
}

void plot(RooRealVar *mass, RooAbsPdf *pdf, RooAbsData *data, string name,vector<string> flashggCats_, int status, double *prob, double chi2FromFit=-1, bool fitWasBlinded=false){
/* Plot single pdf vs data, with chi2 from fit */
    
  // Calculate goodness of fit
  *prob = getGoodnessOfFit(mass,pdf,data,name);
  
  // Use the data directly without conversion to preserve workspace binning
  RooPlot *plot = mass->frame();
  
  // Blinding: escludi solo la SR [100,135] GeV
  mass->setRange("blind_low", mass->getMin(), 320);   // Sideband basso
  mass->setRange("blind_high", 980, mass->getMax()); // Sideband alto
  
  // Fit: usa solo sidebands
  mass->setRange("fit_low", mass->getMin(), 320);     // Low sideband: mass->getMin()-100 GeV
  mass->setRange("fit_high", 980, mass->getMax());   // High sideband: 135-mass->getMax() GeV
  
  // Plot data con blinding SR se richiesto
  if (fitWasBlinded) {
    data->plotOn(plot, MarkerStyle(20), MarkerSize(0.8), CutRange("blind_low,blind_high"));
  } else {
    data->plotOn(plot, MarkerStyle(20), MarkerSize(0.8));
  }

  // Normalizzazione: conta solo eventi nei sidebands
  double Nvis = -1;
  RooArgSet *final_params = pdf->getParameters((const RooArgSet*)(0));
  RooRealVar* bkg_norm = nullptr;
  TIterator* iter = final_params->createIterator();
  RooAbsArg* arg;
  while ((arg = (RooAbsArg*)iter->Next())) {
    if (TString(arg->GetName()).Contains("bkg_norm")) {
      bkg_norm = (RooRealVar*)arg;
    }
  }
  delete iter;
  if (bkg_norm) {
    double fitIntegral = pdf->createIntegral(*mass, NormSet(*mass), Range("blind_low,blind_high"))->getVal();
    double normVal = bkg_norm->getVal();
    Nvis = normVal * (fitIntegral );
    std::cout << "[INFO] Using bkg_norm for Nvis: " << Nvis << " events" << std::endl;
  } else {
    if (fitWasBlinded) {
      Nvis = data->sumEntries("1", "blind_low,blind_high");
      std::cout << "[INFO] Blinded fit: Nvis = " << Nvis << " events in sidebands" << std::endl;
    } else {
      Nvis = data->sumEntries();
      std::cout << "[INFO] Unblinded fit: Nvis = " << Nvis << " events total" << std::endl;
    }
  }

  // Plot PDF su tutto il range con normalizzazione naturale
  std::cout << "[INFO] Plotting function over full mass range with natural normalization" << std::endl;
  pdf->plotOn(plot,LineColor(kCyan),LineWidth(2),Name("total"),Range("blind_low,blind_high"),
                     Normalization(Nvis, RooAbsReal::NumEvent));

 
  
  // Use chi2 from fit if provided, otherwise calculate from plot
  double chi2;
  if (chi2FromFit > 0) {
    chi2 = chi2FromFit;
    cout << "[INFO] Using chi2 from fit: " << chi2 << endl;
  } else {
    // Fallback: calcolare dal grafico (può essere inaffidabile per RooDataHist)
    int np = pdf->getParameters(*data)->getSize();
    chi2 = plot->chiSquare(np);
    cout << "[WARNING] Using chi2 from plot (may be unreliable): " << chi2 << endl;
  }

  TCanvas *canv = new TCanvas();
  // canv->Divide(1,2); // Disabilita il grafico di pull per ora
  // canv->cd(1);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.10);
  // gPad->SetPad(0.01,0.2,0.99,0.99);
  plot->GetXaxis()->SetTitleSize(0.04);
  plot->GetXaxis()->SetTitle("m_{\ell#pi} (GeV)");
  plot->GetYaxis()->SetTitleSize(0.04);;
  //plot->GetYaxis()->SetTitle("Entries");
  plot->GetYaxis()->SetTitleOffset(1.35);

  // PDF già aggiunta sopra, aggiungi solo i parametri
  pdf->paramOn(plot,RooFit::Layout(0.34,0.85,0.89),RooFit::Format("NEA",AutoPrecision(1)));
  plot->getAttText()->SetTextSize(0.025);
  plot->SetMaximum(plot->GetMaximum()*1.4);
  plot->SetTitle("");
  plot->Draw();
  TLatex *lat = new TLatex();
  lat->SetNDC();
  lat->SetTextFont(42);
  lat->SetTextSize(0.034);
  lat->DrawLatex(0.15,0.94,Form("#chi^{2}/ndof = %.3f, Prob = %.2f, Fit Status = %d ",chi2,*prob,status));

  TLatex *lab = new TLatex();
  lab->SetNDC();
  lab->SetTextFont(42);
  lab->SetTextSize(0.034);
  std::string delimiter = "/";
  size_t pos = 0;
  std::string s = name;
  std::string token;
  while ((pos = s.find(delimiter)) != std::string::npos) {
    token = s.substr(0, pos);
    //std::cout << token << std::endl;
    s.erase(0, pos + delimiter.length());
  }
  //delimiter = ".";
  //token = s.substr(0, s.find(delimiter));
  //std::cout << s << std::endl;
  lab->DrawLatex(0.55,0.25,s.c_str());

  /* TODO: Fix pull calculation for RooDataHist
  // pulls
  canv->cd(2);
  gPad->SetLeftMargin(0.15);
  gPad->SetPad(0.01,0.01,0.99,0.2);
  gPad->SetGridy();
  
  // Calculate pulls explicitly using named objects
  RooHist* hpull = plot->pullHist("data", "pdf");
  
  RooPlot *plot2 = mass->frame();
  plot2->GetYaxis()->SetNdivisions(504);
  plot2->GetYaxis()->SetLabelSize(0.17);
  plot2->GetYaxis()->SetTitleSize(0.17);
  plot2->GetYaxis()->SetTitleOffset(0.24);
  //plot2->GetYaxis()->SetRangeUser(-3.0,3.0);
  plot2->SetMinimum(-3.);
  plot2->SetMaximum(3.);
  plot2->GetYaxis()->SetTitle("Pulls") ;
  plot2->GetXaxis()->SetTitle("");
  plot2->GetXaxis()->SetLabelOffset(5);
  plot2->addPlotable(hpull,"P"); 
  plot2->Draw();
  */

  canv->SaveAs(Form("%s.pdf",name.c_str()));
  canv->SaveAs(Form("%s.png",name.c_str()));

  delete canv;
  delete lat;
}
void plot(RooRealVar *mass, RooMultiPdf *pdfs, RooCategory *catIndex, RooAbsData *data, string name, vector<string> flashggCats_, int cat, int bestFitPdf=-1, bool fitWasBlinded=false){
/* Plot MultiPdf vs data */
  
  int color[7] = {kCyan+1,kGreen+2,kAzure+2,kTeal+1,kSpring+2,kCyan+3,kGreen+1};
  TLegend *leg = new TLegend(0.6,0.65,0.95,0.90);
  leg->SetFillColor(0);
  leg->SetLineColor(1);
  RooPlot *plot = mass->frame();

  // Blinding: escludi solo la SR [100,135] GeV
  mass->setRange("blind_low", mass->getMin(), 800);   // Sideband basso
  mass->setRange("blind_high", 850, mass->getMax()); // Sideband alto
  
  // Fit: usa solo sidebands
  mass->setRange("fit_low", mass->getMin(), 800);     // Low sideband: mass->getMin()-100 GeV
  mass->setRange("fit_high", 850, mass->getMax());   // High sideband: 135-mass->getMax() GeV
  
  data->plotOn(plot,MarkerStyle(20),MarkerSize(0.8),CutRange("blind_low,blind_high")); // Plot escludendo la regione di segnale
  
  TCanvas *canv = new TCanvas();
  //TPad *pad1 = new TPad("pad1","pad1",0,0,1,1);
  //pad1->SetBottomMargin(0.18);
  //pad1->Draw();
  //pad1->cd();

  int currentIndex = catIndex->getIndex();
  TObject *datLeg = plot->getObject(int(plot->numItems()-1));
  leg->AddEntry(datLeg,"Data");
  int style=1;
  RooAbsPdf *pdf;
  RooCurve *nomBkgCurve;
  double Nvis;
  if (fitWasBlinded) {
    Nvis = data->sumEntries("1", "blind_low,blind_high");
    std::cout << "[INFO] Blinded fit: Nvis = " << Nvis << " events in sidebands" << std::endl;
  } else {
    Nvis = data->sumEntries();
    std::cout << "[INFO] Unblinded fit: Nvis = " << Nvis << " events total" << std::endl;
  }
  int bestcol= -1;
  for (int icat=0;icat<catIndex->numTypes();icat++){
    int col;
    if (icat<=6) col=color[icat];
    else {col=kBlack; style++;}
    catIndex->setIndex(icat);
    pdfs->getCurrentPdf()->fitTo(*data,RooFit::Minos(0),RooFit::Minimizer("Minuit2","minimize"),RooFit::Strategy(0),RooFit::PrintLevel(-1),RooFit::Optimize(1));  
    if (fitWasBlinded) pdfs->getCurrentPdf()->plotOn(plot,LineColor(kRed),LineWidth(2),NormRange("blind_low,blind_high"), Range(mass->getMin(),mass->getMax()),
                     Normalization(Nvis, RooAbsReal::NumEvent)); // fix blinded to full range norm
    else pdfs->getCurrentPdf()->plotOn(plot,LineColor(col),LineStyle(style)); // Let RooFit handle natural normalization
  //  pdfs->getCurrentPdf()->plotOn(plot,LineColor(col),LineStyle(style),RooFit::Range("full"),RooFit::NormRange("full"));
    TObject *pdfLeg = plot->getObject(int(plot->numItems()-1));
    std::string ext = "";
    if (bestFitPdf==icat) {
      ext=" (Best Fit Pdf) ";
      pdf= pdfs->getCurrentPdf();
      nomBkgCurve = (RooCurve*)plot->getObject(plot->numItems()-1);
      bestcol = col;
    }
    string pdfName = pdfs->getCurrentPdf()->GetName();
    std:: cout << "PDF Name: " << pdfName << std::endl; 
    size_t start = pdfName.find("_") + 1;
    size_t end = pdfName.find("_", start);
    std::string family = pdfName.substr(start, end - start); // Risultato: "Bernstein"
    leg->AddEntry(pdfLeg,Form("%s%s",family.c_str(),ext.c_str()),"L");
  }
  plot->SetTitle(Form("Category 2"));
  plot->SetMaximum(plot->GetMaximum()*1.4);
  plot->GetXaxis()->SetTitle("m_{l\\pi} (GeV)");
  plot->Draw();
  leg->Draw("same");
  CMS_lumi( canv, 0, 0);
  canv->SaveAs(Form("%s.pdf",name.c_str()));
  canv->SaveAs(Form("%s.png",name.c_str()));
  catIndex->setIndex(currentIndex);
  delete canv;
}

void plot(RooRealVar *mass, map<string,RooAbsPdf*> pdfs, RooAbsData *data, string name, vector<string> flashggCats_, int cat, int bestFitPdf=-1, bool fitWasBlinded=true){
/* Plot several Pdfs vs data, without ratio plot, (used for the "truth") */
  
  int color[7] = {kCyan+1,kGreen+2,kAzure+2,kTeal+1,kSpring+2,kCyan+3,kGreen+1};
  TCanvas *canv = new TCanvas();
  TLegend *leg = new TLegend(0.6,0.65,0.88,0.88);
  leg->SetFillColor(0);
  leg->SetLineColor(0);
  RooPlot *plot = mass->frame();

  // Blinding: escludi solo la SR [100,135] GeV
  mass->setRange("blind_low", mass->getMin(), 800);   // Sideband basso
  mass->setRange("blind_high", 850, mass->getMax()); // Sideband alto
  
  // Fit: usa solo sidebands
  mass->setRange("fit_low", mass->getMin(), 800);     // Low sideband: mass->getMin()-100 GeV
  mass->setRange("fit_high", 850, mass->getMax());   // High sideband: 135-mass->getMax() GeV
  
  data->plotOn(plot,MarkerStyle(20),MarkerSize(0.8), CutRange("blind_low,blind_high")); // Plot escludendo la regione di segnale

  TObject *datLeg = plot->getObject(int(plot->numItems()-1));
  if(flashggCats_.size() >0){
  leg->AddEntry(datLeg,Form("Data - %s",flashggCats_[cat].c_str()),"LEP");
  } else {
  leg->AddEntry(datLeg,Form("Data - %d",cat),"LEP");
  }
  int i=0;
  int style=1;
  for (map<string,RooAbsPdf*>::iterator it=pdfs.begin(); it!=pdfs.end(); it++){
    int col;
    if (i<=6) col=color[i];
    else {col=kBlack; style++;}
    double Nvis = -1;
    RooArgSet *final_params = it->second->getParameters((const RooArgSet*)(0));
    RooRealVar* bkg_norm = nullptr;
    TIterator* iter = final_params->createIterator();
    RooAbsArg* arg;
    while ((arg = (RooAbsArg*)iter->Next())) {
      if (TString(arg->GetName()).Contains("bkg_norm")) {
        bkg_norm = (RooRealVar*)arg;
      }
    }
    delete iter;
    bkg_norm = 0;
    if (bkg_norm) {
      Nvis = bkg_norm->getVal();
      std::cout << "[INFO] Using bkg_norm for Nvis: " << Nvis << " events" << std::endl;
    } else {
      if (fitWasBlinded) {
        Nvis = data->sumEntries("1", "blind_low,blind_high");
        std::cout << "[INFO] Blinded fit: Nvis = " << Nvis << " events in sidebands" << std::endl;
      } else {
        Nvis = data->sumEntries();
        std::cout << "[INFO] Unblinded fit: Nvis = " << Nvis << " events total" << std::endl;
      }
    }
    std :: cout << " "  << std::endl;
    std :: cout << " "  << std::endl;
    std :: cout << " "  << std::endl;
    std :: cout << " "  << std::endl;
    std :: cout << "Normalizing events " << Nvis << std::endl;
    std :: cout << " " << std::endl;
    std :: cout <<  " "<< std::endl;
    std :: cout << " " << std::endl;
    //if (fitWasBlinded)
    it->second->plotOn(plot,LineColor(col),LineWidth(2),Range(mass->getMin(),mass->getMax()),
                     Normalization(Nvis, RooAbsReal::NumEvent)); // fix blinded to full range norm
    //else it->second->plotOn(plot,LineColor(col),LineStyle(style)); // Let RooFit handle natural normalization
     // Let RooFit handle natural normalization
    TObject *pdfLeg = plot->getObject(int(plot->numItems()-1));
    std::string ext = "";
    if (bestFitPdf==i) ext=" (Best Fit Pdf) ";
    leg->AddEntry(pdfLeg,Form("%s%s",it->first.c_str(),ext.c_str()),"L");
    i++;
  }
  plot->SetMaximum(plot->GetMaximum()*1.4);
  if (cat >= 0 && cat < (int)flashggCats_.size()) {
    plot->SetTitle(Form(" %s",flashggCats_[cat].c_str()));
  } else {
    plot->SetTitle(Form("cat%d",cat));
  }
  plot->Draw();
  leg->Draw("same");
  CMS_lumi( canv, 0, 0);
  canv->SaveAs(Form("%s.pdf",name.c_str()));
  canv->SaveAs(Form("%s.png",name.c_str()));
  delete canv;
}

void transferMacros(TFile *inFile, TFile *outFile){
  
  TIter next(inFile->GetListOfKeys());
  TKey *key;
  while ((key = (TKey*)next())){
    if (string(key->ReadObj()->ClassName())=="TMacro") {
      //cout << key->ReadObj()->ClassName() << " : " << key->GetName() << endl;
      TMacro *macro = (TMacro*)inFile->Get(key->GetName());
      outFile->cd();
      macro->Write();
    }
  }
}

int getBestFitFunction(RooMultiPdf *bkg, RooAbsData *data, RooCategory *cat, bool silent=false){
/* Get index of the best fit pdf (minimum NLL, including correction) among functions in the multipdf.
   All fits are performed again. */

  double global_minNll = 1E10;
  int best_index = 0;
  int number_of_indeces = cat->numTypes();
    
  RooArgSet snap,clean;
  RooArgSet *params = bkg->getParameters((const RooArgSet*)0);
  params->remove(*cat);  // pdf_index is RooCategory, removed from parameters
  params->snapshot(snap);
  params->snapshot(clean);
  if (!silent) {
    //params->Print("V");
  }
 
  // Uncomment to try to make converge a failed fit
  //bkg->setDirtyInhibit(1);
  //RooAbsReal *nllm = bkg->createNLL(*data);
  //RooMinimizer minim(*nllm);
  //minim.setStrategy(1);
  
  for (int id=0;id<number_of_indeces;id++){    
    params->assignValueOnly(clean);
    cat->setIndex(id);

    //RooAbsReal *nllm = bkg->getCurrentPdf()->createNLL(*data);

    if (!silent) {
      /*
      std::cout << "BEFORE  MAKING FIT" << std::endl;
      params->Print("V");
      std::cout << "-----------------------" << std::endl;    
      */
    }
    
    //minim.minimize("Minuit2","minimize");
    double minNll=0; //(nllm->getVal())+bkg->getCorrection();
    int fitStatus=1;    
    runFit(bkg->getCurrentPdf(),data,&minNll,&fitStatus,/*max iterations*/7);
    // Add the penalty

    minNll=minNll+bkg->getCorrection();

    if (!silent) {
      /*
      std::cout << "After Minimization ------------------  " <<std::endl;
      std::cout << bkg->getCurrentPdf()->GetName() << " " << minNll <<std::endl;
      bkg->Print("v");
      bkg->getCurrentPdf()->getParameters(*data)->Print("V");
      std::cout << " ------------------------------------  " << std::endl;
  
      params->Print("V");
      */
      std::cout << "[INFO] AFTER FITTING" << std::endl;
      std::cout << "[INFO] Function was " << bkg->getCurrentPdf()->GetName() <<std::endl;
      std::cout << "[INFO] Correction Applied is " << bkg->getCorrection() <<std::endl;
      std::cout << "[INFO] NLL + c = " <<  minNll << std::endl;
      std::cout << "-----------------------" << std::endl;
    }
      
    if (minNll < global_minNll){
      global_minNll = minNll;
      snap.assignValueOnly(*params);
      best_index=id;
    }
  } // end loop over pdf_index 
  cat->setIndex(best_index);
  params->assignValueOnly(snap);

  std::cout << "[INFO] Best fit Function -- " << bkg->getCurrentPdf()->GetName() << " " << cat->getIndex() <<std::endl;
  std::cout << "[INFO] Best fit parameters " << std::endl;
  params->Print("V");
  
  return best_index;
}

// Function to create background PDF with optional turn-on function and Z resonance
RooAbsPdf* createBackgroundWithTurnOn(PdfModelBuilder &pdfsModel, string type, int order, RooRealVar* mass, const char* ext, bool includeTurnOn = true, bool includeZ = false, RooWorkspace* dataWS = nullptr, std::string turnOnType = "Fermi") {
  
  // Get the base background PDF from the PdfModelBuilder using the existing getPdf function
  RooAbsPdf *basePdf = getPdf(pdfsModel, type, order, ext);
  
  if (!basePdf) {
    cerr << "[ERROR] Could not create base PDF of type " << type << " order " << order << endl;
    return nullptr;
  }
  
  // Fix normalization for RooAddPdf if it's one
  if (basePdf->InheritsFrom("RooAddPdf")) {
    RooArgSet normSet(*mass);
    ((RooAddPdf*)basePdf)->fixCoefNormalization(normSet);
  }
  
  RooAbsPdf *currentPdf = basePdf;
  
  // Add turn-on function if requested
  if (includeTurnOn) {  // RE-ENABLE turn-on
    // Create turn-on function (Fermi-Dirac) - constrained parameters for stability
    // Turn-on function: 1/(1 + exp((cutoff - x)/beta))
    // Tight ranges around physically motivated values
    RooRealVar *cutoff = new RooRealVar(Form("%s_turnon_cutoff", ext), "Turn-on cutoff", 50., 43,55 );
    RooRealVar *beta = new RooRealVar(Form("%s_turnon_beta", ext), "Turn-on beta", 6.3, 3, 30);
    
    cout << "[INFO] Adding turn-on function for " << type << " order " << order << endl;
    
    // Alternative turn-on functions to try:
    RooGenericPdf *turnOnPdf = nullptr;
    
    // Option 1: Fermi-Dirac (current default)
    if (turnOnType == "Fermi") {
      turnOnPdf = new RooGenericPdf(Form("%s_turnon", ext), "Fermi-Dirac turn-on", 
                                   "1.0/(1.0 + TMath::Exp((@1-@0)/@2))", 
                                   RooArgList(*mass, *cutoff, *beta));
    }
    // Option 2: Error Function turn-on (smoother)
    else if (turnOnType == "Erf" || type == "PowerLaw" || type == "PowerLawSingle") {
      turnOnPdf = new RooGenericPdf(Form("%s_turnon", ext), "Error function turn-on", 
                                   "0.5*(1.0 + TMath::Erf((@0-@1)/@2))", 
                                   RooArgList(*mass, *cutoff, *beta));
    }
    // Option 3: Double exponential (sharp cut)
    else if (turnOnType == "DExp" || type == "Exponential") {
      turnOnPdf = new RooGenericPdf(Form("%s_turnon", ext), "Double exponential turn-on", 
                                   "1.0 - TMath::Exp(-TMath::Exp((@0-@1)/@2))", 
                                   RooArgList(*mass, *cutoff, *beta));
    }
    // Default: Fermi-Dirac
    else {
      turnOnPdf = new RooGenericPdf(Form("%s_turnon", ext), "Turn-on function", 
                                   "1.0/(1.0 + TMath::Exp((@1-@0)/@2))", 
                                   RooArgList(*mass, *cutoff, *beta));
    }
    
    currentPdf = new RooProdPdf(Form("_%s_%s_with_turnon",type.c_str(), ext), 
                               "Background with turn-on", 
                               RooArgList(*currentPdf, *turnOnPdf));
  }
  
  // Add Z resonance if requested
  if (includeZ && dataWS) {  // RE-ENABLE Z with FIXED approach
    RooAbsPdf *zModel = dataWS->pdf("model_Z_c2");
    if (zModel) {
      cout << "[INFO] Adding Z resonance model for " << type << " order " << order << endl;
      
      // Freeze Z shape parameters but release some key ones for fitting
      RooArgSet *zParams = zModel->getParameters(RooArgSet(*mass));
      TIterator *iter = zParams->createIterator();
      RooRealVar *param;
      while ((param = (RooRealVar*)iter->Next())) {
        // Release Z mass and width for better fitting
          param->setConstant(kTRUE);
          cout << "[INFO] Freezing Z parameter: " << param->GetName() << " = " << param->getVal() << " (fixed)" << endl;
        }
      
      delete iter;
      delete zParams;
      
      // Debug: Check Z PDF normalization
      RooArgSet normSet(*mass);
      double zIntegral = zModel->createIntegral(normSet)->getVal();
      cout << "[DEBUG] Z PDF integral over mass range: " << zIntegral << endl;
      
      // Set up Z normalization - INDEPENDENT NORMALIZATION approach
      
      cout << "[INFO] Adding Z resonance using independent normalizations" << endl;
      
      // Create independent normalizations for background and Z
      RooRealVar *bkg_norm = new RooRealVar(Form("%s_bkg_norm", ext), "Background normalization", 
                                           1400000, 0., 2000000.); // ~data events for background
      
      RooRealVar *z_norm = new RooRealVar(Form("%s_z_norm", ext), "Z normalization", 
                                         16000, 14000, 18000.); // Start with ~3.5% of background, wider range
      
      cout << "[INFO] Using independent normalizations:" << endl;
      cout << "  Background norm: " << bkg_norm->getVal() << " events" << endl;
      cout << "  Z norm: " << z_norm->getVal() << " events" << endl;
      
      // Combine with independent normalizations (no coefficients!)
      currentPdf = new RooAddPdf(Form("%s_with_z", ext),
                                "Background with Z",
                                RooArgList(*currentPdf, *zModel),
                                RooArgList(*bkg_norm, *z_norm));
      
      cout << "[INFO] Z added with independent normalizations" << endl;
    } else {
      cout << "[WARNING] Could not find Z model 'model_Z_c2' in workspace" << endl;
    }
  }
  
  return currentPdf;
}

void runIterativeFits(RooRealVar* mass, RooAbsData* data, PdfModelBuilder& pdfsModel, RooWorkspace* dataWS, const std::string& ext, int order, const std::vector<std::string>& flashggCats_, int cat, const std::string& outDir, bool iterativeMode, const std::string& turnOnType) {
    // Primo fit: [mass->getMin(),78] e >140 senza Z
    mass->setRange("fit_low", mass->getMin(), 800);
    mass->setRange("fit_high", 850, mass->getMax());
    mass->setRange("blind_low", mass->getMin(), 800);
    mass->setRange("blind_high", 850, mass->getMax());
    bool includeZ = false;
    std::cout << "[ITERATIVE FIT] Primo fit: " << ext << " [ " << mass->getMin()<< ",78] e >140 senza Z" << std::endl;
    RooAbsPdf* bkgPdf1 = createBackgroundWithTurnOn(pdfsModel, ext, order, mass, ext.c_str(), true, includeZ, dataWS, turnOnType);
    int fitStatus1 = 0;
    double chi2FromFit1 = -1;
    double thisNll1 = 0.;
    runFit(bkgPdf1, data, &thisNll1, &fitStatus1, 7, &chi2FromFit1, true);
    plot(mass, bkgPdf1, data, outDir+"/iterativeFit1_"+ext, flashggCats_, fitStatus1, &chi2FromFit1, chi2FromFit1,true);
    std::cout << "[ITERATIVE FIT] Primo fit status: " << fitStatus1 << std::endl;

    // Secondo fit: [mass->getMin(),100] e >140 con Z
    mass->setRange("fit_low", mass->getMin(), 800);
    mass->setRange("fit_high", 850, mass->getMax());
    mass->setRange("blind_low", mass->getMin(), 800);
    mass->setRange("blind_high", 850, mass->getMax());
    includeZ = true;
    std::cout << "[ITERATIVE FIT] Secondo fit: " << ext << " [mass->getMin(),100] e >140 con Z" << std::endl;
    RooAbsPdf* bkgPdf2 = createBackgroundWithTurnOn(pdfsModel, ext, order, mass, ext.c_str(), true, includeZ, dataWS, turnOnType);
    int fitStatus2 = 0;
    double chi2FromFit2 = -1;
    double thisNll2 = 0.;
    runFit(bkgPdf2, data, &thisNll2, &fitStatus2, 7, &chi2FromFit2, true);
    plot(mass, bkgPdf2, data, outDir+"/iterativeFit2_"+ext, flashggCats_, fitStatus2, &chi2FromFit2, chi2FromFit2,true);
    std::cout << "[ITERATIVE FIT] Secondo fit status: " << fitStatus2 << std::endl;
}

int main(int argc, char* argv[]){
 
  setTDRStyle();
  writeExtraText = true;       // if extra text
  extraText  = "Preliminary";  // default extra text is "Preliminary"
  lumi_8TeV  = "19.1 fb^{-1}"; // default is "19.7 fb^{-1}"
  lumi_7TeV  = "4.9 fb^{-1}";  // default is "5.1 fb^{-1}"
  lumi_sqrtS = "13.6 TeV";       // used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
  string year_ = "2024";
  //int year_ = 2017;

  string fileName;
  // Default to your specific workspace
  string workspaceFile = "/afs/cern.ch/work/e/elfontan/private/dijetAnalysis_ScoutingRun3/BKGModelling/FittingDijetMass/ws/dijetWS_res.root";
  //string workspaceFile = "/pnfs/psi.ch/cms/trivcat/store/user/gcelotto/ws_dt_fit/workspace_step2.root";
  //string workspaceFile = " /t3home/ratramon/ggHBB/Btob_studies/bb_analysis/scripts/workspace.root";
  int ncats;
  int singleCategory;
  int catOffset;
  string datfile;
  string outDir;
  string outfilename;
  bool is2011=false;
  bool verbose=true;
  bool saveMultiPdf=false;
  bool includeTurnOn=false;  // Disable turn-on function by default for speed
  bool includeZ=false;      // Option to include Z resonance
  bool blindSignalRegion=false; // Option to blind signal region [100-140] GeV in fit
  int rebinFactor=1;        // Rebin factor for histogram
  int isFlashgg_ =0;
  string flashggCatsStr_;
  vector<string> flashggCats_;
  bool isData_ =0;
  bool iterativeFit = false;
  std::string turnOnType = "Fermi"; // Default

  po::options_description desc("Allowed options");
  desc.add_options()
    ("help,h",                                                                                  "Show help")
    ("infilename,i", po::value<string>(&fileName),                                              "In file name")
    ("workspace,w", po::value<string>(&workspaceFile)->default_value("/afs/cern.ch/work/e/elfontan/private/dijetAnalysis_ScoutingRun3/BKGModelling/FittingDijetMass/ws/dijetWS_res.root"), "Workspace file with RooDataHist and Z model")
    ("includeTurnOn", po::value<bool>(&includeTurnOn)->default_value(false),                    "Include turn-on function for low mass threshold")
    ("includeZ", po::value<bool>(&includeZ)->default_value(false),                             "Include Z resonance model from workspace")
    ("blindSignalRegion", po::value<bool>(&blindSignalRegion)->default_value(false),           "Blind signal region [100-140] GeV in fit (not just plots)")
    ("ncats,c", po::value<int>(&ncats)->default_value(5),                                       "Number of categories")
    ("singleCat", po::value<int>(&singleCategory)->default_value(0),                           "Run A single Category (default: category 0)")
    ("datfile,d", po::value<string>(&datfile)->default_value("dat/fTest.dat"),                  "Right results to datfile for BiasStudy")
    ("outDir,D", po::value<string>(&outDir)->default_value("plots/fTest"),                      "Out directory for plots")
    ("saveMultiPdf", po::value<string>(&outfilename),                                           "Save a MultiPdf model with the appropriate pdfs")
    ("runFtestCheckWithToys",                                                                   "When running the F-test, use toys to calculate pvals (and make plots) ")
    ("is2011",                                                                                  "Run 2011 config")
    ("is2012",                                                                                  "Run 2012 config")
    ("isFlashgg",  po::value<int>(&isFlashgg_)->default_value(1),                               "Use Flashgg output ")
    ("isData",  po::value<bool>(&isData_)->default_value(0),                                    "Use Data not MC ")
    ("flashggCats,f", po::value<string>(&flashggCatsStr_)->default_value("UntaggedTag_0,UntaggedTag_1,UntaggedTag_2,UntaggedTag_3,UntaggedTag_4,VBFTag_0,VBFTag_1,VBFTag_2,TTHHadronicTag,TTHLeptonicTag,VHHadronicTag,VHTightTag,VHLooseTag,VHEtTag"),                  "Flashgg category names to consider")
    ("year", po::value<string>(&year_)->default_value("2016"),                                  "Dataset year")
    ("catOffset", po::value<int>(&catOffset)->default_value(0),                                 "Category numbering scheme offset")
    ("mN", po::value<float>(&mN)->default_value(2.75),                                          "Mass of the peak, for center of window")
    ("sigma", po::value<float>(&sigma)->default_value(0.025),                                   "Sigma of the peak, for size of window")
    ("nsigma", po::value<int>(&nsigma)->default_value(10),                                       "Sigma multiplier, for size of window")
    ("rebinFactor,r", po::value<int>()->default_value(1),                                       "Rebin factor to reduce number of bins (1=no rebinning)")
    ("verbose,v",                                                                               "Run with more output")
    ("iterativeFit", po::value<bool>(&iterativeFit)->default_value(false), "Enable iterative fitting: first fit [mass->getMin(),78]+>140 w/o Z, then [mass->getMin(),100]+>140 with Z")
    ("turnOnType", po::value<std::string>(&turnOnType)->default_value("Fermi"), "Type of turn-on function: Fermi, Erf, DExp")
  ;
  po::variables_map vm;
  po::store(po::parse_command_line(argc,argv,desc),vm);
  po::notify(vm);
  if (vm.count("help")) { cout << desc << endl; exit(1); }
  if (vm.count("is2011")) is2011=true;
  saveMultiPdf = vm.count("saveMultiPdf");

  if (vm.count("verbose")) verbose=true;
  if (vm.count("runFtestCheckWithToys")) runFtestCheckWithToys=true;
  if (vm.count("rebinFactor")) rebinFactor = vm["rebinFactor"].as<int>();

  std::cout << "DEBUG mN=" << mN << std::endl; 

  mN_low  = mN - nsigma * sigma;
  mN_high = mN + nsigma * sigma;

  if (!verbose) {
    RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
    RooMsgService::instance().setSilentMode(true);
    gErrorIgnoreLevel=kWarning;
  }
  
  // Speed up numerical integrations globally
 // RooAbsReal::defaultIntegratorConfig()->method1D().setLabel("RooAdaptiveGaussKronrodIntegrator1D");
 // RooAbsReal::defaultIntegratorConfig()->getConfigSection("RooAdaptiveGaussKronrodIntegrator1D").setRealValue("maxSeg", 100);
 // RooAbsReal::defaultIntegratorConfig()->getConfigSection("RooAdaptiveGaussKronrodIntegrator1D").setRealValue("epsAbs", 1e-6);
 // RooAbsReal::defaultIntegratorConfig()->getConfigSection("RooAdaptiveGaussKronrodIntegrator1D").setRealValue("epsRel", 1e-6);
  
  
  // DO NOT override ncats for flashgg - keep single category behavior
  // if (isFlashgg_==1){
  //   ncats= flashggCats_.size();
  // }

  if(verbose) std::cout << "[INFO] SaveMultiPdf? " << saveMultiPdf << std::endl;
  TFile *outputfile;
  RooWorkspace *outputws;

  if (saveMultiPdf){
  outputfile = new TFile(outfilename.c_str(),"RECREATE");
  outputws = new RooWorkspace(); outputws->SetName("multipdf");
  }

  system(Form("mkdir -p %s",outDir.c_str()));
  
  // Load the main data file (if provided)
  TFile *inFile = nullptr;
  RooWorkspace *inWS = nullptr;
  
  if (!fileName.empty()) {
    inFile = TFile::Open(fileName.c_str());
    std::cout<< "Loading file: " << fileName << std::endl;
    if(isFlashgg_){
      if (isData_){
        inWS = (RooWorkspace*)inFile->Get("wS");
      } else {
        inWS = (RooWorkspace*)inFile->Get("wS");
      }
    } else {
      inWS = (RooWorkspace*)inFile->Get("wS");
    }
  }
  
  // Load the workspace with RooDataHist and Z model
  TFile *wsFile = TFile::Open(workspaceFile.c_str());
  if (!wsFile || wsFile->IsZombie()) {
    cerr << "[ERROR] Cannot open workspace file: " << workspaceFile << endl;
    return 1;
  }
  
  RooWorkspace *dataWS = (RooWorkspace*)wsFile->Get("workspace_step2");
  //RooWorkspace *dataWS = (RooWorkspace*)wsFile->Get("w");
  if (!dataWS) {
    // Try alternative workspace names
    dataWS = (RooWorkspace*)wsFile->Get("dijet_ws");
    if (!dataWS) {
      dataWS = (RooWorkspace*)wsFile->Get("dijet_ws");
      if (!dataWS) {
        dataWS = (RooWorkspace*)wsFile->Get("dijet_ws");
        if (!dataWS) {
          cerr << "[ERROR] Cannot find workspace in " << workspaceFile << endl;
          wsFile->ls();
          return 1;
        }
      }
    }
  }
  
  std::cout << "[INFO] Loaded workspace from: " << workspaceFile << std::endl;
  if (verbose) {
    std::cout << "[INFO] Workspace contents:" << std::endl;
    dataWS->Print();
  }
  
  // Use dataWS for RooDataHist and inWS for other things (if available)
  if (!inWS) {
    inWS = dataWS;  // Use dataWS as fallback
  }
  if (verbose){
   std::cout << "[INFO]  inWS open " << inWS << std::endl;
   inWS->Print();
  }
  if (saveMultiPdf){
    transferMacros(inFile,outputfile);

    RooRealVar *intL; 
    RooRealVar *sqrts;

    if (isFlashgg_){
      //intL  = (RooRealVar*)inWS->var("IntLumi");
      intL  = intLumi_;
      sqrts = (RooRealVar*)inWS->var("SqrtS");
      if (!sqrts){ sqrts = new RooRealVar("SqrtS","SqrtS",13); }
    std::cout << "[INFO] got intL and sqrts " << intL << ", " << sqrts << std::endl;


    } else {
      //intL  = (RooRealVar*)inWS->var("IntLumi");
      intL  = intLumi_;
      sqrts = (RooRealVar*)inWS->var("Sqrts");
    }
    outputws->import(*intL);
    outputws->import(*sqrts);
    std::cout << "[INFO] got intL and sqrts " << intL << ", " << sqrts << std::endl;
  }

  // Set up which families of functions you want to test
  vector<string> functionClasses;
  functionClasses.push_back("Bernstein");
  functionClasses.push_back("Exponential");
  // functionClasses.push_back("ExponentialSum");  // Removed: doesn't fit well
  // functionClasses.push_back("Chebychev");       // Removed
  functionClasses.push_back("PowerLaw");         // Regular PowerLaw (not Generic)
 // functionClasses.push_back("PowerLawSingle");   // Alternative PowerLaw implementation
  // functionClasses.push_back("PowerLawGeneric"); // Replaced with regular PowerLaw
  // functionClasses.push_back("PowerLawSum");     // Removed: doesn't fit well
  // functionClasses.push_back("Laurent");         // Removed: similar to polynomials
  // functionClasses.push_back("KeysPdf");      // Needs special setup with setKeysPdfAttributes
//  functionClasses.push_back("Chebychev");
//  functionClasses.push_back("Polynomial");
  map<string,string> namingMap;
  namingMap.insert(pair<string,string>("Bernstein","pol"));
  namingMap.insert(pair<string,string>("Exponential","exp"));
  namingMap.insert(pair<string,string>("ExponentialSum","expsum"));
  namingMap.insert(pair<string,string>("Chebychev","che"));
  // namingMap.insert(pair<string,string>("PowerLawSum","pow"));
  // namingMap.insert(pair<string,string>("Laurent","lau"));
//  namingMap.insert(pair<string,string>("Chebychev","che"));
 // namingMap.insert(pair<string,string>("Polynomial","pol"));

  FILE *resFile ;
  //if  (singleCategory >-1) resFile = fopen(Form("%s/fTestResults_%s.txt",outDir.c_str(),flashggCats_[singleCategory].c_str()),"w");
  resFile = fopen(Form("%s/fTestResults.txt",outDir.c_str()),"w");
  vector<map<string,int> > choices_vec;
  vector<map<string,std::vector<int> > > choices_envelope_vec;
  vector<map<string,RooAbsPdf*> > pdfs_vec;

  PdfModelBuilder pdfsModel;
  RooRealVar *mass = (RooRealVar*)dataWS->var("CMS_dijet_mass");
  //RooRealVar *mass = (RooRealVar*)dataWS->var("dijet_ws"); 
  //RooRealVar *mass = (RooRealVar*)dataWS->var("mjj"); 
  std:: cout << "[INFO] Got mass from ws " << mass << std::endl;
  mass->Print("v");
  std:: cout << "[INFO] Mass range: [" << mass->getMin() << ", " << mass->getMax() << "]" << std::endl;
  
  // Define blinding ranges for fit and plot
  //mass->setRange("blind_low", mass->getMin(), 100);   // Sideband basso
  //mass->setRange("blind_high", 135, mass->getMax());   // Sideband alto
  mass->setRange("blind_low", mass->getMin(), 800);   // Sideband basso
  mass->setRange("blind_high", 850, mass->getMax());   // Sideband alto
  if (blindSignalRegion) {
    std::cout << "[INFO] FAKE Signal region [800-850] GeV will be blinded in fit" << std::endl;
  }
  
  pdfsModel.setObsVar(mass);
  double upperEnvThreshold = 0.1; // upper threshold on prob_ftest to include function in envelope (looser than truth function)
  double minGofThreshold = 0.01;  // minimal goodness of fit to include function in envelope

  fprintf(resFile,"Truth Model & d.o.f & $\\Delta NLL_{N+1}$ & $p(\\chi^{2}>\\chi^{2}_{(N\\rightarrow N+1)})$ \\\\\n");
  fprintf(resFile,"\\hline\n");

  std::string ext = is2011 ? "7TeV" : "8TeV";
        if( isFlashgg_ ){
          if( year_ == "all" ){ ext = "13TeV"; }
          else{ ext = Form("%s_13TeV",year_.c_str()); }
        }

  std::cout << "[INFO] Number of categories to process: " << ncats << std::endl;
  // Start processing categories
  for (int cat=0; cat<ncats; cat++){

    map<string,int> choices;
    map<string,std::vector<int> > choices_envelope;
    map<string,RooAbsPdf*> pdfs;
    map<string,RooAbsPdf*> allPdfs;
    string catname;
    if (isFlashgg_){
      if (cat >= 0 && cat < (int)flashggCats_.size()) {
        catname = Form("%s",flashggCats_[cat].c_str());
      } else {
        std::cerr << "[ERROR] Invalid category index " << cat << " for flashggCats_ (size=" << flashggCats_.size() << ")" << std::endl;
        catname = Form("cat%d",cat);
      }
    } else {
      catname = Form("cat%d",cat);
    }
    
    // Option 1: Use as input an unbinned RooDataSet and bin it
    /*
    RooDataSet *dataFull;
    RooDataSet *dataFull0;
    if (isData_) {
    dataFull = (RooDataSet*)inWS->data(Form("Data_13TeV_%s",catname.c_str()));
    if (verbose) std::cout << "[INFO] opened data for  "  << Form("Data_%s",catname.c_str()) <<" - " << dataFull <<std::endl;
    }
    else 
    {dataFull = (RooDataSet*)inWS->data(Form("data_mass_%s",catname.c_str()));
    if (verbose) std::cout << "[INFO] opened data for  "  << Form("data_mass_%s",catname.c_str()) <<" - " << dataFull <<std::endl;
    }

    // mass->setBins(nBinsForFit); // Binning now taken from workspace
    RooDataSet *data;
    string thisdataBinned_name;

    if ( isFlashgg_){
      thisdataBinned_name =Form("CAT_roohist_data_mass_%s",flashggCats_[cat].c_str());
    } else {
      thisdataBinned_name= Form("CAT_roohist_data_mass_cat%d",cat);
    }
    RooDataHist thisdataBinned(thisdataBinned_name.c_str(),"data",*mass,*dataFull);
    data = (RooDataSet*)&thisdataBinned; 
    // both "data" and "thisdataBinned" are binned (number of bins is given by "mass" bins) 
    // "dataFull" is unbinned
    */

    // Option 2 (equivalente): Usa come input un RooDataHist già binned
    string data_name ="h_data_bkg_catRES";// Form("CAT_roohist_data_mass_%s",flashggCats_[cat].c_str()); 
    //string data_name ="rooHist_data_cat2";// Form("CAT_roohist_data_mass_%s",flashggCats_[cat].c_str()); 
    //string data_name ="data_obs";// Form("CAT_roohist_data_mass_%s",flashggCats_[cat].c_str()); 
    RooDataHist  *data_in       = (RooDataHist*)dataWS->data(data_name.c_str());
    if (!data_in) {
      cerr << "[ERROR] Could not find dataset '" << data_name << "' in workspace" << endl;
      continue; // Salta questa categoria e prova con la successiva
    }
    
    // Ottieni il binning dalla variabile di massa nella workspace (che dovrebbe corrispondere al RooDataHist)
    nBinsForFit = mass->getBinning().numBins();
    cout << "[INFO] Using binning from workspace mass variable: " << nBinsForFit << " bins" << endl;
    
    std::cout << "entries dataset " << data_in->sumEntries() << std::endl;
    
    // Debug: Stampa alcune informazioni sul RooDataHist originale
    cout << "[DEBUG] Original RooDataHist info:" << endl;
    data_in->Print("v");
    
    // Implementa il rebinning se richiesto
    RooAbsData *data;
    if (rebinFactor > 1) {
      cout << "[INFO] Rebinning histogram with factor " << rebinFactor << endl;
      
      // Ottieni la variabile di massa originale e il suo binning
      RooRealVar* originalMass = (RooRealVar*)data_in->get()->first();
      double xmin = originalMass->getMin();
      double xmax = originalMass->getMax();
      int originalBins = originalMass->getBinning().numBins();
      int newBins = originalBins / rebinFactor;
      
      cout << "[DEBUG] Original bins: " << originalBins << ", requested new bins: " << newBins << endl;
      
      // Controlla se il rebinning è fattibile
      if (originalBins % rebinFactor != 0) {
        cout << "[WARNING] Original bins (" << originalBins << ") not divisible by rebinFactor (" << rebinFactor << ")" << endl;
        cout << "[INFO] Adjusting newBins to closest integer: " << newBins << endl;
      }
      
      // Crea una nuova variabile di massa con il binning rebinned
      RooRealVar* rebinnedMass = new RooRealVar(originalMass->GetName(), originalMass->GetTitle(), 
                                               xmin, xmax);
      rebinnedMass->setBins(newBins);
      
      // Converti il RooDataHist originale in TH1, crea la versione corretta rebinned
      TH1* originalHist = data_in->createHistogram("temp_hist", *originalMass);
      
      // Crea un nuovo istogramma con il binning corretto invece di usare Rebin
      TH1F* rebinnedHist = new TH1F("rebinned_hist", "rebinned histogram", 
                                    newBins, xmin, xmax);
      
      // Esegui il rebinning manualmente sommando i bin adiacenti
      double binWidth = (xmax - xmin) / newBins;
      for (int i = 1; i <= newBins; i++) {
        double newBinCenter = xmin + (i - 0.5) * binWidth;
        double binContent = 0;
        double binError2 = 0;
        
        // Trova tutti i bin originali che contribuiscono a questo nuovo bin
        double binLow = xmin + (i - 1) * binWidth;
        double binHigh = xmin + i * binWidth;
        
        for (int j = 1; j <= originalBins; j++) {
          double origBinCenter = originalHist->GetBinCenter(j);
          if (origBinCenter >= binLow && origBinCenter < binHigh) {
            binContent += originalHist->GetBinContent(j);
            binError2 += originalHist->GetBinError(j) * originalHist->GetBinError(j);
          }
        }
        
        rebinnedHist->SetBinContent(i, binContent);
        rebinnedHist->SetBinError(i, sqrt(binError2));
      }
      
      // Crea un nuovo RooDataHist dall'istogramma rebinned manualmente
      RooDataHist* rebinnedData = new RooDataHist("rebinned_data", "rebinned data", 
                                                  RooArgSet(*rebinnedMass), rebinnedHist);
      
      // Aggiorna la variabile di massa per utilizzare la versione rebinned
      mass = rebinnedMass;
      nBinsForFit = newBins;
      data = rebinnedData;
      
      cout << "[INFO] Manual rebinning completed: " << originalBins << " -> " << newBins << " bins" << endl;
      
      delete originalHist;
      delete rebinnedHist;
    } else {
      // Usa direttamente il RooDataHist - cast a RooAbsData per compatibilità
      data = data_in;
    }
    
    // Debug: Stampa alcune informazioni sui dati che useremo
    cout << "[DEBUG] Using data type: " << data->ClassName() << endl;
    
    std::cout << "entries dataset " << data->sumEntries() << std::endl;
    
  


    
    RooArgList storedPdfs("store");

    fprintf(resFile,"\\multicolumn{4}{|c|}{\\textbf{Category %d}} \\\\\n",cat);
    fprintf(resFile,"\\hline\n");

    double MinimimNLLSoFar=1e10;
    int simplebestFitPdfIndex = 0;

    // Standard F-Test to find the truth functions
    for (vector<string>::iterator funcType=functionClasses.begin(); 
        funcType!=functionClasses.end(); funcType++){

      // Enable turn-on for all functions
      bool currentIncludeTurnOn = true;
      if (*funcType == "PowerLaw" || *funcType == "PowerLawSingle") {
        currentIncludeTurnOn = true;   // Enable turn-on for PowerLaw functions
        std::cout << "====> Enabling turn-on for " << *funcType << " (PowerLaw variant)" << std::endl;
      }
      // Enable turn-on for all other functions by default
      if (*funcType == "Bernstein" || *funcType == "Exponential") {
        currentIncludeTurnOn = true;
        std::cout << "====> Enabling turn-on for " << *funcType << std::endl;
      }

      std::cout << "======================================= " << std::endl;
      std::cout << "====> FAMILY " << funcType->c_str() << std::endl;
      std::cout << "======================================= " << std::endl;

      double thisNll=0.; double prevNll=0.; double chi2=0.; double prob=0; 
      int order=1; int prev_order=0; int cache_order=0;

      RooAbsPdf *prev_pdf=NULL;
      RooAbsPdf *cache_pdf=NULL;
      std::vector<int> pdforders;

      std::cout << "===> F-TEST for Truth determination" << std::endl;

      int counter =0;
      while (prob<0.05 && order < 7){ 
        cout << "==> " << *funcType << " " << order << endl;
        RooAbsPdf *bkgPdf = createBackgroundWithTurnOn(pdfsModel,*funcType,order,mass,Form("ftest_pdf_%d_%s",(cat+catOffset),ext.c_str()), currentIncludeTurnOn, includeZ, dataWS, turnOnType);
        
        // Additional check: fix normalization for any RooAddPdf
        if (bkgPdf && bkgPdf->InheritsFrom("RooAddPdf")) {
          RooArgSet normSet(*mass);
          ((RooAddPdf*)bkgPdf)->fixCoefNormalization(normSet);
        }
        
        if (!bkgPdf){
          // assume this order is not allowed
          order++;
        }
        else {

          //bkgPdf->Print();
          int fitStatus = 0;
          double chi2FromFit = -1;
          if (iterativeFit){//inserisci una chiamata a runIterativeFits
            runIterativeFits(mass, data, pdfsModel, dataWS, *funcType, order, flashggCats_, singleCategory, outDir, iterativeFit, turnOnType);
            }else{
            runFit(bkgPdf,data,&thisNll,&fitStatus,/*max iterations*/7,&chi2FromFit,blindSignalRegion);//bkgPdf->fitTo(*data,Save(true),RooFit::Minimizer("Minuit2","minimize"));
            }
          if (fitStatus!=0) std::cout << "[WARNING] Warning -- Fit status for " << bkgPdf->GetName() << " at " << fitStatus <<std::endl;
       
          chi2 = 2.*(prevNll-thisNll);
          if (chi2<0. && order>1) chi2=0.;
          if (prev_pdf!=NULL){
            prob = getProbabilityFtest(chi2,order-prev_order,prev_pdf,bkgPdf,mass,data
                ,Form("%s/Ftest_from_%s%d_cat%d.pdf",outDir.c_str(),funcType->c_str(),order,(cat+catOffset)));
            std::cout << "[INFO] F-test Prob == " << prob << std::endl;
          } else {
            prob = 0;
          }
          double gofProb=0;
          // Always make normal plot first to establish correct normalization
          plot(mass,bkgPdf,data,Form("%s/%s%d_%s",outDir.c_str(),funcType->c_str(),order,catname.c_str()),flashggCats_,fitStatus,&gofProb,chi2FromFit,true);
          
          // Also create component plot when Z is included
          if (includeZ && dataWS) {
            RooAbsPdf *zModel = dataWS->pdf("model_Z_c2");
            if (zModel) {
              plotComponents(mass,bkgPdf,zModel,data,Form("%s/%s%d_%s",outDir.c_str(),funcType->c_str(),order,catname.c_str()),flashggCats_,fitStatus,&gofProb,chi2FromFit,blindSignalRegion);
            }
          }
          cout << "[INFO] function type, order, prevNLL, thisNLL, chi2, prob " << endl;
          cout << "[INFO] " << *funcType << " " << order << " " << prevNll << " " << thisNll << " " << chi2 << " " << prob << endl;
          prevNll=thisNll;
          cache_order=prev_order;
          cache_pdf=prev_pdf;
          prev_order=order;
          prev_pdf=bkgPdf;
          order++;
        }
        counter++;
      } // end condition for performing f-test

      // next line is commented, as we want to save only the final result (that takes into account both GOF and F-test results
      //fprintf(resFile,"%15s & %d & %5.3f & %5.3f \\\\\n",funcType->c_str(),cache_order+1,chi2,prob);
      choices.insert(pair<string,int>(*funcType,cache_order));
      pdfs.insert(pair<string,RooAbsPdf*>(Form("%s%d",funcType->c_str(),cache_order),cache_pdf));

      int truthOrder = cache_order;

      // Now run loop to determine functions inside envelope
      std::cout << "===> F-TEST and GOF for ENVELOPE determination" << std::endl;
      if (saveMultiPdf){
        chi2=0.;
        thisNll=0.;
        prevNll=0.;
        prob=0.;
        order=1;
        prev_order=0;
        cache_order=0;
        std::cout << "[INFO] Upper end Threshold for highest order function " << upperEnvThreshold <<std::endl;

        while (prob<upperEnvThreshold){
          cout << "==> " << *funcType << " " << order << endl;
          RooAbsPdf *bkgPdf = createBackgroundWithTurnOn(pdfsModel,*funcType,order,mass,Form("env_pdf_%d_%s",(cat+catOffset),ext.c_str()), currentIncludeTurnOn, includeZ, dataWS, turnOnType);
          
          // Additional check: fix normalization for any RooAddPdf
          if (bkgPdf && bkgPdf->InheritsFrom("RooAddPdf")) {
            RooArgSet normSet(*mass);
            ((RooAddPdf*)bkgPdf)->fixCoefNormalization(normSet);
          }
          
          // Optionally add Z resonance  
          /*if (includeZ && bkgPdf) {
             cout << "[INFO] Adding Z resonance using independent normalizations" << endl;
      
      // Create independent normalizations for background and Z
             RooRealVar *bkg_norm = new RooRealVar(Form("%s_bkg_norm", ext), "Background normalization", 
                                           1400000, 0., 2000000.); // ~data events for background
      
             RooRealVar *z_norm = new RooRealVar(Form("%s_z_norm", ext), "Z normalization", 
                                         16000, 14000, 18000.); // Start with ~3.5% of background, wider range
      
             cout << "[INFO] Using independent normalizations:" << endl;
             cout << "  Background norm: " << bkg_norm->getVal() << " events" << endl;
             cout << "  Z norm: " << z_norm->getVal() << " events" << endl;
      
      // Combine with independent normalizations (no coefficients!)
      compositePdf = new RooAddPdf(Form("%s_with_z", ext),
                                "Background with Z",
                                RooArgList(*currentPdf, *zModel),
                                RooArgList(*bkg_norm, *z_norm));
      
      cout << "[INFO] Z added with independent normalizations" << endl;
              // Fix normalization set to avoid coefficient ambiguity
      //        RooArgSet normSet(*mass);
        //      compositePdf->fixCoefNormalization(normSet);
              
              bkgPdf = compositePdf;
            }*/
          
          if (!bkgPdf ){
            // assume this order is not allowed
            if (order >6) { std::cout << " [WARNING] could not add ] " << std::endl; break ;}
            order++;
          }
          else {
            // Fit and chi-square calculation is repeated
            //RooFitResult *fitRes;
            int fitStatus=0;
            double chi2FromFit = -1;
            runFit(bkgPdf,data,&thisNll,&fitStatus,/*max iterations*/7,&chi2FromFit,blindSignalRegion);//bkgPdf->fitTo(*data,Save(true),RooFit::Minimizer("Minuit2","minimize"));
            //thisNll = fitRes->minNll();
            if (fitStatus!=0) std::cout << "[WARNING] Warning -- Fit status for " << bkgPdf->GetName() << " at " << fitStatus <<std::endl;
            double myNll = 2.*thisNll;
            chi2 = 2.*(prevNll-thisNll);
            if (chi2<0. && order>1) chi2=0.;
            prob = TMath::Prob(chi2,order-prev_order); 

            cout << "[INFO] function type, order, prevNLL, thisNLL, chi2, prob " << endl;
            cout << "[INFO] " << *funcType << " " << order << " " << prevNll << " " << thisNll << " " << chi2 << " " << prob << endl;
            prevNll=thisNll;
            cache_order=prev_order;
            cache_pdf=prev_pdf;

            // Calculate goodness of fit (will use toys for lowstats)
            double gofProb =0; 
            plot(mass,bkgPdf,data,Form("%s/%s%d_%s",outDir.c_str(),funcType->c_str(),order,catname.c_str()),flashggCats_,fitStatus,&gofProb,chi2FromFit);
            
            // Also create component plot when Z is included
            if (includeZ && dataWS) {
              RooAbsPdf *zModel = dataWS->pdf("model_Z_c2");
              if (zModel) {
                plotComponents(mass,bkgPdf,zModel,data,Form("%s/%s%d_%s",outDir.c_str(),funcType->c_str(),order,catname.c_str()),flashggCats_,fitStatus,&gofProb,chi2FromFit,blindSignalRegion);
              }
            }

            if ((prob < upperEnvThreshold) ) { // Looser requirements for the envelope

              if (gofProb > minGofThreshold || order == truthOrder ) {  // Good looking fit or one of our regular truth functions
              //if (gofProb > minGofThreshold) { // minimal requirement on the goodness of fit

                std::cout << "[INFO] Adding to Envelope " << bkgPdf->GetName() << " "<< gofProb 
                  << " 2xNLL + c is " << myNll + bkgPdf->getVariables()->getSize() <<  std::endl;
                allPdfs.insert(pair<string,RooAbsPdf*>(Form("%s%d",funcType->c_str(),order),bkgPdf));
                storedPdfs.add(*bkgPdf);
                pdforders.push_back(order);

                // Keep track but we shall redo this later
                if ((myNll + bkgPdf->getVariables()->getSize()) < MinimimNLLSoFar) {
                  simplebestFitPdfIndex = storedPdfs.getSize()-1;
                  MinimimNLLSoFar = myNll + bkgPdf->getVariables()->getSize();
                }
            //  }
            }
            prev_order=order;
            prev_pdf=bkgPdf;
            order++;
          }
  }
        } // end while

        fprintf(resFile,"%15s & %d & %5.3f & %5.3f \\\\\n",funcType->c_str(),cache_order+1,chi2,prob);
        choices_envelope.insert(pair<string,std::vector<int> >(*funcType,pdforders));
      }
    } // end loop over families

    fprintf(resFile,"\\hline\n");
    choices_vec.push_back(choices);
    choices_envelope_vec.push_back(choices_envelope);
    pdfs_vec.push_back(pdfs);

    plot(mass,pdfs,data,Form("%s/truths_%s",outDir.c_str(),catname.c_str()),flashggCats_,cat,true);

    if (saveMultiPdf){
      // Put selectedModels into a MultiPdf
      string catindexname;
      string catname;
      if (isFlashgg_){
        catindexname = Form("pdfindex_%s_%s",std::to_string(cat).c_str(),ext.c_str());
        catname = Form("%s",std::to_string(cat).c_str());
      } else {
        catindexname = Form("pdfindex_%d_%s",(cat+catOffset),ext.c_str());
        catname = Form("cat%d",(cat+catOffset));
      }
      RooCategory catIndex(catindexname.c_str(),"c");
      RooRealVar nBackground(Form("CMS_hgg_%s_%s_bkgshape_norm",catname.c_str(),ext.c_str()),"nbkg",data->sumEntries(),0,3*data->sumEntries());
      RooMultiPdf *pdf = new RooMultiPdf(Form("CMS_hgg_%s_%s_bkgshape",catname.c_str(),ext.c_str()),"all pdfs",catIndex,storedPdfs);
       //nBackground.removeRange(); // bug in roofit will break combine until dev branch brought in
      //double check the best pdf!
      int bestFitPdfIndex = getBestFitFunction(pdf,data,&catIndex,!verbose);
     
      catIndex.setIndex(bestFitPdfIndex);
      std::cout << "// ------------------------------------------------------------------------- //" <<std::endl; 
      std::cout << "[INFO] Created MultiPdf " << pdf->GetName() << ", in Category " << cat << " with a total of " << catIndex.numTypes() << " pdfs"<< std::endl;
      storedPdfs.Print();
      std::cout << "[INFO] Best Fit Pdf = " << bestFitPdfIndex << ", " << storedPdfs.at(bestFitPdfIndex)->GetName() << std::endl;

      std::cout << "[INFO] Simple check of index "<< simplebestFitPdfIndex <<std::endl;
      std::cout << "// ------------------------------------------------------------------------- //" <<std::endl;

      // mass->setBins(nBinsForFit); // Binning now taken from workspace
      //RooDataHist dataBinned(Form("roohist_data_mass_%s",catname.c_str()),"data",*mass,*dataFull);

      // Save it (also a binned version of the dataset
      outputws->import(*pdf);
      outputws->import(nBackground);
      outputws->import(catIndex);
      //outputws->import(dataBinned);
      outputws->import(*data);
      plot(mass,pdf,&catIndex,data,Form("%s/multipdf_%s",outDir.c_str(),catname.c_str()),flashggCats_,cat,bestFitPdfIndex);

    } // end if saveMultiPdf

  }// end loop over categories 

  if (saveMultiPdf){
    outputfile->cd();
    outputws->Write();
    outputfile->Close();  
  }

  // Write recommended options to screen and to file
/*  FILE *dfile = fopen(datfile.c_str(),"w");
  cout << "[RESULT] Recommended options based on truth" << endl;

  for (int cat=startingCategory; cat<ncats; cat++){
    cout << "Cat " << cat << endl;
    fprintf(dfile,"cat=%d\n",(cat+catOffset)); 
    for (map<string,int>::iterator it=choices_vec[cat-startingCategory].begin(); it!=choices_vec[cat-startingCategory].end(); it++){
      cout << "\t" << it->first << " - " << it->second << endl;
      fprintf(dfile,"truth=%s:%d:%s%d\n",it->first.c_str(),it->second,namingMap[it->first].c_str(),it->second);
    }
    fprintf(dfile,"\n");
  }


  cout << "[RESULT] Recommended options for envelope" << endl;
  for (int cat=startingCategory; cat<ncats; cat++){
    cout << "Cat " << cat << endl;
    fprintf(dfile,"cat=%d\n",(cat+catOffset)); 
    for (map<string,std::vector<int> >::iterator it=choices_envelope_vec[cat-startingCategory].begin(); it!=choices_envelope_vec[cat-startingCategory].end(); it++){
      std::vector<int> ords = it->second;
      for (std::vector<int>::iterator ordit=ords.begin(); ordit!=ords.end(); ordit++){
        cout << "\t" << it->first << " - " << *ordit << endl;
        fprintf(dfile,"envel=%s:%d:%s%d\n",it->first.c_str(),*ordit,namingMap[it->first].c_str(),*ordit);
      }
    }
    fprintf(dfile,"\n");
  }*/

  inFile->Close();
}
