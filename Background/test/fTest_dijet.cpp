#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <limits>
#include <sstream>
#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <regex>

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
#include "RooEffProd.h"
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
#include "RooHistPdf.h"
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

// --------------------------- Define Helpers ---------------------------------------
struct Interval { double lo, hi; };

static std::string trimcpy(const std::string &s) {
  auto a = s.begin(); while (a!=s.end() && std::isspace(static_cast<unsigned char>(*a))) ++a;
  auto b = s.end();   while (b!=a && std::isspace(static_cast<unsigned char>(*(b-1)))) --b;
  return std::string(a,b);
}
static std::vector<std::string> split_commas(const std::string &s) {
  std::vector<std::string> out; std::string cur; std::stringstream ss(s);
  while (std::getline(ss, cur, ',')) out.push_back(trimcpy(cur));
  return out;
}
static bool parse_endpoint(const std::string &tok, double massMin, double massMax, double &val) {
  std::string t = trimcpy(tok);
  if (t.empty() || t=="min") { val = massMin; return true; }
  if (t=="max") { val = massMax; return true; }
  char *end=nullptr; val = std::strtod(t.c_str(), &end);
  return end && *end=='\0';
}
static std::vector<Interval> mergeIntervals(std::vector<Interval> v) {
  if (v.empty()) return v;
  std::sort(v.begin(), v.end(), [](auto &a, auto &b){ return a.lo<b.lo || (a.lo==b.lo && a.hi<b.hi);});
  std::vector<Interval> m; m.reserve(v.size()); m.push_back(v[0]);
  for (size_t i=1;i<v.size();++i) {
    if (v[i].lo <= m.back().hi) m.back().hi = std::max(m.back().hi, v[i].hi);
    else m.push_back(v[i]);
  }
  return m;
}
static std::vector<Interval> complement(const std::vector<Interval>& excl, double a, double b) {
    std::vector<Interval> out;
    double cur = a;
    for (const auto& iv : excl) {
        double L = std::max(a, iv.lo);
        double H = std::min(b, iv.hi);
        if (H <= a || L >= b) {
            continue;
        }
        if (L > cur) {
            out.push_back({cur, L});
        }
        cur = std::max(cur, H);
    }

    if (cur < b) {
        out.push_back({cur, b});
    }
    return out;
}
static std::vector<Interval> parse_intervals(const std::string &spec, double massMin, double massMax) {
    std::vector<Interval> out;

    if (spec.empty()) {
        return out;
    }

    for (const auto& tok : split_commas(spec)) {
        if (tok.empty()) {
            continue;
        }

        auto dash = tok.find_first_of("-:");
        if (dash == std::string::npos) {
            continue;
        }

        double lo, hi;
        bool ok1 = parse_endpoint(tok.substr(0, dash),
                                  massMin, massMax, lo);
        bool ok2 = parse_endpoint(tok.substr(dash + 1),
                                  massMin, massMax, hi);
        if (ok1 && ok2) {
            if (hi < lo) {
                std::swap(lo, hi);
            }
            out.push_back({lo, hi});
        }
    }

    for (auto& iv : out) {
        iv.lo = std::max(iv.lo, massMin);
        iv.hi = std::min(iv.hi, massMax);
    }

    out.erase(
        std::remove_if(out.begin(), out.end(),
                       [](const Interval& iv) {
                           return iv.hi <= iv.lo;
                       }),
        out.end()
    );
  return mergeIntervals(out);
}

static double nDataInRange(RooAbsData* data){
  return data->sumEntries();
}

static void plotDataAndPdf(RooRealVar* x, RooPlot* fr_crsr,
                           RooAbsData* data, RooAbsPdf* pdf,
                           const char* dataName="data", const char* pdfName="pdf",
                           int lineColor=kCyan, int lineWidth=2)
{
  // Data
  data->plotOn(fr_crsr, RooFit::Name(dataName), MarkerStyle(20), MarkerSize(0.8));

  // PDF normalised to the number of data events in the SAME plotted range
  const double n = nDataInRange(data);

  pdf->plotOn(fr_crsr, RooFit::Name(pdfName),
	      LineColor(lineColor), LineWidth(lineWidth),
	      RooFit::Normalization(n, RooAbsReal::NumEvent)); 
}

// --- Helpers for best fit functions and related legend
static std::string toLower(std::string s){
  std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c){ return std::tolower(c); });
  return s;
}

static std::string prettyPdfLabel(const RooAbsPdf* pdf){
  if(!pdf) return "NULL";
  std::string n = toLower(pdf->GetName());

  std::smatch m;

  if (std::regex_search(n, m, std::regex("(bernstein|bern)([0-9]+)"))) return "Bernstein" + m[2].str();
  if (std::regex_search(n, m, std::regex("(chebychev|cheb)([0-9]+)")))  return "Chebychev" + m[2].str();
  if (std::regex_search(n, m, std::regex("(dijet)([0-9]+)")))            return "Dijet" + m[2].str();
  if (std::regex_search(n, m, std::regex("(exp|exponential)([0-9]+)")))  return "Exponential" + m[2].str();
  if (std::regex_search(n, m, std::regex("(pow|powerlaw)([0-9]+)")))     return "PowerLaw" + m[2].str();
  if (std::regex_search(n, m, std::regex("(lau|laurent)([0-9]+)")))      return "Laurent" + m[2].str();

  return pdf->GetName();
}

struct FitQuality {
  double chi2red = -1; // chi2/ndof
  int ndof = -1;
  int nfloat = -1;
  int nbins = -1;
};

static int countFloatingParams(RooAbsPdf* pdf, RooAbsData* data){
  if(!pdf || !data) return 0;
  std::unique_ptr<RooArgSet> pars(pdf->getParameters(*data));
  std::unique_ptr<RooAbsCollection> floats(pars->selectByAttrib("Constant", kFALSE));
  return floats ? floats->getSize() : 0;
}

// --------------------------- End Helpers ------------------------------------------

// Defining maximum order allowed for the FTest (configurable from CLI)
int gMaxFtestOrder = 5;     // Reduced from 8 to 5 to reduce too high orders 
int gMaxEnvelopeOrder = 6; // Consider reducing from 10 to 6 to keep the envelope under control

// Testing alternative criteria to GoF
const double maxChi2red = 20.0;        // chi2/ndof guardrail; avoids GoF p-value underflow issues
const int    maxPerFamily = 3;        // keep up to N orders per family (set to 999 to keep all)
const int    maxTotalEnvelope = 20;   // cap total PDFs per category to avoid blow-ups

// Mass range configuration: [mN - nsigma*sigma, mN + nsigma*sigma] 
float mN = 2.75;
float sigma = 0.025;
int nsigma = 10;
float mN_low  = 50;
float mN_high = 4500;
//float mN_low  = 300;
//float mN_high = 1000;

// Not configured in main
int nBinsForFit  = (mN_high-mN_low)/2;  // 2 GeV per bin instead of 0.25 GeV
int nBinsForPlot = (mN_high-mN_low)/2;

RooRealVar *intLumi_ = new RooRealVar("IntLumi","hacked int lumi", 1000.);

TRandom3 *RandomGen = new TRandom3();

RooAbsPdf* getPdf(PdfModelBuilder &pdfsModel, string type, int order, const char* ext=""){
  if (type=="Bernstein") return pdfsModel.getBernstein(Form("%s_bern%d",ext,order),order); 
  else if (type=="Exponential") return pdfsModel.getExponentialSingle(Form("%s_exp%d",ext,order),order); 
  else if (type=="ExponentialSum") return pdfsModel.getExponential(Form("%s_expsum%d",ext,order),order); 
  else if (type=="PowerLaw") {
    if (order > 3) {  // Increased from 3 to 4 to allow more functions in the envelope
      std::cout << "[WARNING] PowerLaw order " << order << " may cause numerical instability, skipping..." << std::endl;
      return NULL;
    }
    return pdfsModel.getPowerLawSingle(Form("%s_pow%d",ext,order),order); 
  } 
  else if (type=="PowerLawSingle") return pdfsModel.getPowerLawSingle(Form("%s_powsing%d",ext,order),order); 
  else if (type=="PowerLawSum") return pdfsModel.getPowerLaw(Form("%s_pow%d",ext,order),order); 
  else if (type=="PowerLawGeneric") return pdfsModel.getPowerLawGeneric(Form("%s_powgen%d",ext,order),order); 
  else if (type=="PowerLawSimple") return pdfsModel.getPowerLawSimple(Form("%s_powsimp%d",ext,order),order); 
  else if (type=="PowerLawCutoff") return pdfsModel.getPowerLawCutoff(Form("%s_powcut%d",ext,order),order); 
  else if (type=="ExpLog") return pdfsModel.getExpLog(Form("%s_explog%d",ext,order),order); 
  else if (type=="Dijet") {
    if (order > 4) {
      std::cout << "[WARNING] Dijet order " << order << " may cause numerical instability, skipping..." << std::endl;
      return NULL;
    }
    return pdfsModel.getDijet(Form("%s_dijet%d",ext,order),order); 
  } 
  else if (type=="Laurent") return pdfsModel.getLaurentSeries(Form("%s_lau%d",ext,order),order); 
  else if (type=="Chebychev") return pdfsModel.getChebychev(Form("%s_cheb%d",ext,order),order); 
  else {
    cerr << "[ERROR] -- getPdf() -- type " << type << " not recognised." << endl;
    return NULL;
  }
}


struct Chi2Result { double chi2, chi2red; int ndof; double pval; };

Chi2Result computeChi2ReducedExtended(RooRealVar* mass, RooAbsPdf* pdf_in, RooAbsData* data,
                                      const char* rangeList)
{
  RooDataHist* dh = dynamic_cast<RooDataHist*>(data);
  std::unique_ptr<RooDataHist> dh_owner;
  if (!dh) {
    dh_owner = std::make_unique<RooDataHist>("__chi2_dh","__chi2_dh", RooArgSet(*mass), *data);
    dh = dh_owner.get();
  }

  // Make an extended pdf for chi2 calculation (so "NumEvent" normalization is meaningful)
  RooAbsPdf* pdf_ext = nullptr;
  std::unique_ptr<RooExtendPdf> owned_ext;
  std::unique_ptr<RooRealVar> owned_norm; // keep norm alive

  if (pdf_in->InheritsFrom("RooExtendPdf")) {
    pdf_ext = pdf_in;
  } else {
    // Use expectedEvents if available; otherwise derive a norm from the integral in the range
    double norm_val = pdf_in->expectedEvents(RooArgSet(*mass));
    if (!(norm_val > 0 && std::isfinite(norm_val))) {
      auto frac = std::unique_ptr<RooAbsReal>(
        pdf_in->createIntegral(RooArgSet(*mass),
                               RooFit::NormSet(RooArgSet(*mass)),
                               RooFit::Range(rangeList))
      );
      const double n_in_range = dh->sumEntries(nullptr, rangeList);
      const double f_in_range = frac ? frac->getVal() : 1.0;
      norm_val = (f_in_range > 0) ? (n_in_range / f_in_range) : n_in_range;
    }
    owned_norm = std::make_unique<RooRealVar>("__gof_norm","__gof_norm", norm_val, 0., 1e30);
    owned_norm->setConstant(kTRUE);
    owned_ext = std::make_unique<RooExtendPdf>("__gof_ext","__gof_ext", *pdf_in, *owned_norm);
    pdf_ext = owned_ext.get();
  }

  RooChi2Var chi2var("__chi2","__chi2", *pdf_ext, *dh,
                     RooFit::Range(rangeList),
                     RooFit::DataError(RooAbsData::Poisson));

  const double chi2 = chi2var.getVal();

  // --- bins used: for RooDataHist, numEntries() is the number of histogram bins
  int nbins_used = dh->numEntries();
  if (nbins_used <= 0) nbins_used = 1;

  // --- floating params
  std::unique_ptr<RooArgSet> allPars(pdf_ext->getParameters(*dh));
  std::unique_ptr<RooAbsCollection> floats(allPars->selectByAttrib("Constant", kFALSE));
  const int nfloat = floats ? floats->getSize() : 0;

  int ndof = nbins_used - nfloat;
  if (ndof < 1) ndof = 1;

  const double chi2red = chi2 / ndof;
  const double pval = TMath::Prob(chi2, ndof);

  return {chi2, chi2red, ndof, pval};
}

void runFit(RooRealVar* mass, RooAbsPdf *pdf, RooAbsData *data, double *NLL, int *stat_t, int MaxTries, FitQuality *q){
  int ntries=0;
  RooArgSet *params_test = pdf->getParameters((const RooArgSet*)(0));
  
  data->Print("v");
  params_test->Print("v");

  int stat=1;
  double minnll=1e9;;

    while (stat != 0) {
    if (ntries >= MaxTries) break;

    std::cout << "[INFO] Fitting with full mass range" << std::endl;

    std::unique_ptr<RooFitResult> fitTest(
      pdf->fitTo(*data,
                 RooFit::Save(true),
                 RooFit::Minimizer("Minuit2","minimize"),
                 RooFit::Strategy(0),
                 RooFit::PrintLevel(-1),
                 RooFit::Optimize(true),
                 RooFit::SumW2Error(kFALSE))
    );

    stat   = fitTest->status();
    minnll = fitTest->minNll();

    if (stat != 0) params_test->assignValueOnly(fitTest->randomizePars());
    ntries++;
  }

  *stat_t = stat;
  *NLL    = minnll;

  // --- Protection on fit result:
  //if (!std::isfinite(minnll) || minnll > 1e7) {
  if (!std::isfinite(minnll)) {
    std::cout << "[WARNING] Fit resulted in invalid NLL (" << minnll << "), marking as failed" << std::endl;
    *stat_t = 999;
    *NLL = 1e99;
    if (q) { *q = FitQuality{}; }
    return;
  }

  // Fill FitQuality only if requested and fit converged
  if (q) { *q = FitQuality{}; }
  
  if (q && stat == 0) {
    mass->setRange("fullRange", mass->getMin(), mass->getMax());
    auto r = computeChi2ReducedExtended(mass, pdf, data, "fullRange");
    
    q->chi2red = r.chi2red;
    q->ndof    = r.ndof;
    //q->nfloat  = r.nfloat;
    //q->nbins   = r.nbins;
  }
}        


static int bestIndexInPdfMap(RooRealVar* mass,
                             const std::map<std::string,RooAbsPdf*>& pdfs,
                             RooAbsData* data)
{
  int bestIdx = -1;
  double bestScore = 1e300;
  int i = 0;
  for (auto const& kv : pdfs) {
    RooAbsPdf* pdf = kv.second;
    if(!pdf) { ++i; continue; }
    
    int st=0; double nll=0;
    FitQuality fq;
    runFit(mass, pdf, data, &nll, &st, 7, &fq);
    if(st!=0 || !std::isfinite(nll)) { ++i; continue; }

    int nfloat = countFloatingParams(pdf, data);
    double score = 2.0*nll + nfloat;

    if(score < bestScore) { bestScore = score; bestIdx = i; }
    ++i;
  }
  return bestIdx;
}

double getProbabilityFtest(double chi2, int ndof,RooAbsPdf *pdfNull, RooAbsPdf *pdfTest, RooRealVar *mass, RooAbsData *data, std::string name){
  std::cout << "[DEBUG] F-test calculation: chi2=" << chi2 << ", ndof=" << ndof << ", runFtestCheckWithToys=" << runFtestCheckWithToys << std::endl;
  
  // For F-test: chi2 large = significant improvement = little proability
  double prob_asym = TMath::Prob(chi2,ndof);
  std::cout << "[DEBUG] TMath::Prob result: " << prob_asym << std::endl;
  
  // Fix for overflow: if probability is too little, it is treated as a 0
  if (prob_asym < 1e-10) {
    std::cout << "[DEBUG] Probability underflow, setting to 1e-10" << std::endl;
    prob_asym = 1e-10;
  }
  
  // Fix for range: if chi2 is negative or too low, no improvement
  if (chi2 <= 0) {
    std::cout << "[DEBUG] Chi2 <= 0 (" << chi2 << "), no improvement, setting prob=0.9" << std::endl;
    prob_asym = 0.9; 
  }
  
  std::cout << "[DEBUG] Final F-test probability: " << prob_asym << std::endl;
  
  if (!runFtestCheckWithToys) return prob_asym;
  int ndata = data->sumEntries();
  RooFitResult *fitNullData;
  RooFitResult *fitTestData;
  fitNullData = pdfNull->fitTo(*data,RooFit::Save(1),RooFit::Strategy(0)
    ,RooFit::Minimizer("Minuit2","minimize"),RooFit::PrintLevel(-1),RooFit::Optimize(1));
  fitTestData = pdfTest->fitTo(*data,RooFit::Save(1),RooFit::Strategy(0)
    ,RooFit::Minimizer("Minuit2","minimize"),RooFit::PrintLevel(-1),RooFit::Optimize(1)); 
  RooArgSet *params_null = pdfNull->getParameters((const RooArgSet*)(0));
  RooArgSet preParams_null;
  params_null->snapshot(preParams_null);
  RooArgSet *params_test = pdfTest->getParameters((const RooArgSet*)(0));
  RooArgSet preParams_test;
  params_test->snapshot(preParams_test);
  int ntoys = 100;
  TCanvas *can = new TCanvas("can","can",1000,1000); //EF
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
  for (int itoy = 0 ; itoy < ntoys ; itoy++){
    params_null->assignValueOnly(preParams_null);
    params_test->assignValueOnly(preParams_test);
    RooDataSet *binnedtoy = pdfNull->generate(RooArgSet(*mass),ndata,0,1);
    int stat_n=1;
    int stat_t=1;
    int ntries = 0;
    double nllNull,nllTest;
    int MaxTries = 2;
    while (stat_n!=0){
      if (ntries>=MaxTries) break;
      RooFitResult *fitNull;
      fitNull = pdfNull->fitTo(*binnedtoy,RooFit::Save(1),RooFit::Strategy(0)
                                              ,RooFit::Minimizer("Minuit2","minimize"),RooFit::Minos(0),RooFit::Hesse(0),RooFit::PrintLevel(-1),RooFit::Optimize(1));
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
  }
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
  params_null->assignValueOnly(preParams_null);
  params_test->assignValueOnly(preParams_test);
  delete can; delete stas;
  delete gChi2;
  delete leg;
  delete lat;
  return prob_asym;
}

double getGoodnessOfFit(RooRealVar *mass, RooAbsPdf *mpdf, RooAbsData *data, std::string name, double chi2_from_fit, double *chi2_out=nullptr){
  double prob;
  int ntoys = 50;
  name+="_gofTest.pdf";
  RooRealVar norm("norm","norm",data->sumEntries(),0,10E6);
  RooExtendPdf *pdf = new RooExtendPdf("ext","ext",*mpdf,norm);
  
  // Use chi2 already calculated during fit instead of recalculating
  int np = pdf->getParameters(*data)->getSize();
  double chi2_reduced = chi2_from_fit;
  
  // Get actual number of bins used in the fit
  int actualBins = 0;
  TH1* dataHist = data->createHistogram("temp", *mass);
  if (dataHist) {
    actualBins = dataHist->GetNbinsX();
    delete dataHist;
  } else {
    actualBins = nBinsForFit; // fallback
  }
  int ndof = actualBins - np;
  double chi2_total = chi2_reduced * ndof;
  std::cout << "[DEBUG] GOF: actualBins=" << actualBins << ", np=" << np << ", ndof=" << ndof << std::endl;
  std::cout << "[DEBUG] GOF: chi2_reduced=" << chi2_reduced << ", chi2_total=" << chi2_total << std::endl;
  std::cout << "[INFO] Calculating GOF for pdf " << pdf->GetName() << ", using " <<np << " fitted parameters" <<std::endl;
  if ((double)data->sumEntries()/nBinsForFit < 5 ){
    std::cout << "[INFO] Running toys for GOF test " << std::endl;
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
      RooPlot *plot_t = mass->frame();
      binnedtoy->plotOn(plot_t);
      pdf->fitTo(*binnedtoy,RooFit::Save(1),RooFit::Strategy(0)
                                          ,RooFit::Minimizer("Minuit2","minimize"),RooFit::Minos(0),RooFit::Hesse(0),RooFit::PrintLevel(-1),RooFit::Optimize(1));
      pdf->plotOn(plot_t);
      double chi2_t = plot_t->chiSquare(np);
      if( chi2_t>=chi2_reduced) npass++;
      toy_chi2.push_back(chi2_t*ndof);
      delete plot_t;
    }
    std::cout << "[INFO] complete" << std::endl;
    prob = (double)npass / ntoys;
    TCanvas *can = new TCanvas("can","can",1000,1000); //EF
    double medianChi2 = toy_chi2[(int)(((float)ntoys)/2)];
    double rms = TMath::Sqrt(medianChi2);
    TH1F toyhist(Form("gofTest_%s.pdf",pdf->GetName()),";Chi2;",50,medianChi2-5*rms,medianChi2+5*rms);
    for (std::vector<double>::iterator itx = toy_chi2.begin();itx!=toy_chi2.end();itx++){
      toyhist.Fill((*itx));
    }
    toyhist.Draw();
    TArrow lData(chi2_total,toyhist.GetMaximum(),chi2_total,0);
    lData.SetLineWidth(2);
    lData.Draw();
    can->SaveAs(name.c_str());
    params->assignValueOnly(preParams);
  } else {
    prob = TMath::Prob(chi2_total, ndof);
  }
  std::cout << "[INFO] GOF Chi2 total =  " << chi2_total << std::endl;
  std::cout << "[INFO] GOF Chi2/ndof =  " << chi2_reduced << std::endl;
  std::cout << "[INFO] GOF ndof =  " << ndof << std::endl;
  std::cout << "[INFO] GOF p-value  =  " << prob << std::endl;
  if (chi2_out) *chi2_out = chi2_reduced;  // Return chi2/ndof for display
  delete pdf;
  return prob;
}

void plot(RooRealVar *mass, RooAbsPdf *pdf, RooAbsData *data, string name, vector<string> flashggCats_,
	  int status, double *prob, const FitQuality &fq, bool showPulls=false,RooWorkspace* dataWS = nullptr){

  double chi2tot = fq.chi2red * fq.ndof;
  std::cout << "[DEBUG] chi2red=" << fq.chi2red
	    << " ndof=" << fq.ndof
	    << " chi2tot=" << chi2tot
	    << " p=" << TMath::Prob(chi2tot, fq.ndof)
	    << std::endl;

  double pval = (status==0 && fq.chi2red>0 && fq.ndof>0)
    ? TMath::Prob(fq.chi2red * fq.ndof, fq.ndof)
    : -1.0;

  if (prob) *prob = pval;
  
  RooPlot *plot = mass->frame(mass->getMin(), mass->getMax());

  data->plotOn(plot, RooFit::Name("data"),
	       MarkerStyle(20),
	       MarkerSize(0.8));
  
  const double n = data->sumEntries();
  
  pdf->plotOn(plot, Name("total"),
	      LineColor(kCyan), LineWidth(2),
	      Normalization(n, RooAbsReal::NumEvent));
  
  TCanvas *canv;

  if (showPulls) {                                                      
    canv = new TCanvas("canv","canv",1000,1000);                                                                                                     
    canv->Divide(1,2);                                                                                                       
    canv->cd(1);                                                                             
    gPad->SetLeftMargin(0.15);                                                                                     
    gPad->SetRightMargin(0.1);                                                                                     
    gPad->SetBottomMargin(0.02);                                                                                    
    gPad->SetPad(0.01,0.3,0.99,0.99);                                                                        
    gPad->SetLogy(); // Set log scale for y-axis                                                                
    plot->GetXaxis()->SetTitle("");                                                                           
    plot->GetXaxis()->SetLabelSize(0);                                                                          
    plot->GetYaxis()->SetTitleSize(0.05);                                                                     
    plot->GetYaxis()->SetLabelSize(0.04);                                                             
    plot->GetYaxis()->SetTitleOffset(1.2);                                                         
    plot->GetYaxis()->SetTitle("Events");                                                               
  } else {                                                                                    
    canv = new TCanvas("canv","canv",1000,1000);                                                                                                     
    canv->SetLogy(); // Set log scale for y-axis                                                   
    gPad->SetLeftMargin(0.15);                                                                            
    gPad->SetRightMargin(0.1);                                                                                     
    gPad->SetBottomMargin(0.10);                                                                  
    plot->GetXaxis()->SetTitleSize(0.04);                                                     
    plot->GetXaxis()->SetTitle("m_{jj} [GeV]");                                              
    plot->GetYaxis()->SetTitleSize(0.04);                  
    plot->GetYaxis()->SetTitleOffset(1.35);                 
  }
    
  plot->SetTitle("");
  plot->SetMaximum(plot->GetMaximum()*1.4);
  double minVal = plot->GetMinimum();

  // Better minimum for log scale: find a reasonable minimum based on data
  double dataMinimum = minVal;
  if (minVal <= 0) dataMinimum = 0.01;
  else if (minVal < 1e-3) dataMinimum = 1e-3;
  else dataMinimum = minVal * 0.3;  // 30% below minimum for log scale
  plot->SetMinimum(dataMinimum);

  plot->Draw();
  
  if (!showPulls) {
    pdf->paramOn(plot,RooFit::Layout(0.34,0.85,0.89),RooFit::Format("NEA",AutoPrecision(1)));
    plot->getAttText()->SetTextSize(0.025);
  }

  TLatex *lat = new TLatex();
  lat->SetNDC();
  lat->SetTextFont(42);
  lat->SetTextSize(0.034);

  double p = (status==0 && fq.ndof>0) ? TMath::Prob(chi2tot, fq.ndof) : 1.0;
  double logp = -std::log10(std::max(p, std::numeric_limits<double>::min()));
  
  lat->DrawLatex(0.15,0.94,
		 Form("#chi^{2}/ndof = %.3f, -log_{10}(p)=%.1f, status=%d",
		      fq.chi2red, logp, status));
  
  //double logp = (pval>0) ? -std::log10(pval) : 999.;
  //lat->DrawLatex(0.15,0.94,
  //		 Form("#chi^{2}/ndof = %.3f, -log_{10}(p)=%.1f, Fit Status=%d",
  //		      fq.chi2red, logp, status));  
  //lat->DrawLatex(0.15,0.94,
  //		 Form("#chi^{2}/ndof = %.3f, Prob = %.2f, Fit Status = %d",
  //		      fq.chi2red,
  //		      pval,
  //		      status));   
  
  if (!showPulls) {
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
      s.erase(0, pos + delimiter.length());
    }
    lab->DrawLatex(0.55,0.25,s.c_str());
    delete lab;
  }
  
  if (showPulls) {
    canv->cd(2);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.1);                                                                                     
    gPad->SetBottomMargin(0.4);
    gPad->SetPad(0.01,0.01,0.99,0.3);
    gPad->SetGridy();

    // Build pulls using objects already on the TOP plot:
    RooHist* hpull = plot->pullHist("data", "total", true);
    RooPlot* pullFrame = mass->frame(mass->getMin(), mass->getMax());
    pullFrame->addPlotable(hpull, "P");
    
    pullFrame->SetTitle("");
    pullFrame->GetYaxis()->SetTitle("Pull");
    pullFrame->GetYaxis()->SetNdivisions(504);
    pullFrame->GetYaxis()->SetTitleSize(0.12);
    pullFrame->GetYaxis()->SetLabelSize(0.10);
    pullFrame->GetYaxis()->SetTitleOffset(0.45);
    
    pullFrame->GetXaxis()->SetTitle("m_{jj} [GeV]");
    pullFrame->GetXaxis()->SetTitleSize(0.14);
    pullFrame->GetXaxis()->SetLabelSize(0.12);
    pullFrame->GetXaxis()->SetTitleOffset(1.05);
    
    pullFrame->SetMinimum(-3);
    pullFrame->SetMaximum( 3);
    
    pullFrame->Draw();
  }

  canv->SaveAs(Form("%s.pdf",name.c_str()));
  canv->SaveAs(Form("%s.png",name.c_str()));
  delete canv;
  delete lat;
}

void plot(RooRealVar *mass, RooMultiPdf *pdfs, RooCategory *catIndex, RooAbsData *data, string name, vector<string> flashggCats_, int cat, int bestFitPdf=-1){
  //int color[7] = {kCyan+1,kGreen+2,kAzure+2,kTeal+1,kSpring+2,kCyan+3,kGreen+1};
  int color[10] = {
    kCyan+1,
    kPink+7,
    kAzure+6,
    kMagenta+1,
    kOrange+7,
    kBlue+1,
    kGreen+2,
    kRed+1,
    kViolet+1,
    kBlack
  };

  TLegend *leg = new TLegend(0.6,0.65,0.95,0.90);
  leg->SetFillColor(0);
  leg->SetLineColor(1);

  RooPlot *plot = mass->frame();
  data->plotOn(plot,MarkerStyle(20),MarkerSize(0.8));
  TCanvas *canv = new TCanvas();
  canv->SetLogy(); // Set log scale for y-axis
  int currentIndex = catIndex->getIndex();
  TObject *datLeg = plot->getObject(int(plot->numItems()-1));
  leg->AddEntry(datLeg,"Data");
  int style=1;
  RooAbsPdf *pdf;
  RooCurve *nomBkgCurve;
  double Nvis;
  Nvis = data->sumEntries();
  std::cout << "[INFO] Unblinded fit: Nvis = " << Nvis << " events total" << std::endl;

  for (int icat=0;icat<catIndex->numTypes();icat++){
    int col;
    if (icat<=6) col=color[icat];
    else {col=kBlack; style++;}
    catIndex->setIndex(icat);
    //pdfs->getCurrentPdf()->fitTo(*data,RooFit::Minos(0),RooFit::Minimizer("Minuit2","minimize"),RooFit::Strategy(0),RooFit::PrintLevel(-1),RooFit::Optimize(1));  

    pdfs->getCurrentPdf()->plotOn(plot,LineColor(col),LineStyle(style));
    TObject *pdfLeg = plot->getObject(int(plot->numItems()-1));
    std::string ext = "";

    if (bestFitPdf==icat) {
      ext=" (Best Fit Pdf) ";
      pdf= pdfs->getCurrentPdf();
      nomBkgCurve = (RooCurve*)plot->getObject(plot->numItems()-1);
    }
    
    string pdfName = pdfs->getCurrentPdf()->GetName();
    std:: cout << "PDF Name: " << pdfName << std::endl; 
    
    std::cout << "[DEBUG] Full PDF name: " << pdfName << std::endl;
    
    // Simplified logic to extract the name of the MultiPdf
    //std::string legendName = Form("PDF_%d", icat); // fallback
    std::string legendName = prettyPdfLabel(pdfs->getCurrentPdf());

    std::cout << "[DEBUG] Legend name: " << legendName << std::endl;
    leg->AddEntry(pdfLeg, (legendName + ext).c_str(), "L");
  }
  plot->SetTitle(Form("Category %d",cat));
  plot->SetMaximum(plot->GetMaximum()*1.4);
  double minVal = plot->GetMinimum();
  // Better minimum for log scale: find a reasonable minimum based on data
  double dataMinimum = minVal;
  if (minVal <= 0) dataMinimum = 0.01;
  else if (minVal < 1e-3) dataMinimum = 1e-3;
  else dataMinimum = minVal * 0.3;  // 30% below minimum for log scale
  plot->SetMinimum(dataMinimum);
  plot->GetXaxis()->SetTitle("m_{l\\pi} (GeV)");
  plot->Draw();
  leg->Draw("same");
  CMS_lumi( canv, 0, 0);
  canv->SaveAs(Form("%s.pdf",name.c_str()));
  canv->SaveAs(Form("%s.png",name.c_str()));
  catIndex->setIndex(currentIndex);
  delete canv;
}

void plot(RooRealVar *mass, map<string,RooAbsPdf*> pdfs, RooAbsData *data, string name, vector<string> flashggCats_, int cat, const std::string& catnameDijet, int bestFitPdf=-1){

  //int color[7] = {kCyan+1,kGreen+2,kAzure+2,kTeal+1,kSpring+2,kCyan+3,kGreen+1};
  int color[10] = {
    kCyan+1,
    kPink+7,
    kAzure+6,
    kMagenta+1,
    kOrange+7,
    kBlue+1,
    kGreen+2,
    kRed+1,
    kViolet+1,
    kBlack
  };
  TCanvas *canv = new TCanvas("canv","canv",1000,1000); //EF
  canv->SetLogy();
  canv->SetRightMargin(0.1);                                                                                     

  TLegend *leg = new TLegend(0.15,0.15,0.55,0.35); 
  leg->SetFillColor(0);
  leg->SetFillStyle(0);  
  leg->SetLineColor(0);
  leg->SetBorderSize(0); 
  leg->SetTextSize(0.03); 
  leg->SetMargin(0.15);  
  
  RooPlot *plot = mass->frame(mass->getMin(), mass->getMax());
  
  data->plotOn(plot,MarkerStyle(20),MarkerSize(0.8));
  TObject *datLeg = plot->getObject(int(plot->numItems()-1));
  leg->AddEntry(datLeg,Form("Data - %s",catnameDijet.c_str()),"LEP");

  int i=0;
  int style=1;
  for (map<string,RooAbsPdf*>::iterator it=pdfs.begin(); it!=pdfs.end(); it++){
    int col;
    if (i<=6) col=color[i];
    else {col=kBlack; style++;}
    double Nvis = -1;
    
    // Safety check for null pointer
    if (!it->second) {
      std::cerr << "[ERROR] Null PDF pointer found for " << it->first << std::endl;
      continue;
    }
    
    it->second->plotOn(plot,LineColor(col),LineWidth(2),Range(mass->getMin(),mass->getMax()));
    TObject *pdfLeg = plot->getObject(int(plot->numItems()-1));
    std::string ext = "";
    if (bestFitPdf==i) ext=" (Best Fit Pdf) ";
    
    // Extract function family and order from the PDF name for better legend
    std::string pdfName = it->first;
    std::string familyName = pdfName;
    std::string orderStr = "";
    
    // Extract order from PDF name (e.g., "Exponential3" -> family="Exponential", order="3")
    std::size_t lastDigitPos = pdfName.find_last_not_of("0123456789");
    if (lastDigitPos != std::string::npos && lastDigitPos < pdfName.length() - 1) {
      familyName = pdfName.substr(0, lastDigitPos + 1);
      orderStr = pdfName.substr(lastDigitPos + 1);
    }
    
    // Create better legend entry with order information
    std::string legendName = prettyPdfLabel(it->second);
    leg->AddEntry(pdfLeg, (legendName + ext).c_str(), "L");
    
    i++;
  }
  
  plot->SetMaximum(plot->GetMaximum()*1.4);
  double minVal = plot->GetMinimum();
  // Better minimum for log scale: use 0.1 or minVal/10, whichever is larger
  double logMinimum = std::max(0.1, minVal > 0 ? minVal/10.0 : 0.1);
  plot->SetMinimum(logMinimum);
  plot->SetTitle(Form("cat%d",cat));

  plot->Draw();
  leg->Draw("same");
  CMS_lumi( canv, 0, 0);
  canv->SaveAs(Form("%s.pdf",name.c_str()));
  canv->SaveAs(Form("%s.png",name.c_str()));
  delete canv;
}

void transferMacros(TFile *inFile, TFile *outFile){
  if (!inFile || !outFile) {
    std::cout << "[WARNING] transferMacros: inFile or outFile is null, skipping macro transfer" << std::endl;
    return;
  }
  TIter next(inFile->GetListOfKeys());
  TKey *key;
  while ((key = (TKey*)next())){
    if (string(key->ReadObj()->ClassName())=="TMacro") {
      TMacro *macro = (TMacro*)inFile->Get(key->GetName());
      outFile->cd();
      macro->Write();
    }
  }
}

int getBestFitFunction(RooRealVar* mass, RooMultiPdf *bkg, RooAbsData *data, RooCategory *cat, RooWorkspace* ws, const std::string& snapBase, bool silent=false){
  double best = 1e100;
  int bestIdx = 0;

  for (int id=0; id<cat->numTypes(); ++id){
    cat->setIndex(id);

    double nll=0; int st=1;
    FitQuality fq;           
    runFit(mass, bkg->getCurrentPdf(), data, &nll, &st, 7, &fq);

    double prob = -1;
    if (st==0 && fq.chi2red>0 && fq.ndof>0) {
      double chi2tot = fq.chi2red * fq.ndof;
      prob = TMath::Prob(chi2tot, fq.ndof);
    }

    if (st!=0 || !std::isfinite(nll)) nll = 1e99;
    nll += bkg->getCorrection();

    // Snapshot ONLY the parameters of the current pdf
    std::unique_ptr<RooArgSet> pars(bkg->getCurrentPdf()->getParameters(*data));
    std::string snapName = Form("%s_idx%d", snapBase.c_str(), id);
    ws->saveSnapshot(snapName.c_str(), *pars, true);

    if (!silent){
      std::cout << "[INFO] idx="<<id<<" pdf="<<bkg->getCurrentPdf()->GetName()
                <<" NLL+c="<<nll<<" status="<<st<<"\n";
    }

    if (nll < best){ best = nll; bestIdx = id; }
  }

  cat->setIndex(bestIdx);
  std::string bestSnap = Form("%s_idx%d", snapBase.c_str(), bestIdx);
  ws->loadSnapshot(bestSnap.c_str());

  std::cout << "[INFO] Best fit = idx "<<bestIdx<<" pdf="<<bkg->getCurrentPdf()->GetName()
            <<" NLL+c="<<best<<"\n";
  return bestIdx;
}

// Function to create background PDF with optional turn-on function and Z resonance
RooAbsPdf* createBackgroundWithTurnOn(PdfModelBuilder &pdfsModel, string type, int order, RooRealVar* mass,
				      const char* ext, bool includeTurnOn = true,
				      RooWorkspace* dataWS = nullptr, std::string turnOnType = "Erf") {

  RooAbsPdf *basePdf = getPdf(pdfsModel, type, order, ext);
  if (!basePdf) {
    cerr << "[ERROR] Could not create base PDF of type " << type << " order " << order << endl;
    return nullptr;
  }
  if (basePdf->InheritsFrom("RooAddPdf")) {
    RooArgSet normSet(*mass);
    ((RooAddPdf*)basePdf)->fixCoefNormalization(normSet);
  }

  RooAbsPdf *currentPdf = basePdf;

  if (includeTurnOn) {
    //RooRealVar *cutoff = new RooRealVar(Form("%s_turnon_cutoff", ext), "Turn-on cutoff", 100, 43,250 );
    //RooRealVar *beta = new RooRealVar(Form("%s_turnon_beta", ext), "Turn-on beta", 5, 0.1, 100);

    // IMPORTANT: make turn-on params unique per PDF (family+order), otherwise models share them
    // and the behavior changes when you add/remove other models.
    TString tag = Form("%s_%s%d", ext, type.c_str(), order);

    RooGenericPdf *turnOnPdf = nullptr;
    /*
    // VBF (good), RES/PP (might need tuning)
    // --------------------------------------
    RooRealVar *cutoff = new RooRealVar(Form("turnon_cutoff_%s", tag.Data()),
                                        "Turn-on cutoff", 360., 200.0, 600.0); 
    RooRealVar *beta   = new RooRealVar(Form("turnon_beta_%s", tag.Data()),
                                       "Turn-on beta",   50.0, 1.0, 300.0); 
    */

    // BOO
    // --------------------------------------
    RooRealVar *cutoff = new RooRealVar(Form("turnon_cutoff_%s", tag.Data()),
                                        "Turn-on cutoff", 80., 50.0, 450.0); 
    RooRealVar *beta   = new RooRealVar(Form("turnon_beta_%s", tag.Data()),
                                       "Turn-on beta",   5.0, 1.0, 200.0); 

    cout << "[INFO] Adding turn-on function for " << type << " order " << order << endl;

    RooFormulaVar *eff = nullptr;

    if (turnOnType == "Fermi") {
      eff = new RooFormulaVar(Form("%s_turnon_eff", tag.Data()),
			      "1.0/(1.0+TMath::Exp((@1-@0)/@2))",
			      RooArgList(*mass, *cutoff, *beta));
    } else if (turnOnType == "Erf") {
      auto eff_raw = new RooFormulaVar(
				       Form("%s_turnon_eff_raw", tag.Data()),
				       "eff raw",
				       "0.5*(1.0+TMath::Erf((@0-@1)/(TMath::Sqrt(2.)*@2)))",
				       RooArgList(*mass, *cutoff, *beta)
				       );
      
      // eff = max(eps, min(1-eps, eff_raw))
      const double eps = 1e-9;
      eff = new RooFormulaVar(
			      Form("%s_turnon_eff", tag.Data()),
			      "eff clamped",
			      Form("TMath::Min(%.12g, TMath::Max(%.12g, @0))", 1.0 - eps, eps),
			      RooArgList(*eff_raw)
			      );
    } else if (turnOnType == "DExp") {
      eff = new RooFormulaVar(Form("%s_turnon_eff", tag.Data()),
			      "1.0 - TMath::Exp(-TMath::Exp((@0-@1)/@2))",
			      RooArgList(*mass, *cutoff, *beta));
    } else {
      eff = new RooFormulaVar(Form("%s_turnon_eff", tag.Data()),
			      "1.0/(1.0+TMath::Exp((@1-@0)/@2))",
			      RooArgList(*mass, *cutoff, *beta));
    }

    if (!eff) {
      cerr << "[ERROR] Turn-on efficiency not created (turnOnType=" << turnOnType << ")\n";
      return nullptr;
    }

    currentPdf = new RooEffProd(
				Form("%s_with_turnon", tag.Data()),
				"Background with turn-on efficiency",
				*currentPdf,
				*eff
				);
  }
  
  return currentPdf;
}


int main(int argc, char* argv[]){
  setTDRStyle();
  writeExtraText = true;
  extraText  = "Preliminary";
  lumi_8TeV  = "19.1 fb^{-1}";
  lumi_7TeV  = "4.9 fb^{-1}";
  lumi_sqrtS = "13.6 TeV";

  string year_ = "2024";
  string fileName;
  string workspaceFile = "";
  int ncats;
  int singleCategory;
  int catOffset;
  string datfile;
  string outDir;
  string outfilename;
  bool is2011=false;
  bool verbose=true;
  bool saveMultiPdf=false;
  bool includeTurnOn=false;
  int rebinFactor=1;
  string flashggCatsStr_;
  vector<string> flashggCats_;
  bool isData_ =0;
  std::string turnOnType = "Erf"; // Fermi

  po::options_description desc("Allowed options");
  desc.add_options()
    ("help,h",                                                                                  "Show help")
    ("infilename,i", po::value<string>(&fileName),                                              "In file name")
    ("workspace,w", po::value<string>(&workspaceFile)->default_value(""), "Workspace file with RooDataHist and Z model")
    ("includeTurnOn", po::value<bool>(&includeTurnOn)->default_value(false),                    "Include turn-on function for low mass threshold")
    ("ncats,c", po::value<int>(&ncats)->default_value(5),                                       "Number of categories")
    ("singleCat", po::value<int>(&singleCategory)->default_value(0),                           "Run A single Category (default: category 0)")
    ("datfile,d", po::value<string>(&datfile)->default_value("dat/fTest.dat"),                  "Right results to datfile for BiasStudy")
    ("outDir,D", po::value<string>(&outDir)->default_value("plots/fTest"),                      "Out directory for plots")
    ("saveMultiPdf", po::value<string>(&outfilename),                                           "Save a MultiPdf model with the appropriate pdfs")
    ("runFtestCheckWithToys",                                                                   "When running the F-test, use toys to calculate pvals (and make plots) ")
    ("is2011",                                                                                  "Run 2011 config")
    ("is2012",                                                                                  "Run 2012 config")
    ("isData",  po::value<bool>(&isData_)->default_value(0),                                    "Use Data not MC ")
    ("flashggCats,f", po::value<string>(&flashggCatsStr_)->default_value("UntaggedTag_0,UntaggedTag_1,UntaggedTag_2,UntaggedTag_3,UntaggedTag_4,VBFTag_0,VBFTag_1,VBFTag_2,TTHHadronicTag,TTHLeptonicTag,VHHadronicTag,VHTightTag,VHLooseTag,VHEtTag"),                  "Flashgg category names to consider")
    ("year", po::value<string>(&year_)->default_value("2024"),                                  "Dataset year")
    ("catOffset", po::value<int>(&catOffset)->default_value(0),                                 "Category numbering scheme offset")
    ("mN", po::value<float>(&mN)->default_value(2.75),                                          "Mass of the peak, for center of window")
    ("sigma", po::value<float>(&sigma)->default_value(0.025),                                   "Sigma of the peak, for size of window")
    ("nsigma", po::value<int>(&nsigma)->default_value(10),                                       "Sigma multiplier, for size of window")
    ("rebinFactor,r", po::value<int>()->default_value(1),                                       "Rebin factor to reduce number of bins (1=no rebinning)")
    ("verbose,v",                                                                               "Run with more output")
    ("turnOnType", po::value<std::string>(&turnOnType)->default_value("Erf"), "Type of turn-on function: Fermi, Erf, DExp")
    // New options
    ("cat-name", po::value<std::string>()->default_value("PP") , "Category name")
    ("ws-name", po::value<std::string>()->default_value("") , "Name of the RooWorkspace within the input file")
    ("mass-var", po::value<std::string>()->default_value("CMS_dijet_mass"), "RooRealVar name for dijet mass")
    ("data-name", po::value<std::string>()->default_value("") , "Dataset name (single category)")
    ("data-name-pattern", po::value<std::string>()->default_value("h_data_bkg_cat{CAT}"), "Pattern dataset for each catorgy; {CAT} substituted with the index")
    ("mass-min", po::value<double>()->default_value(std::numeric_limits<double>::quiet_NaN()), "(Optional) Min value of the mass variable")
    ("mass-max", po::value<double>()->default_value(std::numeric_limits<double>::quiet_NaN()), "(Optional) Max value of the mass variable")
    ("max-ftest-order", po::value<int>()->default_value(5), "Max order for F-test (default: 5)")
    ("max-envelope-order", po::value<int>()->default_value(6), "Max order for the envelope (default: 6)")
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

  // Read params for maximum order considered
  // ----------------------------------------
  gMaxFtestOrder = vm["max-ftest-order"].as<int>();
  gMaxEnvelopeOrder = vm["max-envelope-order"].as<int>();
  std::cout << "[INFO] Maximum F-test order: " << gMaxFtestOrder << std::endl;
  std::cout << "[INFO] Maximum envelope order: " << gMaxEnvelopeOrder << std::endl;

  std::cout << "DEBUG mN=" << mN << std::endl; 
  mN_low  = mN - nsigma * sigma;
  mN_high = mN + nsigma * sigma;
  if (!verbose) {
    RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
    RooMsgService::instance().setSilentMode(true);
    gErrorIgnoreLevel=kWarning;
  }

  if(verbose) std::cout << "[INFO] SaveMultiPdf? " << saveMultiPdf << std::endl;
  TFile *outputfile=nullptr;
  RooWorkspace *outputws=nullptr;
  if (saveMultiPdf){
    outputfile = new TFile(outfilename.c_str(),"RECREATE");
    outputws = new RooWorkspace(); outputws->SetName("multipdf");
  }

  // Create output directories
  system(Form("mkdir -p %s", outDir.c_str()));
  system(Form("mkdir -p %s/FTestFunctions", outDir.c_str()));
  system(Form("mkdir -p %s/EnvelopeComponents_Pulls", outDir.c_str()));
 
  // Command line options
  // ----------------------------------------
  auto catnameDijet = vm["cat-name"].as<std::string>();
  auto wsName_opt  = vm["ws-name"].as<std::string>();
  auto massVarName = vm["mass-var"].as<std::string>();
  auto dataName    = vm["data-name"].as<std::string>();
  auto dataPat     = vm["data-name-pattern"].as<std::string>();
  double optMassMin = vm["mass-min"].as<double>();
  double optMassMax = vm["mass-max"].as<double>();

  // Load the main data file (if provided)
  // ----------------------------------------
  TFile *inFile = nullptr;
  RooWorkspace *inWS = nullptr;
  if (!fileName.empty()) {
    inFile = TFile::Open(fileName.c_str());
    std::cout<< "Loading file: " << fileName << std::endl;
    inWS = (RooWorkspace*)inFile->Get("wS");
  }
  // Load the workspace with RooDataHist and Z model
  // -----------------------------------------------
  TFile *wsFile = TFile::Open(workspaceFile.c_str());
  if (!wsFile || wsFile->IsZombie()) {
    cerr << "[ERROR] Cannot open workspace file: " << workspaceFile << endl;
    return 1;
  }
  RooWorkspace *dataWS = nullptr;
  if (!wsName_opt.empty()) {
    dataWS = (RooWorkspace*)wsFile->Get(wsName_opt.c_str());
  } else {
    dataWS = (RooWorkspace*)wsFile->Get("w");
    if (!dataWS) dataWS = (RooWorkspace*)wsFile->Get("cms_hgg_workspace");
    if (!dataWS) dataWS = (RooWorkspace*)wsFile->Get("wS");
    if (!dataWS) dataWS = (RooWorkspace*)wsFile->Get("workspace");
  }
  if (!dataWS) {
    cerr << "[ERROR] Cannot find workspace in " << workspaceFile << endl;
    wsFile->ls();
    return 1;
  }

  std::cout << "[INFO] Loaded workspace from: " << workspaceFile << std::endl;
  if (verbose) {
    std::cout << "[INFO] Workspace contents:" << std::endl;
    dataWS->Print();
  }
  if (!inWS) {
    inWS = dataWS;
  }
  if (verbose){
    std::cout << "------------------------------------------------------" << std::endl;
    std::cout << "[INFO]  inWS open " << inWS << std::endl;
    inWS->Print();
  }
  if (saveMultiPdf){
    if (inFile) {
      transferMacros(inFile,outputfile);
    }
    RooRealVar *intL; 
    intL  = intLumi_;

    RooRealVar* sqrts_ws = dataWS->var("sqrt_s_var");
    if (!sqrts_ws) sqrts_ws = dataWS->var("SqrtS"); // fallback if some WS uses this name

    if (!sqrts_ws) {
      std::cerr << "[FATAL] No sqrt(s) variable found (expected sqrt_s_var or SqrtS)\n";
      dataWS->Print();
      return 1;
    }
    
    std::cout << "------------------------------------------------------" << std::endl;
    std::cout << "[INFO] sqrt(s) var name = " << sqrts_ws->GetName()
	      << " val =" << sqrts_ws->getVal() << std::endl; ;
    std::cout << "------------------------------------------------------" << std::endl;
    sqrts_ws->setConstant(kTRUE);
    
    outputws->import(*intL);
  }

  vector<string> functionClasses;
  functionClasses.push_back("Exponential");
  functionClasses.push_back("Bernstein");
  functionClasses.push_back("Chebychev");
  functionClasses.push_back("Dijet");
  functionClasses.push_back("PowerLaw");
  map<string,string> namingMap;
  namingMap.insert(pair<string,string>("Exponential","exp"));
  namingMap.insert(pair<string,string>("Bernstein","bern"));
  namingMap.insert(pair<string,string>("Chebychev","cheb"));
  namingMap.insert(pair<string,string>("Dijet","dijet"));
  namingMap.insert(pair<string,string>("PowerLaw","plaw"));

  FILE *resFile ;
  resFile = fopen(Form("%s/fTestResults.txt",outDir.c_str()),"w");
  vector<map<string,int> > choices_vec;
  vector<map<string,std::vector<int> > > choices_envelope_vec;
  vector<map<string,RooAbsPdf*> > pdfs_vec;

  PdfModelBuilder pdfsModel;
  RooRealVar *mass = (RooRealVar*)dataWS->var(massVarName.c_str()); 
  std:: cout << "[INFO] Got mass from ws " << mass << std::endl;
  mass->Print("v");
  std:: cout << "[INFO] Mass range: [" << mass->getMin() << ", " << mass->getMax() << "]" << std::endl;
  if (mass && mass->hasMin() && !std::isnan(optMassMin)) mass->setMin(optMassMin);
  if (mass && mass->hasMax() && !std::isnan(optMassMax)) mass->setMax(optMassMax);

  const double mmin = mass->getMin();
  const double mmax = mass->getMax();

  pdfsModel.setObsVar(mass);
  double upperEnvThreshold = 0.1;   // Maximum probability for F-test (it was 1, lowered to 0.1)
  double minGofThreshold = 0.01;     // Minimum GoF considered to exclude bad functions (raised from 0.01 to 0.10)

  fprintf(resFile,"Truth Model & d.o.f & $\\Delta NLL_{N+1}$ & $p(\\chi^{2}>\\chi^{2}_{(N\\rightarrow N+1)})$ \\\n");
  fprintf(resFile,"\\hline\n");

  std::string ext = is2011 ? "7TeV" : "8TeV";
  if( year_ == "all" ){ ext = "13p6TeV"; }
  else{ ext = Form("%s_13p6TeV",year_.c_str()); }

  std::cout << "[INFO] Number of categories to process: " << ncats << std::endl;
  for (int cat=0; cat<ncats; cat++){
    map<string,int> choices;
    map<string,std::vector<int> > choices_envelope;
    map<string,RooAbsPdf*> pdfs; // Map is filled once per family
    map<string,RooAbsPdf*> allPdfs;
    map<string,int> familyMinOrder; // Tracking minimum order for the family in the envelope
    map<string,double> familyBestAIC; // Tracking best AIC for the family in the envelope
    map<string,RooAbsPdf*> familyBestPdf; // Tracking best PDF for the family in the envelope 
    string catname;

    catname = Form("cat%d",cat);

    // Parametric dataset
    // ------------------
    auto datasetNameForCat = [&](int c)->std::string{
      if (!dataName.empty()) return dataName;
      std::string s = dataPat; auto p = s.find("{CAT}"); if (p!=std::string::npos) s.replace(p,5,std::to_string(c));
      return s;
    };
    std::cout << "---------------------- CATEGORY IS = " << cat << std::endl;
    string data_name = datasetNameForCat(cat);
    RooDataHist  *data_in = (RooDataHist*)dataWS->data(data_name.c_str());
    if (!data_in) {
      cerr << "[ERROR] Could not find dataset '" << data_name << "' in workspace" << endl;
      continue;
    }

    // Get initial binning from data, will be updated if rebinning occurs
    // ------------------------------------------------------------------
    nBinsForFit = data_in->numEntries();          // RooDataHist bins
    //nBinsForFit = mass->getBinning().numBins();
    std::cout << "------------------------------------------------------" << std::endl;
    cout << "[INFO] Initial binning from workspace mass variable: " << nBinsForFit << " bins" << endl;
    std::cout << "entries dataset " << data_in->sumEntries() << std::endl;
    cout << "[DEBUG] Original RooDataHist info:" << endl;
    data_in->Print("v");

    RooAbsData *data;
    if (rebinFactor > 1) {
      cout << "[INFO] Rebinning histogram with factor " << rebinFactor << endl;
      RooRealVar* originalMass = (RooRealVar*)data_in->get()->first();
      double xmin = originalMass->getMin();
      double xmax = originalMass->getMax();
      int originalBins = originalMass->getBinning().numBins();
      int newBins = originalBins / rebinFactor;
      cout << "[DEBUG] Original bins: " << originalBins << ", requested new bins: " << newBins << endl;
      if (originalBins % rebinFactor != 0) {
        cout << "[WARNING] Original bins (" << originalBins << ") not divisible by rebinFactor (" << rebinFactor << ")" << endl;
        cout << "[INFO] Adjusting newBins to closest integer: " << newBins << endl;
      }
      RooRealVar* rebinnedMass = new RooRealVar(originalMass->GetName(), originalMass->GetTitle(), 
                                               xmin, xmax);
      rebinnedMass->setBins(newBins);
      TH1* originalHist = data_in->createHistogram("temp_hist", *originalMass);
      TH1F* rebinnedHist = new TH1F("rebinned_hist", "rebinned histogram", 
                                    newBins, xmin, xmax);
      double binWidth = (xmax - xmin) / newBins;
      for (int i = 1; i <= newBins; i++) {
        double newBinCenter = xmin + (i - 0.5) * binWidth;
        double binContent = 0;
        double binError2 = 0;
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
      RooDataHist* rebinnedData = new RooDataHist("rebinned_data", "rebinned data", 
                                                  RooArgSet(*rebinnedMass), rebinnedHist);
      mass = rebinnedMass;
      nBinsForFit = newBins;
      data = rebinnedData;
      cout << "[INFO] Manual rebinning completed: " << originalBins << " -> " << newBins << " bins" << endl;
      delete originalHist;
      delete rebinnedHist;
    } else {
      data = data_in;
    }
    // Final confirmation of binning used for fit
    // ------------------------------------------
    cout << "[INFO] Final binning used for fits: " << nBinsForFit << " bins" << endl;
    cout << "[DEBUG] Using data type: " << data->ClassName() << endl;
    std::cout << "entries dataset " << data->sumEntries() << std::endl;

    RooArgList storedPdfs("store");
    fprintf(resFile,"\\multicolumn{4}{|c|}{\\textbf{Category %d}} \\\n",cat);
    fprintf(resFile,"\\hline\n");

    int totalEnvelopeKeptThisCat = 0;
    double MinimimNLLSoFar=1e10;
    int simplebestFitPdfIndex = 0;
    
    for (vector<string>::iterator funcType=functionClasses.begin(); 
        funcType!=functionClasses.end(); funcType++){
      bool currentIncludeTurnOn = includeTurnOn;
      if ((*funcType == "PowerLaw" || *funcType == "PowerLawSingle") && includeTurnOn) {
        currentIncludeTurnOn = true;
        std::cout << "====> Enabling turn-on for " << *funcType << " (PowerLaw variant)" << std::endl;
      }
      std::cout << "======================================= " << std::endl;
      std::cout << "====> FAMILY " << funcType->c_str() << std::endl;
      std::cout << "======================================= " << std::endl;
      double thisNll=0.; double prevNll=0.; double chi2=0.; double prob=0.01; 
      double prev_chi2_red = 999.0; // Start with bad chi2 to ensure loop starts
      int order=1; int prev_order=0; int cache_order=0;
      RooAbsPdf *prev_pdf=NULL;
      RooAbsPdf *cache_pdf=NULL;
      std::vector<int> pdforders;
      std::cout << "===> F-TEST for Truth determination" << std::endl;

      int counter =0;
      std::cout << "[DEBUG] Initial values: prev_chi2_red=" << prev_chi2_red << ", order=" << order << ", gMaxFtestOrder=" << gMaxFtestOrder << std::endl;

      RooAbsPdf* best_pdf = nullptr;
      int best_order = 0;
      double best_score = 1e300;
      
      for (int order=1; order<=gMaxFtestOrder; ++order) {
	RooAbsPdf *bkgPdf = createBackgroundWithTurnOn(pdfsModel,*funcType,order,mass,Form("ftest_pdf_%d_%s_%s%d",(cat+catOffset),ext.c_str(), funcType->c_str(), order), currentIncludeTurnOn, dataWS, turnOnType);
	
	if(!bkgPdf) continue;
	
	int fitStatus = 0;
	FitQuality fq;
	double thisNll = 0.0;
	runFit(mass, bkgPdf, data, &thisNll, &fitStatus, 7, &fq);
	
	if (fitStatus!=0 || !std::isfinite(thisNll)) continue;
	
	const int nfloat = countFloatingParams(bkgPdf, data);
	const double score = 2.0*thisNll + nfloat;
	double gofP = -1;
	plot(mass, bkgPdf, data,
	     Form("%s/FTestFunctions/ftest_%s%d_%s", outDir.c_str(), funcType->c_str(), order, catname.c_str()),
	     flashggCats_, fitStatus, &gofP, fq, true, dataWS);
	
	std::cout << "[TRUTH] " << *funcType << " order=" << order
		  << " NLL=" << thisNll
		  << " nfloat=" << nfloat
		  << " score(2NLL+nfloat)=" << score
		  << " chi2/ndof=" << fq.chi2red
		  << " gofP=" << gofP
		  << "\n";
	if (gofP > 1e-6 && score < best_score) {
	  best_score = score;
	  best_order = order;
	  best_pdf   = bkgPdf;
	}
      }
            
      cache_order = best_order;
      cache_pdf   = best_pdf;

      choices.insert(pair<string,int>(*funcType,cache_order));
      if (cache_pdf != NULL) {
        pdfs.insert(pair<string,RooAbsPdf*>(Form("%s%d",funcType->c_str(),cache_order),cache_pdf));
      } else {
        std::cout << "[WARNING] cache_pdf is NULL for " << *funcType << ", skipping truth plot entry" << std::endl;
      }
      int truthOrder = cache_order;
      std::cout << "--------------------------------------------------------------------------------------------------" << std::endl;
      std::cout << "===> F-TEST and GOF for ENVELOPE determination" << std::endl;
      if (saveMultiPdf){
        chi2=0.;
        thisNll=0.;
        prevNll=0.;
        order=1;
        prev_order=0;
        cache_order=0;
        
        // Standard CMS envelope criteria:
        // 1. Maximum order: highest order with p^F < 0.1
        // 2. Lower orders: p^F < 0.1 AND GoF > 0.01
        
	std::cout << "--------------------------------------------------------------------------------------------------" << std::endl;
        std::cout << "[INFO] CMS Envelope Selection: Max order with p^F < 0.1, lower orders with p^F < 0.1 AND GoF > 0.01" <<std::endl;

	double prev_nll_env   = 0.0;
	int    prev_nfloat_env = 0;
	bool   have_prev = false;   // we only compare after we have a valid previous fit	
        int max_valid_order = 0;
        double final_f_test_prob = 1.0;
        double final_gofProb = 0.0;

	// Track how many orders we kept for THIS family
	int keptThisFamily = 0;
	
	while (order <= gMaxEnvelopeOrder) {
	  
	  std::cout << "[DEBUG-ENV] Top of loop: funcType=" << *funcType
		    << " order=" << order
		    << " gMaxEnvelopeOrder=" << gMaxEnvelopeOrder
		    << " keptFamily=" << keptThisFamily
		    << " keptTotalCat=" << totalEnvelopeKeptThisCat
		    << std::endl;
	  
	  if (totalEnvelopeKeptThisCat >= maxTotalEnvelope) {
	    std::cout << "[INFO] Reached maxTotalEnvelope=" << maxTotalEnvelope
		      << " for this category, stopping envelope building for remaining families." << std::endl;
	    break;
	  }
	  
	  if (keptThisFamily >= maxPerFamily) {
	    std::cout << "[INFO] Reached maxPerFamily=" << maxPerFamily
		      << " for family " << *funcType << ", stopping this family." << std::endl;
	    break;
	  }
	  
	  RooAbsPdf *bkgPdf = createBackgroundWithTurnOn(
							 pdfsModel, *funcType, order, mass,
							 Form("ftest_pdf_%d_%s_%s%d", (cat+catOffset), ext.c_str(), funcType->c_str(), order),
							 currentIncludeTurnOn, dataWS, turnOnType
							 );
	  
	  if (!bkgPdf) {
	    std::cout << "[WARNING] Could not create PDF for " << *funcType << " order " << order
		      << ", skipping." << std::endl;
	    order++;
	    continue;
	  }
	  
	  if (bkgPdf->InheritsFrom("RooAddPdf")) {
	    RooArgSet normSet(*mass);
	    ((RooAddPdf*)bkgPdf)->fixCoefNormalization(normSet);
	  }
	  
	  // Fit
	  int fitStatus = 0;
	  FitQuality fq;
	  runFit(mass, bkgPdf, data, &thisNll, &fitStatus, 7, &fq);
	  
	  if (fitStatus != 0 || !std::isfinite(thisNll)) {
	    std::cout << "[WARNING] Fit failed for " << bkgPdf->GetName()
		      << " status=" << fitStatus << " skipping this order" << std::endl;
	    order++;
	    continue;
	  }
	  
	  // Floating parameter count
	  const int nfloat_env = countFloatingParams(bkgPdf, data);
	  
	  // ---- F-test vs previous VALID order (correct bookkeeping)
	  double f_test_prob = 0.0;  // order=first valid: auto-accept
	  int dof = 0;
	  double f_stat = 0.0;
	  
	  if (!have_prev) {
	    have_prev = true;
	  } else {
	    f_stat = 2.0 * (prev_nll_env - thisNll);
	    dof    = nfloat_env - prev_nfloat_env;
	    f_test_prob = (f_stat > 0.0 && dof > 0) ? TMath::Prob(f_stat, dof) : 1.0;
	  }
	  
	  std::cout << "[DEBUG-ENV] order=" << order
		    << " nll=" << thisNll
		    << " nfloat=" << nfloat_env
		    << " prev_nfloat=" << prev_nfloat_env
		    << " dof=" << dof
		    << " fstat=" << f_stat
		    << " f_test_prob=" << f_test_prob
		    << " chi2red=" << fq.chi2red
		    << std::endl;
	  
	  // Update "previous" only once, AFTER computing the F-test
	  prev_nll_env    = thisNll;
	  prev_nfloat_env = nfloat_env;
	  
	  // Plot (and also compute chi2/p shown in plot header)
	  double gofProb_dummy = -1; // no longer used for selection
	  plot(mass, bkgPdf, data,
	       Form("%s/EnvelopeComponents_Pulls/envelope_%s%d_%s",
		    outDir.c_str(), funcType->c_str(), order, catnameDijet.c_str()),
	       flashggCats_, fitStatus, &gofProb_dummy, fq, true, dataWS);
	  
	  // ---- Selection: replace GoF(p) with chi2red cut
	  const bool passes_f_test = (!have_prev /*never true here*/ ? true : (f_test_prob < upperEnvThreshold));
	  const bool passes_chi2   = (fq.chi2red > 0.0 && fq.ndof > 0 && fq.chi2red < maxChi2red);
	  
	  // CMS-like stop condition: if improvement no longer significant, stop increasing order
	  // (only after we have a previous to compare)
	  if (have_prev && dof > 0 && !passes_f_test) {
	    std::cout << "[INFO] Order " << order << " fails F-test (p^F=" << f_test_prob
		      << " >= " << upperEnvThreshold << "), stopping this family." << std::endl;
	    break;
	  }
	  
	  // Decide inclusion
	  bool include_in_envelope = false;
	  if (!have_prev /*not possible*/ || order == 1) {
	    // keep order 1 if chi2 is sane
	    include_in_envelope = passes_chi2;
	  } else {
	    // for higher orders: require significant improvement + sane chi2
	    include_in_envelope = passes_f_test && passes_chi2;
	  }
	  
	  if (!passes_chi2) {
	    std::cout << "[INFO] Excluding " << *funcType << " order " << order
		      << " due to chi2/ndof=" << fq.chi2red
		      << " (cut < " << maxChi2red << ")" << std::endl;
	  }
	  
	  if (include_in_envelope) {
	    
	    // Parameter stability check (your existing logic; kept)
	    bool parametersStable = true;
	    
	    RooArgSet* params = bkgPdf->getParameters(*data);
	    std::unique_ptr<TIterator> iter(params->createIterator());
	    RooRealVar* param = nullptr;
	    
	    while ((param = (RooRealVar*)iter->Next())) {
	      if (!param->isConstant() && param->hasError()) {
		double relError = (param->getVal() != 0.0)
		  ? std::fabs(param->getError() / param->getVal())
		  : 999.0;
		
		bool atBoundary = (std::fabs(param->getVal() - param->getMin()) < 1e-6) ||
		  (std::fabs(param->getVal() - param->getMax()) < 1e-6);
		
		bool hugeError = (param->getError() > 50.0) || (relError > 50.0);
		
		if (atBoundary || hugeError) {
		  std::cout << "[WARNING] Parameter " << param->GetName()
			    << " potentially unstable: value=" << param->getVal()
			    << " error=" << param->getError()
			    << " relError=" << relError
			    << " atBoundary=" << atBoundary
			    << std::endl;
		  if (hugeError) parametersStable = false;
		}
	      }
	    }
	    delete params;
	    
	    if (parametersStable) {
	      // Keep it
	      const double myNll = 2.0 * thisNll;
	      std::cout << "[INFO] Adding to Envelope " << bkgPdf->GetName()
			<< " p^F=" << f_test_prob
			<< " chi2red=" << fq.chi2red
			<< " score(2NLL+nfloat)=" << (myNll + nfloat_env)
			<< std::endl;
	      
	      // Store ALL passing (family,order)
	      allPdfs.insert({Form("%s%d", funcType->c_str(), order), bkgPdf});
	      storedPdfs.add(*bkgPdf);
	      pdforders.push_back(order);
	      
	      keptThisFamily++;
	      totalEnvelopeKeptThisCat++;
	      
	      // Track "best overall" for reference / highlighting if you still use it
	      if ((myNll + nfloat_env) < MinimimNLLSoFar) {
		simplebestFitPdfIndex = storedPdfs.getSize() - 1;
		MinimimNLLSoFar = myNll + nfloat_env;
	      }
	      
	      // Keep per-family best (optional backward compatibility)
	      const string familyName = *funcType;
	      const double current_penalty = myNll + nfloat_env;
	      if (familyBestAIC.find(familyName) == familyBestAIC.end() || current_penalty < familyBestAIC[familyName]) {
		familyBestAIC[familyName] = current_penalty;
		familyBestPdf[familyName] = bkgPdf;
		familyMinOrder[familyName] = order;
	      }
	      
	      // bookkeeping for report
	      max_valid_order = order;
	      final_f_test_prob = f_test_prob;
	      final_gofProb = -1.0; // no longer meaningful for selection
	      
	    } else {
	      std::cout << "[INFO] Skipping " << *funcType << " order " << order
			<< " due to parameter instability" << std::endl;
	    }
	  }
	  
	  order++;
	}
	
        fprintf(resFile,"%15s & %d & %5.3f & %5.3f \\\n",funcType->c_str(),max_valid_order,final_f_test_prob,final_gofProb);
        choices_envelope.insert(pair<string,std::vector<int> >(*funcType,pdforders));
      }
    }

    fprintf(resFile,"\\hline\n");
    choices_vec.push_back(choices);
    choices_envelope_vec.push_back(choices_envelope);
    pdfs_vec.push_back(pdfs);
    plot(mass,pdfs,data,Form("%s/truths_%s",outDir.c_str(),catnameDijet.c_str()),flashggCats_,cat,catnameDijet,false);

    // Plot envelope-only functions (familyBestPdf contains one best function per family)
    if (saveMultiPdf && !allPdfs.empty()){
      int bestEnvOnlyIdx = bestIndexInPdfMap(mass, allPdfs, data);
      plot(mass, allPdfs, data,
	   Form("%s/EnvelopeComponents_Pulls/envelope_only_%s", outDir.c_str(), catnameDijet.c_str()),
	   flashggCats_, cat, catnameDijet, bestEnvOnlyIdx);
      std::cout << "[INFO] Generated envelope-only plot with " << allPdfs.size() << " functions" << std::endl;
    }
    
    //if (saveMultiPdf && !familyBestPdf.empty()){
    //  int bestEnvOnlyIdx = bestIndexInPdfMap(mass, familyBestPdf, data);
    //  plot(mass,familyBestPdf,data,Form("%s/EnvelopeComponents_Pulls/envelope_only_%s",outDir.c_str(),catnameDijet.c_str()),flashggCats_,cat,catnameDijet,bestEnvOnlyIdx);
    //  std::cout << "[INFO] Generated envelope-only plot with " << familyBestPdf.size() << " functions" << std::endl;
    //}

    if (saveMultiPdf){
      string catindexname;
      string catname2;
      catindexname = Form("pdfindex_%d_%s",(cat+catOffset),ext.c_str());
      catname2 = Form("cat%d",(cat+catOffset));

      RooCategory catIndex(catindexname.c_str(),"c");
      
      // Build final PDF list using only one function per family
      RooArgList finalStoredPdfs("final_store");
      map<string,RooAbsPdf*> finalPdfs;
      
      // Use envelope best functions (familyBestPdf) instead of F-test functions
      std::cout << "[INFO] =======================================================" << std::endl;
      std::cout << "[INFO] Building MultiPdf workspace from ENVELOPE best functions" << std::endl;
      std::cout << "[INFO] Available envelope functions in 'familyBestPdf' map: " << familyBestPdf.size() << std::endl;
      for (const auto& p : familyBestPdf) {
        std::cout << "[INFO]   - " << p.first << " -> " << (p.second ? p.second->GetName() : "NULL") << std::endl;
      }
      std::cout << "[INFO] =======================================================" << std::endl;

      for (const auto& kv : allPdfs) {
	const std::string& pdfKey = kv.first;   // e.g. Exponential3
	RooAbsPdf* pdf = kv.second;
	if (!pdf) continue;
	
	// Optional extra stability gate at final stage (keep consistent with envelope)
	// If you want strict consistency, you can skip this.
	// (Here: keep everything from envelope loop)
	finalPdfs[pdfKey] = pdf;
	finalStoredPdfs.add(*pdf);
	
	std::cout << "[INFO] Adding to final MultiPdf: " << pdfKey
		  << " (PDF: " << pdf->GetName() << ")" << std::endl;
      }
      
      if (finalStoredPdfs.getSize() == 0) {
	std::cout << "[ERROR] No PDFs found from envelope for category " << cat << "! Cannot create MultiPdf." << std::endl;
      }
      
      if (finalStoredPdfs.getSize() == 0) {
        std::cout << "[ERROR] No PDFs found from envelope for category " << cat << "! Cannot create MultiPdf." << std::endl;
        // Fallback: try to use at least one function from allPdfs map
        if (!allPdfs.empty()) {
          std::cout << "[INFO] Attempting to add first available PDF from allPdfs as fallback..." << std::endl;
          for (const auto& fallbackPdf : allPdfs) {
            string pdfKey = fallbackPdf.first;
            RooAbsPdf* pdf = fallbackPdf.second;
            if (pdf) {
              finalPdfs[pdfKey] = pdf;
              finalStoredPdfs.add(*pdf);
              std::cout << "[INFO] Added fallback PDF: " << pdfKey << std::endl;
              break;
            }
          }
        }
      }
      
      // No external normalization variable needed - each PDF is self-normalized
      std::cout << "[DEBUG] Creating RooMultiPdf with " << finalStoredPdfs.getSize() << " PDFs..." << std::endl;
      RooMultiPdf *pdf = new RooMultiPdf(Form("CMS_hgg_%s_%s_bkgshape",catname2.c_str(),ext.c_str()),"all pdfs",catIndex,finalStoredPdfs);
      std::cout << "[DEBUG] RooMultiPdf created successfully, contains " << pdf->getNumPdfs() << " PDFs" << std::endl;

      // Choose best PDF by actually fitting all candidates (same strategy used elsewhere)
      // This sets catIndex to the best one and restores its parameters.
      int bestFitPdfIndex = getBestFitFunction(mass, pdf, data, &catIndex, outputws, Form("bestfit_cat%d_%s",cat,ext.c_str()), false);
      
      catIndex.setIndex(bestFitPdfIndex);
      std::cout << "// ------------------------------------------------------------------------- //" <<std::endl; 
      std::cout << "[INFO] *** MULTIPDF WORKSPACE CONTENTS (F-test selected functions) ***" << std::endl;
      std::cout << "[INFO] Created MultiPdf " << pdf->GetName() << ", in Category " << cat << " with a total of " << catIndex.numTypes() << " pdfs"<< std::endl;
      std::cout << "[INFO] Using PDF index = " << bestFitPdfIndex << " (no refit applied to preserve stability)" << std::endl;
      std::cout << "[INFO] PDFs in workspace:" << std::endl;
      finalStoredPdfs.Print();
      std::cout << "[INFO] Default selected Pdf = " << bestFitPdfIndex << ", " << finalStoredPdfs.at(bestFitPdfIndex)->GetName() << std::endl;
      std::cout << "[INFO] These are the TRUTH functions from F-test (also in truths_*.pdf plot)" << std::endl;
      std::cout << "// ------------------------------------------------------------------------- //" <<std::endl;
      outputws->import(*pdf);
      outputws->import(catIndex);
      outputws->import(*data);

      // Plot with correct names using finalPdfs
      // ---------------------------------------
      //plot(mass,pdf,&catIndex,data,Form("%s/multipdf_cat%s",outDir.c_str(),catnameDijet.c_str()),flashggCats_,cat,bestFitPdfIndex);
      plot(mass,finalPdfs,data,Form("%s/multipdf_cat%s",outDir.c_str(),catnameDijet.c_str()),flashggCats_,cat,catnameDijet,bestFitPdfIndex);

    } // End saveMultiPdf block
  }  // End category loop
  
  if (saveMultiPdf){
    outputfile->cd();
    outputws->Write();
    outputfile->Close();  
  }
  if (inFile) inFile->Close();
  
  return 0;
}
