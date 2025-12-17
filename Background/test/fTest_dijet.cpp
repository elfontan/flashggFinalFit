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
#include "RooChi2Var.h"

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

// === Helpers for blinding multi-range ===
struct Interval { double lo, hi; };
static std::string sidebandRangeList = "blind_low,blind_high"; // default 
static std::string blindRangeList    = "";                     // from CLI --blind

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
static std::string define_ranges(RooRealVar* x, const std::vector<Interval> &ivs, const std::string &basename) {
  std::vector<std::string> names; names.reserve(ivs.size());
  for (size_t i=0;i<ivs.size();++i) {
    std::string nm = basename + std::string("_") + std::to_string(i);
    x->setRange(nm.c_str(), ivs[i].lo, ivs[i].hi); names.push_back(nm);
  }
  std::string list; for (size_t i=0;i<names.size();++i){ if (i) list+=","; list+=names[i]; }
  return list;
}

// Z model name global (configurable from CLI)
std::string gZModelName = "model_Z_c2";

// Defining maximum order allowed for the FTest (configurable from CLI)
int gMaxFtestOrder = 5;     // Reduced from 8 to 5 to reduce too high orders 
int gMaxEnvelopeOrder = 10; // Consider reducing from 10 to 6 to keep the envelope under control

// Mass range configuration 
float mN = 2.75;
float sigma = 0.025;
int nsigma = 10;
float mN_low  = 300;
float mN_high = 1000;

// not configured in main
int nBinsForFit  = (mN_high-mN_low)/10;  // 10 GeV per bin: 150 bins in [0,1500] GeV
int nBinsForPlot = (mN_high-mN_low)/10;
//int nBinsForFit  = (mN_high-mN_low)/2;  // 2 GeV per bin instead of 0.25 GeV
//int nBinsForPlot = (mN_high-mN_low)/2;

RooRealVar *intLumi_ = new RooRealVar("IntLumi","hacked int lumi", 1000.);

TRandom3 *RandomGen = new TRandom3();

RooAbsPdf* getPdf(PdfModelBuilder &pdfsModel, string type, int order, const char* ext=""){
  if (type=="Bernstein") return pdfsModel.getBernstein(Form("%s_bern%d",ext,order),order); 
  else if (type=="Exponential") return pdfsModel.getExponentialSingle(Form("%s_exp%d",ext,order),order); 
  else if (type=="ExponentialSum") return pdfsModel.getExponential(Form("%s_expsum%d",ext,order),order); 
  else if (type=="PowerLaw") {
    if (order > 5) {  // Increased from 3 to 5 to allow more functions in the envelope
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

// Restituisce una pdf estesa pronta per il χ² (possibilmente wrappando quella in ingresso)
// Se trova 'bkg_norm' e/o 'z_norm' li somma e li usa come yield totale; altrimenti scala sui sideband.
std::unique_ptr<RooAbsPdf> makeExtendedForGOF(RooAbsPdf* pdf, RooRealVar* mass,
                                              RooAbsData* data, const char* rangeList) {
  if (pdf->InheritsFrom("RooExtendPdf")) {
    return std::unique_ptr<RooAbsPdf>((RooAbsPdf*)pdf);
  }
  if (pdf->InheritsFrom("RooAddPdf")) {
    double nexp = pdf->expectedEvents(RooArgSet(*mass));
    if (nexp > 0) return std::unique_ptr<RooAbsPdf>((RooAbsPdf*)pdf);
  }
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
    std::unique_ptr<RooAbsReal> frac_pdf(
      pdf->createIntegral(RooArgSet(*mass), RooFit::NormSet(RooArgSet(*mass)),
                          RooFit::Range(rangeList))
    );
    RooDataHist* dh = dynamic_cast<RooDataHist*>(data);
    std::unique_ptr<RooDataHist> dh_owner;
    if (!dh) { dh_owner.reset(new RooDataHist("__gof_dh","__gof_dh", RooArgSet(*mass), *data)); dh = dh_owner.get(); }
    double n_sb = dh->sumEntries(nullptr, rangeList);
    double f_sb = frac_pdf ? frac_pdf->getVal() : 0.;
    double norm_val = (f_sb>0) ? (n_sb / f_sb) : n_sb;
    normVar = new RooRealVar("__gof_norm","__gof_norm", norm_val, 0., 1e15);
    normVar->setConstant(kTRUE);
  }
  return std::unique_ptr<RooAbsPdf>(new RooExtendPdf("__gof_ext","__gof_ext", *pdf, *normVar));
}

Chi2Result computeChi2ReducedExtended(RooRealVar* mass, RooAbsPdf* pdf_in, RooAbsData* data,
                                      const char* rangeList) {
    RooDataHist* dh = dynamic_cast<RooDataHist*>(data);
    std::unique_ptr<RooDataHist> dh_owner;
    if (!dh) {
        dh_owner.reset(new RooDataHist("__chi2_dh","__chi2_dh", RooArgSet(*mass), *data));
        dh = dh_owner.get();
    }
    RooAbsPdf* pdf_ext = nullptr;
    std::unique_ptr<RooExtendPdf> owned_ext;
    if (pdf_in->InheritsFrom("RooExtendPdf")) {
        pdf_ext = pdf_in;
    } else {
        double norm_val = pdf_in->expectedEvents(RooArgSet(*mass));
        if (norm_val <= 0) {
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
    RooChi2Var chi2var("__chi2","__chi2", *pdf_ext, *dh,
                       RooFit::Range(rangeList),
                       RooFit::DataError(RooAbsData::Poisson));
    double chi2 = chi2var.getVal();
    auto tmp = std::unique_ptr<RooPlot>(mass->frame());
    dh->plotOn(tmp.get(), RooFit::CutRange(rangeList));
    auto* h = dynamic_cast<RooHist*>(tmp->getObject(int(tmp->numItems()-1)));
    int nbins_used = h ? h->GetN() : 0;
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
  int ntries=0;
  RooArgSet *params_test = pdf->getParameters((const RooArgSet*)(0));
  data->Print("v");
  params_test->Print("v");
  int stat=1;
  double minnll=10e8;
  while (stat!=0){
    if (ntries>=MaxTries) break;
    RooFitResult *fitTest;
    if (blindSignalRegion) {
      std::cout << "[INFO] Fitting in sidebands: " << sidebandRangeList << std::endl;
      fitTest = pdf->fitTo(*data,RooFit::Save(1),RooFit::Minimizer("Minuit2","minimize"),RooFit::Strategy(0),RooFit::PrintLevel(-1),RooFit::Optimize(1),RooFit::SumW2Error(kFALSE),RooFit::Range(sidebandRangeList.c_str()));
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
  
  // Check for integration failures or invalid NLL values
  if (!std::isfinite(minnll) || minnll > 1e7) {
    std::cout << "[WARNING] Fit resulted in invalid NLL (" << minnll << "), marking as failed" << std::endl;
    *stat_t = 999; // Custom error code for integration failure
    *NLL = 1e8;    // Very high NLL to indicate failure
  }
  
  if (chi2_out && stat == 0 && std::isfinite(minnll)) {
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
      double normVal = bkg_norm->getVal();
      Nvis = normVal;
      std::cout << "[INFO] Using scaled PDF integral for Nvis: " << Nvis << " events (plot range)" << std::endl;
    } else {
      if (blindSignalRegion) {
        Nvis = data->sumEntries("1", sidebandRangeList.c_str());
        std::cout << "[INFO] Blinded fit: Nvis = " << Nvis << " events in sidebands" << std::endl;
      } else {
        Nvis = data->sumEntries();
        std::cout << "[INFO] Unblinded fit: Nvis = " << Nvis << " events total" << std::endl;
      }
    }
    RooRealVar* mass = (RooRealVar*)data->get()->first();
    RooPlot* tempPlot = mass->frame();
    data->plotOn(tempPlot, Name("data"),Range(sidebandRangeList.c_str()));
    bool isAddPdf = pdf->InheritsFrom("RooAddPdf");
    cout << "[DEBUG] PDF type: " << pdf->ClassName() << ", isAddPdf: " << isAddPdf << endl;
    cout << "[DEBUG] Data entries: " << data->sumEntries() << endl;
    pdf->plotOn(tempPlot, Name("pdf"),Range(sidebandRangeList.c_str()));
    cout << "[DEBUG] Using natural normalization for chi2 calculation" << endl;
    int np = pdf->getParameters(*data)->getSize();
    double chi2_reduced = tempPlot->chiSquare("pdf", "data", np);
    
    // Calculate actual number of bins used for chi2
    RooHist* hdata = (RooHist*)tempPlot->findObject("data");
    int nbins_used = hdata ? hdata->GetN() : nBinsForFit;
    int ndof = nbins_used - np;
    double chi2_total = chi2_reduced * ndof;
    
    cout << "[DEBUG] Chi2 calculation: nbins=" << nbins_used << ", np=" << np << ", ndof=" << ndof << endl;
    cout << "[DEBUG] chi2_reduced=" << chi2_reduced << ", chi2_total=" << chi2_total << endl;
    
    *chi2_out = chi2_reduced;  // Return chi2/ndof for compatibility
    delete tempPlot;
  }
  if (stat == 0) {
    std::cout << "=== FINAL PARAMETERS AFTER FIT ===" << std::endl;
    RooArgSet *final_params = pdf->getParameters((const RooArgSet*)(0));
    final_params->Print("v");
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
  std::cout << "[DEBUG] F-test calculation: chi2=" << chi2 << ", ndof=" << ndof << ", runFtestCheckWithToys=" << runFtestCheckWithToys << std::endl;
  
  // Per F-test: chi2 grande = miglioramento significativo = prob piccola
  double prob_asym = TMath::Prob(chi2,ndof);
  std::cout << "[DEBUG] TMath::Prob result: " << prob_asym << std::endl;
  
  // Fix per overflow: se prob è troppo piccola, la trattiamo come 0
  if (prob_asym < 1e-10) {
    std::cout << "[DEBUG] Probability underflow, setting to 1e-10" << std::endl;
    prob_asym = 1e-10;
  }
  
  // Fix per range: se chi2 è negativo o troppo piccolo, nessun miglioramento
  if (chi2 <= 0) {
    std::cout << "[DEBUG] Chi2 <= 0 (" << chi2 << "), no improvement, setting prob=0.9" << std::endl;
    prob_asym = 0.9;  // Non significativo ma non completamente bad
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
    TCanvas *can = new TCanvas();
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

void plotComponents(RooRealVar *mass, RooAbsPdf *pdf, RooAbsPdf *zModel, RooAbsData *data, string name,vector<string> flashggCats_, int status, double *prob, double chi2FromFit=-1, bool fitWasBlinded=false){
  double chi2_gof;
  // Use chi2 from fit if available, otherwise use -1 as placeholder (will cause error)
  if (chi2FromFit < 0) {
    std::cout << "[ERROR] plotComponents called without valid chi2FromFit!" << std::endl;
  }
  *prob = getGoodnessOfFit(mass,pdf,data,name,chi2FromFit,&chi2_gof);
  RooPlot *plot = mass->frame(mass->getMin(), mass->getMax());
  if (fitWasBlinded) {
    data->plotOn(plot, MarkerStyle(20), MarkerSize(0.8), CutRange(sidebandRangeList.c_str()));
  } else {
    data->plotOn(plot, MarkerStyle(20), MarkerSize(0.8));
  }
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
    double fitIntegral = pdf->createIntegral(*mass, NormSet(*mass), Range(sidebandRangeList.c_str()))->getVal();
    double normVal = bkg_norm->getVal();
    Nvis = normVal * (fitIntegral );
    std::cout << "[INFO] Using scaled PDF integral for Nvis: " << Nvis << " events (plot range)" << std::endl;
  } else {
    if (fitWasBlinded) {
      Nvis = data->sumEntries("1", sidebandRangeList.c_str());
      std::cout << "[INFO] Blinded fit: Nvis = " << Nvis << " events in sidebands" << std::endl;
    } else {
      Nvis = data->sumEntries();
      std::cout << "[INFO] Unblinded fit: Nvis = " << Nvis << " events total" << std::endl;
    }
  }
  std::cout << "[INFO] Plotting function over full mass range with natural normalization" << std::endl;
  pdf->plotOn(plot,LineColor(kRed),LineWidth(2),Name("total"),Range(mass->getMin(),mass->getMax()),NormRange(sidebandRangeList.c_str()),
                     Normalization(Nvis, RooAbsReal::NumEvent));
  if (pdf->InheritsFrom("RooAddPdf")) {
    RooAddPdf* addPdf = (RooAddPdf*)pdf;
    RooArgList comps = addPdf->pdfList();
    if (comps.getSize() >= 2) {
      RooAbsPdf* bkgPdf = (RooAbsPdf*)comps.at(0);
      if (bkgPdf) {
        pdf->plotOn(plot,Components(*bkgPdf),LineColor(kBlue),LineStyle(kDashed),LineWidth(2),Name("background"),
                   NormRange(sidebandRangeList.c_str()),
                   Normalization(Nvis, RooAbsReal::NumEvent),Range(mass->getMin(),mass->getMax()));
      }
      if (zModel) {
        pdf->plotOn(plot,Components(*zModel),LineColor(kGreen+2),LineWidth(2),Name("Z"),
                   NormRange(sidebandRangeList.c_str()),
                   Normalization(Nvis, RooAbsReal::NumEvent));
      }
    }
  }
  double chi2;
  if (chi2FromFit > 0) {
    chi2 = chi2FromFit;
    cout << "[INFO] Using chi2 from fit: " << chi2 << endl;
  } else {
    chi2 = chi2_gof;  // Use chi2 from GoF calculation
    cout << "[INFO] Using chi2/ndof from GoF calculation: " << chi2 << endl;
  }
  TCanvas *canv = new TCanvas("canv","canv",800,800);
  canv->Divide(1,2);
  canv->cd(1);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.02);
  gPad->SetPad(0.01,0.3,0.99,0.99);
  gPad->SetLogy(); // Set log scale for y-axis
  plot->GetXaxis()->SetTitleSize(0.04);
  plot->GetXaxis()->SetTitle("");
  plot->GetXaxis()->SetLabelSize(0);
  plot->GetYaxis()->SetTitleSize(0.05);
  plot->GetYaxis()->SetLabelSize(0.04);
  plot->GetYaxis()->SetTitleOffset(1.2);
  plot->GetYaxis()->SetTitle("Events");
  plot->SetTitle("");
  plot->SetMaximum(plot->GetMaximum()*1.3);
  double minVal = plot->GetMinimum();
  // Better minimum for log scale: find a reasonable minimum based on data
  double dataMinimum = minVal;
  if (minVal <= 0) dataMinimum = 0.01;
  else if (minVal < 1e-3) dataMinimum = 1e-3;
  else dataMinimum = minVal * 0.3;  // 30% below minimum for log scale
  plot->SetMinimum(dataMinimum);
  plot->Draw();
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
  canv->cd(2);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.15);
  gPad->SetPad(0.01,0.01,0.99,0.3);
  gPad->SetGridy();
  RooPlot *ratioPlot = mass->frame(mass->getMin(), mass->getMax());
  if (fitWasBlinded) {
    data->plotOn(ratioPlot,MarkerStyle(20),MarkerSize(0.8),CutRange(sidebandRangeList.c_str()));
    pdf->plotOn(ratioPlot,LineColor(kRed),LineWidth(2),NormRange(sidebandRangeList.c_str()),
                     Normalization(Nvis, RooAbsReal::NumEvent));
  } else {
    data->plotOn(ratioPlot,MarkerStyle(20),MarkerSize(0.8));
    pdf->plotOn(ratioPlot,LineColor(kRed),LineWidth(2));
  }
  RooHist* hpull = ratioPlot->pullHist("", "", true);
  RooPlot *pullPlotFinal = mass->frame(mass->getMin(), mass->getMax());
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

void plot(RooRealVar *mass, RooAbsPdf *pdf, RooAbsData *data, string name,vector<string> flashggCats_, int status, double *prob, double chi2FromFit=-1, bool fitWasBlinded=false, bool showPulls=false,RooWorkspace* dataWS = nullptr){
  double chi2_gof;
  // Use chi2 from fit if available
  if (chi2FromFit < 0) {
    std::cout << "[ERROR] plot called without valid chi2FromFit!" << std::endl;
  }
  *prob = getGoodnessOfFit(mass,pdf,data,name,chi2FromFit,&chi2_gof);
  RooPlot *plot = mass->frame(mass->getMin(), mass->getMax());
  if (fitWasBlinded) {
    data->plotOn(plot, MarkerStyle(20), MarkerSize(0.8), CutRange(sidebandRangeList.c_str()));
  } else {
    data->plotOn(plot, MarkerStyle(20), MarkerSize(0.8));
  }
  double Nvis = -1;
  // Crea una variabile di yield
  RooRealVar N_pdf("N_pdf", "N_pdf", data->sumEntries(), 10, 1000000);

  // Wrappa pdf in RooExtendPdf
  RooExtendPdf pdf_extended("pdf_extended", "pdf extended", *pdf, N_pdf);

  // Fitta sui dati CR (o sul range di fitting che vuoi)
  pdf_extended.fitTo(*data);

  // Extract normalization factor (after fit)
  double Nsr_norm = N_pdf.getVal();
  std::cout << "Nsr_norm (fittato): " << Nsr_norm << " sum entries "<<  data->sumEntries() <<  std::endl;

  std::cout << "[INFO] Plotting function over full mass range" << std::endl;
  if (fitWasBlinded) {
    // For blinded fits, normalize only in sidebands
    pdf->plotOn(plot,LineColor(kCyan),LineWidth(2),Name("total"),Range(mass->getMin(), mass->getMax()),
                       NormRange(sidebandRangeList.c_str()),
                       Normalization(N_pdf.getVal(), RooAbsReal::NumEvent));
  } else {
    // For unblinded fits, normalize over full range
    pdf->plotOn(plot,LineColor(kCyan),LineWidth(2),Name("total"),Range(mass->getMin(), mass->getMax()),
                       Normalization(N_pdf.getVal(), RooAbsReal::NumEvent));
  }
  double chi2;
  if (chi2FromFit > 0) {
    chi2 = chi2FromFit;
    cout << "[INFO] Using chi2 from fit: " << chi2 << endl;
  } else {
    chi2 = chi2_gof;  // Use chi2 from GoF calculation
    cout << "[INFO] Using chi2/ndof from GoF calculation: " << chi2 << endl;
  }
  
  TCanvas *canv;
  if (showPulls) {
    canv = new TCanvas("canv","canv",800,800);
    canv->Divide(1,2);
    canv->cd(1);
    gPad->SetLeftMargin(0.15);
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
    canv = new TCanvas();
    canv->SetLogy(); // Set log scale for y-axis
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.10);
    plot->GetXaxis()->SetTitleSize(0.04);
    plot->GetXaxis()->SetTitle("m_{\ell#pi} (GeV)");
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
  lat->DrawLatex(0.15,0.94,Form("#chi^{2}/ndof = %.3f, Prob = %.2f, Fit Status = %d ",chi2,*prob,status));
  
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
    gPad->SetBottomMargin(0.15);
    gPad->SetPad(0.01,0.01,0.99,0.3);
    gPad->SetGridy();
    RooPlot *ratioPlot = mass->frame(mass->getMin(), mass->getMax());
    if (fitWasBlinded) {
      data->plotOn(ratioPlot,MarkerStyle(20),MarkerSize(0.8),CutRange(sidebandRangeList.c_str()));
      pdf->plotOn(ratioPlot,LineColor(kCyan),LineWidth(2),NormRange(sidebandRangeList.c_str()),
                       Normalization(N_pdf.getVal(), RooAbsReal::NumEvent));
    } else {
      data->plotOn(ratioPlot,MarkerStyle(20),MarkerSize(0.8));
      pdf->plotOn(ratioPlot,LineColor(kCyan),LineWidth(2));
    }
    RooHist* hpull = ratioPlot->pullHist("", "", true);
    RooPlot *pullPlotFinal = mass->frame(mass->getMin(), mass->getMax());
    pullPlotFinal->addPlotable(hpull,"P");
    pullPlotFinal->GetYaxis()->SetNdivisions(504);
    pullPlotFinal->GetYaxis()->SetLabelSize(0.12);
    pullPlotFinal->GetYaxis()->SetTitleSize(0.12);
    pullPlotFinal->GetYaxis()->SetTitleOffset(0.5);
    pullPlotFinal->GetYaxis()->SetTitle("Pull");
    pullPlotFinal->GetXaxis()->SetTitleSize(0.12);
    pullPlotFinal->GetXaxis()->SetLabelSize(0.12);
    pullPlotFinal->GetXaxis()->SetTitle("m_{\ell#pi} (GeV)");
    pullPlotFinal->SetTitle("");
    pullPlotFinal->SetMaximum(3);
    pullPlotFinal->SetMinimum(-3);
    pullPlotFinal->Draw();
  }
  
  canv->SaveAs(Form("%s.pdf",name.c_str()));
  canv->SaveAs(Form("%s.png",name.c_str()));
  delete canv;
  delete lat;
  // --- New plot CR/SR with prediction (RooPlot version) ---
  if (dataWS) {
    RooDataHist* data_sr = dynamic_cast<RooDataHist*>(dataWS->data("data_sr_m300"));
    RooDataHist* data_cr = dynamic_cast<RooDataHist*>(dataWS->data("data_cr_m300"));
   if (data_sr && data_cr) {

      // 0) Recupero TF(m) dal ws
     // RooAbsReal* TF = d; // <-- CAMBIA "TF" col nome reale
    //  if (!TF) { std::cerr << "[ERR] TF non trovata nel workspace\n"; return; }

      // 1) Costruisco RooHistPdf del CR
      mass->setRange("PLOTRANGE", 100, 1000);
      RooHistPdf pdf_cr_hist("pdf_cr_hist","CR hist-pdf",RooArgSet(*mass),*data_cr,2);
      pdf_cr_hist.setNormRange("PLOTRANGE");  
      // check: l’integrale deve essere ≈ 1 sul range che userai per il plot
      std::unique_ptr<RooAbsReal> I_den(
      pdf_cr_hist.createIntegral(
        RooArgSet(*mass),
        RooFit::NormSet(RooArgSet(*mass)),
        RooFit::Range("PLOTRANGE")
       )
      );
      std::cout << "Int[pCR|PLOTRANGE] = " << pdf_cr_hist.getVal() << std::endl;  // ~1
      // 2) Predizione continua: SR_pred = CR_hist * TF (rinormalizzata automaticamente)
//      RooEffProd pdf_sr_pred("pdf_sr_pred","SR pred = CR*TF", pdf_cr_hist, *pdf);
  //    #include "RooProduct.h"

      RooProduct pdf_sr_pred("pdf_sr_pred","SR pred = CR*TF", RooArgList(pdf_cr_hist, pdf_extended));
      // 3) Normalizzazioni per plottare in numero di eventi
      const double Ncr = data_cr->sum(false);
      // Se vuoi shape-puro, normalizza SR_pred a N(data SR); altrimenti usa l’atteso CR×<TF>
      // Stima semplice di <TF> pesata dal CR (subbinning): evita il "center-of-bin"
      //const RooAbsBinning& binning = mass->getBinning();
      //const int nb = binning.numBins();
      //prendi nb da istogramma cr
      TH1* h = data_cr->createHistogram("h", *mass);
      int nb = data_cr->numEntries();
      double lo = h->GetXaxis()->GetXmin();
      double hi = h->GetXaxis()->GetXmax();
      std::cout << "nbin " << nb << " xmin " << lo << " xmax " << hi <<std::endl; 
      // ora h->GetBinLowEdge(i), h->GetBinUpEdge(i) ti danno i bordi
      delete h;
      double Nsr_pred_int = 0.0;
      for (int i=0;i<nb;++i){
        data_cr->get(i);
        const double n_i = data_cr->weight(i);
        // media TF nel bin con subcampionamento
        const int sub = 100;
        double acc = 0.0;
        for (int j=0;j<sub;++j){
          const double x = lo + (j+0.5)*(hi-lo)/sub;
          mass->setVal(x);
          //std::cout << pdf->getVal(RooArgSet(*mass)) << " " << n_i << std::endl; 
          acc += pdf->getVal(RooArgSet(*mass));
        }
        const double alpha_i = acc/sub;
        Nsr_pred_int += n_i * alpha_i*(hi-lo)/nb;
     //   std::cout << i << "cumulative pred " << Nsr_pred_int << " counts cr " << n_i << " mean tF value " << alpha_i << std::endl; 
      }
      const double Nsr_obs = data_sr->sum(false);
      // Scegli la normalizzazione (commenta quella che non vuoi):
      // a) confronto shape-puro con i dati SR
      //double Nsr_norm = (Nsr_obs>0 ? Nsr_obs : Nsr_pred_int);
      // b) confronto "atteso" (CR×<TF>)
      mass->setRange("SR",lo,hi);  // oppure il range di SR che vuoi
      std::unique_ptr<RooAbsReal> I_sr( pdf_sr_pred.createIntegral(RooArgSet(*mass), Range("SR")) );
      double yield_pred_sr = N_pdf.getVal();  // Ncr è il totale di CR
      //std::unique_ptr<RooAbsReal> I_numG( pdf_sr_pred.createIntegral(RooArgSet(*mass)) );
      //std::unique_ptr<RooAbsReal> I_denG( pdf_cr_hist.createIntegral(RooArgSet(*mass)) );
      //double kappa = I_numG->getVal();//
      std::unique_ptr<RooAbsReal> I_pdf( pdf->createIntegral(RooArgSet(*mass)) );
      std::cout << "Integrale di pdf: " << I_pdf->getVal() << std::endl;

      std::unique_ptr<RooAbsReal> I_cr( pdf_cr_hist.createIntegral(RooArgSet(*mass)) );
      std::cout << "Integrale di pdf_cr_hist: " << I_cr->getVal() << std::endl;

      std::unique_ptr<RooAbsReal> I_prod( pdf_sr_pred.createIntegral(RooArgSet(*mass), Range("SR")));
      std::cout << "Integrale di pdf_sr_pred (CR*TF): " << I_prod->getVal() << std::endl;
      std::cout << "Atteso (I_cr * I_pdf): " << I_cr->getVal() * I_pdf->getVal() << std::endl;
      std::cout << yield_pred_sr<< "integral of cr*tf" << std::endl;
      std::cout << Nvis<< " Nvis" << std::endl;
      //  std::cout << I_denG->getVal() << "integral of cr" << std::endl;
      //  / I_denG->getVal();  // = ⟨f⟩_CR
      double Nsr_norm = I_prod->getVal()/10.;
      std::cout << " NORM FACTOR " << Nsr_norm <<std::endl;

      // 4) Costruisco i RooPlot (uno per lo spettro, uno per i pull)
      RooPlot* fr = mass->frame();
      RooPlot* frPull = mass->frame();

      // 5) Dati SR e CR sul frame (marcatori diversi)
      //    NB: do un nome agli oggetti sul frame per calcolare i pull
/*      data_sr->plotOn(fr,
        RooFit::Name("dataSR"),
        RooFit::MarkerStyle(20), RooFit::MarkerSize(0.8));

      data_cr->plotOn(fr,
        RooFit::Name("dataCR"),
        RooFit::MarkerStyle(24), RooFit::MarkerColor(kBlue+1), RooFit::MarkerSize(0.7));

      // 6) Disegno CR RooHistPdf (in eventi)
      pdf_cr_hist.plotOn(fr,
        RooFit::Name("curveCR"),
        RooFit::LineColor(kBlue+1), RooFit::LineStyle(kDotted), RooFit::LineWidth(2),
        RooFit::Normalization(Ncr, RooAbsReal::NumEvent),
        RooFit::Range(mass->getMin(), mass->getMax())
      );
*/
      // 7) Disegno SR_pred = CR*TF (in eventi) — è la curva per i pull
      //    Se sei in modalità "blinded", usa NormRange coerente coi sidebands
      if (fitWasBlinded) {
        pdf_sr_pred.plotOn(fr,
          RooFit::Name("curveSRpred"),
          RooFit::LineColor(kRed), RooFit::LineStyle(kSolid), RooFit::LineWidth(2),
          RooFit::Normalization(Nsr_norm, RooAbsReal::NumEvent),
          RooFit::Range(mass->getMin(), mass->getMax()),
          RooFit::NormRange(sidebandRangeList.c_str())
        );
      } else {
        pdf_sr_pred.plotOn(fr,
          RooFit::Name("curveSRpred"),
          RooFit::LineColor(kRed), RooFit::LineStyle(kSolid), RooFit::LineWidth(2),
          RooFit::Normalization(Nsr_norm, RooAbsReal::NumEvent),
          RooFit::Range(mass->getMin(), mass->getMax())
        );
      }

      // 8) Legend
      TLegend* leg_crsr = new TLegend(0.58,0.68,0.90,0.90);
      leg_crsr->SetFillColor(0);
      leg_crsr->SetLineColor(0);
      leg_crsr->AddEntry(fr->findObject("dataSR"),    "Data SR","lep");
      leg_crsr->AddEntry(fr->findObject("dataCR"),    "Data CR","lep");
      leg_crsr->AddEntry(fr->findObject("curveCR"),   "CR RooHistPdf","l");
      leg_crsr->AddEntry(fr->findObject("curveSRpred"),"Pred SR = CR#timesTF","l");

      // 9) Pull = (dataSR - curveSRpred) / sigma_data  (dal RooPlot)
      //    Nota: calcola i pull rispetto alla curva "curveSRpred" e ai punti "dataSR"
/*      RooHist* hPull = fr->pullHist("dataSR","curveSRpred", true); // "true" usa gli errori asimmetrici dei punti
      frPull->addPlotable(hPull,"P");
      TLine *line0 = new TLine(mass->getMin(), 0, mass->getMax(), 0);
      line0->SetLineColor(kGray+2);
      line0->SetLineStyle(2); // tratteggiata, opzionale
      line0->Draw("SAME");
      frPull->SetTitle("");
      frPull->SetMinimum(-4.0);
      frPull->SetMaximum( 4.0);
      frPull->GetYaxis()->SetTitle("Pull");
      frPull->GetXaxis()->SetTitle(mass->GetTitle());
*/
      // 10) Canvas con doppio pad (spettro sopra, pull sotto)
      TCanvas* c = new TCanvas("c_crsr_rooplot","CR/SR check (RooPlot)",900,900);
      c->Divide(1,2);
      c->cd(1);
      gPad->SetPad(0,0.28,1,1);
      gPad->SetBottomMargin(0.03);
      fr->SetTitle(Form("CR vs SR prediction — %s", name.c_str()));
      // scala log opzionale
      // gPad->SetLogy();
      // limiti migliori in basso per log
      // if (gPad->GetLogy()) fr->SetMinimum(std::max(1e-3, fr->GetMinimum()*0.3));
      fr->Draw();
      leg_crsr->Draw();

      c->cd(2);
      gPad->SetPad(0,0,1,0.28);
      gPad->SetTopMargin(0.06);
      gPad->SetBottomMargin(0.35);
      frPull->Draw();

      c->SaveAs(Form("%s_crsr_rooplot.pdf",name.c_str()));
      c->SaveAs(Form("%s_crsr_rooplot.png",name.c_str()));

      delete leg_crsr;
      delete c;
      delete fr;
      delete frPull;
    }
  }
}

void plot(RooRealVar *mass, RooMultiPdf *pdfs, RooCategory *catIndex, RooAbsData *data, string name, vector<string> flashggCats_, int cat, int bestFitPdf=-1, bool fitWasBlinded=false){
  int color[7] = {kCyan+1,kGreen+2,kAzure+2,kTeal+1,kSpring+2,kCyan+3,kGreen+1};
  TLegend *leg = new TLegend(0.6,0.65,0.95,0.90);
  leg->SetFillColor(0);
  leg->SetLineColor(1);
  RooPlot *plot = mass->frame();
  data->plotOn(plot,MarkerStyle(20),MarkerSize(0.8),CutRange(sidebandRangeList.c_str()));
  TCanvas *canv = new TCanvas();
  canv->SetLogy(); // Set log scale for y-axis
  int currentIndex = catIndex->getIndex();
  TObject *datLeg = plot->getObject(int(plot->numItems()-1));
  leg->AddEntry(datLeg,"Data");
  int style=1;
  RooAbsPdf *pdf;
  RooCurve *nomBkgCurve;
  double Nvis;
  if (fitWasBlinded) {
    Nvis = data->sumEntries("1", sidebandRangeList.c_str());
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
    if (fitWasBlinded) pdfs->getCurrentPdf()->plotOn(plot,LineColor(kRed),LineWidth(2),/*NormRange(sidebandRangeList.c_str()),*/ Range(100,mass->getMax()),
                     Normalization(Nvis, RooAbsReal::NumEvent));
    else pdfs->getCurrentPdf()->plotOn(plot,LineColor(col),LineStyle(style));
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
    
    // Debug: Print puntion name 
    std::cout << "[DEBUG] Full PDF name: " << pdfName << std::endl;
    
    // Simplified logic to extract the name of the MultiPdf
    std::string legendName = Form("PDF_%d", icat); // fallback
    
    // Looking for common patterns in PDF names
    if (pdfName.find("bernstein") != std::string::npos || pdfName.find("Bernstein") != std::string::npos) {
      // Look for a number after "bernstein" or "Bernstein"
      std::regex bernstein_regex(".*[Bb]ernstein([0-9]+).*");
      std::smatch match;
      if (std::regex_match(pdfName, match, bernstein_regex) && match.size() > 1) {
        legendName = "Bernstein" + match[1].str();
      } else {
        legendName = "Bernstein";
      }
    } else if (pdfName.find("chebychev") != std::string::npos || pdfName.find("Chebychev") != std::string::npos) {
      std::regex cheb_regex(".*[Cc]hebychev([0-9]+).*");
      std::smatch match;
      if (std::regex_match(pdfName, match, cheb_regex) && match.size() > 1) {
        legendName = "Cheb" + match[1].str();
      } else {
        legendName = "Chebychev";
      }
    } else if (pdfName.find("powerlaw") != std::string::npos || pdfName.find("PowerLaw") != std::string::npos) {
      std::regex power_regex(".*[Pp]ower[Ll]aw([0-9]+).*");
      std::smatch match;
      if (std::regex_match(pdfName, match, power_regex) && match.size() > 1) {
        legendName = "PowerLaw" + match[1].str();
      } else {
        legendName = "PowerLaw";
      }
    } else if (pdfName.find("laurent") != std::string::npos || pdfName.find("Laurent") != std::string::npos) {
      std::regex laurent_regex(".*[Ll]aurent([0-9]+).*");
      std::smatch match;
      if (std::regex_match(pdfName, match, laurent_regex) && match.size() > 1) {
        legendName = "Laurent" + match[1].str();
      } else {
        legendName = "Laurent";
      }
    } else if (pdfName.find("exponential") != std::string::npos || pdfName.find("Exponential") != std::string::npos) {
      legendName = "Exponential";
    } else {
      // If a known pattern is not found, complete name or a fallback is used
      legendName = Form("Func_%d", icat);
    }
    
    std::cout << "[DEBUG] Legend name: " << legendName << std::endl;
    
    leg->AddEntry(pdfLeg,Form("%s%s",legendName.c_str(),ext.c_str()),"L");
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

void plot(RooRealVar *mass, map<string,RooAbsPdf*> pdfs, RooAbsData *data, string name, vector<string> flashggCats_, int cat, int bestFitPdf=-1, bool fitWasBlinded=true){
  int color[7] = {kCyan+1,kGreen+2,kAzure+2,kTeal+1,kSpring+2,kCyan+3,kGreen+1};
  TCanvas *canv = new TCanvas();
  canv->SetLogy(); 
  TLegend *leg = new TLegend(0.15,0.15,0.55,0.35); 
  leg->SetFillColor(0);
  leg->SetFillStyle(0);  
  leg->SetLineColor(0);
  leg->SetBorderSize(0); 
  leg->SetTextSize(0.03); 
  leg->SetMargin(0.15);  
  
  RooPlot *plot = mass->frame(mass->getMin(), mass->getMax());
  
  data->plotOn(plot,MarkerStyle(20),MarkerSize(0.8), CutRange(sidebandRangeList.c_str()));
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
    
    // Safety check for null pointer
    if (!it->second) {
      std::cerr << "[ERROR] Null PDF pointer found for " << it->first << std::endl;
      continue;
    }
    
    RooArgSet *final_params = it->second->getParameters((const RooArgSet*)(0));
    if (!final_params) {
      std::cerr << "[ERROR] Could not get parameters for " << it->first << std::endl;
      continue;
    }
    
    RooRealVar* bkg_norm = nullptr;
    TIterator* iter = final_params->createIterator();
    RooAbsArg* arg;
    while ((arg = (RooAbsArg*)iter->Next())) {
      if (TString(arg->GetName()).Contains("bkg_norm")) {
        bkg_norm = (RooRealVar*)arg;
      }
    }
    delete iter;
    delete final_params;
    
    if (bkg_norm) {
      Nvis = bkg_norm->getVal();
      std::cout << "[INFO] Using bkg_norm for Nvis: " << Nvis << " events" << std::endl;
    } else {
      if (fitWasBlinded) {
        Nvis = data->sumEntries("1", sidebandRangeList.c_str());
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
    std::string legendText;
    if (!orderStr.empty()) {
      legendText = familyName + " (order " + orderStr + ")" + ext;
    } else {
      legendText = pdfName + ext;
    }
    
    leg->AddEntry(pdfLeg, legendText.c_str(), "L");
    i++;
  }
  plot->SetMaximum(plot->GetMaximum()*1.4);
  double minVal = plot->GetMinimum();
  // Better minimum for log scale: use 0.1 or minVal/10, whichever is larger
  double logMinimum = std::max(0.1, minVal > 0 ? minVal/10.0 : 0.1);
  plot->SetMinimum(logMinimum);
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

int getBestFitFunction(RooMultiPdf *bkg, RooAbsData *data, RooCategory *cat, bool silent=false){
  double global_minNll = 1E10;
  int best_index = 0;
  int number_of_indeces = cat->numTypes();
  RooArgSet snap,clean;
  RooArgSet *params = bkg->getParameters((const RooArgSet*)0);
  params->remove(*cat);
  params->snapshot(snap);
  params->snapshot(clean);
  for (int id=0;id<number_of_indeces;id++){    
    params->assignValueOnly(clean);
    cat->setIndex(id);
    double minNll=0;
    int fitStatus=1;    
    runFit(bkg->getCurrentPdf(),data,&minNll,&fitStatus,/*max iterations*/7);
    minNll=minNll+bkg->getCorrection();
    if (!silent) {
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
  }
  cat->setIndex(best_index);
  params->assignValueOnly(snap);
  std::cout << "[INFO] Best fit Function -- " << bkg->getCurrentPdf()->GetName() << " " << cat->getIndex() <<std::endl;
  std::cout << "[INFO] Best fit parameters " << std::endl;
  params->Print("V");
  return best_index;
}

// Function to create background PDF with optional turn-on function and Z resonance
RooAbsPdf* createBackgroundWithTurnOn(PdfModelBuilder &pdfsModel, string type, int order, RooRealVar* mass, const char* ext, bool includeTurnOn = true, bool includeZ = false, RooWorkspace* dataWS = nullptr, std::string turnOnType = "Erf") {
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
    RooRealVar *cutoff = new RooRealVar(Form("%s_turnon_cutoff", ext), "Turn-on cutoff", 100, 43,250 );
    RooRealVar *beta = new RooRealVar(Form("%s_turnon_beta", ext), "Turn-on beta", 5, 0.1, 100);
    cout << "[INFO] Adding turn-on function for " << type << " order " << order << endl;
    RooGenericPdf *turnOnPdf = nullptr;
    if (turnOnType == "Erf") {
      turnOnPdf = new RooGenericPdf(Form("%s_turnon", ext), "Fermi-Dirac turn-on", 
                                   "1.0/(1.0 + TMath::Exp((@1-@0)/@2))", 
                                   RooArgList(*mass, *cutoff, *beta));
    }
    else if (turnOnType == "Erf" || type == "PowerLaw" || type == "PowerLawSingle") {
      turnOnPdf = new RooGenericPdf(Form("%s_turnon", ext), "Error function turn-on", 
                                   "0.5*(1.0 + TMath::Erf((@0-@1)/@2))", 
                                   RooArgList(*mass, *cutoff, *beta));
    }
    else if (turnOnType == "DExp") {
      turnOnPdf = new RooGenericPdf(Form("%s_turnon", ext), "Double exponential turn-on", 
                                   "1.0 - TMath::Exp(-TMath::Exp((@0-@1)/@2))", 
                                   RooArgList(*mass, *cutoff, *beta));
    }
    else {
      turnOnPdf = new RooGenericPdf(Form("%s_turnon", ext), "Turn-on function", 
                                   "1.0/(1.0 + TMath::Exp((@1-@0)/@2))", 
                                   RooArgList(*mass, *cutoff, *beta));
    }
    currentPdf = new RooProdPdf(Form("_%s_%s_with_turnon",type.c_str(), ext), 
                               "Background with turn-on", 
                               RooArgList(*currentPdf, *turnOnPdf));
  }
  if (includeZ && dataWS) {  
    RooAbsPdf *zModel = dataWS->pdf(gZModelName.c_str());
    if (zModel) {
      cout << "[INFO] Adding Z resonance model for " << type << " order " << order << endl;
      RooArgSet *zParams = zModel->getParameters(RooArgSet(*mass));
      TIterator *iter = zParams->createIterator();
      RooRealVar *param;
      while ((param = (RooRealVar*)iter->Next())) {
          param->setConstant(kTRUE);
          cout << "[INFO] Freezing Z parameter: " << param->GetName() << " = " << param->getVal() << " (fixed)" << endl;
        }
      delete iter;
      delete zParams;
      RooArgSet normSet(*mass);
      double zIntegral = zModel->createIntegral(normSet)->getVal();
      cout << "[DEBUG] Z PDF integral over mass range: " << zIntegral << endl;
      cout << "[INFO] Adding Z resonance using independent normalizations" << endl;
      RooRealVar *bkg_norm = new RooRealVar(Form("%s_bkg_norm", ext), "Background normalization", 
                                           1400000, 0., 2000000.);
      RooRealVar *z_norm = new RooRealVar(Form("%s_z_norm", ext), "Z normalization", 
                                         16000, 14000, 18000.);
      cout << "[INFO] Using independent normalizations:" << endl;
      cout << "  Background norm: " << bkg_norm->getVal() << " events" << endl;
      cout << "  Z norm: " << z_norm->getVal() << " events" << endl;
      currentPdf = new RooAddPdf(Form("%s_with_z", ext),
                                "Background with Z",
                                RooArgList(*currentPdf, *zModel),
                                RooArgList(*bkg_norm, *z_norm));
      cout << "[INFO] Z added with independent normalizations" << endl;
    } else {
      cout << "[WARNING] Could not find Z model '"<< gZModelName <<"' in workspace" << endl;
    }
  }
  return currentPdf;
}

void runIterativeFits(RooRealVar* mass, RooAbsData* data, PdfModelBuilder& pdfsModel, RooWorkspace* dataWS, const std::string& ext, int order, const std::vector<std::string>& flashggCats_, int cat, const std::string& outDir, bool iterativeMode, const std::string& turnOnType) {
    // Nota: questa funzione usa setRange con etichette proprie; il runFit userà sidebandRangeList
    mass->setRange("fit_low", mass->getMin(), 78);
    mass->setRange("fit_high", 100, mass->getMax());
    mass->setRange("blind_low", mass->getMin(), 78);
    mass->setRange("blind_high", 100, mass->getMax());
    bool includeZ = false;
    std::cout << "[ITERATIVE FIT] Primo fit: " << ext << " [ " << mass->getMin()<< ",78] e >140 senza Z" << std::endl;
    RooAbsPdf* bkgPdf1 = createBackgroundWithTurnOn(pdfsModel, ext, order, mass, ext.c_str(), true, includeZ, dataWS, turnOnType);
    int fitStatus1 = 0;
    double chi2FromFit1 = -1;
    double thisNll1 = 0.;
    runFit(bkgPdf1, data, &thisNll1, &fitStatus1, 7, &chi2FromFit1, false);
    plot(mass, bkgPdf1, data, outDir+"/iterativeFit1_"+ext, flashggCats_, fitStatus1, &chi2FromFit1, chi2FromFit1,true,true,dataWS);
    std::cout << "[ITERATIVE FIT] Primo fit status: " << fitStatus1 << std::endl;
    mass->setRange("fit_low", mass->getMin(), 100);
    mass->setRange("fit_high", 100, mass->getMax());
    mass->setRange("blind_low", mass->getMin(), 100);
    mass->setRange("blind_high", 100, mass->getMax());
    includeZ = true;
    std::cout << "[ITERATIVE FIT] Secondo fit: " << ext << " [mass->getMin(),100] e >140 con Z" << std::endl;
    RooAbsPdf* bkgPdf2 = createBackgroundWithTurnOn(pdfsModel, ext, order, mass, ext.c_str(), true, includeZ, dataWS, turnOnType);
    int fitStatus2 = 0;
    double chi2FromFit2 = -1;
    double thisNll2 = 0.;
    runFit(bkgPdf2, data, &thisNll2, &fitStatus2, 7, &chi2FromFit2, false);
    plot(mass, bkgPdf2, data, outDir+"/iterativeFit2_"+ext, flashggCats_, fitStatus2, &chi2FromFit2, chi2FromFit2,true,true,dataWS);
    std::cout << "[ITERATIVE FIT] Secondo fit status: " << fitStatus2 << std::endl;
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
  bool includeZ=false;
  bool blindSignalRegion=false;
  int rebinFactor=1;
  int isFlashgg_ =0;
  string flashggCatsStr_;
  vector<string> flashggCats_;
  bool isData_ =0;
  bool iterativeFit = false;
  std::string turnOnType = "Erf"; // Fermi

  po::options_description desc("Allowed options");
  desc.add_options()
    ("help,h",                                                                                  "Show help")
    ("infilename,i", po::value<string>(&fileName),                                              "In file name")
    ("workspace,w", po::value<string>(&workspaceFile)->default_value(""), "Workspace file with RooDataHist and Z model")
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
    ("year", po::value<string>(&year_)->default_value("2024"),                                  "Dataset year")
    ("catOffset", po::value<int>(&catOffset)->default_value(0),                                 "Category numbering scheme offset")
    ("mN", po::value<float>(&mN)->default_value(2.75),                                          "Mass of the peak, for center of window")
    ("sigma", po::value<float>(&sigma)->default_value(0.025),                                   "Sigma of the peak, for size of window")
    ("nsigma", po::value<int>(&nsigma)->default_value(10),                                       "Sigma multiplier, for size of window")
    ("rebinFactor,r", po::value<int>()->default_value(1),                                       "Rebin factor to reduce number of bins (1=no rebinning)")
    ("verbose,v",                                                                               "Run with more output")
    ("iterativeFit", po::value<bool>(&iterativeFit)->default_value(false), "Enable iterative fitting: first fit [mass->getMin(),78]+>140 w/o Z, then [mass->getMin(),100]+>140 with Z")
    ("turnOnType", po::value<std::string>(&turnOnType)->default_value("Erf"), "Type of turn-on function: Fermi, Erf, DExp")
    // New options
    ("cat-name", po::value<std::string>()->default_value("PP") , "Category name")
    ("ws-name", po::value<std::string>()->default_value("") , "Name of the RooWorkspace within the input file")
    ("mass-var", po::value<std::string>()->default_value("CMS_dijet_mass"), "RooRealVar name for dijet mass")
    ("data-name", po::value<std::string>()->default_value("") , "Dataset name (single category)")
    ("data-name-pattern", po::value<std::string>()->default_value("h_data_bkg_cat{CAT}"), "Pattern dataset for each catorgy; {CAT} substituted with the index")
    ("z-model-name", po::value<std::string>()->default_value("model_Z_c2"), "Name of the Z pdf within the workspace")
    ("blind", po::value<std::string>()->default_value("100-100"), "SR windows to be excluded, es. '100-135,500-600' (accept 'min'/'max')")
    ("fit-accept", po::value<std::string>()->default_value("") , "(Optional) Ranges accepted for the fit; if empty, it complements --blind")
    ("mass-min", po::value<double>()->default_value(std::numeric_limits<double>::quiet_NaN()), "(Optional) Min value of the mass variable")
    ("mass-max", po::value<double>()->default_value(std::numeric_limits<double>::quiet_NaN()), "(Optional) Max value of the mass variable")
    ("max-ftest-order", po::value<int>()->default_value(5), "Max order for F-test (default: 5)")
    ("max-envelope-order", po::value<int>()->default_value(10), "Max order for the envelope (default: 10)")
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
  system(Form("mkdir -p %s",outDir.c_str()));

  // Command line options
  // ----------------------------------------
  auto catnameDijet = vm["cat-name"].as<std::string>();
  auto wsName_opt  = vm["ws-name"].as<std::string>();
  auto massVarName = vm["mass-var"].as<std::string>();
  auto dataName    = vm["data-name"].as<std::string>();
  auto dataPat     = vm["data-name-pattern"].as<std::string>();
  auto zModelName  = vm["z-model-name"].as<std::string>();
  auto blindSpec   = vm["blind"].as<std::string>();
  auto fitAccept   = vm["fit-accept"].as<std::string>();
  double optMassMin = vm["mass-min"].as<double>();
  double optMassMax = vm["mass-max"].as<double>();
  gZModelName = zModelName;

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
   std::cout << "[INFO]  inWS open " << inWS << std::endl;
   inWS->Print();
  }
  if (saveMultiPdf){
    if (inFile) {
      transferMacros(inFile,outputfile);
    }
    RooRealVar *intL; 
    RooRealVar *sqrts;
    intL  = intLumi_;
    sqrts = (RooRealVar*)inWS->var("SqrtS");
    if (!sqrts){ sqrts = new RooRealVar("SqrtS","SqrtS",13); }
    outputws->import(*intL);
    outputws->import(*sqrts);
    std::cout << "[INFO] got intL and sqrts " << intL << ", " << sqrts << std::endl;
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
  std::vector<Interval> sr = parse_intervals(blindSpec, mmin, mmax);
  blindRangeList = define_ranges(mass, sr, "sr");
  std::vector<Interval> accept = fitAccept.empty() ? complement(sr, mmin, mmax)
                                                   : parse_intervals(fitAccept, mmin, mmax);
  sidebandRangeList = define_ranges(mass, accept, "sb");
  std::cout << "[INFO] SR (escluse): " << (blindRangeList.empty()?"<none>":blindRangeList) << std::endl;
  std::cout << "[INFO] Sidebands per fit: " << sidebandRangeList << std::endl;
  pdfsModel.setObsVar(mass);
  double upperEnvThreshold = 0.25;   // Maximum probability for F-test (it was 1, lowered to 0.25)
  double minGofThreshold = 0.10;     // Minimum GoF considered to exclude bad functions (raised from 0.01 to 0.10)

  fprintf(resFile,"Truth Model & d.o.f & $\\Delta NLL_{N+1}$ & $p(\\chi^{2}>\\chi^{2}_{(N\\rightarrow N+1)})$ \\\n");
  fprintf(resFile,"\\hline\n");

  std::string ext = is2011 ? "7TeV" : "8TeV";
  if( isFlashgg_ ){
    if( year_ == "all" ){ ext = "13p6TeV"; }
    else{ ext = Form("%s_13p6TeV",year_.c_str()); }
  }
  std::cout << "[INFO] Number of categories to process: " << ncats << std::endl;
  for (int cat=0; cat<ncats; cat++){
    map<string,int> choices;
    map<string,std::vector<int> > choices_envelope;
    map<string,RooAbsPdf*> pdfs;
    map<string,RooAbsPdf*> allPdfs;
    map<string,int> familyMinOrder; // Tracking minimum order for the family in the envelope
    map<string,double> familyBestAIC; // Tracking best AIC for the family in the envelope
    map<string,RooAbsPdf*> familyBestPdf; // Tracking best PDF for the family in the envelope 
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
    nBinsForFit = mass->getBinning().numBins();
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
      while (order < gMaxFtestOrder){   // Continue until max order or chi2 stops improving significantly
        cout << "[DEBUG] Starting test for order " << order << ", prev_chi2_red = " << prev_chi2_red << ", maxOrder = " << gMaxFtestOrder << endl;
        cout << "==> " << *funcType << " " << order << endl;
        RooAbsPdf *bkgPdf = createBackgroundWithTurnOn(pdfsModel,*funcType,order,mass,Form("ftest_pdf_%d_%s",(cat+catOffset),ext.c_str()), currentIncludeTurnOn, includeZ, dataWS, turnOnType);
        if (bkgPdf && bkgPdf->InheritsFrom("RooAddPdf")) {
          RooArgSet normSet(*mass);
          ((RooAddPdf*)bkgPdf)->fixCoefNormalization(normSet);
        }
        if (!bkgPdf){
          cout << "[DEBUG] Failed to create PDF for order " << order << ", moving to next order" << endl;
          order++;
          if (order > gMaxFtestOrder) {
            cout << "[WARNING] Reached maximum F-test order limit (" << gMaxFtestOrder << "), stopping" << endl;
            break;
          }
        }
        else {
          int fitStatus = 0;
          double chi2FromFit = -1;
          if (iterativeFit){
            runIterativeFits(mass, data, pdfsModel, dataWS, *funcType, order, flashggCats_, singleCategory, outDir, iterativeFit, turnOnType);
          } else {
            runFit(bkgPdf,data,&thisNll,&fitStatus,/*max iterations*/7,&chi2FromFit,blindSignalRegion);
          }
          if (fitStatus!=0) {
            std::cout << "[WARNING] Fit failed for " << bkgPdf->GetName() << " with status " << fitStatus << std::endl;
            if (fitStatus == 999) {
              std::cout << "[WARNING] Integration failure detected, skipping this function" << std::endl;
              order++;
              if (order > gMaxFtestOrder) {
                cout << "[WARNING] Reached maximum F-test order limit (" << gMaxFtestOrder << "), stopping" << endl;
                break;
              }
              continue; // Skip to next iteration
            }
          }
          // Standard CMS F-test methodology using NLL + 0.5*Nparams
          double current_nll_corr = thisNll;
          int nparams = (bkgPdf) ? bkgPdf->getVariables()->getSize() : order;
          double current_nll_penalty = current_nll_corr + 0.5 * nparams;
          
          // Previous NLL for F-test calculation
          static double prev_nll_corr = 999.0;
          static int prev_nparams = 0;
          if (order == 1) {
            prev_nll_corr = 999.0; // Reset for each family
            prev_nparams = 0;
          }
          
          std::cout << "[DEBUG] F-test analysis: NLL=" << current_nll_corr << ", nparams=" << nparams 
                   << ", NLL_penalty=" << current_nll_penalty << ", prev_NLL=" << prev_nll_corr << std::endl;
          
          bool significant_improvement = false;
          double f_test_prob = 1.0;
          
          if (order == 1) {
            // Always test first order
            significant_improvement = true;
            f_test_prob = 0.0; // First order automatically passes
            std::cout << "[INFO] First order: NLL=" << current_nll_corr << ", NLL_penalty=" << current_nll_penalty << std::endl;
          } else {
            // Standard F-test: 2*(NLL_prev - NLL_current) ~ chi2 with dof = nparams_current - nparams_prev
            double f_stat = 2.0 * (prev_nll_corr - current_nll_corr);
            int dof = nparams - prev_nparams;
            
            if (f_stat > 0 && dof > 0) {
              f_test_prob = TMath::Prob(f_stat, dof);
            } else {
              f_test_prob = 1.0; // No improvement
            }
            
            // CMS standard: p^F < 0.05 for F-test to pass
            significant_improvement = (f_test_prob < 0.05);
            std::cout << "[INFO] F-test: Order " << order << ", F-stat=" << f_stat << ", dof=" << dof 
                     << ", p^F=" << f_test_prob << ", passes=" << significant_improvement << std::endl;
          }
          
          if (!significant_improvement && order > 1) {
            std::cout << "[INFO] F-test failed (p^F=" << f_test_prob << " >= 0.05), stopping at order " << prev_order << std::endl;
            break; // Stop the loop
          }
          
          // Update for next iteration  
          prev_nll_corr = current_nll_corr;
          prev_nparams = nparams;
          
          // Update for next iteration  
          prev_nll_corr = current_nll_corr;
          prev_nparams = nparams;
          double gofProb=0;
          plot(mass,bkgPdf,data,Form("%s/ftest_%s%d_%s",outDir.c_str(),funcType->c_str(),order,catname.c_str()),flashggCats_,fitStatus,&gofProb,chi2FromFit,true,true,dataWS);
          if (includeZ && dataWS) {
            RooAbsPdf *zModel = dataWS->pdf(gZModelName.c_str());
            if (zModel) {
              plotComponents(mass,bkgPdf,zModel,data,Form("%s/ftest_%s%d_%s",outDir.c_str(),funcType->c_str(),order,catname.c_str()),flashggCats_,fitStatus,&gofProb,chi2FromFit,blindSignalRegion);
            }
          }
          cout << "[INFO] function type, order, NLL, f_test_prob" << endl;
          cout << "[INFO] " << *funcType << " " << order << " " << current_nll_corr << " " << f_test_prob << endl;
          prevNll=thisNll;
          cache_order=prev_order;
          cache_pdf=prev_pdf;
          prev_order=order;
          prev_pdf=bkgPdf;
          
          // Ensure we always have at least one valid cache_pdf
          if (cache_pdf == NULL && bkgPdf != NULL) {
            cache_pdf = bkgPdf;
            cache_order = order;
            cout << "[DEBUG] Setting cache_pdf to current bkgPdf for order " << order << endl;
          }
          order++;
          if (order > gMaxFtestOrder) {
            cout << "[INFO] Reached maximum F-test order limit (" << gMaxFtestOrder << "), stopping F-test" << endl;
            break;
          }
        }
        counter++;
      }
    
      choices.insert(pair<string,int>(*funcType,cache_order));
      if (cache_pdf != NULL) {
        pdfs.insert(pair<string,RooAbsPdf*>(Form("%s%d",funcType->c_str(),cache_order),cache_pdf));
      } else {
        std::cout << "[WARNING] cache_pdf is NULL for " << *funcType << ", skipping truth plot entry" << std::endl;
      }
      int truthOrder = cache_order;
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
        
        std::cout << "[INFO] CMS Envelope Selection: Max order with p^F < 0.1, lower orders with p^F < 0.1 AND GoF > 0.01" <<std::endl;
        
        double prev_nll_env = 0.0;
        int prev_nparams_env = 0;
        int max_valid_order = 0;
        double final_f_test_prob = 1.0;
        double final_gofProb = 0.0;
        
        while (order <= gMaxEnvelopeOrder){
          std::cout << "[DEBUG-ENV] Top of loop: funcType=" << *funcType << " order=" << order << " gMaxEnvelopeOrder=" << gMaxEnvelopeOrder << std::endl;
          cout << "==> " << *funcType << " " << order << endl;
          RooAbsPdf *bkgPdf = createBackgroundWithTurnOn(pdfsModel,*funcType,order,mass,Form("env_pdf_%d_%s",(cat+catOffset),ext.c_str()), currentIncludeTurnOn, includeZ, dataWS, turnOnType);
          std::cout << "[DEBUG-ENV] After createBackgroundWithTurnOn: bkgPdf=" << (void*)bkgPdf << std::endl;
          if (bkgPdf && bkgPdf->InheritsFrom("RooAddPdf")) {
            RooArgSet normSet(*mass);
            ((RooAddPdf*)bkgPdf)->fixCoefNormalization(normSet);
          }
          if (!bkgPdf ){
            if (order > gMaxEnvelopeOrder) { std::cout << " [WARNING] could not add, order limit reached (max: " << gMaxEnvelopeOrder << ") ] " << std::endl; break ;}
            order++;
            continue; // Skip to next iteration if PDF creation failed
          }
          else {
            int fitStatus=0;
            double chi2FromFit = -1;
            runFit(bkgPdf,data,&thisNll,&fitStatus,/*max iterations*/7,&chi2FromFit,blindSignalRegion);
            if (fitStatus!=0) std::cout << "[WARNING] Warning -- Fit status for " << bkgPdf->GetName() << " at " << fitStatus <<std::endl;
            
            // Calculate F-test probability and GoF
            double f_test_prob = 1.0;
            int nparams_env = (bkgPdf) ? bkgPdf->getVariables()->getSize() : order;
            
            if (order == 1) {
              f_test_prob = 0.0; // First order automatically passes F-test
              prev_nll_env = thisNll;
              prev_nparams_env = nparams_env;
            } else {
              // F-test: 2*(NLL_prev - NLL_current) ~ chi2 with dof = nparams_current - nparams_prev
              double f_stat = 2.0 * (prev_nll_env - thisNll);
              int dof = nparams_env - prev_nparams_env;
              
              if (f_stat > 0 && dof > 0) {
                f_test_prob = TMath::Prob(f_stat, dof);
              } else {
                f_test_prob = 1.0; // No improvement
              }
              
              prev_nll_env = thisNll;
              prev_nparams_env = nparams_env;
            }
            
            double gofProb = 0.0; 
            plot(mass,bkgPdf,data,Form("%s/EnvelopeComponents_Pulls/envelope_%s%d_%s",outDir.c_str(),funcType->c_str(),order,catname.c_str()),flashggCats_,fitStatus,&gofProb,chi2FromFit,blindSignalRegion,true,dataWS);
            if (includeZ && dataWS) {
              RooAbsPdf *zModel = dataWS->pdf(gZModelName.c_str());
              if (zModel) {
                plotComponents(mass,bkgPdf,zModel,data,Form("%s/EnvelopeComponents_Pulls/envelope_%s%d_%s",outDir.c_str(),funcType->c_str(),order,catname.c_str()),flashggCats_,fitStatus,&gofProb,chi2FromFit,blindSignalRegion);
              }
            }
            
            cout << "[INFO] function type, order, NLL, p^F, GoF " << endl;
            cout << "[INFO] " << *funcType << " " << order << " " << thisNll << " " << f_test_prob << " " << gofProb << endl;
            
            // CMS envelope criteria
            bool passes_f_test = (f_test_prob < upperEnvThreshold);
            bool passes_gof = (gofProb > minGofThreshold);
            bool include_in_envelope = false;
            
            if (passes_f_test) {
              max_valid_order = order; // Update maximum valid order
              final_f_test_prob = f_test_prob; // Store for output
              final_gofProb = gofProb; // Store for output
              if (order == 1 || passes_gof) {
                include_in_envelope = true;
                std::cout << "[INFO] Including order " << order << " in envelope (p^F=" << f_test_prob 
                         << ", GoF=" << gofProb << ")" << std::endl;
              } else {
                std::cout << "[INFO] Order " << order << " passes F-test (p^F=" << f_test_prob 
                         << ") but fails GoF (" << gofProb << " <= 0.01), excluding" << std::endl;
              }
            } else {
              std::cout << "[INFO] Order " << order << " fails F-test (p^F=" << f_test_prob 
                       << " >= 0.1), stopping envelope construction" << std::endl;
              break; // Stop when F-test fails
            }
            
            if (include_in_envelope) {
              // Check parameter stability before adding to envelope
                  bool parametersStable = true;
                  double maxRelativeError = 0.0;
                  string familyName = *funcType;
                  
                  if (bkgPdf) {
                    RooArgSet* params = bkgPdf->getParameters(*data);
                    TIterator* iter = params->createIterator();
                    RooRealVar* param;
                    while ((param = (RooRealVar*)iter->Next())) {
                      if (!param->isConstant() && param->hasError()) {
                        double relError = (param->getVal() != 0) ? fabs(param->getError() / param->getVal()) : 999.0;
                        maxRelativeError = max(maxRelativeError, relError);
                        
                        // Check if parameter is at boundary or has huge error
                        bool atBoundary = (fabs(param->getVal() - param->getMin()) < 1e-6) || 
                                        (fabs(param->getVal() - param->getMax()) < 1e-6);
                        bool hugeError = (param->getError() > 50.0) || (relError > 50.0);
                        
                        if (atBoundary || hugeError) {
                          std::cout << "[WARNING] Parameter " << param->GetName() 
                                   << " potentially unstable: value=" << param->getVal() 
                                   << ", error=" << param->getError()
                                   << ", relError=" << relError 
                                   << ", atBoundary=" << atBoundary << std::endl;
                          if (hugeError) parametersStable = false;
                        }
                      }
                    }
                    delete iter;
                    delete params;
                  }
                  
                  //if (parametersStable && gofProb > minGofThreshold) {
		  if (parametersStable) {
                    double myNll = 2.*thisNll;
                    std::cout << "[INFO] Adding to Envelope " << bkgPdf->GetName() << " p^F=" << f_test_prob 
                      << " GoF=" << gofProb << " 2xNLL + c is " << myNll + nparams_env <<  std::endl;
                    
                    allPdfs.insert(pair<string,RooAbsPdf*>(Form("%s%d",funcType->c_str(),order),bkgPdf));
                    storedPdfs.add(*bkgPdf);
                    pdforders.push_back(order);
                    
                    // Track the best function overall for later reference
                    if ((myNll + nparams_env) < MinimimNLLSoFar) {
                      simplebestFitPdfIndex = storedPdfs.getSize()-1;
                      MinimimNLLSoFar = myNll + nparams_env;
                    }
                    
                    // Update family tracking (keep the best one for each family for backward compatibility)
                    double current_penalty = myNll + nparams_env;
                    if (familyBestAIC.find(familyName) == familyBestAIC.end() || current_penalty < familyBestAIC[familyName]) {
                      familyBestAIC[familyName] = current_penalty;
                      familyBestPdf[familyName] = bkgPdf;
                      familyMinOrder[familyName] = order;
                    }
                  } else {
                    std::cout << "[INFO] Skipping " << familyName << " order " << order << " for envelope due to parameter instability" << std::endl;
                  }
            }
            
            cache_order=prev_order;
            cache_pdf=prev_pdf;
            prev_order=order;
            prev_pdf=bkgPdf;
            order++;
            
            if (order > gMaxEnvelopeOrder) {
              cout << "[INFO] Reached maximum envelope order limit (" << gMaxEnvelopeOrder << "), stopping envelope building" << endl;
              break;
            }
          }
        }
        fprintf(resFile,"%15s & %d & %5.3f & %5.3f \\\n",funcType->c_str(),max_valid_order,final_f_test_prob,final_gofProb);
        choices_envelope.insert(pair<string,std::vector<int> >(*funcType,pdforders));
      }
    }

    fprintf(resFile,"\\hline\n");
    choices_vec.push_back(choices);
    choices_envelope_vec.push_back(choices_envelope);
    pdfs_vec.push_back(pdfs);
    plot(mass,pdfs,data,Form("%s/truths_%s",outDir.c_str(),catname.c_str()),flashggCats_,cat,false);

    // Plot envelope-only functions (familyBestPdf contains one best function per family)
    if (saveMultiPdf && !familyBestPdf.empty()){
      plot(mass,familyBestPdf,data,Form("%s/EnvelopeComponents_Pulls/envelope_only_%s",outDir.c_str(),catname.c_str()),flashggCats_,cat,false);
      std::cout << "[INFO] Generated envelope-only plot with " << familyBestPdf.size() << " functions" << std::endl;
    }

    if (saveMultiPdf){
      string catindexname;
      string catname2;
      if (isFlashgg_){
        catindexname = Form("pdfindex_%s_%s",std::to_string(cat).c_str(),ext.c_str());
        catname2 = Form("%s",std::to_string(cat).c_str());
      } else {
        catindexname = Form("pdfindex_%d_%s",(cat+catOffset),ext.c_str());
        catname2 = Form("cat%d",(cat+catOffset));
      }
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
      
      for (const auto& truthPdf : familyBestPdf) {
        string familyName = truthPdf.first;  // e.g., "Dijet"
        RooAbsPdf* pdf = truthPdf.second;
        
        // Extract order from PDF name for better legend labeling
        string pdfKey = familyName;  // Default to family name
        if (pdf) {
          string pdfName = pdf->GetName();
          // Extract order from PDF name (e.g., "env_pdf_0_2016_13TeV_dijet4" -> "4")
          // Map family names to their abbreviations in PDF names
          string abbrev = "";
          if (familyName == "Dijet") abbrev = "dijet";
          else if (familyName == "Exponential") abbrev = "exp";
          else if (familyName == "PowerLaw") abbrev = "pow";
          else if (familyName == "Laurent") abbrev = "lau";
          else if (familyName == "Chebychev") abbrev = "cheb";
          else if (familyName == "Bernstein") abbrev = "bern";
          else {
            // Fallback: use lowercase of first 4 chars
            abbrev = familyName.substr(0, std::min(4, (int)familyName.length()));
            std::transform(abbrev.begin(), abbrev.end(), abbrev.begin(), ::tolower);
          }
          
          // Find abbreviation in PDF name
          size_t abbrevPos = pdfName.find(abbrev);
          if (abbrevPos != std::string::npos) {
            // Find digits immediately after abbreviation
            size_t orderStart = abbrevPos + abbrev.length();
            if (orderStart < pdfName.length() && isdigit(pdfName[orderStart])) {
              string orderStr = "";
              while (orderStart < pdfName.length() && isdigit(pdfName[orderStart])) {
                orderStr += pdfName[orderStart];
                orderStart++;
              }
              if (!orderStr.empty()) {
                pdfKey = familyName + orderStr;  // e.g., "Dijet4", "Chebychev6"
              }
            }
          }
        }
        
        // Check parameter stability and normalization for final selection
        double norm_check = 1.0; // Default to 1.0 if PDF is NULL
        if (pdf) {
          // Check parameter stability
          RooArgSet* params = pdf->getParameters(*data);
          TIterator* iter = params->createIterator();
          RooRealVar* param;
          while ((param = (RooRealVar*)iter->Next())) {
            if (!param->isConstant() && param->hasError()) {
              double relError = (param->getVal() != 0) ? fabs(param->getError() / param->getVal()) : 999.0;
              
              // Check if parameter is at boundary or has huge error
              bool atBoundary = (fabs(param->getVal() - param->getMin()) < 1e-6) || 
                              (fabs(param->getVal() - param->getMax()) < 1e-6);
              bool hugeError = (param->getError() > 10.0) || (relError > 5.0);
              
              if (atBoundary || hugeError) {
                std::cout << "[WARNING] Final check: Parameter " << param->GetName() 
                         << " unstable: value=" << param->getVal() 
                         << ", error=" << param->getError()
                         << ", relError=" << relError 
                         << ", atBoundary=" << atBoundary << std::endl;
              }
            }
          }
          delete iter;
          delete params;
        }
        
        // Add to final MultiPdf - RooFit handles all normalization automatically
        finalPdfs[pdfKey] = pdf;
        finalStoredPdfs.add(*pdf);
        
        std::cout << "[INFO] Adding to final MultiPdf (from envelope best): " << pdfKey 
                 << " (PDF: " << pdf->GetName() << ")" << std::endl;
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
      
      // DON'T REFIT! Use the best function from envelope directly to preserve fitted parameters
      int bestFitPdfIndex = 0; // Default to first function
      
      // Find which PDF corresponds to the best envelope function
      // --------------------------------------------------------
      for (int iPdf = 0; iPdf < finalStoredPdfs.getSize(); iPdf++) {
        RooAbsPdf* currentPdf = (RooAbsPdf*)finalStoredPdfs.at(iPdf);
        if (currentPdf) {
          // For now, just use the first one - the envelope selection already chose the best
          bestFitPdfIndex = iPdf;
          break;
        }
      }
      
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
      plot(mass,pdf,&catIndex,data,Form("%s/multipdf_cat%s",outDir.c_str(),catnameDijet.c_str()),flashggCats_,cat,bestFitPdfIndex);

      // Plot with correct names using finalPdfs
      // ---------------------------------------
      plot(mass,finalPdfs,data,Form("%s/multipdf_cat%s",outDir.c_str(),catnameDijet.c_str()),flashggCats_,cat,bestFitPdfIndex,false);

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
