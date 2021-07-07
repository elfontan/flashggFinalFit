import ROOT
import json
import numpy as np
from scipy.optimize import minimize
import scipy.stats
from collections import OrderedDict as od
from array import array

from backgroundFunctions import *
from fittingTools import *

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class modelBuilder:
  # Constructor
  def __init__(self,_name,_cat,_xvar,_data,_functionFamilies,_nBins,_blindingRegion,_minimizerMethod,verbose=False):
    self.name = _name
    self.cat = _cat
    self.xvar = _xvar
    self.data = _data
    self.functionFamilies = _functionFamilies
    self.nBins = _nBins
    self.blind = True
    self.blindingRegion = _blindingRegion
    self.minimizerMethod = _minimizerMethod
    # Containers to store model objects
    self.params = od()
    self.formulas = od()
    self.functions = od()
    self.pdfs = od()
    self.envelopePdfs = od()
    self.DataHist = None
    self.minNLL = 1e10
    self.bestfitPdfIndex = 0
    self.bestfitSnapshot = None
    self.maxTries = 3
    # Prepare datahist to perform fit to
    self.xvar.setBins(self.nBins)
    self.prepareDataHists()
    # Fit containers
    self.FitParameters = None # RooArgList of parameters which are being fitted
    self.FitPDF = None # Current RooPdf to be fit
    self.Ndof = None
    self.chi2 = None
    self.NLL = None
    self.FitResult = None
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Function toggle blinding region on/off
  def setBlind(self,blind=True):
    self.blind = blind

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Function to build pdf
  def buildPdf(self,funcType,order,ext=""):
    prefix = "%s_%s%g%s"%(self.name,self.functionFamilies[funcType]['name'][1],order,ext)
    return getPdf(self,prefix,funcType,order)

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
  # Function for setting N degrees of freedom
  def setNdof(self,_ndof): self.Ndof = _ndof

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Function to extract param bounds
  def extractXBounds(self):
    XBounds = []
    for i in range(len(self.FitParameters)): XBounds.append((self.FitParameters[i].getMin(),self.FitParameters[i].getMax()))
    return np.asarray(XBounds)

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Function to extract initial param value vector
  def extractX0(self):
    X0 = []
    for i in range(len(self.FitParameters)): X0.append(self.FitParameters[i].getVal())
    return np.asarray(X0)


  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Function to convert data to RooDataHist
  def prepareDataHists(self):
    self.DataHist = ROOT.RooDataHist("%s_hist"%self.data.GetName(),"%s_hist"%self.data.GetName(),ROOT.RooArgSet(self.xvar),self.data)

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Function to get pdf normalisation: match with data in sidebands and extrapolate for pdf in signal region
  def getNorm(self,_pdf):
    # Using numpy
    nPdf_sb, nData_sb = [], []
    nPdf_tot = []
    # Loop over bins
    for i in range(self.DataHist.numEntries()):
      p = self.DataHist.get(i)
      self.xvar.setVal(p.getRealValue(self.xvar.GetName()))
      ndata = self.DataHist.weight()
      npdf = self.FitPDF.getVal( ROOT.RooArgSet(self.xvar) )*self.DataHist.binVolume()

      # If blinding
      if( self.blind ):
        if( self.xvar.getVal() > self.blindingRegion[0] )&( self.xvar.getVal() < self.blindingRegion[1] ): 
          nPdf_tot.append( npdf )
        else:
          nData_sb.append( ndata )
          nPdf_sb.append( npdf )
          nPdf_tot.append( npdf )
      else:
        nData_sb.append( ndata )
        nPdf_sb.append( npdf )
        nPdf_tot.append( npdf )

    # Convert to numpy array
    nPdf_tot = np.asarray(nPdf_tot)
    nPdf_sb = np.asarray(nPdf_sb)
    nData_sb = np.asarray(nData_sb)

    # Re-normalise pdf to have same yield as data: accounts for blinding region
    normFactor = nData_sb.sum()/nPdf_sb.sum()
    nPdf_tot = nPdf_tot*normFactor

    return nPdf_tot.sum()

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  def runFit(self,_pdf,_mode='NLL',_verbose=False):
    # Extract fit variables: remove xvar from fit params
    fv = _pdf.getVariables().Clone()
    fv.remove(self.xvar)
    self.FitParameters = ROOT.RooArgList(fv)
    self.FitPDF = _pdf

    # Create initial vector of parameters and calculate initial Chi2
    if _verbose: print "\n --> (%s) Initialising fit parameters"%_pdf.GetName()
    x0 = self.extractX0()
    xbounds = self.extractXBounds()

    if _mode == 'NLL': self.NLL = self.getNLL()
    elif _mode == "chi2": self.chi2 = self.getChi2()
    else: 
      print " --> [ERROR] Fitting mode (%s) not supported. Use NLL or chi2"%_mode
      sys.exit(1)

    # Print parameter pre-fit values
    if _verbose: self.printFitParameters(title="Pre-fit (%s)"%_pdf.GetName(),_mode=_mode)

    # Run fit
    if _verbose: print " --> (%s) Running fit"%_pdf.GetName()
    if _mode == 'NLL': self.FitResult = minimize(NLL,x0,args=self,bounds=xbounds,method=self.minimizerMethod)
    else: self.FitResult = minimize(Chi2,x0,args=self,bounds=xbounds,method=self.minimizerMethod)

    # Add re-tries: randomize parameter set

    # Extract final fit result
    if _verbose: self.printFitParameters(title="Post-fit (%s)"%_pdf.GetName(),_mode=_mode)
    if _mode == 'NLL': 
      self.NLL = self.getNLL()
      return self.NLL, self.FitResult
    else: 
      self.chi2 = self.getChi2()
      return self.chi2, self.FitResult

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Function to calculate NLL after setting fit parameters
  def getNLL(self,_verbose=False):
    x = self.extractX0()
    self.NLL = NLL(x,self,verbose=_verbose)
    return self.NLL

  # Function to calculate chi2 after setting fit parameters
  def getChi2(self,_verbose=False):
    x = self.extractX0()
    self.chi2 = Chi2(x,self,verbose=_verbose)
    return self.chi2
 
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Function to print fit values
  def printFitParameters(self,title="Fit",_mode="NLL"):
    print " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    print " --> %s parameter values:"%(title)
    # Skip MH
    for i in range(0,len(self.FitParameters)): print "    * %-20s = %.6f   [%.6f,%.6f]"%(self.FitParameters[i].GetName(),self.FitParameters[i].getVal(),self.FitParameters[i].getMin(),self.FitParameters[i].getMax())
    print "    ~~~~~~~~~~~~~~~~"
    if _mode == "NLL":
      print "    * NLL = %.6f, n(dof) = %g"%(self.getNLL(),int(self.Ndof))
      #print "    ~~~~~~~~~~~~~~~~"
      #print "    * [VERBOSE] NLL = %.6f"%(self.getNLL(_verbose=True))
    elif _mode == "chi2":
      print "    * chi2 = %.6f, n(dof) = %g --> chi2/n(dof) = %.3f"%(self.getChi2(),int(self.Ndof),self.getChi2()/int(self.Ndof))
      #print "    ~~~~~~~~~~~~~~~~"
      #print "    * [VERBOSE] chi2 = %.6f"%(self.getChi2(_verbose=True))
    print " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"


  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # TODO: envelope and goodness of fit calculations with new model
  # Function to return prob for fTest: asymptotic or using toys
  def getProbabilityFTest(self,dchi2,ndof):
    prob_asymptotic = ROOT.TMath.Prob(dchi2,ndof)
    return prob_asymptotic
    # TODO: option for extract pval from toys

  # Function to extract goodness-of-fit from pdf to data
  def getGOF(self,_pdfInfo,_verbose=True):
    
    # Convert NLL to chi2Get chi2 value from plotting (reweighted) data and pdf
    chi2 = 2*_pdfInfo['NLL']
    ndof = _pdfInfo['Ndof']
    pval = ROOT.TMath.Prob(chi2,ndof)

    if _verbose:
      print " --> Calculating GOF for pdf: %s"%_pdfInfo['name']
      print "    * Number of degrees of freedom: %g"%ndof
      print "    * Observed chi2 = %.6f"%chi2
      print "    * GOF pval = %.6f"%pval
    
    return pval 
    # TODO: option to extract gof from toys



  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Function to perform fTest: save functions that pass to self.functions
  def fTest(self,_maxOrder=6,_pvalFTest=0.05):

    # Loop over function families
    for ff in functionFamilies:
      prevNLL, prob = 0., 0.
      order, prevOrder = 1, 1

      # Loop over orders in function until pval stopping criterion is reached
      while(( prob < _pvalFTest )&( order < _maxOrder+1 )):

        # Build PDF
        bkgPdf = self.buildPdf(ff,order)

        if not bkgPdf:
          # Order not allowed
          order += 1

        else:
          # Run fit and save info
          print "\n --------------------------------------------------------------------"
          print " --> Running fit for bkg pdf: %s\n"%bkgPdf.GetName()
          NLL, fitStatus = self.runFit(bkgPdf,_mode='NLL',_verbose=True)

          # Add to pdfs container
          self.pdfs[(ff,order)] = od()
          self.pdfs[(ff,order)]['name'] = bkgPdf.GetName()
          self.pdfs[(ff,order)]['pdf'] = bkgPdf
          self.pdfs[(ff,order)]['status'] = fitStatus
          self.pdfs[(ff,order)]['NLL'] = NLL
          self.pdfs[(ff,order)]['family'] = ff
          self.pdfs[(ff,order)]['order'] = order
          # Extract norm of function
          self.pdfs[(ff,order)]['norm'] = self.getNorm(bkgPdf)
          # Save n(d.o.f)
          self.pdfs[(ff,order)]['Ndof'] = self.Ndof

          # If fit status is not zero: fit has failed
          #if self.pdfs[(ff,order)]['status'] != 0:
          #  print " --> [WARNING] fit did NOT converge successfully for pdf: %s"%self.pdfs[(ff,order)]['name']
          #  print "     * Assumes NLL > 1e10 --> will not attempt higher orders"
          
          # Calculate dChi2 between this order and previous          
          if order > 1:
            dchi2 = 2.*(prevNLL-NLL) if prevNLL>NLL else 0.
            prob = self.getProbabilityFTest(dchi2,order-prevOrder)
          else:
            prob = 0

          # Add label to dict if order passes fTest
          print "\n --> (%s, order=%g): Prob( chi2 > chi2(data) ) = %.10f"%(ff,order,prob)
          print " --------------------------------------------------------------------"

          if prob < _pvalFTest: self.pdfs[(ff,order)]['ftest'] = True
          else: self.pdfs[(ff,order)]['ftest'] = False
         
          # Store vals and add one to order
          prevNLL = NLL
          prevOrder = order
          order += 1

  
  # Function to loop over possible pdfs which pass the Ftest and implement goodness of fit criteria
  def goodnessOfFit(self,_gofCriteria=0.01):
    # Loop over pdfs in model and select those passing fTest
    for k, v in self.pdfs.iteritems():
      if v['ftest']:
        gof = self.getGOF(v,_verbose=True)
        if gof > _gofCriteria:
          print " --> Adding pdf to envelope: (%s,%s)\n"%(k[0],k[1])
          self.envelopePdfs[k] = v
        else:
          print " --> Not adding pdf to envelope: (%s,%s)\n"%(k[0],k[1])


  # Function to build envelope of bkg functions using RooMultiPdf class
  def buildEnvelope(self,_extension=""):
    # Check if zero pdfs have satisfied fit criteria
    if len(self.envelopePdfs) == 0:
      print " --> [ERROR] No bkg functions satisfy the fit criteria. Try to increase opt.pvalFTest or decrease opt.gofCriteria"
      sys.exit(1)

    # Create pdf index to label pdf in envelope
    pdfIndexName = "pdfindex_%s%s"%(self.cat,_extension)
    self.pdfIndex = ROOT.RooCategory(pdfIndexName,"c")

    # Loop over pdfs in envelope and store in container
    pdflist = ROOT.RooArgList()
    for k,v in self.envelopePdfs.iteritems(): pdflist.add(v['pdf'])

    # Create RooMultiPdf
    self.multipdf = ROOT.RooMultiPdf("CMS_hgg_%s%s_bkgshape"%(self.cat,_extension),"Envelope",self.pdfIndex,pdflist)


  # Function to calculate the best-fit bkg model function and set at nominal function in envelope
  def getBestfit(self,verbose=True):
    bfIndex = 0
    NLL_min = 1e10
    snap, clean, aset = ROOT.RooArgSet(), ROOT.RooArgSet(), ROOT.RooArgSet()
    params = self.multipdf.getParameters(aset)
    params.remove(self.pdfIndex)
    params.snapshot(clean)
    params.snapshot(snap)
    # Loop over functions in envelope
    for i in range(self.pdfIndex.numTypes()):
      params.assignValueOnly(clean)
      self.pdfIndex.setIndex(i)
      # Run fit
      NLL, fitStatus = self.runFit(self.multipdf)
      # Add penalty term for number of free parameters in model
      NLL += self.multipdf.getCorrection()

      if verbose:
        print "  --> Function: %s"%self.multipdf.getCurrentPdf().GetName()
        print "    * Penalty term: %g"%self.multipdf.getCorrection()
        print "    * NLL + penalty = %.3f"%NLL

      # If NLL<NLL_min set bestfit
      if NLL < NLL_min:
        NLL_min = NLL
        snap.assignValueOnly(params)
        self.bestfitPdfIndex = i

    # Set pdfindex 
    self.pdfIndex.setIndex(self.bestfitPdfIndex)
    params.assignValueOnly(snap)
    self.bestfitSnapshot = params

    if verbose:
      print " --> Best-fit function: %s (index=%g)"%(self.multipdf.getCurrentPdf().GetName(),self.pdfIndex.getIndex())

  # Function to calculate norm of background function
  #def getNorm(self,index=None):
    # If index = None then use pre-calculated best-fit
    
  
    
      
