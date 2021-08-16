# Script to hold definitions of bkg functions
from collections import OrderedDict as od
import ROOT
import math
import ctypes

functionFamilies = od()

# Bernstein polynomial
#functionFamilies['Bernstein'] = od()
#functionFamilies['Bernstein']['name'] = ['Bernstein','bern']

# Exponential
functionFamilies['Exponential'] = od()
functionFamilies['Exponential']['name'] = ['Exponential','exp']

# Power law
functionFamilies['PowerLaw'] = od()
functionFamilies['PowerLaw']['name'] = ['PowerLaw','pow']

# Laurent series
functionFamilies['Laurent'] = od()
functionFamilies['Laurent']['name'] = ['Laurent','lau']


# Define function families
def getPdf(model,prefix,funcType,order):
 
  # For bernstein polynomial...
  if funcType == "Bernstein":
    paramList = ROOT.RooArgList()
    for i in range(order):
      pname = "%s_p%g"%(prefix,i)
      p = ROOT.RooRealVar(pname,pname,0.1*(i+1),-15.,15.)
      f = ROOT.RooFormulaVar("%s_sq"%pname,"%s_sq"%pname,"@0*@0",ROOT.RooArgList(p))
      # Add params (and functions) to model
      model.params[pname] = p
      model.formulas["%s_sq"%pname] = f
      paramList.add(model.formulas["%s_sq"%pname])
    # Make Bernsteins
    return ROOT.RooBernsteinFast(order)(prefix,prefix,model.xvar,paramList)
      
  # For Exponential
  if funcType == "Exponential":
    # Only odd orders allowed
    if order%2==0: return False
    else:
      nFracs = (order-1)/2
      nExps = order-nFracs
      assert(nFracs==nExps-1)
      fracs = ROOT.RooArgList()
      exps = ROOT.RooArgList()
      for i in range(1,nFracs+1):
        fname = "%s_frac%g"%(prefix,i)
        model.params[fname] = ROOT.RooRealVar(fname,fname,0.9-float(i-1)*1./nFracs,0.0001,0.9999)
        fracs.add(model.params[fname])
      for i in range(1,nExps+1):
        pname = "%s_p%g"%(prefix,i)
        funcname = "%s_e%g"%(prefix,i)
        model.params[pname] = ROOT.RooRealVar(pname,pname,max(-1.,-0.04*(i+1)),-1.,0.)
        model.functions[funcname] = ROOT.RooExponential(funcname,funcname,model.xvar,model.params[pname])
        exps.add(model.functions[funcname])
      # Add up exponentials
      return ROOT.RooAddPdf(prefix,prefix,exps,fracs,True)

  # For Power Law
  elif funcType == "PowerLaw":
    # Only odd orders allowed
    if order%2==0: return False
    else:
      nFracs = (order-1)/2
      nPows = order-nFracs
      assert(nFracs==nPows-1)
      fracs = ROOT.RooArgList()
      pows = ROOT.RooArgList()
      for i in range(1,nFracs+1):
        fname = "%s_frac%g"%(prefix,i)
        model.params[fname] = ROOT.RooRealVar(fname,fname,0.9-float(i-1)*1./nFracs,0.0001,0.9999)
        fracs.add(model.params[fname])
      for i in range(1,nPows+1):
        pname = "%s_p%g"%(prefix,i)
        funcname = "%s_e%g"%(prefix,i)
        model.params[pname] = ROOT.RooRealVar(pname,pname,max(-9.,-1.*(i+1)),-9.,1.)
        model.functions[funcname] = ROOT.RooPower(funcname,funcname,model.xvar,model.params[pname])
        pows.add(model.functions[funcname])
      # Add up powers
      return ROOT.RooAddPdf(prefix,prefix,pows,fracs,True)

  # For Laurent series
  elif funcType == "Laurent":
    nLower = int(math.ceil(order/2.)) 
    nHigher = order-nLower
    pows = ROOT.RooArgList()
    plist = ROOT.RooArgList()
    # 0-th order
    funcname = "%s_pow0"%prefix
    model.functions[funcname] = ROOT.RooPower(funcname,funcname,model.xvar,ROOT.RooFit.RooConst(-4.))
    pows.add(model.functions[funcname])
    # Even terms
    for i in range(1,nLower+1):
      pname = "%s_l%g"%(prefix,i)
      model.params[pname] = ROOT.RooRealVar(pname,pname,0.25/order,0.000001,0.999999)
      plist.add(model.params[pname])
      funcname = "%s_powl%g"%(prefix,i)
      model.functions[funcname] = ROOT.RooPower(funcname,funcname,model.xvar,ROOT.RooFit.RooConst(-4.-i))
      pows.add(model.functions[funcname])
    # Odd terms
    for i in range(1,nHigher+1):
      pname = "%s_h%g"%(prefix,i)
      model.params[pname] = ROOT.RooRealVar(pname,pname,0.25/order,0.000001,0.999999)
      plist.add(model.params[pname])
      funcname = "%s_powh%g"%(prefix,i)
      model.functions[funcname] = ROOT.RooPower(funcname,funcname,model.xvar,ROOT.RooFit.RooConst(-4.+i))
      pows.add(model.functions[funcname])
    # Add up terms
    return ROOT.RooAddPdf(prefix,prefix,pows,plist,True)

