# Script to hold definitions of bkg functions
from collections import OrderedDict as od
import ROOT

functionFamilies = od()

# Bernstein polynomial
#functionFamilies['Bernstein'] = od()
#functionFamilies['Bernstein']['name'] = ['Bernstein','bern']

# Exponential
functionFamilies['Exponential'] = od()
functionFamilies['Exponential']['name'] = ['Exponential','exp']

# Power law
#functionFamilies['PowerLaw'] = od()
#functionFamilies['PowerLaw']['name'] = ['PowerLaw','pow']

# Laurent series
#functionFamilies['Laurent'] = od()
#functionFamilies['Laurent']['name'] = ['Laurent','lau']


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
    if order==1: return ROOT.RooBernsteinFast<1>(prefix,prefix,model.xvar,paramList)
    elif order==2: return ROOT.RooBernsteinFast<2>(prefix,prefix,model.xvar,paramList)
    elif order==3: return ROOT.RooBernsteinFast<3>(prefix,prefix,model.xvar,paramList)
    elif order==4: return ROOT.RooBernsteinFast<4>(prefix,prefix,model.xvar,paramList)
    elif order==5: return ROOT.RooBernsteinFast<5>(prefix,prefix,model.xvar,paramList)
    elif order==6: return ROOT.RooBernsteinFast<6>(prefix,prefix,model.xvar,paramList)
    else: return False
      
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
        ename = "%s_e%g"%(prefix,i)
        model.params[pname] = ROOT.RooRealVar(pname,pname,max(-1.,-0.04*(i+1)),-1.,0.)
        model.functions[ename] = ROOT.RooExponential(ename,ename,model.xvar,model.params[pname])
        exps.add(model.functions[ename])
      # Add up exponentials
      return ROOT.RooAddPdf(prefix,prefix,exps,fracs,True)
        
        

