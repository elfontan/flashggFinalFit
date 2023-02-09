import ROOT
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def plot(mhs, vals, spline_name):
  plt.scatter(mhs, vals)
  plt.xlabel(r"$m_H$")
  plt.ylabel(spline_name)

  plt.plot([mhs[0], mhs[-1]], [vals[0], vals[-1]])

  plt.savefig("%s.png"%spline_name)
  plt.clf()

f = ROOT.TFile(sys.argv[1],"read")
w = f.Get("wsig_13TeV")

splines = ["frac_g0_GG2H_2018_UntaggedTag_0_13TeV","frac_g1_GG2H_2018_UntaggedTag_0_13TeV", "dm_g0_GG2H_2018_UntaggedTag_0_13TeV", "dm_g1_GG2H_2018_UntaggedTag_0_13TeV", "dm_g2_GG2H_2018_UntaggedTag_0_13TeV", "sigma_fit_g0_GG2H_2018_UntaggedTag_0_13TeV", "sigma_fit_g1_GG2H_2018_UntaggedTag_0_13TeV", "sigma_fit_g2_GG2H_2018_UntaggedTag_0_13TeV", "fea_GG2H_2018_UntaggedTag_0_13TeV"]

#mhs = [5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70]
mhs = [5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70]
#mhs = [15, 20, 25, 30, 35, 40, 45, 50, 60, 65, 70]
#import numpy as np
#mhs = np.linspace(10, 65, 100)

for spline in splines:
  vals = []
  for mh in mhs:
    w.var("MH").setVal(mh)
    vals.append(w.function(spline).getVal())
  print(vals)
  plot(mhs, vals, spline)





