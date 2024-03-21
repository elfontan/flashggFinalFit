import ROOT
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
#import mplhep as hep
import argparse
import sys
import os

ROOT.gROOT.SetBatch()
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetOptStat(0)

argparser = argparse.ArgumentParser(description='Parser used for non default arguments', formatter_class=argparse.ArgumentDefaultsHelpFormatter, add_help=True)
argparser.add_argument('--o', dest='o', default='./', help='Output dir')
argparser.add_argument('--c', dest='c', default='cat0', help='Category')
argparser.add_argument('--w', dest='w', default='cat0', help='Workspace')

args = argparser.parse_args()
outdir = args.o
cat = args.c
wFile = args.w


def plot(mhs, vals, spline_name):
  plt.scatter(mhs, vals)
  plt.style.use('classic')
  plt.xlabel(r"$m_{\gamma\gamma}$ [GeV]", fontsize=18)
  plt.ylabel(spline_name, fontsize=14)

  if (cat == "cat0"):
    plt.plot([mhs[0], mhs[-1]], [vals[0], vals[-1]], color="royalblue", label='Category 0')
  elif (cat == "cat1"):
    plt.plot([mhs[0], mhs[-1]], [vals[0], vals[-1]], color="red", label='Category 1')

  #plt.style.use(hep.style.ROOT)
  #hep.cms.label("Preliminary", data=False)
  plt.title("CMS Simulation Preliminary", fontsize=22, loc='left')
  plt.title("13 TeV", fontsize=22, loc='right')
  #plt.title("CMS Simulation Preliminary                13 TeV", fontsize=14, loc='right')
  plt.legend(loc='lower right', frameon=False, fontsize=18)
  #plt.legend(loc='upper left', frameon=False, fontsize=18)
  #plt.legend(loc='upper center', frameon=False, fontsize=18)
  plt.savefig(outdir+"/%s.png"%spline_name)
  plt.clf()

#f = ROOT.TFile(sys.argv[1],"read")
f = ROOT.TFile(wFile,"read")
w = f.Get("wsig_13TeV")

if (cat == "cat0"):
  splines = ["mean_dcb_GG2H_2018_UntaggedTag_0_13TeV", "sigma_dcb_GG2H_2018_UntaggedTag_0_13TeV", 
             "a1_dcb_GG2H_2018_UntaggedTag_0_13TeV", "a2_dcb_GG2H_2018_UntaggedTag_0_13TeV", 
             "n1_dcb_GG2H_2018_UntaggedTag_0_13TeV", "n2_dcb_GG2H_2018_UntaggedTag_0_13TeV", 
             "sigma_gaus_GG2H_2018_UntaggedTag_0_13TeV", "frac_GG2H_2018_UntaggedTag_0_13TeV" 
             #, "dm_dcb_GG2H_2018_UntaggedTag_0_13TeV", "fea_GG2H_2018_UntaggedTag_0_13TeV"
                ]
elif (cat == "cat1"):
  splines = ["mean_dcb_GG2H_2018_UntaggedTag_1_13TeV", "sigma_dcb_GG2H_2018_UntaggedTag_1_13TeV", 
             "a1_dcb_GG2H_2018_UntaggedTag_1_13TeV", "a2_dcb_GG2H_2018_UntaggedTag_1_13TeV", 
             "n1_dcb_GG2H_2018_UntaggedTag_1_13TeV", "n2_dcb_GG2H_2018_UntaggedTag_1_13TeV", 
             "sigma_gaus_GG2H_2018_UntaggedTag_1_13TeV", "frac_GG2H_2018_UntaggedTag_1_13TeV" 
             #, "dm1_dcb_GG2H_2018_UntaggedTag_1_13TeV", "fea_GG2H_2018_UntaggedTag_1_13TeV"
           ]

mhs = [10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70]

for spline in splines:
  vals = []
  for mh in mhs:
    w.var("MH").setVal(mh)
    vals.append(w.function(spline).getVal())
  print(vals)
  plot(mhs, vals, spline)
