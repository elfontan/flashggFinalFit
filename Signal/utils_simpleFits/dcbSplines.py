import ROOT
import sys
import argparse
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def plot(mhs, vals, spline_name):
    plt.scatter(mhs, vals)
    plt.xlabel(r"$m_{\gamma\gamma}$")
    plt.ylabel(spline_name)

    plt.plot([mhs[0], mhs[-1]], [vals[0], vals[-1]])

    plt.savefig("%s.png" % spline_name)
    plt.clf()

parser = argparse.ArgumentParser(description='Plot splines for different categories.')
parser.add_argument('input_file', type=str, help='Path to the input ROOT file')
parser.add_argument('--category', type=int, choices=[0, 1], default=0, help='Category: 0 or 1')
args = parser.parse_args()

f = ROOT.TFile(args.input_file, "read")
w = f.Get("wsig_13TeV")

if args.category == 0:
    splines = ["mean_dcb_GG2H_2018_UntaggedTag_0_13TeV", "sigma_dcb_GG2H_2018_UntaggedTag_0_13TeV", 
               "a1_dcb_GG2H_2018_UntaggedTag_0_13TeV", "a2_dcb_GG2H_2018_UntaggedTag_0_13TeV", 
               "n1_dcb_GG2H_2018_UntaggedTag_0_13TeV", "n2_dcb_GG2H_2018_UntaggedTag_0_13TeV", 
               "dm_dcb_GG2H_2018_UntaggedTag_0_13TeV", "sigma_gaus_GG2H_2018_UntaggedTag_0_13TeV", 
               "frac_GG2H_2018_UntaggedTag_0_13TeV", #"fea_GG2H_2018_UntaggedTag_0_13TeV"
           ]
elif args.category == 1:
    splines = ["mean_dcb_GG2H_2018_UntaggedTag_1_13TeV", "sigma_dcb_GG2H_2018_UntaggedTag_1_13TeV", 
               "a1_dcb_GG2H_2018_UntaggedTag_1_13TeV", "a2_dcb_GG2H_2018_UntaggedTag_1_13TeV", 
               "n1_dcb_GG2H_2018_UntaggedTag_1_13TeV", "n2_dcb_GG2H_2018_UntaggedTag_1_13TeV", 
               "dm_dcb_GG2H_2018_UntaggedTag_1_13TeV", "sigma_gaus_GG2H_2018_UntaggedTag_1_13TeV", 
               "frac_GG2H_2018_UntaggedTag_1_13TeV", #"fea_GG2H_2018_UntaggedTag_1_13TeV"
           ]
else:
    print("Invalid category. Choose either 0 or 1.")
    sys.exit(1)

mhs = [10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70]

for spline in splines:
    vals = []
    for mh in mhs:
        w.var("MH").setVal(mh)
        vals.append(w.function(spline).getVal())
    print(vals)
    plot(mhs, vals, spline)
