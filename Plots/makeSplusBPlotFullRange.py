import ROOT
import json
import math
import pandas
import numpy as np
import re
import array
from ROOT import kBlue, kRed, kCyan, kGray

def SetErrors(binEnd,shift,scale,h1Epre,h1Erpre,h2Epre,h2Erpre,h1E,h1Er,h2E,h2Er):
  for i in range(0,binEnd):
    h1E.SetPoint(i,h1Epre.GetPointX(i+shift),scale*h1Epre.GetPointY(i+shift))
    h1E.SetPointError(i,h1Epre.GetErrorXlow(i+shift),h1Epre.GetErrorXhigh(i+shift),scale*h1Epre.GetErrorYlow(i+shift),scale*h1Epre.GetErrorYhigh(i+shift))
    h1Er.SetPoint(i,h1Erpre.GetPointX(i+shift),scale*h1Erpre.GetPointY(i+shift))
    h1Er.SetPointError(i,h1Erpre.GetErrorXlow(i+shift),h1Erpre.GetErrorXhigh(i+shift),scale*h1Erpre.GetErrorYlow(i+shift),scale*h1Erpre.GetErrorYhigh(i+shift))
    h2E.SetPoint(i,h2Epre.GetPointX(i+shift),scale*h2Epre.GetPointY(i+shift))
    h2E.SetPointError(i,h2Epre.GetErrorXlow(i+shift),h2Epre.GetErrorXhigh(i+shift),scale*h2Epre.GetErrorYlow(i+shift),scale*h2Epre.GetErrorYhigh(i+shift))
    h2Er.SetPoint(i,h2Erpre.GetPointX(i+shift),scale*h2Erpre.GetPointY(i+shift))
    h2Er.SetPointError(i,h2Erpre.GetErrorXlow(i+shift),h2Erpre.GetErrorXhigh(i+shift),scale*h2Erpre.GetErrorYlow(i+shift),scale*h2Erpre.GetErrorYhigh(i+shift))

#mass 12
hD1f = ROOT.TFile("./SplusBModels_AllData_cat0/BOnly/cat0_CMS_hgg_mass_12_hD.root","READ")
hD1 = hD1f.Get("h_data_cat0__CMS_hgg_mass")

hS1f = ROOT.TFile("./SplusBModels_AllData_cat0/BOnly/cat0_CMS_hgg_mass_12_hSB.root","READ")
hS1 = hS1f.Get("h_sb_pdfNBins_cat0__CMS_hgg_mass")
hSr1f = ROOT.TFile("./SplusBModels_AllData_cat0/BOnly/cat0_CMS_hgg_mass_12_hSr.root","READ")
hSr1 = hSr1f.Get("h_sb_pdfNBins_cat0__CMS_hgg_mass")

hB1f = ROOT.TFile("./SplusBModels_AllData_cat0/BOnly/cat0_CMS_hgg_mass_12_hB.root","READ")
hB1 = hB1f.Get("h_b_pdfNBins_cat0__CMS_hgg_mass")
hBr1f = ROOT.TFile("./SplusBModels_AllData_cat0/BOnly/cat0_CMS_hgg_mass_12_hBr.root","READ")
hBr1 = hBr1f.Get("h_b_pdfNBins_cat0__CMS_hgg_mass")

h1E1f = ROOT.TFile("./SplusBModels_AllData_cat0/BOnly/cat0_CMS_hgg_mass_12_1sig.root","READ")
h1E1pre = h1E1f.Get("gr_1sig")
h1E1 = ROOT.TGraphAsymmErrors()
h1Er1f = ROOT.TFile("./SplusBModels_AllData_cat0/BOnly/cat0_CMS_hgg_mass_12_1sig_r.root","READ")
h1Er1pre = h1Er1f.Get("gr_1sig_r")
h1Er1 = ROOT.TGraphAsymmErrors()

h2E1f = ROOT.TFile("./SplusBModels_AllData_cat0/BOnly/cat0_CMS_hgg_mass_12_2sig.root","READ")
h2E1pre = h2E1f.Get("gr_2sig")
h2E1 = ROOT.TGraphAsymmErrors()
h2Er1f = ROOT.TFile("./SplusBModels_AllData_cat0/BOnly/cat0_CMS_hgg_mass_12_2sig_r.root","READ")
h2Er1pre = h2Er1f.Get("gr_2sig_r")
h2Er1 = ROOT.TGraphAsymmErrors()

hDr1 = ROOT.TH1F("hDr1","hDr1",70,9.0,16.0)
hDr1.Add(hB1)
hDr1.Scale(-0.025)
hDr1.Add(hD1,1)

for h in [hD1,hDr1,hB1,hBr1,hS1,hSr1]:
  h.GetXaxis().SetRangeUser(10.0,14.5)
  h.Scale(1.0/h.GetBinWidth(0))

SetErrors(46,10,10,h1E1pre,h1Er1pre,h2E1pre,h2Er1pre,h1E1,h1Er1,h2E1,h2Er1)

#mass 33
hD2f = ROOT.TFile("./SplusBModels_AllData_cat0/BOnly/cat0_CMS_hgg_mass_33_hD.root","READ")
hD2 = hD2f.Get("h_data_cat0__CMS_hgg_mass")

hS2f = ROOT.TFile("./SplusBModels_AllData_cat0/BOnly/cat0_CMS_hgg_mass_33_hSB.root","READ")
hS2 = hS2f.Get("h_sb_pdfNBins_cat0__CMS_hgg_mass")
hSr2f = ROOT.TFile("./SplusBModels_AllData_cat0/BOnly/cat0_CMS_hgg_mass_33_hSr.root","READ")
hSr2 = hSr2f.Get("h_sb_pdfNBins_cat0__CMS_hgg_mass")

hB2f = ROOT.TFile("./SplusBModels_AllData_cat0/BOnly/cat0_CMS_hgg_mass_33_hB.root","READ")
hB2 = hB2f.Get("h_b_pdfNBins_cat0__CMS_hgg_mass")
hBr2f = ROOT.TFile("./SplusBModels_AllData_cat0/BOnly/cat0_CMS_hgg_mass_33_hBr.root","READ")
hBr2 = hBr2f.Get("h_b_pdfNBins_cat0__CMS_hgg_mass")

h1E2f = ROOT.TFile("./SplusBModels_AllData_cat0/BOnly/cat0_CMS_hgg_mass_33_1sig.root","READ")
h1E2pre = h1E2f.Get("gr_1sig")
h1E2 = ROOT.TGraphAsymmErrors()
h1Er2f = ROOT.TFile("./SplusBModels_AllData_cat0/BOnly/cat0_CMS_hgg_mass_33_1sig_r.root","READ")
h1Er2pre = h1Er2f.Get("gr_1sig_r")
h1Er2 = ROOT.TGraphAsymmErrors()

h2E2f = ROOT.TFile("./SplusBModels_AllData_cat0/BOnly/cat0_CMS_hgg_mass_33_2sig.root","READ")
h2E2pre = h2E2f.Get("gr_2sig")
h2E2 = ROOT.TGraphAsymmErrors()
h2Er2f = ROOT.TFile("./SplusBModels_AllData_cat0/BOnly/cat0_CMS_hgg_mass_33_2sig_r.root","READ")
h2Er2pre = h2Er2f.Get("gr_2sig_r")
h2Er2 = ROOT.TGraphAsymmErrors()

hDr2 = ROOT.TH1F("hDr2","hDr2",155,13.0,44.0)
hDr2.Add(hB2)
hDr2.Scale(-0.025)
hDr2.Add(hD2,1)

for h in [hD2,hDr2,hB2,hBr2,hS2,hSr2]:
  h.GetXaxis().SetRangeUser(14.5,40.0)
  h.Scale(1.0/h.GetBinWidth(0))

SetErrors(129,7,5,h1E2pre,h1Er2pre,h2E2pre,h2Er2pre,h1E2,h1Er2,h2E2,h2Er2)


#mass 50
hD3f = ROOT.TFile("./SplusBModels_AllData_cat0/BOnly/cat0_CMS_hgg_mass_50_hD.root","READ")
hD3 = hD3f.Get("h_data_cat0__CMS_hgg_mass")

hS3f = ROOT.TFile("./SplusBModels_AllData_cat0/BOnly/cat0_CMS_hgg_mass_50_hSB.root","READ")
hS3 = hS3f.Get("h_sb_pdfNBins_cat0__CMS_hgg_mass")
hSr3f = ROOT.TFile("./SplusBModels_AllData_cat0/BOnly/cat0_CMS_hgg_mass_50_hSr.root","READ")
hSr3 = hSr3f.Get("h_sb_pdfNBins_cat0__CMS_hgg_mass")

hB3f = ROOT.TFile("./SplusBModels_AllData_cat0/BOnly/cat0_CMS_hgg_mass_50_hB.root","READ")
hB3 = hB3f.Get("h_b_pdfNBins_cat0__CMS_hgg_mass")
hBr3f = ROOT.TFile("./SplusBModels_AllData_cat0/BOnly/cat0_CMS_hgg_mass_50_hBr.root","READ")
hBr3 = hBr3f.Get("h_b_pdfNBins_cat0__CMS_hgg_mass")

h1E3f = ROOT.TFile("./SplusBModels_AllData_cat0/BOnly/cat0_CMS_hgg_mass_50_1sig.root","READ")
h1E3pre = h1E3f.Get("gr_1sig")
h1E3 = ROOT.TGraphAsymmErrors()
h1Er3f = ROOT.TFile("./SplusBModels_AllData_cat0/BOnly/cat0_CMS_hgg_mass_50_1sig_r.root","READ")
h1Er3pre = h1Er3f.Get("gr_1sig_r")
h1Er3 = ROOT.TGraphAsymmErrors()

h2E3f = ROOT.TFile("./SplusBModels_AllData_cat0/BOnly/cat0_CMS_hgg_mass_50_2sig.root","READ")
h2E3pre = h2E3f.Get("gr_2sig")
h2E3 = ROOT.TGraphAsymmErrors()
h2Er3f = ROOT.TFile("./SplusBModels_AllData_cat0/BOnly/cat0_CMS_hgg_mass_50_2sig_r.root","READ")
h2Er3pre = h2Er3f.Get("gr_2sig_r")
h2Er3 = ROOT.TGraphAsymmErrors()

hDr3 = ROOT.TH1F("hDr3","hDr3",60,36.0,60.0)
hDr3.Add(hB3)
hDr3.Scale(-0.025)
hDr3.Add(hD3,1)

for h in [hD3,hDr3,hB3,hBr3,hS3,hSr3]:
  h.GetXaxis().SetRangeUser(40.0,54.5)
  h.Scale(1.0/h.GetBinWidth(0))

SetErrors(37,9,2.5,h1E3pre,h1Er3pre,h2E3pre,h2Er3pre,h1E3,h1Er3,h2E3,h2Er3)

#mass 68
hD4f = ROOT.TFile("./SplusBModels_AllData_cat0/BOnly/cat0_CMS_hgg_mass_68_hD.root","READ")
hD4 = hD4f.Get("h_data_cat0__CMS_hgg_mass")

hS4f = ROOT.TFile("./SplusBModels_AllData_cat0/BOnly/cat0_CMS_hgg_mass_68_hSB.root","READ")
hS4 = hS4f.Get("h_sb_pdfNBins_cat0__CMS_hgg_mass")
hSr4f = ROOT.TFile("./SplusBModels_AllData_cat0/BOnly/cat0_CMS_hgg_mass_68_hSr.root","READ")
hSr4 = hSr4f.Get("h_sb_pdfNBins_cat0__CMS_hgg_mass")

hB4f = ROOT.TFile("./SplusBModels_AllData_cat0/BOnly/cat0_CMS_hgg_mass_68_hB.root","READ")
hB4 = hB4f.Get("h_b_pdfNBins_cat0__CMS_hgg_mass")
hBr4f = ROOT.TFile("./SplusBModels_AllData_cat0/BOnly/cat0_CMS_hgg_mass_68_hBr.root","READ")
hBr4 = hBr4f.Get("h_b_pdfNBins_cat0__CMS_hgg_mass")

h1E4f = ROOT.TFile("./SplusBModels_AllData_cat0/BOnly/cat0_CMS_hgg_mass_68_1sig.root","READ")
h1E4pre = h1E4f.Get("gr_1sig")
h1E4 = ROOT.TGraphAsymmErrors()
h1Er4f = ROOT.TFile("./SplusBModels_AllData_cat0/BOnly/cat0_CMS_hgg_mass_68_1sig_r.root","READ")
h1Er4pre = h1Er4f.Get("gr_1sig_r")
h1Er4 = ROOT.TGraphAsymmErrors()

h2E4f = ROOT.TFile("./SplusBModels_AllData_cat0/BOnly/cat0_CMS_hgg_mass_68_2sig.root","READ")
h2E4pre = h2E4f.Get("gr_2sig")
h2E4 = ROOT.TGraphAsymmErrors()
h2Er4f = ROOT.TFile("./SplusBModels_AllData_cat0/BOnly/cat0_CMS_hgg_mass_68_2sig_r.root","READ")
h2Er4pre = h2Er4f.Get("gr_2sig_r")
h2Er4 = ROOT.TGraphAsymmErrors()

hDr4 = ROOT.TH1F("hDr4","hDr4",56,49.0,77.0)
hDr4.Add(hB4)
hDr4.Scale(-0.025)
hDr4.Add(hD4,1)

for h in [hD4,hDr4,hB4,hBr4,hS4,hSr4]:
  h.GetXaxis().SetRangeUser(54.5,70.0)
  h.Scale(1.0/h.GetBinWidth(0))

SetErrors(33,10,2,h1E4pre,h1Er4pre,h2E4pre,h2Er4pre,h1E4,h1Er4,h2E4,h2Er4)

hDf = ROOT.TFile("./SplusBModels_AllData_cat0/BOnly/data_UntaggedTag_0_FullRange.root","READ")
hD = hDf.Get("h_data__CMS_hgg_mass")
hDr = hD.Clone("hDr")

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

canv = ROOT.TCanvas("canv","canv",1800,1200)
#canv = ROOT.TCanvas("canv","canv",1800,1200)
pad1 = ROOT.TPad("pad1","pad1",0,0.25,1,1)
pad1.SetTickx()
pad1.SetTicky()
pad1.SetBottomMargin(0.17)
pad1.SetLeftMargin(0.12)
pad1.SetRightMargin(0.04)
pad1.Draw()

pad2 = ROOT.TPad("pad2","pad2",0,0,1,0.35)
pad2.SetTickx()
pad2.SetTicky()
pad2.SetTopMargin(0.01)
pad2.SetBottomMargin(0.30)
pad2.SetLeftMargin(0.12)
pad2.SetRightMargin(0.04)
pad2.Draw()
padSizeRatio = 0.75/0.35

# Axis options 
ROOT.TGaxis.SetMaxDigits(4)
ROOT.TGaxis.SetExponentOffset(-0.05,0.00,"y")

edges=array.array('f',[10.0, 10.1, 10.2, 10.3, 10.4, 10.5, 10.6, 10.7, 10.8, 10.9, 11.0, 11.1, 11.2, 11.3, 11.4, 11.5, 11.6, 11.7, 11.8, 11.9, 12.0, 12.1, 12.2, 12.3, 12.4, 12.5, 12.6, 12.7, 12.8, 12.9, 13.0, 13.1, 13.2, 13.3, 13.4, 13.5, 13.6, 13.7, 13.8, 13.9, 14.0, 14.1, 14.2, 14.3, 14.4, 14.6, 14.8, 15.0, 15.2, 15.4, 15.6, 15.8, 16.0, 16.2, 16.4, 16.6, 16.8, 17.0, 17.2, 17.4, 17.6, 17.8, 18.0, 18.2, 18.4, 18.6, 18.8, 19.0, 19.2, 19.4, 19.6, 19.8, 20.0, 20.2, 20.4, 20.6, 20.8, 21.0, 21.2, 21.4, 21.6, 21.8, 22.0, 22.2, 22.4, 22.6, 22.8, 23.0, 23.2, 23.4, 23.6, 23.8, 24.0, 24.2, 24.4, 24.6, 24.8, 25.0, 25.2, 25.4, 25.6, 25.8, 26.0, 26.2, 26.4, 26.6, 26.8, 27.0, 27.2, 27.4, 27.6, 27.8, 28.0, 28.2, 28.4, 28.6, 28.8, 29.0, 29.2, 29.4, 29.6, 29.8, 30.0, 30.2, 30.4, 30.6, 30.8, 31.0, 31.2, 31.4, 31.6, 31.8, 32.0, 32.2, 32.4, 32.6, 32.8, 33.0, 33.2, 33.4, 33.6, 33.8, 34.0, 34.2, 34.4, 34.6, 34.8, 35.0, 35.2, 35.4, 35.6, 35.8, 36.0, 36.2, 36.4, 36.6, 36.8, 37.0, 37.2, 37.4, 37.6, 37.8, 38.0, 38.2, 38.4, 38.6, 38.8, 39.0, 39.2, 39.4, 39.6, 39.8, 40.0, 40.4, 40.8, 41.2, 41.6, 42.0, 42.4, 42.8, 43.2, 43.6, 44.0, 44.4, 44.8, 45.2, 45.6, 46.0, 46.4, 46.8, 47.2, 47.6, 48.0, 48.4, 48.8, 49.2, 49.6, 50.0, 50.4, 50.8, 51.2, 51.6, 52.0, 52.4, 52.8, 53.2, 53.6, 54.0, 54.5, 55.0, 55.5, 56.0, 56.5, 57.0, 57.5, 58.0, 58.5, 59.0, 59.5, 60.0, 60.5, 61.0, 61.5, 62.0, 62.5, 63.0, 63.5, 64.0, 64.5, 65.0, 65.5, 66.0, 66.5, 67.0, 67.5, 68.0, 68.5, 69.0, 69.5, 70.0])

# Nominal plot
pad1.cd()
h_axes = ROOT.TH1F("h_axes","h_axes",len(edges)-1,edges)
h_axes.SetMaximum(5500.0)
h_axes.SetMinimum(0.)
h_axes.SetTitle("")
#h_axes.GetXaxis().SetRangeUser(10.0,70.0)
h_axes.GetXaxis().SetTitle("")
h_axes.GetXaxis().SetLabelSize(0)
h_axes.GetYaxis().SetTitleSize(0.07)
h_axes.GetYaxis().SetTitleOffset(0.86)
h_axes.GetYaxis().SetLabelSize(0.060)
h_axes.GetYaxis().SetLabelOffset(0.011)
h_axes.GetYaxis().SetTitle("< Events / GeV >")
h_axes.Draw()

lines_masses = [14.5, 40, 54.5]
line1 = ROOT.TLine(14.5, 0.0, 14.5, 3250.0)
line1.SetLineColor(kRed-7)
line1.SetLineStyle(7)
line1.SetLineWidth(2)
line1.Draw("same")
line2 = ROOT.TLine(40.0, 0.0, 40.0, 2700.0)
line2.SetLineColor(kRed-7)
line2.SetLineStyle(7)
line2.SetLineWidth(2)
line2.Draw("same")
line3 = ROOT.TLine(54.5, 0.0, 54.5, 2700.0)
line3.SetLineColor(kRed-7)
line3.SetLineStyle(7)
line3.SetLineWidth(2)
line3.Draw("same")

# Add bands
col_sig1 = ROOT.TColor.GetColor("#85D1FB")
col_sig2 = ROOT.TColor.GetColor("#FFDF7F")
for h2E in [h2E1,h2E2,h2E3,h2E4]:
  h2E.SetFillColor(col_sig2)
  h2E.SetFillStyle(1001)
  h2E.Draw("LE3SAME")
for h1E in [h1E1,h1E2,h1E3,h1E4]:
  h1E.SetFillColor(col_sig1)
  h1E.SetFillStyle(1001)
  h1E.Draw("LE3SAME")
# Set pdf style
#for hS in [hS1,hS2,hS3,hS4]:
#  hS.SetLineWidth(2)
#  hS.SetLineColor(kBlue)
#  hS.Scale(0.025)
#  hS.Draw("Hist same ][")
#  print hS.GetMaximum()
for hB in [hB1,hB2,hB3,hB4]:
  hB.SetLineWidth(1)
  hB.SetLineColor(kBlue)
  hB.Scale(0.025)
  hB.SetLineStyle(1) #EF
  hB.Draw("Hist same ][")
  print(hB.GetMaximum())
# Set data style
for hD in [hD1,hD2,hD3,hD4]:
  hD.SetMarkerStyle(20)
  hD.SetMarkerColor(1)
  hD.SetLineColor(1)
  hD.Draw("Same PE")

latex = ROOT.TLatex()
latex.SetTextFont(42)
latex.SetTextAlign(22)

window_labels = [("W1", 12.2), ("W2", 28.), ("W3", 47.), ("W4", 62.2)]
for label, pos in window_labels:
  latex.SetTextColor(kRed-7)
  latex.SetTextSize(0.055)
  latex.DrawLatex(pos, 350.0, label)

# Add legend
leg = ROOT.TLegend(0.56,0.53,0.91,0.86)
leg.SetFillStyle(1001)
leg.SetFillColor(0)
leg.SetLineColor(0)
leg.SetTextSize(0.06)
#leg.SetHeader("#scale[1.25]{H #rightarrow #gamma#gamma}", "C")
hD1_leg = hD1.Clone("hD1_leg")
hD1_leg.SetMarkerSize(2) 
leg.AddEntry(hD1_leg,"Data","ep")
#leg.AddEntry(hS1,"S+B fit","l")
leg.AddEntry(hB1,"B component of S+B fit","l")
leg.AddEntry(h1E1,"#pm1 #sigma (B-only)","F")
leg.AddEntry(h2E1,"#pm2 #kern[-0.26]{#sigma} (B-only)","F")
leg.Draw("Same")

# Add TLatex to plot
lat0 = ROOT.TLatex()
lat0.SetTextFont(42)
lat0.SetTextAlign(11)
lat0.SetNDC()
lat0.SetTextSize(0.06)
lat0.DrawLatex(0.78,0.92,"54.4 fb^{-1} (13 TeV)")
#lat0.DrawLatex(0.16,0.70,"#scale[0.75]{H #rightarrow #gamma#gamma}")

lat1 = ROOT.TLatex()
lat1.SetTextFont(61)
lat1.SetTextAlign(11)
lat1.SetNDC()
lat1.SetTextSize(0.08)
lat1.DrawLatex(0.15,0.81,"CMS")

lat2 = ROOT.TLatex()
lat2.SetTextFont(52)
lat2.SetTextAlign(11)
lat2.SetNDC()
lat2.SetTextSize(0.06)
lat2.DrawLatex(0.24,0.81,"")
#lat2.DrawLatex(0.24,0.81,"Preliminary")
#lat2.DrawLatex(0.23,0.81,"Work in progress")


pad1.Update()

# Ratio plot
pad2.cd()
h_axes_ratio = h_axes.Clone()
#h_axes_ratio.Reset()
h_axes_ratio.SetMaximum(500.0)
h_axes_ratio.SetMinimum(-500.0)
h_axes_ratio.SetTitle("")
h_axes_ratio.GetXaxis().SetRangeUser(10.0,70.0)
h_axes_ratio.GetXaxis().SetTitleSize(0.07*padSizeRatio)
h_axes_ratio.GetXaxis().SetTitleOffset(0.92)
h_axes_ratio.GetXaxis().SetLabelSize(0.06*padSizeRatio)
h_axes_ratio.GetXaxis().SetLabelOffset(0.020)
h_axes_ratio.GetXaxis().SetTickLength(0.03*padSizeRatio)
h_axes_ratio.GetXaxis().SetTitle("m_{#gamma#gamma} [GeV]")
h_axes_ratio.GetYaxis().SetTitleSize(0.07*padSizeRatio)
h_axes_ratio.GetYaxis().SetTitleOffset(0.40)
h_axes_ratio.GetYaxis().SetLabelSize(0.06*padSizeRatio)
h_axes_ratio.GetYaxis().SetLabelOffset(0.011)
h_axes_ratio.GetYaxis().SetTitle("Data - Bkg.")
h_axes_ratio.GetYaxis().SetNdivisions(305,True)
#h_axes_ratio.GetYaxis().SetNdivisions(205,True)
#h_axes_ratio.GetYaxis().ChangeLabel(5,-1,0,-1,-1,-1,"")

h_axes_ratio.Draw()

line1r = ROOT.TLine(14.5, -500.0, 14.5, 500.0)
line1r.SetLineColor(kRed-7)
line1r.SetLineStyle(7)
line1r.SetLineWidth(2)
line1r.Draw("same")
line2r = ROOT.TLine(40.0, -500.0, 40.0, 500.0)
line2r.SetLineColor(kRed-7)
line2r.SetLineStyle(7)
line2r.SetLineWidth(2)
line2r.Draw("same")
line3r = ROOT.TLine(54.5, -500.0, 54.5, 500.0)
line3r.SetLineColor(kRed-7)
line3r.SetLineStyle(7)
line3r.SetLineWidth(2)
line3r.Draw("same")

# Draw bands 
for h2Er in [h2Er1,h2Er2,h2Er3,h2Er4]:
  h2Er.SetFillColor(col_sig2)
  h2Er.SetFillStyle(1001)
  h2Er.Draw("LE3SAME")
for h1Er in [h1Er1,h1Er2,h1Er3,h1Er4]:
  h1Er.SetFillColor(col_sig1)
  h1Er.SetFillStyle(1001)
  h1Er.Draw("LE3SAME")
# Set pdf style
#for hSr in [hSr1,hSr2,hSr3,hSr4]:
#  hSr.SetLineWidth(2)
#  hSr.SetLineColor(kBlue)
#  hSr.Draw("Hist same ][")
for hBr in [hBr1,hBr2,hBr3,hBr4]:
  hBr.SetLineWidth(1)
  hBr.SetLineStyle(1)
  hBr.SetLineColor(kBlue)
  hBr.Draw("Hist same ][")
# Set data style
for hDr in [hDr1,hDr2,hDr3,hDr4]:
  hDr.SetMarkerStyle(20)
  hDr.SetMarkerColor(1)
  hDr.SetLineColor(1)
  hDr.Draw("Same PE")

#for label, pos in window_labels:
#  latex.SetTextSize(0.06)
#  latex.DrawLatex(pos, -425.0, label)

# Add TLatex to ratio plot
lat3 = ROOT.TLatex()
lat3.SetTextFont(42)
lat3.SetTextAlign(33)
lat3.SetNDC(1)
lat3.SetTextSize(0.060*padSizeRatio)
#lat3.DrawLatex(0.94,0.45,"B component subtracted from data")

pad2.Update()

# Save canvas
canv.Update()
canv.SaveAs("/eos/user/e/elfontan/www/Hgg_veryLowMass_Paper/SplusBFits/2_bkgModeling_CWR.C")
canv.SaveAs("/eos/user/e/elfontan/www/Hgg_veryLowMass_Paper/SplusBFits/2_bkgModeling_CWR.png")
canv.SaveAs("/eos/user/e/elfontan/www/Hgg_veryLowMass_Paper/SplusBFits/2_bkgModeling_CWR.pdf")
