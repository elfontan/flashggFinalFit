#!/usr/bin/env python3

import argparse
import os
from array import array

import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot parametric signal model shapes from a saved RooWorkspace."
    )
    parser.add_argument(
        "--input",
        type=str,
        required=True,
        help="Input ROOT file containing the RooWorkspace.",
    )
    parser.add_argument(
        "--workspace",
        type=str,
        default=None,
        help="Workspace name. If omitted, auto-detect the first RooWorkspace in the file.",
    )
    parser.add_argument(
        "--pdf",
        type=str,
        default=None,
        help="PDF name. If omitted, auto-detect the first pdf ending in '_pdf'.",
    )
    parser.add_argument(
        "--xvar",
        type=str,
        default="mjj",
        help="Observable name in the workspace.",
    )
    parser.add_argument(
        "--mhvar",
        type=str,
        default="MH",
        help="Mass-hypothesis variable name in the workspace.",
    )
    parser.add_argument(
        "--outdir",
        type=str,
        default="plots_parametric_workspace_scan",
        help="Output directory.",
    )
    parser.add_argument(
        "--outname",
        type=str,
        default="parametric_workspace_scan",
        help="Output file base name.",
    )
    parser.add_argument(
        "--mhMin",
        type=float,
        default=400.0,
        help="Minimum MH value to plot.",
    )
    parser.add_argument(
        "--mhMax",
        type=float,
        default=2000.0,
        help="Maximum MH value to plot.",
    )
    parser.add_argument(
        "--mhStep",
        type=float,
        default=200.0,
        help="Step in MH.",
    )
    parser.add_argument(
        "--xmin",
        type=float,
        default=None,
        help="Optional x-axis minimum override.",
    )
    parser.add_argument(
        "--xmax",
        type=float,
        default=3500.0,
        help="x-axis maximum override. Default: 3500 GeV.",
    )
    parser.add_argument(
        "--nbins",
        type=int,
        default=800,
        help="Number of sampling bins.",
    )
    parser.add_argument(
        "--label",
        type=str,
        default="Parametric VBF signal model",
        help="Label drawn on canvas.",
    )
    parser.add_argument(
        "--noUnitNorm",
        action="store_true",
        help="Do not normalize each sampled curve to unit area.",
    )
    return parser.parse_args()


def set_style():
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)
    ROOT.gStyle.SetTitleBorderSize(0)
    ROOT.gStyle.SetLegendBorderSize(0)
    ROOT.gStyle.SetFrameLineWidth(2)
    ROOT.gStyle.SetLineWidth(2)
    ROOT.gStyle.SetPadLeftMargin(0.14)
    ROOT.gStyle.SetPadRightMargin(0.05)
    ROOT.gStyle.SetPadBottomMargin(0.13)
    ROOT.gStyle.SetPadTopMargin(0.08)
    ROOT.gStyle.SetLabelSize(0.04, "XYZ")
    ROOT.gStyle.SetTitleSize(0.05, "XYZ")
    ROOT.gStyle.SetTitleOffset(1.05, "X")
    ROOT.gStyle.SetTitleOffset(1.35, "Y")


def draw_cms_label():
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextFont(62)
    latex.SetTextSize(0.050)
    latex.DrawLatex(0.145, 0.935, "CMS")
    latex.SetTextFont(52)
    latex.SetTextSize(0.038)
    latex.DrawLatex(0.235, 0.935, "Simulation Preliminary")
    latex.SetTextFont(42)
    latex.SetTextSize(0.038)
    latex.SetTextAlign(31)
    latex.DrawLatex(0.95, 0.935, "2024 (13.6 TeV)")


def find_workspace(root_file, requested_name=None):
    if requested_name:
        ws = root_file.Get(requested_name)
        if not ws:
            raise RuntimeError(f"Workspace '{requested_name}' not found in {root_file.GetName()}")
        return ws

    for key in root_file.GetListOfKeys():
        obj = key.ReadObj()
        if obj.InheritsFrom("RooWorkspace"):
            return obj
    raise RuntimeError(f"No RooWorkspace found in {root_file.GetName()}")


def find_pdf(ws, requested_name=None):
    if requested_name:
        pdf = ws.pdf(requested_name)
        if not pdf:
            raise RuntimeError(f"PDF '{requested_name}' not found in workspace '{ws.GetName()}'")
        return pdf

    all_pdfs = ws.allPdfs()
    iterator = all_pdfs.createIterator()
    chosen = None
    while True:
        obj = iterator.Next()
        if not obj:
            break
        if obj.GetName().endswith("_pdf"):
            chosen = obj
            break
        if chosen is None:
            chosen = obj
    if not chosen:
        raise RuntimeError(f"No RooAbsPdf found in workspace '{ws.GetName()}'")
    return chosen


def build_mh_values(mh_min, mh_max, mh_step):
    values = []
    current = mh_min
    while current <= mh_max + 1e-6:
        values.append(float(current))
        current += mh_step
    return values


def sample_pdf(pdf, xvar, mhvar, mh_value, nbins, xmin, xmax, unit_norm=True):
    mhvar.setVal(float(mh_value))
    hist = ROOT.TH1D(f"h_mh_{int(round(mh_value))}", "", nbins, xmin, xmax)
    hist.Sumw2()
    norm_set = ROOT.RooArgSet(xvar)
    for ibin in range(1, nbins + 1):
        xvar.setVal(hist.GetBinCenter(ibin))
        yval = max(0.0, float(pdf.getVal(norm_set)))
        hist.SetBinContent(ibin, yval)
    if unit_norm and hist.Integral() > 0.0:
        hist.Scale(1.0 / hist.Integral("width"))
    return hist


def color_sequence():
    return [
        ROOT.kAzure + 1,
        ROOT.kOrange + 7,
        ROOT.kGreen + 2,
        ROOT.kMagenta + 1,
        ROOT.kRed + 1,
        ROOT.kBlue + 2,
        ROOT.kTeal + 3,
        ROOT.kPink + 7,
        ROOT.kSpring + 5,
        ROOT.kViolet + 7,
        ROOT.kCyan + 2,
        ROOT.kOrange - 3,
        ROOT.kGreen - 2,
        ROOT.kRed - 7,
        ROOT.kBlue - 7,
        ROOT.kMagenta - 7,
        ROOT.kGray + 2,
    ]


def assign_line_color(index):
    base_colors = color_sequence()
    return base_colors[index % len(base_colors)]


def main():
    args = parse_args()
    set_style()
    os.makedirs(args.outdir, exist_ok=True)

    root_file = ROOT.TFile.Open(args.input)
    if not root_file or root_file.IsZombie():
        raise RuntimeError(f"Could not open input file: {args.input}")

    ws = find_workspace(root_file, args.workspace)
    pdf = find_pdf(ws, args.pdf)
    xvar = ws.var(args.xvar)
    mhvar = ws.var(args.mhvar)

    if not xvar:
        raise RuntimeError(f"Observable '{args.xvar}' not found in workspace '{ws.GetName()}'")
    if not mhvar:
        raise RuntimeError(f"Mass variable '{args.mhvar}' not found in workspace '{ws.GetName()}'")

    xmin = args.xmin if args.xmin is not None else xvar.getMin()
    xmax = args.xmax if args.xmax is not None else xvar.getMax()
    mh_values = build_mh_values(args.mhMin, args.mhMax, args.mhStep)

    hists = []
    ymax = 0.0
    for idx, mh_value in enumerate(mh_values):
        hist = sample_pdf(
            pdf,
            xvar,
            mhvar,
            mh_value,
            args.nbins,
            xmin,
            xmax,
            unit_norm=(not args.noUnitNorm),
        )
        hist.SetDirectory(0)
        hist.SetLineColor(assign_line_color(idx))
        hist.SetLineWidth(3)
        hist.SetFillStyle(0)
        ymax = max(ymax, hist.GetMaximum())
        hists.append((mh_value, hist))

    canvas = ROOT.TCanvas("c_workspace_scan", "", 1200, 1000)
    canvas.SetLeftMargin(0.14)
    canvas.SetBottomMargin(0.13)

    frame_hist = ROOT.TH1D("frame_hist", "", 1, xmin, xmax)
    frame_hist.SetDirectory(0)
    frame_hist.SetTitle("")
    frame_hist.GetXaxis().SetTitle("m_{jj} [GeV]")
    frame_hist.GetYaxis().SetTitle("Arbitrary units" if args.noUnitNorm else "a.u. / normalized")
    frame_hist.SetMinimum(0.0)
    frame_hist.SetMaximum(1.25 * ymax if ymax > 0.0 else 1.0)
    frame_hist.Draw("AXIS")

    first = True
    for _, hist in hists:
        if first:
            hist.Draw("HIST SAME")
            first = False
        else:
            hist.Draw("HIST SAME")

    legend = ROOT.TLegend(0.62, 0.48, 0.90, 0.88)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.030)
    for mh_value, hist in hists:
        legend.AddEntry(hist, f"m_{{H}} = {int(round(mh_value))} GeV", "l")
    legend.Draw()

    text = ROOT.TLatex()
    text.SetNDC()
    text.SetTextFont(42)
    text.SetTextSize(0.032)
    text.DrawLatex(0.16, 0.84, args.label)

    draw_cms_label()

    outbase = os.path.join(args.outdir, args.outname)
    canvas.SaveAs(outbase + ".png")
    canvas.SaveAs(outbase + ".pdf")

    root_file.Close()


if __name__ == "__main__":
    main()
