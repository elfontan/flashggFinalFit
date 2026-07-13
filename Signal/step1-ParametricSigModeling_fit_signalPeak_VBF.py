#!/usr/bin/env python3

"""
Robust parametric fit for dijet VBF signal peak histograms.

Strategy:
1. Fit each mass point independently with a Double Crystal Ball.
2. Fit each DCB parameter as a polynomial function of the signal mass.
3. Build a parametric RooCrystalBall model from those polynomial trends.
"""

import argparse
import json
import math
import os
import re
from array import array
from collections import OrderedDict as od

import numpy as np
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

DEFAULT_PARAM_ORDERS = {}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Robust parametric fit of VBF dijet signal histograms with a Double Crystal Ball."
    )
    parser.add_argument(
        "--indir",
        type=str,
        default="/afs/cern.ch/work/e/elfontan/private/dijetAnalysis_ScoutingRun3/MC_Samples/SignalOptimisationStudies/VBF-Cat/VBF-dijetMass-Histos_ForFIT",
        help="Directory containing vbf-m<MASS>-<ALGO>.root files.",
    )
    parser.add_argument(
        "--outdir",
        type=str,
        default="plots_VBFSignalFits_parametric",
        help="Output directory.",
    )
    parser.add_argument(
        "--hist",
        type=str,
        default="h_widejet_peak_centralFirst_final_5GeV",
        help="Histogram name to fit.",
    )
    parser.add_argument(
        "--algo",
        type=str,
        default=None,
        help="Optional algorithm token filter for legacy files named vbf-m<MASS>-<ALGO>.root.",
    )
    parser.add_argument(
        "--recursive",
        action="store_true",
        help="Scan subdirectories of --indir as well.",
    )
    parser.add_argument(
        "--massPoints",
        type=str,
        default=None,
        help="Comma-separated subset of masses to use, e.g. 300,500,750,1000,2000,3000.",
    )
    parser.add_argument(
        "--ext",
        type=str,
        default="parametric",
        help="Extension tag used in output names.",
    )
    parser.add_argument(
        "--mhPolyOrder",
        type=int,
        default=2,
        help="Default polynomial order in mass for the DCB parameters.",
    )
    parser.add_argument(
        "--paramOrders",
        type=str,
        default=None,
        help="Optional comma-separated per-parameter overrides, e.g. nL:0,sigmaL:2,mean:2.",
    )
    parser.add_argument(
        "--fitRangeMin",
        type=float,
        default=0.100,
        help="Global fit-range minimum in m_jj.",
    )
    parser.add_argument(
        "--fitRangeMax",
        type=float,
        default=3300.0,
        help="Global fit-range maximum in m_jj.",
    )
    parser.add_argument(
        "--fitWindowLowScale",
        type=float,
        default=0.25,
        help="Per-mass fit window low edge as a fraction of the signal mass.",
    )
    parser.add_argument(
        "--fitWindowHighScale",
        type=float,
        default=1.75,
        help="Per-mass fit window high edge as a fraction of the signal mass.",
    )
    parser.add_argument(
        "--massWindowOverrides",
        type=str,
        default=None,
        help="Optional comma-separated per-mass fit-window overrides as mass:low:high, e.g. 300:0.40:1.45,500:0.45:1.55.",
    )
    parser.add_argument(
        "--drawWindowLowScale",
        type=float,
        default=0.20,
        help="Per-mass display window low edge as a fraction of the signal mass.",
    )
    parser.add_argument(
        "--drawWindowHighScale",
        type=float,
        default=2.10,
        help="Per-mass display window high edge as a fraction of the signal mass.",
    )
    parser.add_argument(
        "--modelMassMax",
        type=float,
        default=3200.0,
        help="Upper edge of the parametric MH range.",
    )
    parser.add_argument(
        "--referenceMass",
        type=float,
        default=750.0,
        help="Reference mass used to center the polynomial parameterization.",
    )
    parser.add_argument(
        "--plotAsHist",
        action="store_true",
        help="Sample the fitted pdf into a TH1 instead of drawing a smooth RooFit curve.",
    )
    parser.add_argument(
        "--plotMarginScale",
        type=float,
        default=7.5,
        help="Plot half-width in units of max(sigmaL,sigmaR).",
    )
    parser.add_argument(
        "--plotMinHalfWidth",
        type=float,
        default=300.0,
        help="Minimum half-width of the displayed mass window in GeV.",
    )
    parser.add_argument(
        "--useWidePlotWindow",
        action="store_true",
        help="Use the full configured draw window for the single-mass and overlay plots.",
    )
    parser.add_argument(
        "--saveWorkspace",
        action="store_true",
        help="Save a RooWorkspace containing the parametric DCB model.",
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


PARAM_TREND_CONFIG = {
    "mean": {"mode": "linear", "floor": None},
    "sigmaL": {"mode": "log_shifted", "floor": 1.0},
    "sigmaR": {"mode": "log_shifted", "floor": 1.0},
    "alphaL": {"mode": "log_shifted", "floor": 0.05},
    "nL": {"mode": "log_shifted", "floor": 0.5},
    "alphaR": {"mode": "log_shifted", "floor": 0.05},
    "nR": {"mode": "log_shifted", "floor": 0.5},
}


def parse_signal_file(path, relpath):
    match = re.search(r"vbf-m(?P<mass>\d+)-(?P<algo>[^/]+)\.root$", relpath)
    if match:
        parent_dir = os.path.dirname(relpath)
        return {
            "path": path,
            "mass": int(match.group("mass")),
            "algo": match.group("algo"),
            "tag": os.path.basename(parent_dir) if parent_dir not in ("", ".") else "base",
        }
    match = re.search(r"vbf-m(?P<mass>\d+)\.root$", relpath)
    if not match:
        return None
    parent_dir = os.path.dirname(relpath)
    return {
        "path": path,
        "mass": int(match.group("mass")),
        "algo": "default",
        "tag": os.path.basename(parent_dir) if parent_dir not in ("", ".") else "base",
    }


def discover_files(indir, recursive, algo_filter, selected_masses=None):
    file_infos = []
    if recursive:
        for root_dir, _, filenames in os.walk(indir):
            for filename in filenames:
                if not filename.endswith(".root"):
                    continue
                full_path = os.path.join(root_dir, filename)
                relpath = os.path.relpath(full_path, indir)
                info = parse_signal_file(full_path, relpath)
                if info is None:
                    continue
                if algo_filter and info["algo"] != algo_filter:
                    continue
                if selected_masses and info["mass"] not in selected_masses:
                    continue
                file_infos.append(info)
    else:
        for filename in os.listdir(indir):
            if not filename.endswith(".root"):
                continue
            full_path = os.path.join(indir, filename)
            info = parse_signal_file(full_path, filename)
            if info is None:
                continue
            if algo_filter and info["algo"] != algo_filter:
                continue
            if selected_masses and info["mass"] not in selected_masses:
                continue
            file_infos.append(info)
    return sorted(file_infos, key=lambda item: (item["tag"], item["algo"], item["mass"]))


def load_histogram(info, hist_name):
    root_file = ROOT.TFile.Open(info["path"])
    if not root_file or root_file.IsZombie():
        raise RuntimeError(f"Could not open file: {info['path']}")
    hist = root_file.Get(hist_name)
    if not hist:
        root_file.Close()
        raise RuntimeError(f"Missing histogram {hist_name} in {info['path']}")
    hist = hist.Clone(f"{hist.GetName()}_{info['mass']}_{info['algo']}_{info['tag']}")
    hist.SetDirectory(0)
    root_file.Close()
    if hist.Integral() <= 0:
        raise RuntimeError(f"Empty histogram in {info['path']}")
    return hist


def sanitized_name(token):
    return re.sub(r"[^A-Za-z0-9_]+", "_", token)


def parse_mass_window_overrides(spec):
    overrides = {}
    if not spec:
        return overrides
    for token in spec.split(","):
        token = token.strip()
        if not token:
            continue
        parts = token.split(":")
        if len(parts) != 3:
            raise ValueError(
                f"Invalid mass-window override '{token}'. Expected mass:low:high."
            )
        mass = int(parts[0])
        low = float(parts[1])
        high = float(parts[2])
        if high <= low:
            raise ValueError(
                f"Invalid mass-window override '{token}': high must be greater than low."
            )
        overrides[mass] = (low, high)
    return overrides


def build_datahist(hist, xvar, name):
    return ROOT.RooDataHist(name, name, ROOT.RooArgList(xvar), ROOT.RooFit.Import(hist))


def get_single_fit_window(mass, hist_xmin, hist_xmax, global_xmin, global_xmax, low_scale=0.30, high_scale=2.10):
    return (
        max(hist_xmin, global_xmin, low_scale * mass),
        min(hist_xmax, global_xmax, high_scale * mass),
    )


def get_single_draw_range(mass, hist_xmin, hist_xmax, global_xmin, global_xmax, low_scale=0.20, high_scale=2.10):
    return (
        max(hist_xmin, global_xmin, low_scale * mass),
        min(hist_xmax, global_xmax, high_scale * mass),
    )


def determine_plot_window(mean, sigma_left, sigma_right, xmin, xmax, margin_scale=6.0, min_half_width=220.0):
    width = max(sigma_left, sigma_right)
    margin = max(float(min_half_width), float(margin_scale) * width)
    low = max(xmin, mean - margin)
    high = min(xmax, mean + margin)
    if high <= low:
        low, high = xmin, xmax
    return low, high


def dcb_shape_value(xval, mean, sigma_left, sigma_right, alpha_left, n_left, alpha_right, n_right):
    alpha_left = max(0.05, float(alpha_left))
    alpha_right = max(0.05, float(alpha_right))
    n_left = max(0.5, float(n_left))
    n_right = max(0.5, float(n_right))
    sigma_left = max(1.0, float(sigma_left))
    sigma_right = max(1.0, float(sigma_right))

    if xval < mean:
        tval = (xval - mean) / sigma_left
        if tval < -alpha_left:
            aval = (n_left / alpha_left) ** n_left * math.exp(-0.5 * alpha_left * alpha_left)
            bval = n_left / alpha_left - alpha_left
            return aval * (bval - tval) ** (-n_left)
        return math.exp(-0.5 * tval * tval)

    tval = (xval - mean) / sigma_right
    if tval > alpha_right:
        aval = (n_right / alpha_right) ** n_right * math.exp(-0.5 * alpha_right * alpha_right)
        bval = n_right / alpha_right - alpha_right
        return aval * (bval + tval) ** (-n_right)
    return math.exp(-0.5 * tval * tval)


def make_dcb_histogram(params, template_hist, name, xmin, xmax):
    bin_width = template_hist.GetXaxis().GetBinWidth(1)
    nbins = max(1, int(round((xmax - xmin) / bin_width)))
    hist = ROOT.TH1D(name, "", nbins, xmin, xmax)
    hist.Sumw2()
    for ibin in range(1, nbins + 1):
        xval = hist.GetBinCenter(ibin)
        val = dcb_shape_value(
            xval,
            params["mean"],
            params["sigmaL"],
            params["sigmaR"],
            params["alphaL"],
            params["nL"],
            params["alphaR"],
            params["nR"],
        )
        hist.SetBinContent(ibin, val * hist.GetBinWidth(ibin))
    if hist.Integral() > 0.0 and template_hist.Integral() > 0.0:
        hist.Scale(template_hist.Integral() / hist.Integral())
    return hist


def compute_single_fit_gof(hist, params, fit_low, fit_high, n_float_params=7):
    model_hist = make_dcb_histogram(
        params,
        hist,
        f"{hist.GetName()}_gof_model",
        hist.GetXaxis().GetXmin(),
        hist.GetXaxis().GetXmax(),
    )

    first_bin = hist.FindBin(fit_low + 1e-6)
    last_bin = hist.FindBin(fit_high - 1e-6)
    n_bins = max(0, last_bin - first_bin + 1)
    ndof = max(0, n_bins - n_float_params)

    pearson = 0.0
    baker_cousins = 0.0

    for ibin in range(first_bin, last_bin + 1):
        obs = float(hist.GetBinContent(ibin))
        exp = float(model_hist.GetBinContent(ibin))
        if exp > 0.0:
            pearson += (obs - exp) ** 2 / exp
            if obs > 0.0:
                baker_cousins += 2.0 * (exp - obs + obs * math.log(obs / exp))
            else:
                baker_cousins += 2.0 * exp

    return {
        "pearson_chi2": float(pearson),
        "baker_cousins_chi2": float(baker_cousins),
        "ndof": int(ndof),
        "n_bins_fit": int(n_bins),
    }


def fit_single_mass(
    hist,
    mass,
    fit_xmin,
    fit_xmax,
    fit_scale_low,
    fit_scale_high,
    draw_scale_low,
    draw_scale_high,
):
    hist_xmin = hist.GetXaxis().GetXmin()
    hist_xmax = hist.GetXaxis().GetXmax()
    fit_low, fit_high = get_single_fit_window(
        mass,
        hist_xmin,
        hist_xmax,
        fit_xmin,
        fit_xmax,
        fit_scale_low,
        fit_scale_high,
    )
    draw_low, draw_high = get_single_draw_range(
        mass,
        hist_xmin,
        hist_xmax,
        fit_xmin,
        fit_xmax,
        draw_scale_low,
        draw_scale_high,
    )

    if fit_high <= fit_low:
        raise RuntimeError(
            f"Invalid fit window for mass {mass}: [{fit_low}, {fit_high}] from global range "
            f"[{fit_xmin}, {fit_xmax}] and scales [{fit_scale_low}, {fit_scale_high}]"
        )
    if draw_high <= draw_low:
        draw_low, draw_high = fit_low, fit_high

    xvar = ROOT.RooRealVar(f"mjj_{mass}", "m_{jj}", fit_low, fit_high)
    datahist = build_datahist(hist, xvar, f"datahist_m{mass}")

    mean = ROOT.RooRealVar(f"mean_{mass}", "mean", float(mass), 0.75 * mass, 1.25 * mass)
    sigma_left = ROOT.RooRealVar(
        f"sigmaL_{mass}", "sigmaL", max(10.0, 0.08 * mass), 10.0, 600.0
    )
    sigma_right = ROOT.RooRealVar(
        f"sigmaR_{mass}", "sigmaR", max(12.0, 0.12 * mass), 1.0, 500.0
    )
    alpha_left = ROOT.RooRealVar(f"alphaL_{mass}", "alphaL", 5.0, 0.05, 30.0)
    n_left = ROOT.RooRealVar(f"nL_{mass}", "nL", 1.0, 0.5, 50.0)
    alpha_right = ROOT.RooRealVar(f"alphaR_{mass}", "alphaR", 1.2, 0.03, 7.0)
    n_right = ROOT.RooRealVar(f"nR_{mass}", "nR", 2.0, 0.1, 50.0)

    model = ROOT.RooCrystalBall(
        f"dcb_{mass}",
        f"dcb_{mass}",
        xvar,
        mean,
        sigma_left,
        sigma_right,
        alpha_left,
        n_left,
        alpha_right,
        n_right,
    )

    fitres = model.fitTo(
        datahist,
        ROOT.RooFit.Save(True),
        ROOT.RooFit.PrintLevel(-1),
        ROOT.RooFit.Strategy(1),
    )

    temp_result = {
        "mean": float(mean.getVal()),
        "sigmaL": float(sigma_left.getVal()),
        "sigmaR": float(sigma_right.getVal()),
        "alphaL": float(alpha_left.getVal()),
        "nL": float(n_left.getVal()),
        "alphaR": float(alpha_right.getVal()),
        "nR": float(n_right.getVal()),
    }
    gof = compute_single_fit_gof(hist, temp_result, fit_low, fit_high)

    result = {
        "mass": int(mass),
        "fitRange": [float(fit_low), float(fit_high)],
        "drawRange": [float(draw_low), float(draw_high)],
        "mean": float(mean.getVal()),
        "mean_err": float(mean.getError()),
        "sigmaL": float(sigma_left.getVal()),
        "sigmaL_err": float(sigma_left.getError()),
        "sigmaR": float(sigma_right.getVal()),
        "sigmaR_err": float(sigma_right.getError()),
        "alphaL": float(alpha_left.getVal()),
        "alphaL_err": float(alpha_left.getError()),
        "nL": float(n_left.getVal()),
        "nL_err": float(n_left.getError()),
        "alphaR": float(alpha_right.getVal()),
        "alphaR_err": float(alpha_right.getError()),
        "nR": float(n_right.getVal()),
        "nR_err": float(n_right.getError()),
        "status": int(fitres.status()),
        "covQual": int(fitres.covQual()),
        "pearson_chi2": gof["pearson_chi2"],
        "baker_cousins_chi2": gof["baker_cousins_chi2"],
        "chi2_ndof": gof["ndof"],
        "n_bins_fit": gof["n_bins_fit"],
    }
    return {
        "hist": hist,
        "xvar": xvar,
        "datahist": datahist,
        "model": model,
        "fitres": fitres,
        "result": result,
    }


def draw_single_mass_fit(fit_info, outdir, tag, plot_range=None, plot_as_hist=True, plot_margin_scale=6.0, plot_min_half_width=220.0):
    mass = fit_info["result"]["mass"]
    hist = fit_info["hist"]
    xvar = fit_info["xvar"]
    datahist = fit_info["datahist"]
    model = fit_info["model"]
    result = fit_info["result"]

    if plot_range is not None:
        plot_low, plot_high = plot_range
    elif "drawRange" in result:
        plot_low, plot_high = result["drawRange"]
    else:
        plot_low, plot_high = determine_plot_window(
            result["mean"],
            result["sigmaL"],
            result["sigmaR"],
            hist.GetXaxis().GetXmin(),
            hist.GetXaxis().GetXmax(),
            margin_scale=plot_margin_scale,
            min_half_width=plot_min_half_width,
        )

    canvas = ROOT.TCanvas(f"c_single_m{mass}_{tag}", "", 1200, 1000)
    canvas.SetLeftMargin(0.14)
    canvas.SetBottomMargin(0.13)

    if plot_as_hist:
        data_plot = hist.Clone(f"{hist.GetName()}_singlemass_plot")
        data_plot.SetDirectory(0)
        data_plot.SetMarkerStyle(20)
        data_plot.SetMarkerSize(0.9)
        data_plot.SetMarkerColor(ROOT.kBlack)
        data_plot.SetLineColor(ROOT.kBlack)
        data_plot.GetXaxis().SetRangeUser(plot_low, plot_high)
        data_plot.GetXaxis().SetTitle("m_{jj} [GeV]")
        data_plot.GetYaxis().SetTitle("Events / 5 GeV")

        pdf_plot = make_dcb_histogram(
            result,
            hist,
            f"h_pdf_single_m{mass}_{tag}",
            plot_low,
            plot_high,
        )
        pdf_plot.SetLineColor(ROOT.kAzure + 2)
        pdf_plot.SetLineWidth(3)
        pdf_plot.SetFillStyle(0)

        data_plot.SetTitle("")
        ymax = 1.35 * max(data_plot.GetMaximum(), pdf_plot.GetMaximum())
        data_plot.SetMaximum(ymax if ymax > 0.0 else 1.0)
        data_plot.Draw("E1")
        pdf_plot.Draw("HIST SAME")
        data_obj = data_plot
        pdf_obj = pdf_plot
    else:
        model.removeStringAttribute("fitrange")
        frame = xvar.frame(ROOT.RooFit.Range(plot_low, plot_high))
        datahist.plotOn(
            frame,
            ROOT.RooFit.Name("data"),
            ROOT.RooFit.MarkerStyle(20),
            ROOT.RooFit.MarkerSize(0.9),
            ROOT.RooFit.LineColor(ROOT.kBlack),
            ROOT.RooFit.CutRange("fitRange"),
        )
        model.plotOn(
            frame,
            ROOT.RooFit.Name("model"),
            ROOT.RooFit.LineColor(ROOT.kAzure + 2),
            ROOT.RooFit.LineWidth(3),
            ROOT.RooFit.Range(plot_low, plot_high),
            ROOT.RooFit.NormRange("fitRange"),
        )
        frame.SetTitle("")
        frame.GetXaxis().SetTitle("m_{jj} [GeV]")
        frame.GetYaxis().SetTitle("Events / 5 GeV")
        frame.Draw()
        data_obj = frame.findObject("data")
        pdf_obj = frame.findObject("model")

    legend = ROOT.TLegend(0.58, 0.72, 0.90, 0.88)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.030)
    legend.SetHeader(f"Single-mass DCB fit, M={mass} GeV", "C")
    legend.AddEntry(data_obj, "Input histogram", "lep")
    legend.AddEntry(pdf_obj, "Double Crystal Ball", "l")
    legend.Draw()

    text = ROOT.TLatex()
    text.SetNDC()
    text.SetTextFont(42)
    text.SetTextSize(0.028)
    text.DrawLatex(0.66, 0.66, f"#mu = {result['mean']:.2f} #pm {result['mean_err']:.2f} GeV")
    text.DrawLatex(0.66, 0.62, f"#sigma_{{L}} = {result['sigmaL']:.2f} GeV")
    text.DrawLatex(0.66, 0.58, f"#sigma_{{R}} = {result['sigmaR']:.2f} GeV")
    if result["chi2_ndof"] > 0:
        chi2_over_ndof = result["baker_cousins_chi2"] / float(result["chi2_ndof"])
        text.DrawLatex(
            0.66,
            0.54,
            f"#chi^{{2}}/ndof = {chi2_over_ndof:.2f}",
        )
        text.DrawLatex(0.66, 0.50, f"status = {result['status']}, covQual = {result['covQual']}")
    else:
        text.DrawLatex(0.66, 0.54, f"status = {result['status']}, covQual = {result['covQual']}")
    draw_cms_label()

    outbase = os.path.join(outdir, f"{tag}_m{mass}_singleMassFit")
    canvas.SaveAs(outbase + ".png")
    canvas.SaveAs(outbase + ".pdf")


def safe_err(value):
    if value is None or not math.isfinite(value) or value <= 0.0:
        return 1.0
    return float(value)


def parse_param_order_overrides(raw_text):
    overrides = dict(DEFAULT_PARAM_ORDERS)
    if not raw_text:
        return overrides
    for token in raw_text.split(","):
        token = token.strip()
        if not token:
            continue
        if ":" not in token:
            raise RuntimeError(f"Invalid --paramOrders entry '{token}'. Expected name:order.")
        name, order_text = token.split(":", 1)
        name = name.strip()
        if name not in PARAM_TREND_CONFIG:
            raise RuntimeError(f"Unknown parameter '{name}' in --paramOrders.")
        overrides[name] = int(order_text.strip())
    return overrides


def get_parameter_order(name, default_order, overrides):
    return int(overrides.get(name, default_order))


def transform_parameter_value(name, value):
    config = PARAM_TREND_CONFIG[name]
    if config["mode"] == "linear":
        return float(value)
    floor = float(config["floor"])
    return math.log(max(float(value) - floor, 1e-6))


def inverse_transform_parameter_value(name, value):
    config = PARAM_TREND_CONFIG[name]
    if config["mode"] == "linear":
        return float(value)
    floor = float(config["floor"])
    return floor + math.exp(float(value))


def evaluate_trend_value(trend_fit, x_value):
    raw_value = evaluate_polynomial(trend_fit["coeffs"], x_value)
    return inverse_transform_parameter_value(trend_fit["name"], raw_value)


def fit_parameter_trend(name, masses, values, errors, order, reference_mass):
    max_order = min(order, len(masses) - 1)
    x = np.asarray([mass - reference_mass for mass in masses], dtype=float)
    y_physical = np.asarray(values, dtype=float)
    y = np.asarray([transform_parameter_value(name, value) for value in values], dtype=float)
    w = np.asarray([1.0 / safe_err(err) for err in errors], dtype=float)
    coeffs_desc = np.polyfit(x, y, deg=max_order, w=w)
    coeffs = [float(value) for value in coeffs_desc[::-1]]
    fitted = [evaluate_trend_value({"name": name, "coeffs": coeffs}, xv) for xv in x]
    chi2 = 0.0
    for yval, yfit, err in zip(y_physical, fitted, errors):
        err_use = safe_err(err)
        chi2 += ((yval - yfit) / err_use) ** 2
    ndof = max(0, len(masses) - len(coeffs))
    return coeffs, fitted, chi2, ndof


def evaluate_polynomial(coeffs, x):
    total = 0.0
    for power, coeff in enumerate(coeffs):
        total += coeff * (x ** power)
    return total


def make_graph_with_errors(masses, values, errors, color):
    graph = ROOT.TGraphErrors(
        len(masses),
        array("d", [float(mass) for mass in masses]),
        array("d", [float(value) for value in values]),
        array("d", [0.0] * len(masses)),
        array("d", [float(err) for err in errors]),
    )
    graph.SetLineColor(color)
    graph.SetMarkerColor(color)
    graph.SetMarkerStyle(20)
    graph.SetMarkerSize(1.1)
    graph.SetLineWidth(3)
    return graph


def make_trend_curve(trend_fit, reference_mass, xmin, xmax, color):
    npoints = 400
    xs = []
    ys = []
    for idx in range(npoints):
        xval = xmin + (xmax - xmin) * idx / float(npoints - 1)
        xs.append(xval)
        ys.append(evaluate_trend_value(trend_fit, xval - reference_mass))
    graph = ROOT.TGraph(npoints, array("d", xs), array("d", ys))
    graph.SetLineColor(color)
    graph.SetLineWidth(3)
    return graph


def save_trend_plots(fit_results, trend_fits, masses, outdir, tag, reference_mass, model_mass_max):
    param_map = od(
        [
            ("mean", "#mu [GeV]"),
            ("sigmaL", "#sigma_{L} [GeV]"),
            ("sigmaR", "#sigma_{R} [GeV]"),
            ("alphaL", "#alpha_{L}"),
            ("nL", "n_{L}"),
            ("alphaR", "#alpha_{R}"),
            ("nR", "n_{R}"),
        ]
    )

    colors = {
        "mean": ROOT.kAzure + 1,
        "sigmaL": ROOT.kOrange + 7,
        "sigmaR": ROOT.kOrange + 2,
        "alphaL": ROOT.kGreen + 2,
        "nL": ROOT.kRed + 1,
        "alphaR": ROOT.kMagenta + 1,
        "nR": ROOT.kBlue + 2,
    }

    parametric_values = od()
    for mass in masses:
        parametric_values[str(mass)] = {}
        for name in param_map:
            parametric_values[str(mass)][name] = evaluate_trend_value(trend_fits[name], mass - reference_mass)

    with open(os.path.join(outdir, f"{tag}_parametric_parameters.json"), "w") as handle:
        json.dump(parametric_values, handle, indent=2, sort_keys=True)

    for name, ytitle in param_map.items():
        values = [fit_results[mass][name] for mass in masses]
        errors = [fit_results[mass][f"{name}_err"] for mass in masses]
        color = colors[name]
        graph = make_graph_with_errors(masses, values, errors, color)
        curve = make_trend_curve(trend_fit=trend_fits[name], reference_mass=reference_mass, xmin=min(masses), xmax=model_mass_max, color=color)

        canvas = ROOT.TCanvas(f"c_{tag}_{name}", "", 1000, 800)
        canvas.SetLeftMargin(0.14)
        canvas.SetBottomMargin(0.13)
        graph.SetTitle(f";Signal mass [GeV];{ytitle}")
        curve_min = min(curve.GetY()[i] for i in range(curve.GetN()))
        curve_max = max(curve.GetY()[i] for i in range(curve.GetN()))
        point_vals = []
        for val, err in zip(values, errors):
            err_use = safe_err(err)
            point_vals.extend([val - err_use, val + err_use])
        y_min = min(point_vals + [curve_min])
        y_max = max(point_vals + [curve_max])
        span = max(1e-6, y_max - y_min)
        lower_pad = 0.18
        upper_pad = 0.22
        if name == "alphaR":
            lower_pad = 0.28
            upper_pad = 0.35
        elif name == "sigmaL":
            lower_pad = 0.24
            upper_pad = 0.32
        graph.SetMinimum(y_min - lower_pad * span)
        graph.SetMaximum(y_max + upper_pad * span)
        graph.Draw("AP")
        curve.Draw("L SAME")

        legend = ROOT.TLegend(0.58, 0.74, 0.90, 0.88)
        legend.SetBorderSize(0)
        legend.SetFillStyle(0)
        legend.SetTextSize(0.030)
        legend.AddEntry(graph, "Single-mass fit", "lep")
        trend_label = f"pol{trend_fits[name]['order']}"
        if trend_fits[name]["mode"] != "linear":
            trend_label += " in log-space"
        legend.AddEntry(curve, trend_label, "l")
        legend.Draw()

        draw_cms_label()

        outbase = os.path.join(outdir, f"{tag}_{name}_vs_mass")
        canvas.SaveAs(outbase + ".png")
        canvas.SaveAs(outbase + ".pdf")


def build_parametric_model(trend_fits, fit_xmin, fit_xmax, mass_min, model_mass_max, reference_mass, tag):
    xvar = ROOT.RooRealVar("mjj", "m_{jj}", fit_xmin, fit_xmax)
    MH = ROOT.RooRealVar("MH", "MH", float(mass_min), float(model_mass_max))
    dMH = ROOT.RooFormulaVar("dMH", f"@0-{reference_mass:.6f}", ROOT.RooArgList(MH))

    keep = {"coeffs": [], "poly": {}, "safe": {}}

    def make_poly(name):
        coeff_list = ROOT.RooArgList(f"{name}_coeffs")
        for idx, coeff in enumerate(trend_fits[name]["coeffs"]):
            var = ROOT.RooRealVar(f"{name}_p{idx}", f"{name}_p{idx}", float(coeff))
            var.setConstant(True)
            keep["coeffs"].append(var)
            coeff_list.add(var)
        poly = ROOT.RooPolyVar(name, name, dMH, coeff_list)
        keep["poly"][name] = poly
        return poly

    def make_physical(name):
        raw = make_poly(name)
        config = PARAM_TREND_CONFIG[name]
        if config["mode"] == "linear":
            return raw
        floor = float(config["floor"])
        transformed = ROOT.RooFormulaVar(
            f"{name}_phys",
            f"exp(@0)+{floor:.12g}",
            ROOT.RooArgList(raw),
        )
        keep["safe"][name] = transformed
        return transformed

    mean = make_physical("mean")
    sigma_left = make_physical("sigmaL")
    sigma_right = make_physical("sigmaR")
    alpha_left = make_physical("alphaL")
    n_left = make_physical("nL")
    alpha_right = make_physical("alphaR")
    n_right = make_physical("nR")

    pdf = ROOT.RooCrystalBall(
        f"{tag}_pdf",
        f"{tag}_pdf",
        xvar,
        mean,
        sigma_left,
        sigma_right,
        alpha_left,
        n_left,
        alpha_right,
        n_right,
    )
    return {
        "xvar": xvar,
        "MH": MH,
        "dMH": dMH,
        "pdf": pdf,
        "mean": mean,
        "keep": keep,
    }


def evaluate_parametric_point(trend_fits, reference_mass, mass):
    dm = float(mass) - float(reference_mass)
    values = {}
    for name in ["mean", "sigmaL", "sigmaR", "alphaL", "nL", "alphaR", "nR"]:
        values[name] = evaluate_trend_value(trend_fits[name], dm)
    return values


def draw_parametric_overlay(model_info, trend_fits, reference_mass, hist, mass, outdir, tag, plot_range=None, plot_as_hist=True, plot_margin_scale=6.0, plot_min_half_width=220.0):
    params = evaluate_parametric_point(trend_fits, reference_mass, mass)
    mean_val = params["mean"]
    sigma_left_val = params["sigmaL"]
    sigma_right_val = params["sigmaR"]
    if plot_range is None:
        plot_low, plot_high = determine_plot_window(
            mean_val,
            sigma_left_val,
            sigma_right_val,
            hist.GetXaxis().GetXmin(),
            hist.GetXaxis().GetXmax(),
            margin_scale=plot_margin_scale,
            min_half_width=plot_min_half_width,
        )
    else:
        plot_low, plot_high = plot_range
    gof = compute_single_fit_gof(hist, params, plot_low, plot_high)

    canvas = ROOT.TCanvas(f"c_param_{mass}_{tag}", "", 1200, 1000)
    canvas.SetLeftMargin(0.14)
    canvas.SetBottomMargin(0.13)

    if plot_as_hist:
        data_plot = hist.Clone(f"{hist.GetName()}_param_plot_{mass}")
        data_plot.SetDirectory(0)
        data_plot.SetMarkerStyle(20)
        data_plot.SetMarkerSize(0.9)
        data_plot.SetMarkerColor(ROOT.kBlack)
        data_plot.SetLineColor(ROOT.kBlack)
        data_plot.GetXaxis().SetRangeUser(plot_low, plot_high)
        data_plot.GetXaxis().SetTitle("m_{jj} [GeV]")
        data_plot.GetYaxis().SetTitle("Events / 5 GeV")

        pdf_plot = make_dcb_histogram(
            params,
            hist,
            f"h_pdf_param_m{mass}_{tag}",
            plot_low,
            plot_high,
        )
        pdf_plot.SetLineColor(ROOT.kOrange + 7)
        pdf_plot.SetLineWidth(3)
        pdf_plot.SetFillStyle(0)

        data_plot.SetTitle("")
        ymax = 1.35 * max(data_plot.GetMaximum(), pdf_plot.GetMaximum())
        data_plot.SetMaximum(ymax if ymax > 0.0 else 1.0)
        data_plot.Draw("E1")
        pdf_plot.Draw("HIST SAME")
        data_obj = data_plot
        pdf_obj = pdf_plot
    else:
        xvar = model_info["xvar"]
        MH = model_info["MH"]
        pdf = model_info["pdf"]
        MH.setVal(float(mass))
        datahist = build_datahist(hist, xvar, f"param_datahist_{tag}_{mass}")
        frame = xvar.frame(ROOT.RooFit.Range(plot_low, plot_high))
        datahist.plotOn(
            frame,
            ROOT.RooFit.Name("data"),
            ROOT.RooFit.MarkerStyle(20),
            ROOT.RooFit.MarkerSize(0.9),
            ROOT.RooFit.LineColor(ROOT.kBlack),
        )
        pdf.plotOn(
            frame,
            ROOT.RooFit.Name("model"),
            ROOT.RooFit.LineColor(ROOT.kOrange + 7),
            ROOT.RooFit.LineWidth(3),
            ROOT.RooFit.Range(plot_low, plot_high),
        )
        frame.SetTitle("")
        frame.GetXaxis().SetTitle("m_{jj} [GeV]")
        frame.GetYaxis().SetTitle("Events / 5 GeV")
        frame.Draw()
        data_obj = frame.findObject("data")
        pdf_obj = frame.findObject("model")

    legend = ROOT.TLegend(0.58, 0.72, 0.90, 0.88)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.030)
    legend.SetHeader(f"Parametric DCB, M={mass} GeV", "C")
    legend.AddEntry(data_obj, "Input histogram", "lep")
    legend.AddEntry(pdf_obj, "Parametric DCB", "l")
    legend.Draw()

    text = ROOT.TLatex()
    text.SetNDC()
    text.SetTextFont(42)
    text.SetTextSize(0.028)
    text.DrawLatex(0.66, 0.66, f"#mu = {mean_val:.2f} GeV")
    text.DrawLatex(0.66, 0.62, f"#sigma_{{L}} = {sigma_left_val:.2f} GeV")
    text.DrawLatex(0.66, 0.58, f"#sigma_{{R}} = {sigma_right_val:.2f} GeV")
    if gof["ndof"] > 0:
        chi2_over_ndof = gof["baker_cousins_chi2"] / float(gof["ndof"])
        text.DrawLatex(0.66, 0.54, f"#chi^{{2}}/ndof = {chi2_over_ndof:.2f}")
    draw_cms_label()

    outbase = os.path.join(outdir, f"{tag}_m{mass}_parametricFit")
    canvas.SaveAs(outbase + ".png")
    canvas.SaveAs(outbase + ".pdf")


def save_workspace(model_info, outdir, tag):
    ws = ROOT.RooWorkspace(f"{tag}_ws", f"{tag}_ws")
    getattr(ws, "import")(model_info["xvar"])
    getattr(ws, "import")(model_info["MH"])
    getattr(ws, "import")(model_info["pdf"])
    for poly in model_info["keep"]["poly"].values():
        getattr(ws, "import")(poly)
    for func in model_info["keep"]["safe"].values():
        getattr(ws, "import")(func)
    fout = ROOT.TFile(os.path.join(outdir, f"{tag}_parametric_workspace.root"), "RECREATE")
    ws.Write()
    fout.Close()


def main():
    args = parse_args()
    set_style()
    os.makedirs(args.outdir, exist_ok=True)
    param_order_overrides = parse_param_order_overrides(args.paramOrders)

    selected_masses = None
    if args.massPoints:
        selected_masses = {int(token) for token in args.massPoints.split(",") if token}
    mass_window_overrides = parse_mass_window_overrides(args.massWindowOverrides)

    file_infos = discover_files(args.indir, args.recursive, args.algo, selected_masses)
    if not file_infos:
        raise RuntimeError(f"No matching ROOT files found in {args.indir}")

    fit_xmin = float(args.fitRangeMin)
    fit_xmax = float(args.fitRangeMax)
    hists = od()
    masses = []

    for info in file_infos:
        hist = load_histogram(info, args.hist)
        hists[info["mass"]] = hist
        masses.append(info["mass"])

    masses = sorted(set(masses))
    mass_min = min(masses)
    model_mass_max = max(float(args.modelMassMax), float(max(masses)))
    reference_mass = float(args.referenceMass)
    if reference_mass < mass_min or reference_mass > max(masses):
        print(
            f"[WARNING] reference mass {reference_mass:.1f} GeV is outside the fitted mass range "
            f"[{mass_min}, {max(masses)}]."
        )

    fit_results = od()
    plot_ranges = od()
    per_mass_summary = od()
    for mass in masses:
        fit_scale_low = args.fitWindowLowScale
        fit_scale_high = args.fitWindowHighScale
        if mass in mass_window_overrides:
            fit_scale_low, fit_scale_high = mass_window_overrides[mass]
        fit_info = fit_single_mass(
            hist=hists[mass],
            mass=mass,
            fit_xmin=fit_xmin,
            fit_xmax=fit_xmax,
            fit_scale_low=fit_scale_low,
            fit_scale_high=fit_scale_high,
            draw_scale_low=args.drawWindowLowScale,
            draw_scale_high=args.drawWindowHighScale,
        )
        fit_results[mass] = fit_info["result"]
        per_mass_summary[str(mass)] = fit_info["result"]
        if args.useWidePlotWindow:
            plot_ranges[mass] = tuple(fit_results[mass]["drawRange"])
        else:
            plot_ranges[mass] = determine_plot_window(
                fit_results[mass]["mean"],
                fit_results[mass]["sigmaL"],
                fit_results[mass]["sigmaR"],
                hists[mass].GetXaxis().GetXmin(),
                hists[mass].GetXaxis().GetXmax(),
                margin_scale=args.plotMarginScale,
                min_half_width=args.plotMinHalfWidth,
            )
        draw_single_mass_fit(
            fit_info,
            args.outdir,
            args.ext,
            plot_range=plot_ranges[mass],
            plot_as_hist=True,
            plot_margin_scale=args.plotMarginScale,
            plot_min_half_width=args.plotMinHalfWidth,
        )

    with open(os.path.join(args.outdir, f"{args.ext}_single_mass_fit_results.json"), "w") as handle:
        json.dump(per_mass_summary, handle, indent=2, sort_keys=True)

    param_names = ["mean", "sigmaL", "sigmaR", "alphaL", "nL", "alphaR", "nR"]
    trend_fits = od()
    for name in param_names:
        values = [fit_results[mass][name] for mass in masses]
        errors = [fit_results[mass][f"{name}_err"] for mass in masses]
        order_to_use = get_parameter_order(name, args.mhPolyOrder, param_order_overrides)
        coeffs, fitted, chi2, ndof = fit_parameter_trend(
            name=name,
            masses=masses,
            values=values,
            errors=errors,
            order=order_to_use,
            reference_mass=reference_mass,
        )
        trend_fits[name] = {
            "name": name,
            "coeffs": coeffs,
            "fitted": fitted,
            "chi2": float(chi2),
            "ndof": int(ndof),
            "order": int(order_to_use),
            "mode": PARAM_TREND_CONFIG[name]["mode"],
            "floor": PARAM_TREND_CONFIG[name]["floor"],
        }

    with open(os.path.join(args.outdir, f"{args.ext}_parametric_coefficients.json"), "w") as handle:
        json.dump(trend_fits, handle, indent=2, sort_keys=True)

    save_trend_plots(
        fit_results=fit_results,
        trend_fits=trend_fits,
        masses=masses,
        outdir=args.outdir,
        tag=args.ext,
        reference_mass=reference_mass,
        model_mass_max=model_mass_max,
    )

    model_info = build_parametric_model(
        trend_fits=trend_fits,
        fit_xmin=fit_xmin,
        fit_xmax=fit_xmax,
        mass_min=mass_min,
        model_mass_max=model_mass_max,
        reference_mass=reference_mass,
        tag=sanitized_name(args.ext),
    )

    for mass in masses:
        draw_parametric_overlay(
            model_info,
            trend_fits,
            reference_mass,
            hists[mass],
            mass,
            args.outdir,
            args.ext,
            plot_range=plot_ranges[mass],
            plot_as_hist=True,
            plot_margin_scale=args.plotMarginScale,
            plot_min_half_width=args.plotMinHalfWidth,
        )

    if args.saveWorkspace:
        save_workspace(model_info, args.outdir, args.ext)

    print("\n[INFO] Robust parametric fit complete")
    print(f"       masses            = {','.join(str(mass) for mass in masses)}")
    print(f"       fit range         = [{fit_xmin:.1f}, {fit_xmax:.1f}] GeV")
    print(
        "       fit window scales = "
        f"[{args.fitWindowLowScale:.2f}, {args.fitWindowHighScale:.2f}] x mass"
    )
    if mass_window_overrides:
        override_text = ",".join(
            f"{mass}:{low:.2f}:{high:.2f}"
            for mass, (low, high) in sorted(mass_window_overrides.items())
        )
        print(f"       mass overrides    = {override_text}")
    print(
        "       draw window scales = "
        f"[{args.drawWindowLowScale:.2f}, {args.drawWindowHighScale:.2f}] x mass"
    )
    print(f"       wide plot window  = {args.useWidePlotWindow}")
    print(f"       parametric MH max = {model_mass_max:.1f} GeV")
    print(f"       reference mass    = {reference_mass:.1f} GeV")
    print(
        "       order overrides   = "
        + ",".join(f"{name}:{order}" for name, order in sorted(param_order_overrides.items()))
    )
    print(f"       output dir        = {args.outdir}")


if __name__ == "__main__":
    main()
