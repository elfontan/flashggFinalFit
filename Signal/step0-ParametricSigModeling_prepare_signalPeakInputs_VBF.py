#!/usr/bin/env python3

"""
Prepare fit-ready VBF signal peak ROOT files from central-widejet study trees.

This script reads current reco-study outputs such as
`CentralVBFHTo2B_M500_centralWideJetVBF_scoutNano_<SELECTION>.root`, builds
the wide-jet full-selection mass histogram from the `Events` tree, and writes
one output file per mass with a selection-transparent naming:

  vbf-m<MASS>.root

The produced ROOT files contain the histogram needed by the parametric fit,
with contents derived directly from the current tree-based outputs.
"""

import argparse
import glob
import os
import re

import ROOT

ROOT.gROOT.SetBatch(True)

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DEFAULT_INPUT_DIR = os.path.join(SCRIPT_DIR, "Output-VBFSelection")
DEFAULT_OUTPUT_DIR = os.path.join(DEFAULT_INPUT_DIR, "PreparedFitInputs")
DEFAULT_SELECTION = "centPt30-Deta1p3_fwdPt24_VBF-Deta4-Mjj500"


def parse_args():
    parser = argparse.ArgumentParser(
        description="Prepare fit-ready VBF signal mass histograms from central-widejet study trees."
    )
    parser.add_argument(
        "--inputs",
        nargs="+",
        default=None,
        help="Explicit input ROOT files or glob patterns.",
    )
    parser.add_argument(
        "--selection",
        type=str,
        default=DEFAULT_SELECTION,
        help="Selection tag in the input filename, e.g. centPt30-Deta1p3_fwdPt24_VBF-Deta4-Mjj500.",
    )
    parser.add_argument(
        "--inputDir",
        type=str,
        default=DEFAULT_INPUT_DIR,
        help="Directory searched when --selection is used.",
    )
    parser.add_argument(
        "--samplePrefix",
        type=str,
        default="CentralVBFHTo2B",
        help="Filename prefix used when auto-discovering inputs from --selection.",
    )
    parser.add_argument(
        "--tree",
        type=str,
        default="Events",
        help="Input tree name.",
    )
    parser.add_argument(
        "--outdir",
        type=str,
        default=DEFAULT_OUTPUT_DIR,
        help="Directory where fit-ready ROOT files are written.",
    )
    parser.add_argument(
        "--algoName",
        type=str,
        default=None,
        help="Deprecated compatibility option. Ignored in the default naming scheme.",
    )
    parser.add_argument(
        "--binWidth",
        type=float,
        default=5.0,
        help="Histogram bin width in GeV.",
    )
    parser.add_argument(
        "--xmaxFactor",
        type=float,
        default=2.1,
        help="Histogram upper edge relative to the signal mass.",
    )
    return parser.parse_args()


def input_patterns(args):
    if args.inputs:
        return args.inputs
    return [
        os.path.join(
            args.inputDir,
            f"{args.samplePrefix}_M*_centralWideJetVBF_scoutNano_{args.selection}.root",
        )
    ]


def expand_inputs(patterns):
    files = []
    for pattern in patterns:
        matches = sorted(glob.glob(pattern))
        if matches:
            files.extend(matches)
        elif os.path.exists(pattern):
            files.append(pattern)
    files = sorted(set(files))
    if not files:
        raise RuntimeError("No input files found.")
    return files


def infer_mass(path):
    match = re.search(r"_M(\d+)_", os.path.basename(path))
    return int(match.group(1)) if match else -1


def infer_selection_tag(path):
    base = os.path.basename(path)
    match = re.search(
        r"_centralWideJetVBF_scoutNano_(.+)\.root$",
        base,
    )
    return match.group(1) if match else "unknown"


def make_hist(tree, name, expr, selection, nbins, xmin, xmax, xtitle):
    hist = ROOT.TH1F(name, f";{xtitle};Events / 5 GeV", nbins, xmin, xmax)
    hist.Sumw2()
    tree.Draw(f"{expr}>>{name}", selection, "goff")
    hist.SetDirectory(0)
    return hist


def write_metadata(root_file, source_path, selection_tag):
    source = ROOT.TNamed("source_file", source_path)
    selection = ROOT.TNamed("selection_tag", selection_tag)
    source.Write()
    selection.Write()


def process_file(path, args):
    mass = infer_mass(path)
    if mass < 0:
        raise RuntimeError(f"Could not infer mass from input file: {path}")

    selection_tag = infer_selection_tag(path)
    root_file = ROOT.TFile.Open(path)
    if not root_file or root_file.IsZombie():
        raise RuntimeError(f"Could not open ROOT file: {path}")

    tree = root_file.Get(args.tree)
    if tree is None:
        root_file.Close()
        raise RuntimeError(f"Could not find tree '{args.tree}' in {path}")

    xmax = max(700.0, args.xmaxFactor * float(mass))
    nbins = max(1, int(round(xmax / args.binWidth)))

    histograms = {
        "h_widejet_peak_centralFirst_final_5GeV": make_hist(
            tree,
            f"h_widejet_peak_centralFirst_final_5GeV_m{mass}",
            "widejet_pair_mass",
            "pass_trigger_baseline == 1 && isVBF == 1",
            nbins,
            0.0,
            xmax,
            "Wide-jet dijet mass [GeV]",
        ),
    }

    root_file.Close()

    empty_histograms = [name for name, hist in histograms.items() if hist.Integral() <= 0]
    if "h_widejet_peak_centralFirst_final_5GeV" in empty_histograms:
        raise RuntimeError(f"Central VBF histogram is empty for {path}")

    os.makedirs(args.outdir, exist_ok=True)
    outpath = os.path.join(args.outdir, f"vbf-m{mass}.root")
    out_file = ROOT.TFile(outpath, "RECREATE")
    for name, hist in histograms.items():
        hist.SetName(name)
        hist.Write()
    write_metadata(out_file, path, selection_tag)
    out_file.Close()

    return {
        "input": path,
        "output": outpath,
        "mass": mass,
        "selection": selection_tag,
        "entries": int(histograms["h_widejet_peak_centralFirst_final_5GeV"].Integral()),
    }


def main():
    args = parse_args()
    files = expand_inputs(input_patterns(args))
    results = [process_file(path, args) for path in files]
    results.sort(key=lambda row: row["mass"])

    print(f"{'Mass':>6s} {'Entries':>10s}  Output")
    for row in results:
        print(f"{row['mass']:6d} {row['entries']:10d}  {row['output']}")

    print(f"\n[INFO] Wrote {len(results)} fit-ready ROOT files to: {args.outdir}")


if __name__ == "__main__":
    main()
