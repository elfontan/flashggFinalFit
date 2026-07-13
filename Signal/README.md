# Parametric Signal Modeling for Run 3 Dijet VBF studies

The streamlined workflow used to build the parametric signal model starting from a certain set of VBF-selection reco outputs is documented here.

The guiding idea is simple:

- copy the ROOT files for one chosen VBF selection into `Output-VBFSelection/`;
- prepare fit inputs with standard naming for each reference mass point used to build the model;
- run the parametric fit on those prepared files.

## Inputs

For a given VBF working point, copy the corresponding reco-study outputs into:

`Output-VBFSelection/`

For example, if you want to consider 300-500-750-1000-2000-3000 GeV mass hypotheses for the set of selections `centPt30-Deta1p3_fwdPt24_VBF-Deta4-Mjj500`, you will have a set of files like this:

- `CentralVBFHTo2B_M300_centralWideJetVBF_scoutNano_centPt30-Deta1p3_fwdPt24_VBF-Deta4-Mjj500.root`
- `CentralVBFHTo2B_M500_centralWideJetVBF_scoutNano_centPt30-Deta1p3_fwdPt24_VBF-Deta4-Mjj500.root`
- `CentralVBFHTo2B_M750_centralWideJetVBF_scoutNano_centPt30-Deta1p3_fwdPt24_VBF-Deta4-Mjj500.root`
- `CentralVBFHTo2B_M1000_centralWideJetVBF_scoutNano_centPt30-Deta1p3_fwdPt24_VBF-Deta4-Mjj500.root`
- `CentralVBFHTo2B_M2000_centralWideJetVBF_scoutNano_centPt30-Deta1p3_fwdPt24_VBF-Deta4-Mjj500.root`
- `CentralVBFHTo2B_M3000_centralWideJetVBF_scoutNano_centPt30-Deta1p3_fwdPt24_VBF-Deta4-Mjj500.root`

The only physics content used for the signal model is the wide-jet invariant mass built after full selection:

- `widejet_pair_mass`
- with selection `pass_trigger_baseline == 1 && isVBF == 1`

No alternative test histograms are part of the nominal workflow here.

## Step 0: prepare fit-ready ROOT files

Run:

```bash
python3 step0-ParametricSigModeling_prepare_signalPeakInputs_VBF.py
```

By default this script:

- searches in `Output-VBFSelection/`;
- assumes the current selection
  `centPt30-Deta1p3_fwdPt24_VBF-Deta4-Mjj500`;
- reads the `Events` tree;
- builds the wide-jet full-selection dijet mass histogram;
- writes one ROOT file per mass into
  `Output-VBFSelection/PreparedFitInputs/`.

The output naming is selection-transparent:

- `vbf-m300.root`
- `vbf-m500.root`
- `vbf-m750.root`
- `vbf-m1000.root`
- `vbf-m2000.root`
- `vbf-m3000.root`

Each output file contains:

- `h_widejet_peak_centralFirst_final_5GeV`

This is the histogram used as input to the parametric fit.

If a different working point is copied into `Output-VBFSelection/`, run:

```bash
python3 step0-ParametricSigModeling_prepare_signalPeakInputs_VBF.py \
  --selection centPt30-Deta1p3_fwdPt21_VBF-Deta4p5-Mjj450
```

## Step 1: run the parametric fit

Once the prepared inputs are available, run the script `step1-ParametricSigModeling_fit_signalPeak_VBF.py` with options like:

* Inputs directory, e.g. `--indir Output-VBFSelection/PreparedFitInputs`
* Name of the histo: `--hist h_widejet_peak_centralFirst_final_5GeV`
* Output directory
* Order of the polynomial function to interpolate the parameters, e.g. `--mhPolyOrder 1`
   * Specific orders for each parameters can be specified, e.g.: `--paramOrders mean:1,sigmaL:3,sigmaR:2,alphaL:2,nL:3,alphaR:3,nR:3`
* Reference mass taken for the parametrisation (arbitrary, usually a central value is convenient), e.g.: `--referenceMass 750``
* List of masses, e.g. `--massPoints 500,750,1000,2000,3000`
```

EXAMPLE: best configuration for now (mass > 500 GeV)
```
python3 step1-ParametricSigModeling_fit_signalPeak_VBF.py
  --indir Output-VBFSelection/PreparedFitInputs \
  --hist h_widejet_peak_centralFirst_final_5GeV \
  --outdir /eos/user/e/elfontan/www/dijetAnaRun3/SIGModelling/VBF-CAT/PARAMETRIC-SIGModel/MASS-Above500 \
  --drawWindowHighScale 2.0 --useWidePlotWindow --fitRangeMax 5000  \
  --massPoints 500,750,1000,2000,3000 \
  --mhPolyOrder 1 \
  --referenceMass 750 \
  --paramOrders mean:1,sigmaL:2,sigmaR:2,alphaL:1,alphaR:1,nL:2,nR:1 --saveWorkspace	
```

EXAMPLE: Best configuration for now (mass < 500 GeV)
```
python3 step1-ParametricSigModeling_fit_signalPeak_VBF.py
  --indir Output-VBFSelection/PreparedFitInputs \
  --hist h_widejet_peak_centralFirst_final_5GeV \
  --outdir /eos/user/e/elfontan/www/dijetAnaRun3/SIGModelling/VBF-CAT/PARAMETRIC-SIGModel/MASS-Below500 \
  --drawWindowHighScale 2.0 --useWidePlotWindow --fitRangeMax 1000  \
  --massPoints 200,300,500 \
  --mhPolyOrder 2 \
  --referenceMass 300 \
  --paramOrders mean:1,sigmaL:2,sigmaR:1,alphaL:2,alphaR:2,nL:1,nR: \
  --massWindowOverrides 500:0.7:1.3 --saveWorkspace


Useful knobs in step 1:

- `--fitWindowLowScale` and `--fitWindowHighScale`
  control the actual single-mass fit window;
- `--massWindowOverrides`
  lets specific masses use a dedicated fit window without changing all the
  others; use the syntax `mass:low:high`, for example
  `300:0.45:1.45,500:0.45:1.55`;
- `--drawWindowLowScale` and `--drawWindowHighScale`
  control the nominal visible mass range;
- `--useWidePlotWindow`
  forces the plots to use that full configured visible range instead of the
  tighter automatic view around the fitted peak;
- `--fitRangeMax`
  sets the global upper edge allowed for the fit and display windows.


## Step 2: scan the saved workspace

If the parametric fit is run with `--saveWorkspace`, the resulting
RooWorkspace can be scanned to visualize the interpolated signal shapes across
mass hypotheses.

Example usage:

```bash
python3 step2-ParametricSigModeling_workspaceScan.py \
  --input /eos/user/e/elfontan/www/dijetAnaRun3/SIGModelling/VBF-CAT/PARAMETRIC-SIGModel/PolOrder1-massRef750-tailTuned_noM300/parametric_parametric_workspace.root \
  --outdir /eos/user/e/elfontan/www/dijetAnaRun3/SIGModelling/VBF-CAT/PARAMETRIC-SIGModel/ \
  --outname workspace_scan_400to2000 \
  --mhMin 400 \
  --mhMax 2000 \
  --mhStep 200 \
  --label "Parametric VBF signal model"
```

python3 step2-ParametricSigModeling_workspaceScan.py   --input /eos/user/e/elfontan/www/dijetAnaRun3/SIGModelling/VBF-CAT/PARAMETRIC-SIGModel/MASS-Above500/parametric_parametric_workspace.root   --outdir /eos/user/e/elfontan/www/dijetAnaRun3/SIGModelling/VBF-CAT/PARAMETRIC-SIGModel/   --outname workspace_scan_500to2900   --mhMin 500   --mhMax 3000   --mhStep 200   --label "Parametric VBF signal model"

python3 step2-ParametricSigModeling_workspaceScan.py   --input /eos/user/e/elfontan/www/dijetAnaRun3/SIGModelling/VBF-CAT/PARAMETRIC-SIGModel/MASS-Below500/parametric_parametric_workspace.root   --outdir /eos/user/e/elfontan/www/dijetAnaRun3/SIGModelling/VBF-CAT/PARAMETRIC-SIGModel/   --outname workspace_scan_300to500   --mhMin 300   --mhMax 500   --mhStep 500  --xmax 1000  --label "Parametric VBF signal model"


Useful options for step 2:

- `--workspace`
  explicitly choose the RooWorkspace name if more than one is stored;
- `--pdf`
  explicitly choose the PDF name if needed;
- `--xmin` and `--xmax`
  override the plotted dijet-mass range;
- `--noUnitNorm`
  draw the sampled curves without normalizing each one to unit area.

To make Step 2 available, run Step 1 with `--saveWorkspace` so that the
parametric RooWorkspace is written to the output directory.

## Directory layout

Inside `Signal/` the intended structure is:

- `Output-VBFSelection/`
  raw reco-study ROOT files for one selected VBF working point;
- `Output-VBFSelection/PreparedFitInputs/`
  fit-ready ROOT files produced by step 0;
- `plots_VBFSignalFits_parametric/`
  outputs of the parametric fit.

## Workflow summary

1. choose one VBF working point;
2. copy the corresponding `CentralVBFHTo2B_M*.root` files into
   `Output-VBFSelection/`;
3. run step 0;
4. run the parametric fitter.

When testing a new selection, the workflow stays the same;
only the copied input files and, if needed, the `--selection` argument of step 0 need to be updated.
