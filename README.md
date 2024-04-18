# Final Fits (lite)

**NOTE** If using Final Fits with the output of HiggsDNA, then please switch to [this](https://github.com/cms-analysis/flashggFinalFit/tree/dev_higgsdnafinalfit) branch

Welcome to the new Final Fits package. Here lies a a series of scripts which are used to run the final stages of the CMS Hgg analysis: signal modelling, background modelling, datacard creation, final statistical interpretation and final result plots.

Slides from the flashgg tutorial series can be found [here](https://indico.cern.ch/event/963619/contributions/4112177/attachments/2151275/3627204/finalfits_tutorial_201126.pdf)

## Download and setup instructions

Setup in lxplus7. As explained in the [Combine page](https://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/latest/#within-cmssw-recommended-for-cms-users), CMSSW_11_3_X runs on slc7, which can be setup using apptainer.
At first, install Combine as per the documentation.
```
cmsrel CMSSW_11_3_4
cd CMSSW_11_3_4/src
cmsenv
git fetch origin
git checkout v9.2.1
scramv1 b clean; scramv1 b # always make a clean build
```

Then, install flashggFinalFit packages:
```
git clone -b EF_lowMassAnalysis git@github.com:elfontan/flashggFinalFit.git
cd flashggFinalFit/
```

In every new shell run the following to add `tools/commonTools` and `tools/commonObjects` to your `${PYTHONPATH}`:
```
cmsenv
source setup.sh
```

## Contents
The Finals Fits package contains several subfolders which are used for the following steps:

* Create the Signal Model (see `Signal` dir)
* Create the Background Model (see `Background` dir)
* Generate a Datacard (see `Datacard` dir)
* Running fits with combine (see `Combine` dir)
* Scripts to produce plots (see `Plots` dir)

The signal modelling, background modelling and datacard creation can be ran in parallel. Of course the final fits (`Combine`) requires the output of these three steps. In addition, the scripts in the `Trees2WS` dir are a series of lightweight scripts for converting standard ROOT trees into a RooWorkspace that can be read by the Final Fits package.

Finally, the objects and tools which are common to all subfolders are defined in the `tools` directory. If your input workspaces differ from the flashgg output workspace structure, then you may need to change the options here.

Each of the relevant folders are documented with specific `README.md` files. Some (temporary) instructions can be found in this [google docs](https://docs.google.com/document/d/1NwUrPvOZ2bByaHNqt_Fr6oYcP7icpbw1mPlw_3lHhEE/edit)
