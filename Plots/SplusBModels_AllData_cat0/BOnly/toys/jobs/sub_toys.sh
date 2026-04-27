#!/bin/bash
ulimit -s unlimited
set -e
cd /afs/cern.ch/work/a/atsatsos/ULLowmassFGG/CMSSW_11_3_4/src/flashggFinalFit/Plots/SplusBModels_AllData_cat0/toys
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scramv1 runtime -sh`

itoy=$1

#Generate command
echo Generating...
combine /afs/cern.ch/work/a/atsatsos/ULLowmassFGG/CMSSW_11_3_4/src/flashggFinalFit/Plots/SplusBModels_AllData_cat0/higgsCombine_initialSnapshot.MultiDimFit.mH12.root -m 12.000 -M GenerateOnly --expectSignal=0 --saveWorkspace --saveToys --toysFrequentist --bypassFrequentistFit -t 1 --setParameters r=0.000 -s -1 -n _${itoy}_gen_step --snapshotName MultiDimFit

#Fit command
echo Fitting...
mv higgsCombine_${itoy}_gen_step*.root gen_${itoy}.root
combine /afs/cern.ch/work/a/atsatsos/ULLowmassFGG/CMSSW_11_3_4/src/flashggFinalFit/Plots/SplusBModels_AllData_cat0/higgsCombine_initialSnapshot.MultiDimFit.mH12.root -t 1 --toysFile=gen_${itoy}.root -m 12.000 -M MultiDimFit -P r --floatOtherPOIs=1 --expectSignal=0 --saveWorkspace --toysFrequentist --bypassFrequentistFit --setParameters r=0.000 -n _${itoy}_fit_step --cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_hideConstants --X-rtd MINIMIZER_multiMin_maskConstraints --X-rtd MINIMIZER_multiMin_maskChannels=2

#Throw command
echo Throwing...
mv higgsCombine_${itoy}_fit_step*.root fit_${itoy}.root
combine fit_${itoy}.root -m 12.000 --snapshotName MultiDimFit -M GenerateOnly --expectSignal=0 --saveToys --saveWorkspace --toysFrequentist --bypassFrequentistFit -t -1 -s -1 -n _${itoy}_throw_step --setParameters r=0

mv higgsCombine_${itoy}_throw_step*.root toy_${itoy}.root
rm gen_${itoy}.root fit_${itoy}.root
