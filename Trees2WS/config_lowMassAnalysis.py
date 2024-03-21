# Input config file for running trees2ws

trees2wsCfg = {

    # Name of RooDirectory storing input tree
    'inputTreeDir':'tagsDumper/trees',
    
    # Variables to be added to dataframe: use wildcard * for common strings
    'mainVars':["CMS_hgg_mass","weight", "vtxdz"], # Vars to add to nominal RooDatasets
    'dataVars':["CMS_hgg_mass","weight"], # Vars for data workspace (trees2ws_data.py script)
    'stxsVar':'', # Var for STXS splitting: if using option doSTXSSplitting
    'notagVars':["weight"], # Vars to add to NOTAG RooDataset
    'systematicsVars':["CMS_hgg_mass","weight"], # Variables to add to sytematic RooDataHists
    'theoryWeightContainers':{}, # Theory weights to add to nominal + NOTAG RooDatasets, value corresponds to number of weights (0-N)
    
    # List of systematics: use string YEAR for year-dependent systematics
    'systematics':[],
    
    # Analysis categories: python list of cats or use 'auto' to extract from input tree
    'cats':'auto'
}
