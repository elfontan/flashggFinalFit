# Config file: options for signal fitting

backgroundScriptCfg = {
  
  # Setup
  'inputWSFile':'/eos/user/a/atsatsos/ULFlashGG_Files/UL18_Data_Lowmassxml_v1/BlindedTree/ws/blinded_tree_merged.root', # location of 'allData.root' file
  'cats':'auto', # auto: automatically inferred from input ws
  'ext':'lite', # extension to add to output directory
  'year':'2018', # Use merged when merging all years in category (for plots)

  # Job submission options
  'batch':'local', # [condor,SGE,IC,local]
  'queue':'microcentury' # for condor e.g. microcentury
  
}
