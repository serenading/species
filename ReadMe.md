#### Features summaries files to use . 

**20201218_184325**: feature summaries have 3016 features (without dorsal-ventral features) with short, matlab-friendly names.  
**20210112_105808**: feature summaries have 3016 features (without dorsal-ventral features) with short, matlab-friendly names; filtered with worm size criteria. Use for most purposes.  
**2021xxxxxx(TBC)**: feature summaries have 3016 features (without dorsal-ventral features) with short, matlab-friendly names; filtered with worm size criteria; generated for various windows around bluelight stimulation. Use for analysing bluelight videos.

# Useful scripts . 

`generateFeatSummary`: First script to run. Combines Tierpsy tables with metadata table to generate a combined *FeatureTable* used for downstream analysis.  
`doPCA`: Script performs PCA analysis, with options to specify which features, which strains, and which light conditions to use.  
`classifySpecies`: Script uses supervised machine learning algorithms to train classifiers for a specified variable based on extracted Tierpsy features. It also has the option to apply sequantial feature selection to identify top features to use for classification.  
`plotFeatAcrossights`: Script plots selected timeseries features across the three light conditions.  
`plotTraj`: Script visualises full frame trajectories from an entire camera view (4x4 wells).  

# Useful functions . 

`findMatchingFileInd`: Function takes prestim file index from featureTable and returns bluelight and poststim file indices, and the well name.  
`getBluelightFeatWindows`: Function generates 3x2 windows for each of the light condition in seconds based on standard Hydra bluelight stimulation conditions.  

# Under development .  

`plotFeatForBluelight`: Script plots selected features across the blue light condition videos.  
`analyzeBlueLightSensitivity`: Script analyses bluelight sensitivity by looking at significant feature changes, including motion state.  