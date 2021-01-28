%% This script uses three files to generate a features summary table from Hydra data.
% author: serenading. Jan 2021

clear
close all

%% Set parameters
% set which feature extraction timestamp to use
extractStamp = '20210112_105808'; % 20201218_184325 for standard feature extraction, 20210112_105808 for filtered data, 20210119_073010 filtered data with bluelight windows.
% set project directory
projectDir = '/Users/sding/OneDrive - Imperial College London/species';
% have wells been annotated?
annotatedWells = true;
% which parts of the filenames_summary_tierpsy full file path should be used to reconstruct the imgstore name?
filenameparts = [6,7];

%% Generate and save featureTable files
if strcmp(extractStamp,'20210119_073010') % '20210119_073010' are for multiple bluelight windows so loop through the function a few times for each window
    windows = [0:8];
    for windowCtr = 1:numel(windows)
        extractStamp = ['20210119_073010_window_' num2str(windows(windowCtr))];
        featureTable = getFeatureTable(extractStamp,projectDir,annotatedWells,filenameparts);
        writetable(featureTable,[projectDir '/Results/fullFeaturesTable_' extractStamp '.csv']);
    end
else
    featureTable = getFeatureTable(extractStamp,projectDir,annotatedWells,filenameparts);
    writetable(featureTable,[projectDir '/Results/fullFeaturesTable_' extractStamp '.csv']);
end