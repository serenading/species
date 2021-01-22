%% Script analyses bluelight sensitivity by looking at significant feature changes between prestim and bluelight videos.
% author: @serenading. Jan 2020

clear
close all

addpath('../AggScreening/')
addpath('../AggScreening/auxiliary/')

%% Set parameters
extractStamp = '20201218_184325'; % 20201218_184325 for standard feature extraction, 20210112_105808 for filtered data
feats = {'motion_mode_forward_fraction','motion_mode_paused_fraction','motion_mode_backward_fraction'};
strain = 'N2'; % 'N2','CB4856','MY23','QX1410','VX34','NIC58','JU1373'; 
n_subsample = 100; % number of replicates per strain to include. Set to NaN to use all files
n_nonFeatVar = 33; % the first n columns of the feature table that do not contain features. =33
resultsDir = '/Users/sding/OneDrive - Imperial College London/species/Results';
bonCorr = true; % apply bonferroni correction for multiple comparisons
n_windows = 2; % two windows: one for prestim, one for bluelight.

%% Find files for the specified strain
% load features table
featureTable = readtable([resultsDir '/fullFeaturesTable_' extractStamp '.csv']);
% get light condition
light_condition = getLightcondition(featureTable);
% filter for prestim files for the specified strain
fileInd = find(strcmp(featureTable.strain_name,strain) & strcmp(light_condition,'prestim'));
% subsample a few files for plotting
if ~isnan(n_subsample)
    fileInd = datasample(fileInd,n_subsample,'Replace',false);
else
    n_subsample = numel(fileInd);
end

%% Go through each file, find matching prestim and bluelight condition file indices

% Preallocate matrix to hold file indices: n x 2 matrix, each column is one time window (prestim vs. bluelight)
allFileInd = NaN(n_subsample,n_windows); 

% Go through each file
for fileCtr = 1:n_subsample
    % get prestim file index
    allFileInd(fileCtr,1) = fileInd(fileCtr);
    % get matching bluelight file index
    [bluelightfileIdx,~,well] = findMatchingFileInd(allFileInd(fileCtr,1),featureTable);
    % some bluelight file index may be empty because worms have all crawled
    % off leaving no extracted features files. Keep the index as NaN in
    % this case for trimming later
    if ~isempty(bluelightfileIdx) 
        allFileInd(fileCtr,2) = bluelightfileIdx;
    end
end

% Remove files that have NaN index in any window
[trimmedFileInd,n_filesDropped] = dropNaNFiles(allFileInd);

%% Extract features for the corresponding time windows
% prellocate 
featVals = NaN(size(trimmedFileInd,1),numel(feats),n_windows); 

% extract values
for windowCtr = 1:n_windows 
    featVals(:,:,windowCtr) = featureTable{trimmedFileInd(:,windowCtr),feats};
end
    
% remove experiments with NaN feature values in any window
[featVals,n_filesDropped,featValsCopy] = dropNaNVals(featVals);

%% Perform statistical analysis