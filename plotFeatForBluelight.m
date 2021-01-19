%% Script plots selected features across the blue light condition videos. 
% author: @serenading. Jan 2021.

clear
close all

addpath('../AggScreening/')

%% Set parameters

feats = {'',''};
strain = 'N2'; % 'N2','CB4856','MY23','QX1410','VX34','NIC58','JU1373'; 
n_subsample = 1; % number of replicates per strain to include.
n_nonFeatVar = 33; % the first n columns of the feature table that do not contain features. =33
bluelightInterval = [60,70; 160,170; 260,270]; 
windowDuration = 10; % 10 (default based on Ziwei's time window analysis). Time window in seconds to retain for feature analysis.
resultsDir = '/Volumes/Ashur DT2/species/Results/';

%% Get feature extraction windows around 
[prelight, bluelight, postlight] = getBluelightFeatWindows(bluelightInterval,windowDuration);
% Enter the resultant windows (in seconds) into Tierpsy feature summariser for feature
% window extraction (takes about a minute per file on my macpro)
% [50:60,65:75,75:85,150:160,165:175,175:185,250:260,265:275,275:285]

%% Find bluelight files for specified strain
% load features table
featureTable = readtable('/Users/sding/OneDrive - Imperial College London/species/Results/fullFeaturesTable_20201218_184325.csv');
% get light condition
light_condition = getLightcondition(featureTable);
% filter for prestim files for the specified strain
strainInd = find(strcmp(featureTable.strain_name,strain) & strcmp(light_condition,'bluelight'));
% subsample a few files for plotting
strainInd = datasample(strainInd,n_subsample,'Replace',false);

