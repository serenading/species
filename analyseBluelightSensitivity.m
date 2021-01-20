%% Script analyses bluelight sensitivity by looking at significant feature changes, including motion state
% author: @serenading. Jan 2020

clear
close all

addpath('../AggScreening/')

%% Set parameters
feats = {'motion_mode_forward'}; % cell array containing feature names as strings.
strain = 'N2'; % 'N2','CB4856','MY23','QX1410','VX34','NIC58','JU1373'; 
n_subsample = 1; % number of replicates per strain to include.
n_nonFeatVar = 33; % the first n columns of the feature table that do not contain features. =33
frameRate = 25;
smoothWindow = frameRate; % window to smooth feature over, in frames (data acquired at 25fps)
resultsDir = '/Volumes/Ashur DT2/species/Results/';

%% Find files for specified strain
% load features table
featureTable = readtable('/Users/sding/OneDrive - Imperial College London/species/Results/fullFeaturesTable_20201218_184325.csv');
% get light condition
light_condition = getLightcondition(featureTable);
% filter for prestim files for the specified strain
strainInd = find(strcmp(featureTable.strain_name,strain) & strcmp(light_condition,'prestim'));
% subsample a few files for plotting
strainInd = datasample(strainInd,n_subsample,'Replace',false);