%% Script plots selected timeseries features across the three light conditions
% author: @serenading. Jan 2021.

clear
close all

addpath('../AggScreening/')

% TODO: Currently plotting one experimental replicate at a time
% (n_subsample = 1). Add ways to average across several replicates with
% confidence interval shading

%% Set parameters
feats = {'motion_mode'};%_forward_fraction','motion_mode_paused_fraction','motion_mode_backward_fraction'}; % {'angular_velocity_head_tip','angular_velocity_head_base','speed_midbody'}; % cell array containing feature names as strings.
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

%% Go through each file, find matching light condition files, and get timeseries feature
% go through each file
for fileCtr = 1:numel(n_subsample)
    % get prestim file index
    prestimfileIdx = strainInd(fileCtr);
    % get matching bluelight/poststim file indices
    [bluelightfileIdx,poststimfileIdx,well] = findMatchingFileInd(prestimfileIdx,featureTable);
    % get full filenames
    prestimfilename = [resultsDir,featureTable.filename{prestimfileIdx},'/metadata_featuresN.hdf5'];
    bluelightfilename = [resultsDir,featureTable.filename{bluelightfileIdx},'/metadata_featuresN.hdf5'];
    poststimfilename = [resultsDir,featureTable.filename{poststimfileIdx},'/metadata_featuresN.hdf5'];
    % load features files
    features_prestim = h5read(prestimfilename,'/timeseries_data');
    features_bluelight = h5read(bluelightfilename,'/timeseries_data');
    features_poststim = h5read(poststimfilename,'/timeseries_data');
    % get logical index for the relevant well 
    wellsLogInd_prestim = strcmp(cellstr(features_prestim.well_name'),well);
    wellsLogInd_bluelight = strcmp(cellstr(features_bluelight.well_name'),well); 
    wellsLogInd_poststim = strcmp(cellstr(features_poststim.well_name'),well); 
    
    %% plot features
    featfigure = figure; hold on
    title(['Sample timeseries feature for ' strain])
    % go through each feature (rows of subplot)
    for featCtr = 1:numel(feats)
        feat = feats{featCtr};
        % prestim
        ax1 = subplot(numel(feats),3,(featCtr-1)*3+1);
        x = features_prestim.timestamp(wellsLogInd_prestim);
        y = features_prestim.(feat)(wellsLogInd_prestim);
        plot(x,smoothdata(y,'movmedian',smoothWindow),'.');
        title('prestim'); xlabel('frames'); ylabel(feat,'Interpreter','none'); xlim([0 8000])
        % bluelight
        ax2 = subplot(numel(feats),3,(featCtr-1)*3+2);
        x = features_bluelight.timestamp(wellsLogInd_bluelight);
        y = features_bluelight.(feat)(wellsLogInd_bluelight);
        plot(x,smoothdata(y,'movmedian',smoothWindow),'.');
        title('bluelight'); xlabel('frames'); ylabel(feat,'Interpreter','none'); xlim([0 8000])
        % poststim
        ax3 = subplot(numel(feats),3,(featCtr-1)*3+3);
        x = features_poststim.timestamp(wellsLogInd_poststim);
        y = features_poststim.(feat)(wellsLogInd_poststim);
        plot(x,smoothdata(y,'movmedian',smoothWindow),'.');
        title('poststim'); xlabel('frames'); ylabel(feat,'Interpreter','none'); xlim([0 8000])
        % format plot
        linkaxes([ax1 ax2 ax3],'xy')
    end
end

