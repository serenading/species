%% Script plots selected features across the blue light condition videos.
% author: @serenading. Jan 2021.

clear
close all

addpath('../AggScreening/')
addpath('../AggScreening/auxiliary/')

%% Set parameters

extractStamp = '20210119_073010'; % '20210119_073010' feature summaries have multiple bluelight windows.
feats = {'motion_mode_forward_fraction','motion_mode_paused_fraction','motion_mode_backward_fraction'};
strains = {'N2','CB4856','MY23','QX1410','VX34','NIC58','JU1373'}; % 'N2','CB4856','MY23','QX1410','VX34','NIC58','JU1373';
n_subsample = NaN; % number of replicates per strain to sample. Set to NaN to use all files
bluelightInterval = [60,70; 160,170; 260,270];
windowDuration = 10; % 10 (default based on Ziwei's time window analysis). Time window in seconds to retain for feature analysis.
windownames = [0:8]; % windows are named 0 through 8
resultsDir = '/Users/sding/OneDrive - Imperial College London/species/Results';

%% Get feature extraction windows around
[midpointAllwindows, prelight, bluelight, postlight] = getBluelightFeatWindows(bluelightInterval,windowDuration);
n_windows = numel(windownames);
assert(numel(midpointAllwindows) == n_windows)
% Enter the resultant windows (in seconds) into Tierpsy feature summariser for feature
% window extraction (takes about a minute per file on my macpro)
% [50:60,65:75,75:85,150:160,165:175,175:185,250:260,265:275,275:285]

%% Go through each strain
for strainCtr = 2:numel(strains)
    strain = strains{strainCtr};

    %% Find file indices (window0) for the specified strain
    % load features table
    featureTable = readtable([resultsDir '/fullFeaturesTable_' extractStamp '_window_0.csv']);
    % get light condition
    light_condition = getLightcondition(featureTable);
    % filter for bluelight files for the specified strain
    fileInd = find(strcmp(featureTable.strain_name,strain) & strcmp(light_condition,'bluelight'));
    % subsample a few files for plotting
    if ~isnan(n_subsample)
        fileInd = datasample(fileInd,n_subsample,'Replace',false);
    end
    n_sample = numel(fileInd);

    %%  Preallocate matrix to hold file indices: n x 9 matrix, each column is one time window (9 windows total)
    allFileInd = NaN(n_sample,n_windows); 
    allFileInd(:,1) = fileInd;
    
    %% Get imgstorenames and wellnames for the subsampled files 
    % (we will use these and bluelight condition to find the same files
    % across all nine windows)
    imgstorenames = featureTable.filename(fileInd);
    wellnames = featureTable.well_name(fileInd);
    
    %% Find matching bluelight files for the subsequent windows 
    % go through each subsequent window
    for windowCtr = 2:9 
        % get window name as output by Tierpsy
        window = windownames(windowCtr);
        % load features table for that window
        featureTable = readtable([resultsDir '/fullFeaturesTable_' extractStamp '_window_' num2str(window) '.csv']);
        % get light condition
        light_condition = getLightcondition(featureTable);
        % get file indices that match the files chosen for window_0, as the
        % indices are not the same across all nine windows
        for sampleCtr = 1:n_sample
            fileIdx = find(strcmp(featureTable.filename,imgstorenames{sampleCtr}) &...
                strcmp(featureTable.well_name,wellnames{sampleCtr}) &...
                strcmp(light_condition,'bluelight'));
            if ~isempty(fileIdx)
                allFileInd(sampleCtr,windowCtr) = fileIdx;
            end
        end
    end
    
    %% Remove files that have NaN index in any window
    [trimmedFileInd,n_filesDropped] = dropNaNFiles(allFileInd);
    
    %% Extract features for the corresponding time windows
    % preallocate
    featVals = NaN(size(trimmedFileInd,1),numel(feats),n_windows);
    
    % go through each time window
    for windowCtr = 1:n_windows
        % load features table
        window = windownames(windowCtr);
        featureTable = readtable([resultsDir '/fullFeaturesTable_' extractStamp '_window_' num2str(window) '.csv']);
        % extract feature values
        featVals(:,:,windowCtr) = featureTable{trimmedFileInd(:,windowCtr),feats};
    end
    
    % remove experiments with NaN feature values in any window
    [featVals,n_filesDropped,featValsCopy] = dropNaNVals(featVals);
    
    %% Plot features
    figure; hold on
    colors = {'m','c','k','r'};
    
    % Go through each feature and plot as shaded error bar
    for featCtr = 1:numel(feats)
        feat = feats{featCtr};
        vals = squeeze(featVals(:,featCtr,:));
        shadedErrorBar(midpointAllwindows,vals,{@nanmean,@(x) nanstd(x)/sqrt(numel(x))},...
            'lineprops',{[colors{featCtr} '-o'],'markerfacecolor',colors{featCtr}},'transparent',1); hold on
    end
    
    % Add bluelight interval shading
    x = horzcat(bluelightInterval,fliplr(bluelightInterval))';
    yL = horzcat(get(gca,'YLim'));
    y = [yL(1),yL(1),yL(2),yL(2)]';
    y = horzcat(y,y,y);
    patch(x,y,'blue','EdgeColor','none','FaceAlpha',0.3)
    
    % Format
    legend(horzcat(feats,{'bluelight pulse'}) ,'Interpreter','none')
    title([strain sprintf('\n') 'n = ' num2str(size(featVals,1))])
    xlim([30 360])
    xlabel('time (s)')
    
    % Save
    savefig([resultsDir '/motion_mode/' strain '.fig'])
end