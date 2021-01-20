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
n_subsample = 175; % number of replicates per strain to sample.
n_nonFeatVar = 33; % the first n columns of the feature table that do not contain features. =33
bluelightInterval = [60,70; 160,170; 260,270];
windowDuration = 10; % 10 (default based on Ziwei's time window analysis). Time window in seconds to retain for feature analysis.
windownames = [0:8]; % windows are named 0 through 8
resultsDir = '/Users/sding/OneDrive - Imperial College London/species/Results';
trimNaNFiles = true;

%% Get feature extraction windows around
[midpointAllwindows, prelight, bluelight, postlight] = getBluelightFeatWindows(bluelightInterval,windowDuration);
% Enter the resultant windows (in seconds) into Tierpsy feature summariser for feature
% window extraction (takes about a minute per file on my macpro)
% [50:60,65:75,75:85,150:160,165:175,175:185,250:260,265:275,275:285]

%% Go through each strain
for strainCtr = 1:numel(strains)
    strain = strains{strainCtr};
    
    %% Preallocate n_sample x n_features x 9 windows array to hold feature values
    featVals = NaN(n_subsample,numel(feats),numel(windownames));
    
    %% Extract feature values from bluelight files for the specified strain for the first window
    % load features table from the first window
    featureTable = readtable([resultsDir '/fullFeaturesTable_' extractStamp '_window_0.csv']);
    % get light condition
    light_condition = getLightcondition(featureTable);
    % filter for bluelight files for the specified strain
    fileInd = find(strcmp(featureTable.strain_name,strain) & strcmp(light_condition,'bluelight'));
    % subsample a few files for plotting
    fileInd = datasample(fileInd,n_subsample,'Replace',false);
    % get imgstorenames and wellnames for the subsampled files - we will use
    % these and bluelight condition to find the same files across all nine
    % windows
    imgstorenames = featureTable.filename(fileInd);
    wellnames = featureTable.well_name(fileInd);
    % extract feature values
    for featCtr = 1:numel(feats)
        feat = feats{featCtr};
        featVals(:,featCtr,1) = featureTable.(feat)(fileInd);
    end
    
    %% Find matching bluelight files for the subsequent windows and extract feature values
    % go through each subsequent window
    for windowCtr = 2:9 % go through the second to the final files one by one
        window = windownames(windowCtr);
        % load features table from the first window
        featureTable = readtable([resultsDir '/fullFeaturesTable_' extractStamp '_window_' num2str(window) '.csv']);
        % get light condition
        light_condition = getLightcondition(featureTable);
        % reset fileInd variable
        fileInd = NaN(n_subsample,1);
        % get file indices that match the files chosen for window_0, as the
        % indices may not necessarily be the same across all nine windows (they
        % are not)
        for sampleCtr = 1:n_subsample
            fileIdx = find(strcmp(featureTable.filename,imgstorenames{sampleCtr}) &...
                strcmp(featureTable.well_name,wellnames{sampleCtr}) &...
                strcmp(light_condition,'bluelight'));
            if ~isempty(fileIdx)
                fileInd(sampleCtr) = fileIdx;
            end
        end
        % trim off empty fileInd
        fileInd = fileInd(~isnan(fileInd));
        % extract feature values
        for featCtr = 1:numel(feats)
            feat = feats{featCtr};
            featVals(1:numel(fileInd),featCtr,windowCtr) = featureTable.(feat)(fileInd);
        end
    end
    
    %% Remove files that have NaN value in any window
    if trimNaNFiles
        featValsKeep = featVals;
        featValsSmall = squeeze(featVals(:,1,:));
        [rows2drop,~] = find(isnan(featValsSmall));
        rowLogInd = true(size(featVals,1),1);
        rowLogInd(unique(rows2drop)) = false;
        featVals = featVals(rowLogInd,:,:);
    end
    
    %% Plot features
    figure; hold on
    colors = {'c','m','k'};
    
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
    xlim([30 400])
    xlabel('time (s)')
    
    % Save
    savefig([resultsDir '/motion_state/' strain '.fig'])
end