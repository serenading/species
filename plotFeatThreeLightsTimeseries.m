%% Script reads featuresN.hdf5 files and plots selected timeseries features across the three light conditions
% (as opposed to using extracted features values in each light window as in plotFeatThreeLights.m).
% author: @serenading. Jan 2021.

clear
close all

addpath('../AggScreening/')
addpath('../AggScreening/auxiliary/')

% TODO: make functions for calculating reversal frequencies etc. from
% extracted time series data

%% Set parameters
extractStamp = '20201218_184325'; % 20201218_184325 for standard feature extraction, 20210112_105808 for filtered data
feats = {'speed_midbody','angular_velocity_head_tip','motion_mode'};%_forward_fraction','motion_mode_paused_fraction','motion_mode_backward_fraction'}; % {'angular_velocity_head_tip','angular_velocity_head_base','speed_midbody'}; % cell array containing feature names as strings.
lightConditions = {'prestim','bluelight','poststim'}; % {'prestim','bluelight','poststim'}
strains = {'N2'}; % 'N2','CB4856','MY23','QX1410','VX34','NIC58','JU1373';
n_subsample = 20; % number of replicates per strain to sample.
frameRate = 25;
smoothWindow = 30; % window to smooth feature over, in seconds
resultsDir = '/Volumes/Ashur DT2/species/Results/';
n_frames = 6*60*frameRate + 10; % maximum duration per video is 6 min. Plus 10 extra frames because the recordings don't always finish at the precise end frame

%% Load features table
featureTable = readtable(['/Users/sding/OneDrive - Imperial College London/species/Results/fullFeaturesTable_' extractStamp '.csv']);

%% Load saved matching file indices and feature sets
% Load matching file indices across three light conditions
load(['matchingFileInd/threelight_' extractStamp '.mat'],'allFileInd','allFileIndWindows');

%% Go through each strain to get file indices
for strainCtr = 1:numel(strains)
    
    % Get matching file indices for this strain
    strain = strains{strainCtr};
    disp(['Extracting feature values for ' strain  ' ...'])
    strainAllFileInd = allFileInd.(strain);
    
    % If specified, subsample a few files for plotting
    if ~isnan(n_subsample)
        allRowInd = 1:size(strainAllFileInd,1);
        keepRowInd= datasample(allRowInd,n_subsample,'Replace',false);
        strainAllFileInd = strainAllFileInd(keepRowInd,:);
    end
    
    % Remove files that have NaN index in any window
    [trimmedFileInd,~] = dropNaNFiles(strainAllFileInd);
    
    %% Go through each feature to pre-allocate variable
    for featCtr = 1:numel(feats)
        
        % Get feature name
        feat = feats{featCtr};
        
        % Preallocate n-files by n-frames by n-lights variable to hold feature values for this feature
        featVals.(feat) = NaN(size(trimmedFileInd,1),n_frames,numel(lightConditions));
    end
    
    %% Go through each light condition
    for lightCtr = 1:numel(lightConditions)
        
        % Get each light condition
        light = lightConditions{lightCtr};
        lightColIdx = find(strcmp(allFileIndWindows,light));
        
        % Get file indices for this light condition
        fileInd = trimmedFileInd(:,lightColIdx);
        
        %% Go through each file
        for fileCtr = 1:numel(fileInd)
            fileIdx = fileInd(fileCtr);
            disp(['Extracting features from file ' num2str(fileCtr) ' out of ' num2str(numel(fileInd)) ' for ' light ' condition...'])
            
            % Get full file name and well number
            filename = [resultsDir,featureTable.filename{fileIdx},'/metadata_featuresN.hdf5'];
            well = featureTable.well_name(fileIdx);
            
            % Load time series data
            features = h5read(filename,'/timeseries_data');
            wellLogInd = strcmp(cellstr(features.well_name'),well);
            
            %% Go through each feature
            for featCtr = 1:numel(feats)
                
                % Get feature name
                feat = feats{featCtr};
                
                % Get time series values for this well
                x = features.timestamp(wellLogInd)+1; % add 1 for adjust for python indexing
                y = features.(feat)(wellLogInd);
                
                % Filter out NaN feature values
                x = x(~isnan(y));
                y = y(~isnan(y));
                
                % Accumulate feature sum for each frame
                featSumByFrame = accumarray(x,y,[n_frames,1]);
                
                % Accumulate count for each frame
                frameSumByFrame = accumarray(x,ones(size(x)),[n_frames,1]);
                
                % Get average feature value per frame
                featValByFrame = featSumByFrame./frameSumByFrame;
                
                % Write to preallocated variable
                featVals.(feat)(fileCtr,:,lightCtr) = featValByFrame;
            end
        end
    end
    
    %% Plot features
    featfigure = figure; hold on
    bluelightSubplotPos = [];
    
    % Go through each feature (rows of subplot)
    for featCtr = 1:numel(feats)
        feat = feats{featCtr};
        
        % Go through each light condition (column of subplot)
        for lightCtr = 1:numel(lightConditions)
            
            % Set subplot
            subplotPos = (featCtr-1)*numel(lightConditions)+lightCtr;
            ax(lightCtr) = subplot(numel(feats), numel(lightConditions),subplotPos); hold on
            
            % Get xy values
            y = squeeze(featVals.(feat)(:,:,lightCtr));
            x = 1:size(y,2);
            
            % Smooth y values over a time window
            y = smoothdata(y,'movmedian',smoothWindow*frameRate);
            
            % Plot
            shadedErrorBar(x,y,{@nanmean,@(x) nanstd(x)/sqrt(numel(x))},'lineProps','-k','transparent',1); hold on
            
            % Format
            title(lightConditions{lightCtr}); xlabel('frames'); ylabel(feat,'Interpreter','none')
            if strcmp(lightConditions{lightCtr},'bluelight')
                xlim([0 9000])
                bluelightSubplotPos = [bluelightSubplotPos,subplotPos];
            else
                xlim([0 7500])
            end
        end
        
        % Link y axes
        linkaxes(ax,'y')
        
        % Add bluelight pulse intervals
        if ~isempty(bluelightSubplotPos)
            % Set blue light interval frames for blue light patch x
            bluelightInterval = [60,70; 160,170; 260,270]*frameRate; % in frames
            x = horzcat(bluelightInterval,fliplr(bluelightInterval))';
            % Get ylim for linked axes to get blue light patch y
            thisFeatAxes = findall(ax(lightCtr),'type','axes');
            yL = get(thisFeatAxes,'YLim');
            y = [yL(1),yL(1),yL(2),yL(2)]';
            y = horzcat(y,y,y);
            % Add blue light patch objects to each of the bluelight subplots
            for subplotCtr = 1:numel(bluelightSubplotPos)
                subplot(numel(feats), numel(lightConditions),bluelightSubplotPos(subplotCtr)); hold on
                patch(x,y,'blue','EdgeColor','none','FaceAlpha',0.3)
            end
        end
    end
end