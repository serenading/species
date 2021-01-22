%% Script plots selected features across the three light condition videos.
% Script creates a 3x3 plot encompassing each strain in a subplot. 
% author: @serenading. Jan 2021.

clear
close all

addpath('../AggScreening/')
addpath('../AggScreening/auxiliary/')

% TODO: find matching indices too slow at the moment. see if there's away to
% avoid using a loop

%% Set parameters
extractStamp = '20201218_184325'; % 20201218_184325 for standard feature extraction, 20210112_105808 for filtered data
feats = {'speed_midbody_norm_90th','speed_midbody_norm_50th','speed_midbody_norm_10th'};
featCollectionName = 'midbody_speed_norm';
strains = {'N2','CB4856','MY23','QX1410','VX34','NIC58','JU1373'}; % 'N2','CB4856','MY23','QX1410','VX34','NIC58','JU1373';
lightInterval = [0,5*60; 5*60,11*60; 11*60,16*60]; % 5 min prestim, 6 min bluelight, 5 min poststim
n_subsample = NaN; % number of replicates per strain to sample. Set to NaN to use all files
n_nonFeatVar = 33; % the first n columns of the feature table that do not contain features. =33
resultsDir = '/Users/sding/OneDrive - Imperial College London/species/Results';

%% Get feature extraction windows
midpointAllwindows = sum(lightInterval,2)/2;
n_windows = numel(midpointAllwindows);

%% load features table
featureTable = readtable([resultsDir '/fullFeaturesTable_' extractStamp '.csv']);
% get light condition
light_condition = getLightcondition(featureTable);
% get species name
species_names = getSpeciesnames(featureTable);

%% Initialise figure
figure; hold on % this will be a nx3 plot, with column 1 being elegans, 2 being briggsae, 3 being tropicalis.
colors = {'m','c','k','r'};
subplotPos = [1;4;7]; % initialise position for subplot figures. each row is a different species
subplotActualPos = [];

%% Go through each strain
for strainCtr = 1:numel(strains)
    
    %% Get strain name
    strain = strains{strainCtr};
    disp(['Generating plot for ' strain  ' ...'])
    
    %% Find files for specified strain
    % filter for prestim files for the specified strain
    fileInd = find(strcmp(featureTable.strain_name,strain) & strcmp(light_condition,'prestim'));
    % subsample a few files for plotting
    if ~isnan(n_subsample)
        fileInd = datasample(fileInd,n_subsample,'Replace',false);
    end
    n_sample = numel(fileInd);
    
    %%  Preallocate matrix to hold file indices: n x 9 matrix, each column is one time window (9 windows total)
    allFileInd = NaN(n_sample,n_windows);
    allFileInd(:,1) = fileInd;
    
    %% Go through each file, find matching light condition files, and get timeseries feature
    % go through each file
    for sampleCtr = 1:n_sample
        % get prestim file index
        prestimfileIdx = allFileInd(sampleCtr,1);
        % get matching bluelight/poststim file indices
        [bluelightfileIdx,poststimfileIdx,well] = findMatchingFileInd(prestimfileIdx,featureTable);
        % recording bluelight/poststim file indices
        if ~isempty(bluelightfileIdx)
            allFileInd(sampleCtr,2) = bluelightfileIdx;
        end
        if ~isempty(poststimfileIdx)
            allFileInd(sampleCtr,3) = poststimfileIdx;
        end
    end
    
    %% Remove files that have NaN index in any window
    [trimmedFileInd,n_filesDropped] = dropNaNFiles(allFileInd);
    
    %% Extract features for the corresponding time windows
    % preallocate
    featVals = NaN(size(trimmedFileInd,1),numel(feats),n_windows);
    
    % go through each time window
    for windowCtr = 1:n_windows
        % extract feature values
        featVals(:,:,windowCtr) = featureTable{trimmedFileInd(:,windowCtr),feats};
    end
    
    % remove experiments with NaN feature values in any window
    [featVals,n_filesDropped,featValsCopy] = dropNaNVals(featVals);
    
    %% Plot features
    % Get species name for strain
    species_name = species_names(trimmedFileInd(1,1));
    
    % Find subplot location
    if strcmp(species_name,'elegans')
        subplotCol = 1;
    elseif strcmp(species_name,'briggsae')
        subplotCol = 2;
    elseif strcmp(species_name,'tropicalis')
        subplotCol = 3;
    else
        subplotCol = NaN;
        warning(['Invalid species name: ' speciesname])
    end
    if isscalar(subplotCol)
        % set subplot location
        ax(strainCtr) = subplot(3,3,subplotPos(subplotCol)); hold on
        % keep track of actual locations used for adding time window annotations later
        subplotActualPos = [subplotActualPos,subplotPos(subplotCol)];
        % update location counter
        subplotPos(subplotCol) = subplotPos(subplotCol)+1;
    end
        
    % Go through each feature and plot as shaded error bar
    for featCtr = 1:numel(feats)
        feat = feats{featCtr};
        vals = squeeze(featVals(:,featCtr,:));
        shadedErrorBar(midpointAllwindows,vals,{@nanmean,@(x) nanstd(x)/sqrt(numel(x))},...
            'lineprops',{[colors{featCtr} '-o'],'markerfacecolor',colors{featCtr}},'transparent',1); hold on
    end

    % Format
    title([strain sprintf('\n') 'n = ' num2str(size(featVals,1))])
    xlim([min(lightInterval(:)),max(lightInterval(:))])
    xlabel('time (s)')
    
end

%% Format plot
% link axis for all subplots
allAxes = findall(0,'type','axes');
linkaxes(allAxes,'xy')
% add legend
subplot(3,3,1)
legend(horzcat(feats, {'light condition'}) ,'Interpreter','none')
% add light interval annotations to each subplot
x = horzcat(lightInterval,fliplr(lightInterval))';
yL = get(allAxes,'YLim');
yL = yL{1};
y = [yL(1),yL(1),yL(2),yL(2)]';
y = horzcat(y,y,y);
for subplotCtr = 1:numel(subplotActualPos)
    subplot(3,3,subplotActualPos(subplotCtr))
    patch(x,y,[0 0 0],'EdgeColor','k','FaceAlpha',0)
    text(midpointAllwindows(1), yL(2)*0.95, 'prestim','HorizontalAlignment','center')
    text(midpointAllwindows(2), yL(2)*0.95, 'bluelight','HorizontalAlignment','center')
    text(midpointAllwindows(3), yL(2)*0.95, 'poststim','HorizontalAlignment','center')
end
% save figure
savefig([resultsDir '/threelight/' featCollectionName '_allstrains.fig'])