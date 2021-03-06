%% Script plots selected features across the blue light condition videos.
% author: @serenading. Jan 2021.

clear
close all

addpath('../AggScreening/')
addpath('../AggScreening/auxiliary/')

%% Set parameters

extractStamp = '20210119_073010'; % '20210119_073010' feature summaries have multiple bluelight windows.
featSetName = 'motion_mode_fraction';
strains = {'N2','CB4856','MY23','QX1410','VX34','NIC58','JU1373'}; % 'N2','CB4856','MY23','QX1410','VX34','NIC58','JU1373';
n_subsample = NaN; % number of replicates per strain to sample. Set to NaN to use all files
bluelightInterval = [60,70; 160,170; 260,270];
windowDuration = 10; % 10 (default based on Ziwei's time window analysis). Time window in seconds to retain for feature analysis.
windownames = [0:8]; % windows are named 0 through 8
resultsDir = '/Users/sding/OneDrive - Imperial College London/species/Results';

%% Get feature extraction windows
[midpointAllwindows, prelight, bluelight, postlight] = getBluelightFeatWindows(bluelightInterval,windowDuration);
n_windows = numel(windownames);
assert(numel(midpointAllwindows) == n_windows)
% Enter the resultant windows (in seconds) into Tierpsy feature summariser for feature
% window extraction (takes about a minute per file on my macpro)
% [50:60,65:75,75:85,150:160,165:175,175:185,250:260,265:275,275:285]

%% Load saved matching file indices and feature sets
% load matching file indices across three light conditions
load(['matchingFileInd/bluelight_' extractStamp '.mat'],'allFileInd');
% get feature sets
load('featureSet/features.mat','features');
feats = features.(featSetName);

%% Initialise figure
figure; hold on % this will be a 3xn plot, with row 1 being elegans, 2 being briggsae, 3 being tropicalis.
colors = {'m','c','k','r'};
subplotPos = [1;4;7]; % initialise position for subplot figures. each row is a different species
subplotActualPos = [];

%% Go through each strain
for strainCtr = 1:numel(strains)
    
    % Get matching file indices for this strain
    strain = strains{strainCtr};
    disp(['Getting file indices for ' strain  ' ...'])
    strainAllFileInd = allFileInd.(strain);

    % If specified, subsample a few files for plotting
    if ~isnan(n_subsample)
        allRowInd = 1:size(strainAllFileInd,1);
        keepRowInd= datasample(allRowInd,n_subsample,'Replace',false);
        strainAllFileInd = strainAllFileInd(keepRowInd,:);
    end

    % Remove files that have NaN index in any window
    [trimmedFileInd.(strain),~] = dropNaNFiles(strainAllFileInd);
    
    % Preallocate variable for holding feature values
    featVals.(strain) = NaN(size(trimmedFileInd.(strain),1),numel(feats),n_windows);
    imagingDates.(strain) = NaN(size(trimmedFileInd.(strain),1),n_windows);
end
    
%% Go through each time window for feature extraction
% (loading the featureTables is the slow step, so loop through strains
% inside this section to avoid having to load the tables for each strain).

for windowCtr = 1:n_windows
    % Get window name
    window = windownames(windowCtr);
    disp(['Extracting features from time window ' num2str(windowCtr) ' out of ' num2str(n_windows) '...'])
    
    % Read featureTable for this window
    featureTable = readtable([resultsDir '/fullFeaturesTable_' extractStamp '_window_' num2str(window) '.csv']);
    
    % Go through each strain to extract features for this time window
    for strainCtr = 1:numel(strains)
        strain = strains{strainCtr};
        featVals.(strain)(:,:,windowCtr) = featureTable{trimmedFileInd.(strain)(:,windowCtr),feats};
        imagingDates.(strain)(:,windowCtr) = featureTable.date_yyyymmdd(trimmedFileInd.(strain)(:,windowCtr));
    end
end

%% Go through each strain

for strainCtr = 1:numel(strains)
    strain = strains{strainCtr};
    disp(['Generating plot for ' strain  ' ...'])

    %% Remove experiments with NaN feature values in any window
    [featVals.(strain),rowLogInd,~] = dropNaNVals(featVals.(strain));
    imagingDates.(strain) = imagingDates.(strain)(rowLogInd,:);
    
     %% Plot features
     
    % Get species name for strain using the currently loaded featureTable
    species_names = getSpeciesnames(featureTable);
    species_name = species_names{trimmedFileInd.(strain)(1,windowCtr)};
    
    % Find subplot location
    if strcmp(species_name,'elegans')
        subplotRow = 1;
    elseif strcmp(species_name,'briggsae')
        subplotRow = 2;
    elseif strcmp(species_name,'tropicalis')
        subplotRow = 3;
    else
        warning(['Invalid species name: ' species_name])
    end
    if isscalar(subplotRow)
        % set subplot location
        subplot(3,3,subplotPos(subplotRow)); hold on
        % keep track of actual locations used for adding time window annotations later
        subplotActualPos = [subplotActualPos,subplotPos(subplotRow)];
        % update location counter
        subplotPos(subplotRow) = subplotPos(subplotRow)+1;
    end
        
    % Go through each feature and plot as shaded error bar
    for featCtr = 1:numel(feats)
        feat = feats{featCtr};
        vals = squeeze(featVals.(strain)(:,featCtr,:));
        shadedErrorBar(midpointAllwindows,vals,{@nanmean,@(x) nanstd(x)/sqrt(numel(x))},...
            'lineprops',{[colors{featCtr} '-o'],'markerfacecolor',colors{featCtr}},'transparent',1); hold on
    end

    % Format
    title([strain newline 'n = ' num2str(size(featVals.(strain),1))])
    xlim([30 300])
    xlabel('time (s)')
    
end

%% Format plot

% link axis for all subplots
allAxes = findall(0,'type','axes');
linkaxes(allAxes,'xy')

% add light interval annotations to each subplot
x = horzcat(bluelightInterval,fliplr(bluelightInterval))';
yL = get(allAxes,'YLim');
yL = yL{1};
y = [yL(1),yL(1),yL(2),yL(2)]';
y = horzcat(y,y,y);
for subplotCtr = 1:numel(subplotActualPos)
    subplot(3,3,subplotActualPos(subplotCtr)); hold on
    patch(x,y,'blue','EdgeColor','none','FaceAlpha',0.3)
end

% add legend
subplot(3,3,1); hold on
legend(horzcat(feats,{'bluelight pulse'}) ,'Interpreter','none') 

% save figure
%savefig([resultsDir '/bluelight/' featSetName '_allstrains_' extractStamp '.fig'])