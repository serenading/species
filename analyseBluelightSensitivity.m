%% Script analyses bluelight sensitivity by looking at significant feature changes between prestim and bluelight videos.
% Options to specify window size as big to use the 5 min prestim vs. 6 min bluelight videos, 
% or specify window size as small to use the 10 sec prelight vs. 10 sec bluelight pulse (first pulse only)
% author: @serenading. Jan 2020

clear
close all

%% Set parameters
windowSize = 'small'; % 'big' to use prestim/bluelight level feature extraction windows (5-6 minutes), 'small' to use 10 second windows surrounding lught pulses the bluelight movies
featSetName = 'Tierpsy_all'; % 'Tierpsy_256','Tierpsy_all','motion_mode','midbody_speed','angular_velocity'
strains = {'N2','CB4856','MY23','QX1410','VX34','NIC58','JU1373'}; % 'N2','CB4856','MY23','QX1410','VX34','NIC58','JU1373';
n_subsample = NaN; % number of replicates per strain to sample. Set to NaN to use all files
resultsDir = '/Users/sding/OneDrive - Imperial College London/species/Results';
bonCorr = true; % apply bonferroni correction for multiple comparisons

%% load data
if strcmp(windowSize,'big')
    windows = [1,2]; % 1 for prestim (5 min), 2 for bluelight (6 min)
    extractStamp = '20210112_105808'; % 20201218_184325 for standard feature extraction, 20210112_105808 for filtered data
    % load features table
    featureTable = readtable([resultsDir '/fullFeaturesTable_' extractStamp '.csv']);
    % load matching file indices across three light conditions
    load(['matchingFileInd/threelight_' extractStamp '.mat'],'allFileInd');
elseif strcmp(windowSize,'small')
    windows = [0,1]; % 0 for first prelight window (10 s), 1 for first bluelight window (10 s)
    extractStamp = '20210119_073010'; % '20210119_073010' feature summaries have multiple bluelight windows.
    % load matching file indices across blue light pulses
    load(['matchingFileInd/bluelight_' extractStamp '.mat'],'allFileInd');
else
    warning('Please specify windowSize as big or small.')
end

% get feature sets
load('featureSet/features.mat','features');
feats = features.(featSetName);

% set alpha for t-test p value
alpha = 0.05;
if bonCorr
    alpha = alpha/numel(feats); % correct for multiple comparison
end

%% Go through each strain
for strainCtr = 1:numel(strains)
    
    %% Get matching file indices for this strain
    strain = strains{strainCtr};
    disp(strain)
    strainAllFileInd = allFileInd.(strain);
      
    %% Initialise variable for containing results
    sensitiveFeat_n.(strain) = 0;
    sensitiveFeatInd.(strain) = [];
    
    %% If specified, subsample a few files for plotting
    if isscalar(n_subsample)
        allRowInd = 1:size(strainAllFileInd,1);
        keepRowInd= datasample(allRowInd,n_subsample,'Replace',false);
        strainAllFileInd = strainAllFileInd(keepRowInd,:);
    end
    
    %% Remove files that have NaN index in any window
    [trimmedFileInd.(strain),~] = dropNaNFiles(strainAllFileInd);
    
    %% Extract features for the corresponding time windows
    % preallocate
    n_windows = numel(windows);
    featVals.(strain) = NaN(size(trimmedFileInd.(strain),1),numel(feats),n_windows);
end
    
%% Go through each time window for feature extraction
% (loading the featureTables is the slow step, so loop through strains
% inside this section to avoid having to load the tables for each strain).

disp(newline)
for windowCtr = 1:n_windows
    % Get window name
    window = windows(windowCtr);
    disp(['Extracting features from time window ' num2str(windowCtr) ' out of ' num2str(n_windows) '...'])
    
    if strcmp(windowSize,'small')
        % Read featureTable for this window
        featureTable = readtable([resultsDir '/fullFeaturesTable_' extractStamp '_window_' num2str(window) '.csv']);
    end
    
    % Go through each strain to extract features for this time window
    for strainCtr = 1:numel(strains)
        strain = strains{strainCtr};
        featVals.(strain)(:,:,windowCtr) = featureTable{trimmedFileInd.(strain)(:,windowCtr),feats};
    end
end

%% Go through each strain

disp(newline)
for strainCtr = 1:numel(strains)
    strain = strains{strainCtr};
    disp(['Analysing blue light sensitivity for ' strain  ' ...'])
    
    %% Remove experiments with NaN feature values in any window
    [featVals.(strain),~,~] = dropNaNVals(featVals.(strain));
    
    %% Perform statistical analysis: t-test correcting for multiple comparisons
    % conduct t-test for each feature
    for featCtr = 1:numel(feats)
        prestimFeat = featVals.(strain)(:,featCtr,1);
        bluelightFeat = featVals.(strain)(:,featCtr,2);
        [h,p] = ttest2(prestimFeat,bluelightFeat,'Alpha',alpha);
        % record sensitive feature
        if h==1
            sensitiveFeat_n.(strain) = sensitiveFeat_n.(strain)+1;
            sensitiveFeatInd.(strain) = [sensitiveFeatInd.(strain), featCtr];
        end
    end
end

%% Display results
disp(newline)
for strainCtr = 1:numel(strains)
    strain = strains{strainCtr};
    sigPercent = round(sensitiveFeat_n.(strain)/numel(feats)*100);
    disp([strain ' has ' num2str(sensitiveFeat_n.(strain)) ' features out of ' num2str(numel(feats))  ' (' num2str(sigPercent) '%) that show significant changes upon bluelight stimulation.'])
end

%% Save results
save([resultsDir '/bluelightSensitivity/' featSetName '_windowSize_' windowSize '_' extractStamp '.mat'],'sensitiveFeat_n','sensitiveFeatInd','feats')