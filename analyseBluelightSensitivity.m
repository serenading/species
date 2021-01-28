%% Script analyses bluelight sensitivity by looking at significant feature changes between prestim and bluelight videos.
% author: @serenading. Jan 2020

clear
close all

%% Set parameters
extractStamp = '20210112_105808'; % 20201218_184325 for standard feature extraction, 20210112_105808 for filtered data
featSetName = 'Tierpsy_all'; % 'Tierpsy_256','Tierpsy_all','motion_mode','midbody_speed','angular_velocity'
strains = {'N2','CB4856','MY23','QX1410','VX34','NIC58','JU1373'}; % 'N2','CB4856','MY23','QX1410','VX34','NIC58','JU1373';
n_subsample = NaN; % number of replicates per strain to sample. Set to NaN to use all files
resultsDir = '/Users/sding/OneDrive - Imperial College London/species/Results';
windows = [1,2]; % 1 for prestim, 2 for bluelight
bonCorr = true; % apply bonferroni correction for multiple comparisons

%% load data
% load features table
featureTable = readtable([resultsDir '/fullFeaturesTable_' extractStamp '.csv']);
% load matching file indices across three light conditions
load(['matchingFileInd/threelight_' extractStamp '.mat'],'allFileInd');
% get feature sets
load('featureSet/features.mat','features');
feats = features.(featSetName);

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
    if ~isnan(n_subsample)
        allRowInd = 1:size(strainAllFileInd,1);
        keepRowInd= datasample(allRowInd,n_subsample,'Replace',false);
        strainAllFileInd = strainAllFileInd(keepRowInd,:);
    end
    
    %% Remove files that have NaN index in any window
    [trimmedFileInd,~] = dropNaNFiles(strainAllFileInd);
    
    %% Extract features for the corresponding time windows
    % preallocate
    n_windows = numel(windows);
    featVals = NaN(size(trimmedFileInd,1),numel(feats),n_windows);
    
    % go through each time window
    for windowCtr = 1:n_windows
        % extract feature values
        window = windows(windowCtr);
        featVals(:,:,windowCtr) = featureTable{trimmedFileInd(:,window),feats};
    end
    
    % remove experiments with NaN feature values in any window
    [featVals,~,~] = dropNaNVals(featVals);
    
    %% Perform statistical analysis: t-test correcting for multiple comparisons
    % set alpha for t-test p value
    alpha = 0.05;
    if bonCorr
        alpha = alpha/numel(feats); % correct for multiple comparison
    end
    % conduct t-test for each feature
    for featCtr = 1:numel(feats)
        prestimFeat = featVals(:,featCtr,1);
        bluelightFeat = featVals(:,featCtr,2);
        [h,p] = ttest2(prestimFeat,bluelightFeat,'Alpha',alpha);
        % record sensitive feature
        if h==1
            sensitiveFeat_n.(strain) = sensitiveFeat_n.(strain)+1;
            sensitiveFeatInd.(strain) = [sensitiveFeatInd.(strain), featCtr];
        end
    end
end

%% Display results
for strainCtr = 1:numel(strains)
    strain = strains{strainCtr};
    sigPercent = round(sensitiveFeat_n.(strain)/numel(feats)*100);
    disp([strain ' has ' num2str(sensitiveFeat_n.(strain)) ' features out of ' num2str(numel(feats))  ' (' num2str(sigPercent) '%) that show significant changes upon bluelight stimulation.'])
end

%% Save results
save([resultsDir '/bluelightSensitivity/' featSetName '_' extractStamp '.mat'],'sensitiveFeat_n','sensitiveFeatInd','feats')