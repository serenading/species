close all
clear

%% Script plots worm length time series data (useful for detectionn of defecation events every 50-60 seconds)
% Author: @serenading. Feb 2021.

extractStamp = '20201218_184325'; % 20201218_184325 for standard feature extraction, 20210112_105808 for filtered data
strain = 'JU1373';
n_subsample = 30;
feature = 'length';
frameRate = 25;
minTrajDuration = 150; % mininum trajectory in seconds
smoothWindow = 1; % smooth feature over the window, in seconds 

%% Load feature table
featureTable = readtable(['/Users/sding/OneDrive - Imperial College London/species/Results/fullFeaturesTable_' extractStamp '.csv'],'Delimiter',',','PreserveVariableNames',true);

%% Find file indices for the strain
fileInd = find(strcmp(featureTable.strain_name,strain));

% subsample as specified
if ~isnan(n_subsample)
    allRowInd = 1:size(fileInd,1);
    keepRowInd= datasample(allRowInd,n_subsample,'Replace',false);
    fileInd = fileInd(keepRowInd,:);
else
    n_subsample = numel(fileInd);
end

%% Load time series data

wormsChecked = 0;
wormsPlotted = 0;

% get file name
for fileCtr = 1:n_subsample
    fileIdx = fileInd(fileCtr);
    filename = ['/Volumes/Ashur DT2/species/Results/' featureTable.filename{fileIdx},'/metadata_featuresN.hdf5'];
    well = featureTable.well_name{fileIdx};
    
    % load timeseries data 
    features = h5read(filename,'/timeseries_data');
    
    % get logical index for the well
    well_name = cellstr(features.well_name');
    wellLogInd = strcmp(well_name,well);
    
    % go through each unique worm in this well
    uniqueWorms = unique(features.worm_index(wellLogInd));
    for wormCtr = 1:numel(uniqueWorms)
        
        % keep track of worms checked
        wormsChecked = wormsChecked +1;
        
        % find worm name and logical indices for this well
        worm = uniqueWorms(wormCtr);
        wormLogInd = wellLogInd & features.worm_index == worm;
        
        % only plot worms with trajectories over the specified minimum
        % trajectory duration
        if nnz(wormLogInd)>= minTrajDuration*frameRate 
            
            % keep count of worms plotted
            wormsPlotted = wormsPlotted+1;
            
            % get feature value
            y = features.(feature)(wormLogInd);
            x = features.timestamp(wormLogInd);
            
            % smooth y value
            y = smoothdata(y,'movmedian',smoothWindow*frameRate);
            
            % align x values at 0 and change x value into seconds
            x = (x-min(x))/frameRate;
            
            %% plot
            figure; hold on
            plot(x,y)
            
            % format
            xlabel('time (s)')
            ylabel(feature,'Interpreter','none')
            legend(num2str(worm))
            title([featureTable.filename{fileIdx} '_' well],'Interpreter','none')
        end  
    end
end

%% Display worm count
disp([num2str(wormsChecked) ' worm trajectories checked in ' num2str(fileCtr) ' wells. Out of these, ' num2str(wormsPlotted) ' trajectories met the min traj criteria and are plotted.'])