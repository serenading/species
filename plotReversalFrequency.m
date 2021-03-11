clear
close all

addpath('../AggScreening/')
addpath('../AggScreening/auxiliary/')

%% Set parameters
extractStamp = '20201218_184325'; % 20201218_184325 for standard feature extraction, 20210112_105808 for filtered data
lightConditions = {'prestim','bluelight','poststim'}; % {'prestim','bluelight','poststim'}
strains = {'N2','CB4856','MY23','QX1410','VX34','NIC58','JU1373'};
n_subsample = 20; % number of replicates per strain to sample.
frameRate = 25;
analysisResultsDir = '/Users/sding/OneDrive - Imperial College London/species/Results/';
trackingResultsDir = '/Volumes/Ashur DT2/species/Results/';
eventWindowSize = 20; % event window in seconds. reversal frequency are calculated for this window size using a sliding window.
smoothWindowSize = 3; % window in seconds to smooth calculated frequency data over.

%% Initialise figure
figure; hold on
bluelightPos = NaN(1,numel(strains)); % keep track of which subplots are bluelight figures

%% Check to see if calculated frequencies for this eventWindowSize does not already exist
if exist([analysisResultsDir 'threelight/reversal_frequency_windowSize' num2str(eventWindowSize) '.mat']) ==2
    load([analysisResultsDir 'threelight/reversal_frequency_windowSize' num2str(eventWindowSize) '.mat'],'featVals')
else
    %% Load features table
    featureTable = readtable([analysisResultsDir 'fullFeaturesTable_' extractStamp '.csv']);
    
    %% Load matching file indices across three light conditions
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
        
        % Preallocate n-files by n-frames by n-lights variable to hold feature values for this feature
        n_frames = 6*60*frameRate + 10; % maximum duration per video is 6 min. Plus 10 extra frames because the recordings don't always finish at the precise end frame
        featVals.(strain) = NaN(size(trimmedFileInd,1),n_frames,numel(lightConditions));
        
        %% Go through each light condition
        for lightCtr = 1:numel(lightConditions)
            
            % Get each light condition
            light = lightConditions{lightCtr};
            lightColIdx = find(strcmp(allFileIndWindows,light));
            
            % Get file indices for this light condition
            fileInd = trimmedFileInd(:,lightColIdx);
            n_files = numel(fileInd);
            
            %% Go through each file
            for fileCtr = 1:n_files
                fileIdx = fileInd(fileCtr);
                disp(['Extracting features from file ' num2str(fileCtr) ' out of ' num2str(n_files) ' for ' light ' condition...'])
                
                % Get full file name and well number
                filename = [trackingResultsDir,featureTable.filename{fileIdx},'/metadata_featuresN.hdf5'];
                well = featureTable.well_name(fileIdx);
                
                % Load time series data
                features = h5read(filename,'/timeseries_data');
                wellLogInd = strcmp(cellstr(features.well_name'),well);
                eventFreq = calculateReversalFrequency(features,wellLogInd,eventWindowSize,frameRate);
                
                % Write to preallocated variable
                featVals.(strain)(fileCtr,1:numel(eventFreq),lightCtr) = eventFreq;
            end
        end
    end
    save([analysisResultsDir 'threelight/reversal_frequency_windowSize' num2str(eventWindowSize) '.mat'],'featVals');
end

%% Plot
for strainCtr = 1:numel(strains)
    strain = strains{strainCtr};
    
    % Go through each light condition (column of subplot)
    for lightCtr = 1:numel(lightConditions)
        
        % Set subplot
        subplotPos = (strainCtr-1)*numel(lightConditions)+lightCtr;
        subplot(numel(strains), numel(lightConditions),subplotPos); hold on
        
        % Get xy values
        y = squeeze(featVals.(strain)(:,:,lightCtr));
        x = 1:size(y,2);
        
        % Smooth y values over a time window
        if ~isnan(smoothWindowSize)
            y = smoothdata(y,'movmedian',smoothWindowSize*frameRate);
        end
        
        % Plot
        shadedErrorBar(x,y,{@nanmean,@(x) nanstd(x)/sqrt(numel(x))},'lineProps','-k','transparent',1); hold on
        
        % Format
        n_files = size(featVals.(strain),1);
        title([strain ', ' lightConditions{lightCtr} newline ' ( n = ' num2str(n_files) ' wells )']); xlabel('frames'); ylabel('reversal frequency','Interpreter','none')
        if strcmp(lightConditions{lightCtr},'bluelight')
            xlim([0 9000])
            bluelightPos(strainCtr) = subplotPos;
        else
            xlim([0 7500])
        end
    end
end

%% Format plot

% Link y axes
ax = findall(0,'type','axes');
linkaxes(ax,'y')

% Add bluelight pulse intervals
if ~isempty(bluelightPos)
    % Set blue light interval frames for blue light patch x
    bluelightInterval = [60,70; 160,170; 260,270]*frameRate; % in frames
    x = horzcat(bluelightInterval,fliplr(bluelightInterval))';
    % Get ylim for linked axes to get blue light patch y
    yL = get(ax,'YLim');
    yL = yL{1};
    y = [yL(1),yL(1),yL(2),yL(2)]';
    y = horzcat(y,y,y);
    % Add blue light patch objects to each of the bluelight subplots
    for subplotCtr = 1:numel(bluelightPos)
        subplot(numel(strains), numel(lightConditions),bluelightPos(subplotCtr)); hold on
        patch(x,y,'blue','EdgeColor','none','FaceAlpha',0.3)
    end
end

%% Save
savefig([analysisResultsDir 'threelight/reversal_frequency_' num2str(eventWindowSize) 's_allstrains_' extractStamp '.fig']);