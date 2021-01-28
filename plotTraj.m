%% Script visualises full frame trajectories from an entire camera view (4x4 wells).
% author: @serenading. Jan 2021.

clear
close all

%% Set parameters
strain = 'N2'; % 'N2','CB4856','MY23','QX1410','VX34','NIC58','JU1373';
light = 'poststim'; % 'prestim', 'bluelight','poststim'. Specify the light condition to plot.
n_subsample = 2; % Specify how full camera FOV per strain to show.
resultsDir = '/Volumes/Ashur DT2/species/Results/';
plotFirstFrame = true;
frameRate = 25;

%% Find files for specified strain
% load features table
featureTable = readtable('/Users/sding/OneDrive - Imperial College London/species/Results/fullFeaturesTable_20201218_184325.csv');
% get light condition
light_condition = getLightcondition(featureTable);
% filter for prestim files for the specified strain
strainInd = find(strcmp(featureTable.strain_name,strain) & strcmp(light_condition,'prestim'));
% subsample a few files for plotting
strainInd = datasample(strainInd,n_subsample,'Replace',false);

%% Go through each file to 
for fileCtr = 1:n_subsample
    % get prestim file index
    prestimfileIdx = strainInd(fileCtr);
    % get matching bluelight/poststim file indices
    [bluelightfileIdx,poststimfileIdx,~] = findMatchingFileInd(prestimfileIdx,featureTable);
    % get full skeletons filename
    if strcmp(light,'prestim')
        imgstorename = featureTable.filename{prestimfileIdx};
    elseif strcmp(light,'bluelight')
        imgstorename = featureTable.filename{bluelightfileIdx};
    elseif strcmp(light,'poststim')
        imgstorename = featureTable.filename{poststimfileIdx};
    end
    filename = [resultsDir,imgstorename,'/metadata_skeletons.hdf5'];
    %% Plot trajectories from the entire camera and annotate wells with wellname and strain name
    trajFig = figure; hold on
    % load trajectories data from skeletons file
    trajData = h5read(filename,'/trajectories_data');
    % load first full frame
    if plotFirstFrame
        fullImages = h5read(strrep(strrep(filename,'Results','MaskedVideos'),'_skeletons.hdf5','.hdf5'),'/full_data');
        firstFrame = fullImages(:,:,1)'; % transpose image to adjust for python/MATLAB differences
        firstFrame = cat(3,firstFrame,firstFrame,firstFrame); % convert to RGB
        imshow(firstFrame); hold on
    end
    % plot trajectories
    pointsize = 5;
    timeCmap = double(trajData.frame_number)/frameRate/60; % color according to time progression
    scatter(trajData.coord_x,trajData.coord_y,pointsize,timeCmap,'filled');
    % add well boundaries
    fov_wells = h5read(strrep(filename,'_skeletons','_featuresN'),'/fov_wells');
    wellnames = cellstr(fov_wells.well_name');
    for wellCtr = 1:numel(fov_wells.x_max)
        wellname = wellnames{wellCtr};
        % get five-point well coordinates
        well_coords = [fov_wells.x_min(wellCtr),fov_wells.y_min(wellCtr); fov_wells.x_min(wellCtr),fov_wells.y_max(wellCtr);...
            fov_wells.x_max(wellCtr),fov_wells.y_max(wellCtr); fov_wells.x_max(wellCtr),fov_wells.y_min(wellCtr);
            fov_wells.x_min(wellCtr),fov_wells.y_min(wellCtr)];
        plot(well_coords(:,1),well_coords(:,2),'r-')
        % add well name
        text(double(fov_wells.x_min(wellCtr)+50),double(fov_wells.y_max(wellCtr))-50,wellname,'Color','red')
        % add strain name. Note: if the well is marked as bad then the strainname won't show
        wellIdx = find(cellfun(@(x) strcmp(x,imgstorename), featureTable.filename) &...
            cellfun(@(x) strcmp(x,wellname), featureTable.well_name));  
        assert(numel(wellIdx) <= 1, ['There should only be 1 well for this imgstorename and this well name, but ' num2str(numel(wellIdx)) ' are found.'])
        strainname = featureTable.strain_name(wellIdx);
        text(double(fov_wells.x_min(wellCtr)+50),double(fov_wells.y_max(wellCtr))-100,strainname,'Color','red')
    end
    axis equal
    colorbar
    title(imgstorename,'Interpreter','none')
end