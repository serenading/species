function [bluelightfileIdx,poststimfileIdx,well] = findMatchingFileInd(prestimfileIdx,featureTable)

%% INPUTS:
% prestim file index
% featureTable

%% OUTPUTS:
% file indices for matching bluelight and poststim files

%% FUNCTION

% Get well name
well = featureTable.well_name{prestimfileIdx};

% Get imgstorenames
prestimimgstore = featureTable.filename{prestimfileIdx}; % get prestim imgstore name
prestimsplit = strsplit(prestimimgstore,'.'); % split prestim imgstorename
cameraname = prestimsplit{2};
prestim = prestimsplit{1}(1:end-6); % remove timestamp
bluelight = strrep(prestim,'prestim','bluelight');
poststim = strrep(prestim,'prestim','poststim');

% Get file indices
bluelightfileIdx = find(cellfun(@(x) contains(x,bluelight),featureTable.filename) &...
    cellfun(@(x) contains(x,cameraname),featureTable.filename) &...
    cellfun(@(x) strcmp(x,well),featureTable.well_name));

poststimfileIdx = find(cellfun(@(x) contains(x,poststim),featureTable.filename) &...
    cellfun(@(x) contains(x,cameraname),featureTable.filename) &...
    cellfun(@(x) strcmp(x,well),featureTable.well_name));

% check that only one index is returned for each light condition
assert(numel(bluelightfileIdx) == 1, ['There should only be one bluelight file index, but ' num2str(numel(bluelightfileIdx)) ' are found.'])
assert(numel(poststimfileIdx) == 1, ['There should only be one poststim file index, but ' num2str(numel(poststimfileIdx)) ' are found.'])

end