%% Function generates allFileInd variable to hold matching file indices across three light conditions for each strain,
% so they can be saved and loaded for use later, as this matching step is rather time consuming. 

function allFileInd = generateThreeLightMatchingIndices(extractStamp,strains)

%% INPUTS: 
% extractStamp: '20201218_184325' for standard feature extraction, '20210112_105808' for filtered data
% strains: cell array containing strain names as strings. % 'N2','CB4856','MY23','QX1410','VX34','NIC58','JU1373';

%% OUTPUT:
% allFileInd: struct indexable by strain name. allFileInd.strain is a
% n_file by 3 matrix where column 1 is prestim file index,  column 2 is
% bluelight file index, and column 3 is poststim file index from the
% featureTable.

%% FUNCTION

%% load features table
featureTable = readtable(['/Users/sding/OneDrive - Imperial College London/species/Results/fullFeaturesTable_' extractStamp '.csv']);
% get light condition
light_condition = getLightcondition(featureTable);

for strainCtr = 1:numel(strains)
    %% Get strain name
    strain = strains{strainCtr};
    disp(['Getting allFileInd for ' strain  ' ...'])
    
    %% Find files for specified strain
    % filter for prestim files for the specified strain
    fileInd = find(strcmp(featureTable.strain_name,strain) & strcmp(light_condition,'prestim'));
    n_sample = numel(fileInd);
    
    %%  Preallocate matrix to hold file indices: n_sample x n_window
    allFileInd.(strain) = NaN(n_sample,3);
    allFileInd.(strain)(:,1) = fileInd;
    
     %% Go through each file, find matching light condition files, and get timeseries feature
    % go through each file
    for sampleCtr = 1:n_sample
        % get prestim file index
        prestimfileIdx = allFileInd.(strain)(sampleCtr,1);
        % get matching bluelight/poststim file indices
        [bluelightfileIdx,poststimfileIdx,~] = findMatchingFileInd(prestimfileIdx,featureTable);
        % recording bluelight/poststim file indices
        if ~isempty(bluelightfileIdx)
            allFileInd.(strain)(sampleCtr,2) = bluelightfileIdx;
        end
        if ~isempty(poststimfileIdx)
            allFileInd.(strain)(sampleCtr,3) = poststimfileIdx;
        end
    end
end

save(['matchingFileInd/threelight_'  extractStamp '.mat'],'allFileInd')

end