%% Function generates a featureTable from Hydra data by combining 
% 1. Compiled Hydra metadata file; 2. Filenames file generated by Tierpsy; 3. Features summary file generated by Tierpsy.
% author: serenading. Jan 2021

function featureTable = getFeatureTable(extractStamp,projectDir,annotatedWells)

%% INPUT:
% extractStamp: yyyymmdd_hhmmss timestamp attached to the tierpsy
% summariser output files.
% projectDir: path to the directory where the tierpsy summariser outputs
% and metadata are saved.
% annotatedWells: are wells annotated? true or false.

%% OUTPUT:
% featureTable: combined table, filtered using well annotator and Tierpsy outputs

%% Import features data and combine with metadata

% load features matrix, correspondong filenames, and metadata
tierpsyFeatureTable = readtable([projectDir '/features_summary_tierpsy_plate_' extractStamp '.csv'],'Delimiter',',');%,'preserveVariableNames',true);
tierpsyFileTable = readtable([projectDir '/filenames_summary_tierpsy_plate_' extractStamp '.csv'],'Delimiter',',','CommentStyle','#');%,'preserveVariableNames',true);
if annotatedWells
    metadataTable = readtable([projectDir '/wells_updated_metadata.csv'],'Delimiter',',','preserveVariableNames',true);
else
    metadataTable = readtable([projectDir '/metadata.csv'],'Delimiter',',','preserveVariableNames',true);
end

% rename metadata column heads to match Tierpsy output
metadataTable.Properties.VariableNames{'imgstore_name'} = 'filename';
metadataTable.Properties.VariableNames{'worm_strain'} = 'strain_name';

% reformat the filenames to match the imgstore name
filenames = tierpsyFileTable.filename;
splitnames = cellfun(@(x) strsplit(x,'/'), filenames,'UniformOutput', false);
reconstructnames = cellfun(@(x) strcat(x{6},'/',x{7}), splitnames,'UniformOutput',false);
tierpsyFileTable.filename = reconstructnames;

%% Join tables

% join the Tierpsy tables to match filenames with file_id. Required in case
% features were not extracted for any files.
combinedTierpsyTable = outerjoin(tierpsyFileTable,tierpsyFeatureTable,'MergeKeys',true);

% finally, join tables to get strain names for each set of features
featureTable = outerjoin(metadataTable,combinedTierpsyTable,'MergeKeys',true);

%% Filter table for good data 
% get row logical index for valid files
% (note: this does not give a balanced set of files where good data for all three light conditions are present)
if annotatedWells
    rowLogInd_wellAnnotator = strcmp(featureTable.is_bad_well,'False') & featureTable.well_label==1;
else
    rowLogInd_wellAnnotator = true(size(featureTable.is_good));
end
rowLogInd_Tierpsy = strcmp(featureTable.is_good,'True') & featureTable.is_good_well==1;
rowLogInd = rowLogInd_wellAnnotator & rowLogInd_Tierpsy; % ~75% data is good

% trim featureTable down to those with valid files
featureTable = featureTable(rowLogInd,:);

%% Display QC messages
% 
n_totalFiles = numel(rowLogInd); 
% 
n_bad_wellAnnotator = n_totalFiles - nnz(rowLogInd_wellAnnotator); 
disp([num2str(n_bad_wellAnnotator) ' out of ' num2str(n_totalFiles) ' files (across all light conditions) are dropped because they are manually marked as bad wells.'])
%
n_bad_Tierpsy = n_totalFiles - nnz(rowLogInd_Tierpsy); 
disp([num2str(n_bad_Tierpsy) ' out of ' num2str(n_totalFiles) ' files (across all light conditions) are dropped because Tierpsy determined they are bad.'])
%
n_goodFiles = nnz(rowLogInd);
disp(['After filtering, a total of ' num2str(n_goodFiles) ' out of ' num2str(n_totalFiles) ' files are retained for the featureTable.'])
%
n_prestimFiles = nnz(contains(featureTable.filename,'prestim'));
n_bluelightFiles = nnz(contains(featureTable.filename,'bluelight'));
n_poststimFiles = nnz(contains(featureTable.filename,'poststim'));
disp(['Out of these  ' num2str(n_goodFiles) ' good files, ' num2str(n_prestimFiles) ' are prestim, ' num2str(n_bluelightFiles) ' are bluelight, and ' num2str(n_poststimFiles) ' are poststim.'])

end