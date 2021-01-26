%% Function filters features table by specified strain and feature requirements before PCA.

function [featureTable, classLabels,filenames] = filterFeatureTable(featureTable,classVar,strains,featSetName)

%% Inputs:
% featureTable: file-by-feature featureTable
% classVar: String or cell array of strings containing the name of the variable to classify for. Each string must be a variable field of the featureTable.
% strains: Cell array containing strains to yse for analysis.
% featSetName: 'Tierpsy_256' or 'Tierpsy_all' 

%% Outputs: 
% featureTable: file-by-feature featureTable with non-feature metadata variables removed.
% classLabels: file x 1 class labels to feed into the classification task. If multiple multiple variables are specied for classVar then classLabels is a struct where each field contains the labels for that variable.

%% FUNCTION: 

%% process rows (observations) first

% retain strains as specified
strainLogInd = ismember(featureTable.strain_name,strains);
featureTable = featureTable(strainLogInd,:);
% get strain classification labels
if ischar(classVar)
    classLabels = featureTable.(classVar); 
else
    for varCtr = 1:numel(classVar)
        var = classVar{varCtr};
        classLabels.(var) = featureTable.(var);
    end
end

%% process columns (features) second

% load feature sets 
load('featureSet/features.mat','features');
featNames = features.(featSetName);
% trim table down to specified feature sets
featureTable = featureTable(:,featNames);

end