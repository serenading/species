%% Function removes incomplete feature values from the featVal matrix.
% An experiment is considered incomplete where features are empty for any
% of the sample windows (i.e. due to worms having crawled off leading to no
% extracted Tierpsy features).

function [featVals,n_filesDropped,featValsCopy] = dropNaNVals(featVals)

%% INPUT:
% featVals: n_sample x n_feature x n_windows array holding feature values

%% OUTPUTS:
% featVals: trimmed feature values array with incomplete experiments removed.
% n_filesDropped: scalar showing how many experiments are dropped.
% featValsCopy: a copy of the feature values matrix prior to processing.

%% FUNCTION:

% make copy of unprocessed featVals
featValsCopy= featVals;
% subsample the first feature out of n features
featValsSmall = squeeze(featVals(:,1,:));
% find row index that contains NaN index
[rows2drop,~] = find(isnan(featValsSmall));
% generate row logical index for rows to keep
rowLogInd = true(size(featVals,1),1);
rowLogInd(unique(rows2drop)) = false;
% trim down file indices by dropping rows with NaN file index
featVals = featVals(rowLogInd,:,:);
% keep count of the number of files dropped
n_filesDropped = numel(unique(rows2drop));
assert(n_filesDropped == size(featValsCopy,1) - size(featVals,1));
% display message
disp([num2str(n_filesDropped) ' files out of ' num2str(numel(rowLogInd)) ' dropped due to NaN feature values in at least one of the time windows.'])

end