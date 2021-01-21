%% Function removes incomplete experiments from the file indices matrix.
% An experiment is considered incomplete where features are empty for any
% of the sample windows (i.e. due to worms having crawled off leading to no
% extracted Tierpsy features). 

function [trimmedFileInd,n_filesDropped] = dropNaNFiles(allFileInd)

%% INPUT:
% allFileInd: n_sample x n_windows array holding file indices

%% OUTPUTS:
% trimmedFileInd: trimmed file indices array with incomplete experiments removed. 
% n_filesDropped: scalar showing how many experiments are dropped 

%% FUNCTION:

% find row index that contains NaN index
[rows2drop,~] = find(isnan(allFileInd));
% generate row logical index for rows to keep
rowLogInd = true(size(allFileInd,1),1);
rowLogInd(unique(rows2drop)) = false;
% trim down file indices by dropping rows with NaN file index
trimmedFileInd = allFileInd(rowLogInd,:,:);
% keep count of the number of files dropped
n_filesDropped = numel(unique(rows2drop));
assert(n_filesDropped == size(allFileInd,1) - size(trimmedFileInd,1));

end