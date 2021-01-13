function featureTable = subsampleData(featureTable,n_subsample)

strains = unique(featureTable.strain_name);
rowLogInd = false(size(featureTable,1));

for strainCtr = 1:numel(strains)
    strain = strains(strainCtr);
    strainInd = find(strcmp(featureTable.strain_name,strain));
    strainInd2keep = datasample(strainInd,n_subsample,'Replace',false);
    rowLogInd(strainInd2keep) = true;
end

featureTable = featureTable(rowLogInd,:);