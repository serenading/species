function light_condition = getLightcondition(featureTable)

filenames = featureTable.filename;
light_condition = repmat({''},size(filenames));
light_condition(find(cellfun(@(x) contains(x,'prestim'),filenames))) = {'prestim'};
light_condition(find(cellfun(@(x) contains(x,'bluelight'),filenames))) = {'bluelight'};
light_condition(find(cellfun(@(x) contains(x,'poststim'),filenames))) = {'poststim'};

end