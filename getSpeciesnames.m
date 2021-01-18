function species_names = getSpeciesnames(featureTable)

strainnames = featureTable.strain_name;
species_names = repmat({''},size(strainnames));
elegansLogInd = cellfun(@(x) strcmp(x,'N2'),strainnames) | cellfun(@(x) strcmp(x,'CB4856'),strainnames) | cellfun(@(x) strcmp(x,'MY23'),strainnames);
briggsaeLogInd = cellfun(@(x) strcmp(x,'QX1410'),strainnames) | cellfun(@(x) strcmp(x,'VX34'),strainnames);
tropicalisLogInd = cellfun(@(x) strcmp(x,'NIC58'),strainnames) | cellfun(@(x) strcmp(x,'JU1373'),strainnames);
species_names(elegansLogInd) = {'elegans'};
species_names(briggsaeLogInd) = {'briggsae'};
species_names(tropicalisLogInd) = {'tropicalis'};

end