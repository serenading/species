clear
close all

% load strains
strains = {'N2','CB4856','MY23','QX1410','VX34','NIC58','JU1373'}; 
% set path to where the extracted features summary files are located. 
resultsDir = '/Users/sding/OneDrive - Imperial College London/species/Results';

% get three lights matching file indices
extractStamp = '20201218_184325';
[allFileInd,allFileIndWindows] = getMatchingIndicesThreeLights(extractStamp,strains,resultsDir);
save(['matchingFileInd/threelight_'  extractStamp '.mat'],'allFileInd','allFileIndWindows')

% get three lights matching file indices for filtered data
extractStamp = '20210112_105808';
[allFileInd,allFileIndWindows] = getMatchingIndicesThreeLights(extractStamp,strains,resultsDir);
save(['matchingFileInd/threelight_'  extractStamp '.mat'],'allFileInd','allFileIndWindows')

% get bluelight matching file indices
extractStamp = '20210119_073010';
windownames = [0:8];
[allFileInd,allFileIndWindows] = getMatchingIndicesBlueLight(extractStamp,strains,windownames,resultsDir);
save(['matchingFileInd/bluelight_'  extractStamp '.mat'],'allFileInd','allFileIndWindows')