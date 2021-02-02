%% Script performs PCA analysis, with options to specify which features, which strains, and which light conditions to use. 
% Script generates individual plots for each strain and combined plots for
% all strains, light conditions, experimental dates, and plate batches. 
% author: @serenading Jan 2021

clear
close all

%% Set analysis parameters
extractStamp = '20201218_184325'; % 20201218_184325 for standard feature extraction, 20210112_105808 for filtered data
lightConditions = {'prestim','bluelight','poststim'}; % 'prestim','bluelight','poststim'
classVar = {'strain_name','light_condition','date_yyyymmdd','date_plates_poured_yyyymmdd'}; % variable names to retain for plotting

% Set filering parameters
strains = {'N2','CB4856','MY23','QX1410','VX34','NIC58','JU1373'}; 
featSetName = 'Tierpsy_256'; % 'Tierpsy_256' or 'Tierpsy_all'

n_subsample = NaN; % number of replicates per strain to include. Set to NaN to include all samples
n_skeletons_range = [50 22500]; % n_skeleton range to use for retaining the well. 25fps x 60s/min x 5 min x 3 worms = 22500 skeletons.

removeOutlier = true; % option to remove outlier coefficients from each PC
n_std2include = 3; % define outliers as being outside mean+/- n standard deviations. Use 3 for conservative outlier removal retaining 99.9% data. Use 10 to remove only extreme outliers.

%% Load featureTable
% load features table
featureTable = readtable(['/Users/sding/OneDrive - Imperial College London/species/Results/fullFeaturesTable_' extractStamp '.csv']);
% get light condition
light_condition = getLightcondition(featureTable);
featureTable.light_condition = light_condition;
% get feature sets
load('featureSet/features.mat','features');
feats = features.(featSetName);

%% Filter featureTable
% Filter for selected light condition
lightLogInd = ismember(light_condition,lightConditions);
featureTable = featureTable(lightLogInd,:);
% Filter for strains
strainLogInd = ismember(featureTable.strain_name,strains);
featureTable = featureTable(strainLogInd,:);
% filter out files with too few skeletons
skelLogInd = featureTable.n_skeletons > n_skeletons_range(1) & featureTable.n_skeletons < n_skeletons_range(2);
featureTable = featureTable(skelLogInd,:);
% subsample data if specified
if isscalar(n_subsample)
    featureTable = subsampleData(featureTable,n_subsample);
end
% filter featureTable based on specified strain and features
[featureTable, classLabels] = filterFeatureTable(featureTable,classVar,strains,featSetName);

% add retained labels to the end of featureTable
for varCtr = 1:numel(classVar)
    var = classVar{varCtr};
    featureTable.(var) = classLabels.(var);
end

%% Extract featureMat
featureMat = table2array(featureTable(:,feats));
n_strains = numel(unique(featureTable.strain_name));

%% Analyze features with PCA
% pre-process feature matrix for PCA
[featureMat,dropLogInd] = preprocessFeatMat(featureMat);
n_feats = size(featureMat,2);
% do pca
[pc, score, ~, ~, explained] = pca(featureMat);

%% Set outliers from each PC to zero
if removeOutlier
    disp(['Removing outliers that fall outside of mean +/- ' num2str(n_std2include) ' deviations range.'])
    for PCCtr = 1:3
        coeff = score(:,PCCtr); % get coefficients for this PC
        coeffRange = [mean(coeff)-n_std2include*std(coeff), mean(coeff)+n_std2include*std(coeff)]; % get acceptable range of PC coefficients
        outlierLogInd = coeff < coeffRange(1) | coeff > coeffRange(2); % find outlier indices
        coeff(outlierLogInd) = 0; % set outlier coefficient to 0
        score(:,PCCtr) = coeff; % write back into scores
        disp([num2str(nnz(outlierLogInd)) ' outliers removed from PC ' num2str(PCCtr)])
    end
end

%% Get strain logical index (range of n per strain: [222 414])
N2LogInd = strcmp(featureTable.strain_name,'N2');
CB4856LogInd = strcmp(featureTable.strain_name,'CB4856');
MY23LogInd = strcmp(featureTable.strain_name,'MY23');
QX1410LogInd = strcmp(featureTable.strain_name,'QX1410');
VX34LogInd = strcmp(featureTable.strain_name,'VX34');
NIC58LogInd = strcmp(featureTable.strain_name,'NIC58');
JU1373LogInd = strcmp(featureTable.strain_name,'JU1373');

%% Plot first two PCs and colour by strain

% % PC 1 and 2
% figure; hold on
% %
% plot(score(N2LogInd,1),score(N2LogInd,2),'b*')
% plot(score(CB4856LogInd,1),score(CB4856LogInd,2),'bo')
% plot(score(MY23LogInd,1),score(MY23LogInd,2),'b.')
% plot(score(QX1410LogInd,1),score(QX1410LogInd,2),'m*')
% plot(score(VX34LogInd,1),score(VX34LogInd,2),'mo')
% plot(score(NIC58LogInd,1),score(NIC58LogInd,2),'g*')
% plot(score(JU1373LogInd,1),score(JU1373LogInd,2),'go')
% %
% xlabel(['PC1 (' num2str(round(explained(1))) ')%'])
% ylabel(['PC2 (' num2str(round(explained(2))) ')%'])
% legend({'N2','CB4856','MY23','QX1410','VX34','NIC58','JU1373'})
% title(['PCA plot with ' num2str(n_strains) ' strains and ' num2str(n_feats) ' features'])
% 
% % PC 1 and 3
% figure; hold on
% %
% plot(score(N2LogInd,1),score(N2LogInd,3),'b*')
% plot(score(CB4856LogInd,1),score(CB4856LogInd,3),'bo')
% plot(score(MY23LogInd,1),score(MY23LogInd,3),'b.')
% plot(score(QX1410LogInd,1),score(QX1410LogInd,3),'m*')
% plot(score(VX34LogInd,1),score(VX34LogInd,3),'mo')
% plot(score(NIC58LogInd,1),score(NIC58LogInd,3),'g*')
% plot(score(JU1373LogInd,1),score(JU1373LogInd,3),'go')
% %
% xlabel(['PC1 (' num2str(round(explained(1))) ')%'])
% ylabel(['PC3 (' num2str(round(explained(2))) ')%'])
% legend({'N2','CB4856','MY23','QX1410','VX34','NIC58','JU1373'})
% title(['PCA plot with ' num2str(n_strains) ' strains and ' num2str(n_feats) ' features'])
% 
% % PC 2 and 3
% figure; hold on
% %
% plot(score(N2LogInd,2),score(N2LogInd,3),'b*')
% plot(score(CB4856LogInd,2),score(CB4856LogInd,3),'bo')
% plot(score(MY23LogInd,2),score(MY23LogInd,3),'b.')
% plot(score(QX1410LogInd,2),score(QX1410LogInd,3),'m*')
% plot(score(VX34LogInd,2),score(VX34LogInd,3),'mo')
% plot(score(NIC58LogInd,2),score(NIC58LogInd,3),'g*')
% plot(score(JU1373LogInd,2),score(JU1373LogInd,3),'go')
% %
% xlabel(['PC2 (' num2str(round(explained(1))) ')%'])
% ylabel(['PC3 (' num2str(round(explained(2))) ')%'])
% legend({'N2','CB4856','MY23','QX1410','VX34','NIC58','JU1373'})
% title(['PCA plot with ' num2str(n_strains) ' strains and ' num2str(n_feats) ' features'])

%% 3D plot of the first three PCs and colour by strain
%
figure;
%
scatter3(score(N2LogInd,1),score(N2LogInd,2),score(N2LogInd,3),'b*','LineWidth',1)
hold on
scatter3(score(CB4856LogInd,1),score(CB4856LogInd,2),score(CB4856LogInd,3),'bo','LineWidth',1)
scatter3(score(MY23LogInd,1),score(MY23LogInd,2),score(MY23LogInd,3),'b.','LineWidth',1)
scatter3(score(QX1410LogInd,1),score(QX1410LogInd,2),score(QX1410LogInd,3),'m*','LineWidth',1)
scatter3(score(VX34LogInd,1),score(VX34LogInd,2),score(VX34LogInd,3),'mo','LineWidth',1)
scatter3(score(NIC58LogInd,1),score(NIC58LogInd,2),score(NIC58LogInd,3),'g*','LineWidth',1)
scatter3(score(JU1373LogInd,1),score(JU1373LogInd,2),score(JU1373LogInd,3),'go','LineWidth',1)

%
xlabel(['PC1 (' num2str(round(explained(1))) ')%'])
ylabel(['PC2 (' num2str(round(explained(2))) ')%'])
zlabel(['PC3 (' num2str(round(explained(3))) ')%'])
%
legend(strains)
title(['PCA plot with ' num2str(n_strains) ' strains and ' num2str(n_feats) ' features'])
set(gca,'fontsize',15)

%% Get logical index for each strain 
strainLogInd = false(numel(featureTable.strain_name),numel(strains));
for strainCtr = 1:numel(strains)
    strain = strains{strainCtr};
    strainLogInd(:,strainCtr) = strcmp(featureTable.strain_name,strain);
end

%% 3D plot of the first three PCs showing one strain at a time
%
for strainCtr = 1:numel(strains)
    strain = strains{strainCtr};
    figure;
    % plot the strain of interest in red
    scatter3(score(strainLogInd(:,strainCtr),1),score(strainLogInd(:,strainCtr),2),score(strainLogInd(:,strainCtr),3),...
        7,'red','filled')
    hold on
    % plot the rest of the strains in gray
    scatter3(score(~strainLogInd(:,strainCtr),1),score(~strainLogInd(:,strainCtr),2),score(~strainLogInd(:,strainCtr),3),...
        7,[0.7,0.7,0.7],'filled')
    xlabel(['PC1 (' num2str(round(explained(1))) ')%'])
    ylabel(['PC2 (' num2str(round(explained(2))) ')%'])
    zlabel(['PC3 (' num2str(round(explained(3))) ')%'])
    %
    legend(strain)
    title(['PCA plot with ' num2str(n_feats) ' features'])
    set(gca,'fontsize',15)
end

%% Get logical index for each light condition 
lightLogInd = false(numel(featureTable.light_condition),numel(lightConditions));
for lightCtr = 1:numel(lightConditions)
    light = lightConditions{lightCtr};
    lightLogInd(:,lightCtr) = strcmp(featureTable.light_condition,light);
end

%% 3D plot of the first three PCs and colour by light condition
%
if numel(lightConditions)>1
    figure;
    plotcolors = {'m','c','k'};
    lightCtr = 1;
    scatter3(score(lightLogInd(:,lightCtr),1),score(lightLogInd(:,lightCtr),2),score(lightLogInd(:,lightCtr),3),[plotcolors{lightCtr} 'o'],'LineWidth',1)
    hold on
    for lightCtr = 2:numel(lightConditions)
        scatter3(score(lightLogInd(:,lightCtr),1),score(lightLogInd(:,lightCtr),2),score(lightLogInd(:,lightCtr),3),[plotcolors{lightCtr} 'o'],'LineWidth',1)
    end
    
    xlabel(['PC1 (' num2str(round(explained(1))) ')%'])
    ylabel(['PC2 (' num2str(round(explained(2))) ')%'])
    zlabel(['PC3 (' num2str(round(explained(3))) ')%'])
    %
    legend(lightConditions)
    title(['PCA plot with ' num2str(numel(lightConditions)) ' light conditions and ' num2str(n_feats) ' features'])
    set(gca,'fontsize',15)
end

%% Get logical index for each experimental date
featureTable.date_yyyymmdd = string(featureTable.date_yyyymmdd); % turn date class into string
dates = unique(featureTable.date_yyyymmdd);
dateLogInd = false(numel(featureTable.date_yyyymmdd),numel(dates));
for dateCtr = 1:numel(dates)
    date = dates(dateCtr);
    dateLogInd(:,dateCtr) = strcmp(featureTable.date_yyyymmdd,date);
end

%% 3D plot of the first three PCs and colour by experimental date
%
if numel(dates)>1
    figure;
    plotcolors = {'m','c','k'};
    dateCtr = 1;
    scatter3(score(dateLogInd(:,dateCtr),1),score(dateLogInd(:,dateCtr),2),score(dateLogInd(:,dateCtr),3),[plotcolors{dateCtr} 'o'],'LineWidth',1)
    hold on
    for dateCtr = 2:numel(dates)
        scatter3(score(dateLogInd(:,dateCtr),1),score(dateLogInd(:,dateCtr),2),score(dateLogInd(:,dateCtr),3),[plotcolors{dateCtr} 'o'],'LineWidth',1)
    end
    
    xlabel(['PC1 (' num2str(round(explained(1))) ')%'])
    ylabel(['PC2 (' num2str(round(explained(2))) ')%'])
    zlabel(['PC3 (' num2str(round(explained(3))) ')%'])
    %
    legend(dates)
    title(['PCA plot with ' num2str(numel(dates)) ' experimental dates and ' num2str(n_feats) ' features'])
    set(gca,'fontsize',15)
end

%% Get logical index for each plate batch
featureTable.date_plates_poured_yyyymmdd = string(featureTable.date_plates_poured_yyyymmdd); % turn date class into string
batches = unique(featureTable.date_plates_poured_yyyymmdd);
batchLogInd = false(numel(featureTable.date_plates_poured_yyyymmdd),numel(batches));
for batchCtr = 1:numel(batches)
    batch = batches(batchCtr);
    batchLogInd(:,batchCtr) = strcmp(featureTable.date_plates_poured_yyyymmdd,batch);
end

%% 3D plot of the first three PCs and colour by plate batch
%
if numel(batches)>1
    figure;
    plotcolors = {'m','c','k'};
    batchCtr = 1;
    scatter3(score(batchLogInd(:,batchCtr),1),score(batchLogInd(:,batchCtr),2),score(batchLogInd(:,batchCtr),3),[plotcolors{batchCtr} 'o'],'LineWidth',1)
    hold on
    for batchCtr = 2:numel(batches)
        scatter3(score(batchLogInd(:,batchCtr),1),score(batchLogInd(:,batchCtr),2),score(batchLogInd(:,batchCtr),3),[plotcolors{batchCtr} 'o'],'LineWidth',1)
    end
    
    xlabel(['PC1 (' num2str(round(explained(1))) ')%'])
    ylabel(['PC2 (' num2str(round(explained(2))) ')%'])
    zlabel(['PC3 (' num2str(round(explained(3))) ')%'])
    %
    legend(batches)
    title(['PCA plot with ' num2str(numel(batches)) ' plate batches and ' num2str(n_feats) ' features'])
    set(gca,'fontsize',15)
end


%% Plot variance explained as a function of number of PCs
%
figure; hold on
%
subplot(1,2,1)
plot(cumsum(explained),'o','LineWidth',2); 
xlim([0 200])
ylim([0 100]); xlabel('Number of PCs'); ylabel('Variance explained (%)');
if strcmp(featSetName,'Tierpsy_256')
    xlim([0 80])
end
%
set(gca,'fontsize',15)
% 
subplot(1,2,2)
plot(explained,'o','LineWidth',2)
xlim([0 200])
xlabel('Number of PCs'); ylabel('Additional variance explained (%)');
%
set(gca,'fontsize',15)

%% See what's inside the first PC
% [feat,featInd] = sort(pc(:,1)); % PC1 
% featureTable.Properties.VariableNames(featInd)'

%% Clustergram
% rowLabels = featureTable.strain_name;
% colLabels = featureTable.Properties.VariableNames(1:end-numel(classVar));
% colLabels = colLabels(~dropLogInd);
% cgObj = clustergram(featureMat,'RowLabels',rowLabels,'ColumnLabels',colLabels,...
%     'Colormap',redbluecmap,'ShowDendrogram','on','OptimalLeafOrder',true)

%% Results
% Good separation at the species level using PC3 with either PC1 or 2 
% elegans and tropicalis separate better from each other than they do from
% briggsae
% PC2 seems good at separating reference (*) and divergent (o) strains
% No clear separation based on light condition or experimental date, but
% plate batch does seem to have an effect..