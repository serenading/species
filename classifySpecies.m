%% Script uses supervised machine learning algorithms to train classifiers for a specified variable
% based on extracted Tierpsy features. It also has the option to apply
% sequantial feature selection to identify top features to use for
% classification

% author: @serenading. Jan 2021.

% TODO: integrate the filterFeatureTable function from agg. repo. Need to
% add Tierpsy full feature list (abbreviated) first and do away with
% n_nonFeatVar.


clear
close all

addpath('../AggScreening/')
addpath('../AggScreening/strainsList/')

%% Specify which variable to classify for
classVar = 'light_condition'; % 'species_name','strain_name','light_condition'

%% Set analysis parameters
extractStamp = '20210112_105808'; % 20201218_184325 for standard feature extraction, 20210112_105808 for size filtered data
n_nonFeatVar = 33; % the first n columns of the feature table that do not contain features. =33
lightCondition = ''; % Leave empty '' to use all. 'prestim','bluelight','poststim'

% Set filering parameters
strains2keep = {'N2','CB4856','MY23','QX1410','VX34','NIC58','JU1373'}; 
strains2drop = {}; % {'VX34','NIC58'} Cell array containing strains to drop from analysis.

feats2keep = {'Tierpsy_256'}; % Use all features if left empty. {'Tierpsy_256'} or {'feat1','feat2'}. Cell array containing features to use for analysis. 
feats2drop = {}; % {'path'};

n_subsample = NaN; % number of replicates per strain to include. Set to NaN to include all samples
n_skeletons_range = [50 22500]; % n_skeleton range to use for retaining the well. 25fps x 60s/min x 5 min x 3 worms = 22500 skeletons.

% Select which tasks to perform
performSequentialFeatureSelection = false;
trainClassifier = true;

%% Optional parameters

% SFS parameters. Note: Currently using linear discriminant analysis for sfs. Otherwise redefine classf function.
if performSequentialFeatureSelection
    % Generate additional diagnostic plots for SFS
    plotSFS = true;
    if isempty(feats2keep)
        % Use top n sorted features with lowest p-values to make sfs run faster.
        % Also if n_feat too high for n_obs then some classifiers do not work.
        % Likely needs tweaking based on input observations.
        n_sortedFeatsInput = 500;
    else
        n_sortedFeatsInput = NaN;
    end
    n_topFeatsOutput = 30;
end

% Classification parameters. How many replicates to use for hold out and cross validation
holdoutRatio = 1/4; 
crossVal_k = 5;

%% Load and process featureTable

% Load featuresTable and set classification variables according to the task
featureTable = readtable(['/Users/sding/OneDrive - Imperial College London/species/Results/fullFeaturesTable_' extractStamp '.csv']);
species_names = getSpeciesnames(featureTable); % Get species information to featureTable
light_condition = getLightcondition(featureTable);

% Filter out poststim videos for light condition classification
dropLogInd = strcmp(light_condition,'bluelight');
featureTable = featureTable(~dropLogInd,:);
light_condition = light_condition(~dropLogInd);

% Get classLabels
if strcmp(classVar,'species_name')
    classLabels = species_names;
elseif strcmp(classVar,'strain_name')
    classLabels = featureTable.strain_name;
elseif strcmp(classVar,'light_condition')
    classLabels = light_condition;
end

% trim down featuresTable to contain only features
if isempty(feats2keep)
    featureTable = featureTable(:,n_nonFeatVar+1:end);
elseif strcmp(feats2keep,'Tierpsy_256')
    feats2keep = table2cell(readtable('strainsList/Tierpsy_256_short.csv','PreserveVariableNames',true,'ReadVariableNames',false));
    featureTable = featureTable(:,feats2keep);
else
    featureTable = featureTable(:,feats2keep);
end

%% Pre-process features matrix and turn back into table

% Split table into matrix and featNames
featureMat = table2array(featureTable);
featNames = featureTable.Properties.VariableNames;
% Preprocess feature matrix: drop zero standard deviation, NaN's, z-normalise, etc.
[featureMat,dropLogInd] = preprocessFeatMat(featureMat);
featNames = featNames(~dropLogInd);
% Put the table back together
featureTable = array2table(featureMat,'VariableNames',featNames);

%% Holdout partition to separate training and test set

% Generate test/training partition
holdoutCVP = cvpartition(classLabels,'holdout',holdoutRatio);

% Extract training and test datasets
dataTrain = featureMat(holdoutCVP.training,:);
grpTrain = classLabels(holdoutCVP.training,:);
dataTest = featureMat(holdoutCVP.test,:);
grpTest = classLabels(holdoutCVP.test,:);

%% Sequential feature selection
% Following instructions from https://uk.mathworks.com/help/stats/examples/selecting-features-for-classifying-high-dimensional-data.html

if performSequentialFeatureSelection
    
    % ecdf plot of anova test values
    grpTrainInd = grp2idx(grpTrain);
    % One-way anova for multiple class labels
        p = NaN(1,size(dataTrain,2));
        for featCtr = 1:size(dataTrain,2)
            p(featCtr) = anova1(dataTrain(:,featCtr),grpTrainInd,'off');
        end
    if plotSFS
        figure; ecdf(p);
        xlabel('P value');
        ylabel('CDF value')
    end
    
    % Sort features
    [~,featureIdxSortbyP] = sort(p,2);
    
    % Define function for sfs
    
    % subspace discriminant model - TODO: needs to write this as a proper function to enable multiple lines
    %     Mdl = fitcensemble(xtrain,ytrain,'Method','subspace'); % trained object using subspace discriminant model.
    %     classf = @(xtrain,ytrain,xtest,ytest) ...
    %         sum(~strcmp(ytest,fitcensemble(xtrain,ytrain,'Method','subspace').predictFcn(xtest)));
    
    % Linear discriminant model
    classf = @(xtrain,ytrain,xtest,ytest) ...
        sum(~strcmp(ytest,classify(xtest,xtrain,ytrain,'linear')));
    
    % Sequential feature selection with k-fold cross-validation
    kfoldCVP = cvpartition(grpTrain,'kfold',crossVal_k);
    if ~isnan(n_sortedFeatsInput)
        featureIdxSortbyP = featureIdxSortbyP(1:n_sortedFeatsInput);
    end
    [fsLocal,historyCV] = sequentialfs(classf,dataTrain(:,featureIdxSortbyP),grpTrain,'cv',kfoldCVP,'nfeatures',n_topFeatsOutput);
    featCols = featureIdxSortbyP(fsLocal);
    keyFeats = featureTable.Properties.VariableNames(featCols)'
    
    if plotSFS
        % Evaluate the performance of the selected model with the features
%         testMCELocal = crossval(classf,featureMat(:,featureIdxSortbyP(fsLocal)),classLabels,'partition',...
%             holdoutCVP)/holdoutCVP.TestSize;
        figure;
        
        % Plot of the cross-validation MCE as a function of the number of top features
        subplot(1,2,1)
        plot(historyCV.Crit,'o');
        xlabel('Number of Features');
        ylabel('Cross-validation Misclassification Error');
        ylim([0 1])
        
        % Plot the derivative of the cross-validation MCE as a function of the number of top features
        dCrit = [historyCV.Crit 0]-[1 historyCV.Crit];
        subplot(1,2,2)
        plot(dCrit,'o');
        xlim([0 numel(historyCV.Crit)])
        xlabel('Number of Features');
        ylabel('Change in CV MSE');
        
        %         % plot of resubstitution MCE values on the training set
        %         % (i.e., without performing cross-validation during the feature selection procedure)
        %         % as a function of the number of top features:
        %         [fsResub,historyResub] = sequentialfs(classf,dataTrain(:,featureIdxSortbyP),grpTrain,...
        %             'cv','resubstitution','nfeatures',n_topFeatsOutput);
        %         figure; plot(1:n_topFeatsOutput, historyCV.Crit,'bo',1:n_topFeatsOutput, historyResub.Crit,'r^');
        %         xlabel('Number of Features');
        %         ylabel('MCE');
        %         legend({[num2str(crossVal_k) '-fold CV MCE'], 'Resubstitution MCE'},'location','NE');
    end
    
    % Slice out key features from the full featureTable
    featCols = [];
    for keyFeatCtr = 1:numel(keyFeats)
        keyFeat = keyFeats{keyFeatCtr};
        keyFeat = string(keyFeat);
        featCols = [featCols, find(strcmp(featureTable.Properties.VariableNames,keyFeat))];
    end
    featureTable = featureTable(:,featCols);
    disp([num2str(numel(featCols)) ' key features retained.'])
end

%% Train classifier

if trainClassifier
    
    % Add classLabels to featureTable for classification
    featureTable.(classVar) = classLabels;
    
    % Divide featureTable into training and test sets
    trainFeatureTable = featureTable(holdoutCVP.training,:);
    testFeatureTable = featureTable(holdoutCVP.test,:);
    
    % Perform model selection using classificationLearner GUI
    disp(['Please use GUI to train model with trainFeatureTable. Use ' num2str(crossVal_k) '-fold cross validation, select the best model, and export to workspace.'])
    classificationLearner
    
    % Wait until trained model has been selected and saved
    while exist('trainedModel')~=1
        pause
    end
    
    %% Classification accuracy on unseen hold-out test set
    yfit = trainedModel.predictFcn(testFeatureTable);
    accuracy = nnz(strcmp(yfit,testFeatureTable.(classVar)))/holdoutCVP.TestSize;
    disp(['Test accuracy from trained model is ' num2str(accuracy) ' on unseen data.'])
    
end

%% Results
% For classifying species:
% - linear discriminant, quadratic SVM, and subspace discriminant ensemble models
% work well (> 90% cross-validation and test accuracies)
% - SFS shows 4 features are good for achieving high
% accuracy
%     {'curvature_midbody_w_forward_abs_50th'                   }
%     {'d_curvature_mean_hips_w_backward_abs_90th'              }
%     {'curvature_std_midbody_abs_50th'                         }
%     {'curvature_std_neck_w_forward_abs_50th'                  }

% For classifying strains:
% - linear discriminant, quadratic SVM, and subspace discriminant ensemble models
% work well (> 80% cross-validation and test accuracies)
% - SFS shows 8 features are good for achieving high
% accuracy
%     {'minor_axis_w_forward_10th'                             }
%     {'curvature_std_neck_w_forward_abs_50th'                 }
%     {'d_rel_to_neck_radial_vel_head_tip_w_forward_90th'      }
%     {'curvature_midbody_w_forward_abs_50th'                  }
%     {'ang_vel_tail_base_w_forward_abs_50th'                  }
%     {'d_curvature_mean_midbody_abs_90th'                     }
%     {'minor_axis_10th'                                       }
%     {'rel_to_neck_radial_vel_head_tip_w_forward_90th'        }

% For classifying three light conditions: prestim + bluelight + poststim:
% - quadratic SVM and ensemble subspace discriminant models give ~50% cross-validation and test accuracies.
% Accuracies are higher for bluelight than for pre- and post-stim. Other
% methods ~40-50% cross-validation accuracies. Random guess 33%. 
% - SFS shows 4 features before accuracy plateaus at 50%
%     {'path_coverage_head'                                  }
%     {'ang_vel_head_tip_abs_50th'                           }
%     {'speed_w_forward_IQR'                                 }
%     {'ang_vel_head_base_w_forward_abs_50th'                }

% For classifying two light conditions: prestim + bluelight:
% - quadratic SVM and ensemble subspace discriminant models give ~65-70% cross-validation and test accuracies.
% - SFS shows 7 features before accuracy plateaus at 70%
%     {'path_coverage_head'                                 }
%     {'speed_w_forward_90th'                               }
%     {'speed_head_tip_w_forward_10th'                      }
%     {'path_transit_time_midbody_95th'                     }
%     {'minor_axis_w_forward_10th'                          }
%     {'d_ang_vel_midbody_w_backward_abs_IQR'               }
%     {'d_ang_vel_head_tip_abs_10th'                        }
% - results seem to suggest some strains are not bluelight sensitive to give very high accuracies between light conditions

% For classifying two light conditions: prestim + poststim:
% - quadratic SVM and ensemble subspace discriminant models give ~57-60%
% cross-validation and test accuracies. Random guess gives 50%.