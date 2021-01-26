%% Script specifies useful feature sets for analysis.
% Re-run script after modifications to update the saved feature set file
% Author: @serenading. Jan 2021.

featSetName = 'Tierpsy_256';
features.(featSetName) = table2cell(readtable('featureSet/Tierpsy_256_short.csv',"ReadVariableNames",false))';

featSetName = 'Tierpsy_all';
featureTable = readtable('/Users/sding/OneDrive - Imperial College London/species/Results/fullFeaturesTable_20210112_105808.csv');
features.(featSetName) = featureTable.Properties.VariableNames(34:end)';

featSetName = 'angular_velocity';
features.(featSetName) = {'ang_vel_head_tip_abs_90th','ang_vel_head_tip_abs_50th','ang_vel_head_base_abs_90th','ang_vel_head_base_abs_50th'};

featSetName = 'head_curvature';
features.(featSetName) = {'curvature_head_abs_90th','curvature_head_abs_50th','curvature_head_abs_10th'};

featSetName = 'head_curvature_norm';
features.(featSetName) = {'curvature_head_norm_abs_90th','curvature_head_norm_abs_50th','curvature_head_norm_abs_10th'};

featSetName = 'midbody_curvature';
features.(featSetName) = {'curvature_midbody_abs_90th','curvature_midbody_abs_50th','curvature_midbody_abs_10th'};

featSetName = 'midbody_curvature_norm';
features.(featSetName) = {'curvature_midbody_norm_abs_90th','curvature_midbody_norm_abs_50th','curvature_midbody_norm_abs_10th'};

featSetName = 'midbody_speed';
features.(featSetName) = {'speed_midbody_90th','speed_midbody_50th','speed_midbody_10th'};

featSetName = 'midbody_speed_norm';
features.(featSetName) = {'speed_midbody_norm_90th','speed_midbody_norm_50th','speed_midbody_norm_10th'};

featSetName = 'motion_mode';
features.(featSetName) = {'motion_mode_forward_fraction','motion_mode_paused_fraction','motion_mode_backward_fraction'};

featSetName = 'path_coverage';
features.(featSetName) = {'path_coverage_head','path_coverage_body','path_coverage_midbody','path_coverage_tail'};

featSetName = 'path_coverage_norm';
features.(featSetName) = {'path_coverage_head_norm','path_coverage_body_norm','path_coverage_midbody_norm','path_coverage_tail_norm'};

save('featureSet/features.mat','features');