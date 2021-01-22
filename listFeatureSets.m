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

save('featureSet/features.mat','features');