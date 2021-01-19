%% Function generates 3x2 windows for each of the light condition in seconds based on standard Hydra bluelight stimulation conditions.

function [prelight, bluelight, postlight] = getBluelightFeatWindows(bluelightInterval,windowDuration)

%% INPUTS: 
% bluelightInterval: [60,70; 160,170; 260,270] (default). nx2 matrix, each row specifying the
% start and end time (in seconds) of the bluelight stimulation window. 
% windowDuration: 10 (default). scalar specifying the time window to retain for feature
% analysis.

%% OUTPUTS:
% prelight: nx2 matrix, each row representing the start and end time (in
% seconds) of the 10 seconds preceding each bluelight interval.
% bluelight: nx2 matrix, each row representing the start and end time (in
% seconds) of the 10 seconds flanking each end of the bluelight interval.
% postlight: nx2 matrix, each row representing the start and end time (in
% seconds) of the 10 seconds immediately following the bluelight windows. 

%% FUNCTION

if nargin<2
    windowDuration = 10; % default bluelight interval is 10 seconds
    if nargin <1
         bluelightInterval: [60,70; 160,170; 260,270]; % default bluelight intervals are 10 second intervals starting at 60, 160, 260s into the bluelight video.
    end
end

% pre-allocate window sizes
prelight = NaN(size(bluelightInterval));
bluelight = NaN(size(bluelightInterval));
postlight = NaN(size(bluelightInterval));

% get light condition windows in seconds
prelight(:,2) = bluelightInterval(:,1);
prelight(:,1) = prelight(:,2) - windowDuration;
bluelight(:,1) = bluelightInterval(:,2) - windowDuration/2;
bluelight(:,2) = bluelightInterval(:,2) + windowDuration/2;
postlight(:,1) = bluelight(:,2);
postlight(:,2) = postlight(:,1) + windowDuration;

end