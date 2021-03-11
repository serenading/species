%% Function calculates reversal frequency over a rolling window from Hydra timeseries data
% Author: @serenading. Feb 2021.

function eventFreq = calculateReversalFrequency(features,wellLogInd,windowSize,frameRate)

%% INPUTS:
% features: struct as read by features = h5read(filename_featuresN.hdf5,'/timeseries_data').
% wellLogInd: logical index for the multi-plate well in consideration.
% windowSize: scalar of rolling window size, in seconds, for event frequency calculation.
% frameRate: scalar of frames per second. Usually 25. 

%% OUTPUT:
% eventFreq: 1 x n_frames vector showing reversal frequency per worm per second at each frame.

%% FUNCTION:

% Get time series values for this well
x = features.timestamp(wellLogInd)+1; % add 1 for adjust for python indexing
y = features.motion_mode(wellLogInd); % get motion mode flags
worms = double(features.worm_index(wellLogInd));

% Estimate number of worms per frame
n_worms_estimate = get_n_worms_estimate(x);
%disp([num2str(n_worms_estimate) ' worms estimated to be present for this well.'])

% Offset motion mode flag by 1 for edge detection
y2 = NaN(size(y));
y2(2:end) = y(1:end-1);
ydiff = y2-y;
eventEdgeLogInd = false(size(y));
eventEdgeLogInd(find(ydiff == -2| ydiff == -1)) = true; % worms reverse either by switching from forward (1) to reverse (-1), or from paused (0) to reverse (-1).
eventCountTotal = nnz(eventEdgeLogInd);
%disp([num2str(eventCountTotal) ' reversal events detected.'])

% Offset worm index flag by 1 to detect changing trajectories
worms2 = NaN(size(worms));
worms2(2:end) = worms(1:end-1);
wormsdiff = worms2-worms;
wormsEdgeLogInd = false(size(worms));
wormsEdgeLogInd(find(wormsdiff ~=0)) = true;

% Check if the edge is where worm index switches and remove edge flag if so
eventEdgeLogInd(wormsEdgeLogInd) = false;
eventCountRemoved = eventCountTotal - nnz(eventEdgeLogInd);
%disp([num2str(eventCountRemoved) ' detected reversal events removed, because they are due to the change of trajectories.'])

% Get event count by frame
eventCount = uint8(eventEdgeLogInd);
eventCountByFrame = accumarray(x,eventCount);
% wormCountByFrame = accumarray(x,worms); % this gives very high values for some reason. why?

% Sum event count over specified rolling window
time = windowSize*frameRate; % get window size time in frames
eventMovSum = movsum(eventCountByFrame,time,'omitnan','Endpoints','shrink');

% Calculate reversal frequency
eventFreq = eventMovSum/n_worms_estimate/windowSize; % get frequency of event per worm per second

end