function n_worms_estimate = get_n_worms_estimate(frame_numbers)

%% Function gets an estimate of the number of worms using the frame_number vector.
% Function is adapted from tierpsy/features/tierpsy_features/helper.py
% Author: @serenading. Feb 2021. 

percentile = 99;

n_per_frame = groupcounts(frame_numbers);

if numel(n_per_frame)>0
    n_worms_estimate = prctile(n_per_frame,percentile);
else
    n_worms_estimate = 0;
end

end