function [boxcar_matrix, condition_names, t] = boxcar_from_para(para_file, sr)

% [unique_conds,~,unique_condition_indices] = unique(P.conds);
% 
% n_conds = length(unique_conds);

try
P = read_para(para_file);
catch
    keyboard
end
condition_names = P.conds;  

n_condition_onsets = length(P.onsets);
pad_duration_in_seconds = 50;
n_samples = ceil((P.onsets(end) + P.durs(end) + pad_duration_in_seconds) * sr);

boxcar_matrix = zeros(n_samples, n_condition_onsets);

sampling_period = 1/sr;
t = (0:n_samples-1)'/sr;
sampling_windows = [t-sampling_period/2, t+sampling_period/2];

condition_onsets = P.onsets;
condition_offsets = P.onsets+P.durs;

for i = 1:n_samples
    
    first_point_in_bounds = sampling_windows(i,1) > condition_onsets ...
        & sampling_windows(i,1) < condition_offsets;
    second_point_in_bounds = sampling_windows(i,2) > condition_onsets ...
        & sampling_windows(i,2) < condition_offsets;
    
    boxcar_matrix(i, first_point_in_bounds & second_point_in_bounds) = 1;
    
    xi = first_point_in_bounds & ~second_point_in_bounds;
    boxcar_matrix(i, xi) = (condition_offsets(xi) - sampling_windows(i,1))/sampling_period;

    xi = ~first_point_in_bounds & second_point_in_bounds;
    boxcar_matrix(i, xi) = (sampling_windows(i,2) - condition_onsets(xi))/sampling_period;
    
end