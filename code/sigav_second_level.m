function sigav_second_level(...
    data_matrix_files, para_files, condition_names_file, ...
    MAT_file_second_level, MAT_files_first_level, varargin)

% Combines results from several runs. Individual run analyses computed by
% sigav.m
% 
% 2016-08-25 - Created, Sam NH

% optional arguments and defaults
I.overwrite = false;
I.onset_delay = 5;
I.offset_delay = 1;
I.remove_run_offsets = true;
I = parse_optInputs_keyvalue(varargin, I);

assert(length(data_matrix_files) == length(para_files));
assert(length(data_matrix_files) == length(MAT_files_first_level));

% first level analysis
n_runs = length(data_matrix_files);
for i = 1:n_runs
    
    % individual run analysis
    if ~exist(MAT_files_first_level{i}, 'file') || I.overwrite
        sigav(data_matrix_files{i}, para_files{i}, ...
            condition_names_file, MAT_files_first_level{i}, ...
            I.onset_delay, I.offset_delay);
    end
    
end

% load results form individual runs
if exist(MAT_file_second_level, 'file') && ~I.overwrite
    return;
end

for i = 1:n_runs
    load(MAT_files_first_level{i}, ...
        'null_response', 'condition_responses', 'mean_signal');
    
    % initialize
    if i == 1
        mean_signal_allruns = nan([size(mean_signal), n_runs]); %#ok<NODEF>
        null_response_allruns = nan([size(null_response), n_runs]); %#ok<NODEF>
        condition_responses_allruns = nan([size(condition_responses), n_runs]); %#ok<NODEF>
    end
    
    % error check
    try
        assert(sum(~isnan(condition_responses(:)))>0);
    catch
        keyboard;
    end
    
    % assign
    mean_signal_allruns(:,:,i) = mean_signal;
    null_response_allruns(:,:,i) = null_response;
    condition_responses_allruns(:,:,i) = condition_responses;
    clear null_response condition responses mean_signal;
end

% remove differences in mean signal across runs
if I.remove_run_offsets
    
    run_offsets = mean_signal_allruns ...
        - repmat( mean( mean_signal_allruns, 3 ), [1,1,n_runs]);
    
    % correct for baseline differences
    null_response_allruns = null_response_allruns - run_offsets;
    condition_responses_allruns = condition_responses_allruns ...
        - repmat(run_offsets, [size(condition_responses,1), 1, 1]);
    
end

% average across runs
condition_responses = nanmean(condition_responses_allruns, 3);
null_response = mean(null_response_allruns, 3);

% convert to % signal chance
X = repmat(null_response, size(condition_responses,1), 1);
psc = 100 * (condition_responses - X) ./ X; %#ok<NASGU>
clear X;

% save results
load(condition_names_file, 'condition_names');
mean_signal = mean( mean_signal_allruns, 3 ); %#ok<NASGU>
save(MAT_file_second_level, 'psc', 'condition_responses', ...
    'null_response', 'mean_signal', 'condition_names', '-v7.3');

