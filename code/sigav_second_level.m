function [matfile_second_level, matfile_first_level] = sigav_second_level(...
    data_matrix_files, para_files, ...
    condition_names_file, output_directory, onset_delay, ...
    offset_delay, remove_run_offsets, overwrite)

% Combines results from several runs. Individual run analyses computed by
% sigav.m
% 
% 2016-08-25 - Created, Sam NH

% default output directory
if nargin < 4
    output_directory = pwd;
    fprintf('No output directory specified\nSaving to working directory\m')
end

if nargin < 5
    onset_delay = 5;
end

if nargin < 6
    offset_delay = 1;
end

if nargin < 7
    remove_run_offsets = true;
end

if nargin < 8
    overwrite = false;
end

n_runs = length(data_matrix_files);
matfile_first_level = cell(1, n_runs);
for i = 1:n_runs
    
    % individual run analysis
    matfile_first_level{i} = [output_directory '/r' num2str(i) '.mat'];
    if ~exist(matfile_first_level{i}, 'file') || overwrite
        sigav(data_matrix_files{i}, para_files{i}, ...
            condition_names_file, matfile_first_level{i}, ...
            onset_delay, offset_delay);
    end
    
end

% load results form individual runs
matfile_second_level = [output_directory '/allruns.mat'];
if ~exist(matfile_second_level, 'file') || overwrite
    
    for i = 1:n_runs
        load(matfile_first_level{i}, ...
            'null_response', 'condition_responses', 'mean_signal');
        
        % initialize
        if i == 1
            mean_signal_allruns = nan([size(mean_signal), n_runs]);
            null_response_allruns = nan([size(null_response), n_runs]); %#ok<NODEF>
            condition_responses_allruns = nan([size(condition_responses), n_runs]); %#ok<NODEF>
        end
        
        % error check
        assert(sum(~isnan(condition_responses(:)))>0);
        
        % assign
        mean_signal_allruns(:,:,i) = mean_signal;
        null_response_allruns(:,:,i) = null_response;
        condition_responses_allruns(:,:,i) = condition_responses;
        clear null_response condition responses mean_signal;
    end
    
    % remove differences in mean signal across runs
    if remove_run_offsets
        
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
    psc = 100 * (condition_responses - X) ./ X;
    clear X;
    
    % save results
    load(condition_names_file, 'condition_names');
    mean_signal = mean( mean_signal_allruns, 3 );
    save(matfile_second_level, 'psc', 'condition_responses', ...
        'null_response', 'mean_signal', 'condition_names', '-v7.3');
end

