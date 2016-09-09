function MAT_file = sigav(data_matrix_file, para_file, ...
    condition_names_file, MAT_file, varargin)

% Calculates percent signal change for each condition by signal averaging a
% fixed number of time-points after stimulus onset

% optional arguments and defaults
I.onset_delay = 5;
I.offset_delay = 1;
I = parse_optInputs_keyvalue(varargin, I);

%% Format data matrix

% load data matrix
% -> time x voxel
load(data_matrix_file, 'data_matrix', 'TR');
Y = data_matrix;
clear data_matrix;

% remove voxels with NaN values
voxels_without_NaN = all(~isnan(Y));
Y = Y(:,voxels_without_NaN);
n_voxels_without_NaN = sum(voxels_without_NaN);

% average signal for each voxel
% -> 1 x voxels
mean_signal = mean(Y,1); 

%% Average response to each event

% timing information about each event
P = read_para(para_file);

n_trials = length(P.onsets);
n_TR = size(Y,1);
t = (0:n_TR-1)*TR;
response = nan(n_trials, n_voxels_without_NaN);
for i = 1:n_trials
    xi = t >= P.onsets(i) + I.onset_delay ...
        & t <= P.onsets(i) + P.durs(i) + I.offset_delay;
    response(i,:) = mean(Y(xi,:),1);
end
clear xi t nTR n_trials;

%% Average across events from the same condition and convert to PSC

% mean response to null
xi = strcmp('NULL', P.conds);
assert(sum(xi) > 0);
null_response = mean(response(xi,:),1);

% mean response to all other conditions
load(condition_names_file, 'condition_names');
n_conditions = length(condition_names); %#ok<USENS>
condition_responses = nan(n_conditions, n_voxels_without_NaN);
for i = 1:n_conditions
    xi = strcmp(condition_names{i}, P.conds);
    if ~isempty(xi)
        condition_responses(i,:) = mean(response(xi,:),1);
    end
end

X = repmat(null_response, n_conditions, 1);
psc = 100*(condition_responses - X) ./ X;
clear X;

%% Fill in NaN entries and save

% fill in NaN entries
psc = fillin_NaN_voxels(psc, voxels_without_NaN, 2); %#ok<NASGU>
condition_responses = fillin_NaN_voxels(condition_responses, voxels_without_NaN, 2); %#ok<NASGU>
null_response = fillin_NaN_voxels(null_response, voxels_without_NaN, 2); %#ok<NASGU>
mean_signal = fillin_NaN_voxels(mean_signal, voxels_without_NaN, 2); %#ok<NASGU>

% save
save(MAT_file, 'psc', 'condition_responses', 'null_response', ...
     'mean_signal', 'condition_names', '-v7.3');



