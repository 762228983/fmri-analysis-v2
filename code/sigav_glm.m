function MAT_file = sigav_glm(data_matrix_file, para_file, parameter_file, ...
    MAT_file, varargin)

% A GLM analysis is applied to signal averaged responses. 
% 
% Calculates percent signal change for a set of events/stimuli by signal
% averaging a fixed number of time-points after event onset. A set of predictors
% is then regressed against the signal averaged values, and the beta weights
% from this analysis are contrasted. Stats are computed using vanilla OLS
% equations, assuming independent, Gaussian errors, and a permutation test is
% used to compute stats by shuffling the order of conditions.
% 
% 2016-09-09: Last edited, Sam NH

% optional arguments and defaults
I.onset_delay = 5;
I.offset_delay = 1;
I = parse_optInputs_keyvalue(varargin, I);

% analysis parameters
P = load(parameter_file);

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
% 1 x voxels
mean_signal = mean(Y,1); 

%% Average response to each event

% timing information about each event
T = read_para(para_file);

n_trials = length(T.onsets);
n_TR = size(Y,1);
t = (0:n_TR-1)*TR;
response = nan(n_trials, n_voxels_without_NaN);
for i = 1:n_trials
    xi = t >= T.onsets(i) + I.onset_delay ...
        & t <= T.onsets(i) + T.durs(i) + I.offset_delay;
    response(i,:) = mean(Y(xi,:),1);
end
clear xi t nTR n_trials;

%% Measure null and convert to psc

% mean response to null
% -> 1 x voxels
xi = strcmp('NULL', T.conds);
assert(sum(xi) > 0);
null_response = mean(response(xi,:),1);

% convert to psc
% -> trials x voxels
psc = response - repmat(null_response, n_trials, 1) ...
    ./ repmat(null_response, n_trials, 1);

%% Weights

% weights applied to each condition for each regressor
% -> condition x regressor
n_regressors = length(P.regressor_names);
W = zeros(n_trials,n_regressors);
for i = 1:n_trials
    xi = strcmp(T.conds{i}, P.condition_names);
    if any(xi)
        W(i,:) = P.regressor_weights(xi,:);
    end
end

% weights with one regressor per condition
n_conditions = length(P.condition_names);
W_one_per_condition = zeros(n_trials,n_conditions);
for i = 1:n_trials
    xi = strcmp(T.conds{i}, P.condition_names);
    W_one_per_condition(i,xi) = 1;
end

%% Regression analyses

% regress weighted event matrix
% -> contrasts x voxels
[beta_contrast, logP_ols, contrast_variance, df, residual] = ...
    regress_stats_ols( psc, W, P.contrast_weights, 'demean', false); %#ok<ASGLU>

% separate beta weight for regressor
% redundant with previous analysis if contrast matrix is the identity
% -> regressors x voxels
beta_one_per_regressor = ...
    regress_stats_ols(psc, W, eye(n_regressors), 'demean', false);

% separate beta weight for each condition
% -> condition x voxels
beta_one_per_condition = ...
    regress_stats_ols(psc, W_one_per_condition, eye(n_conditions), 'demean', false);

%% Fill in NaN entries and save

% fill in NaN entries
psc = fillin_NaN_voxels(psc, voxels_without_NaN, 2); %#ok<NASGU>
null_response = fillin_NaN_voxels(null_response, voxels_without_NaN, 2); %#ok<NASGU>
mean_signal = fillin_NaN_voxels(mean_signal, voxels_without_NaN, 2); %#ok<NASGU>
beta_contrast = fillin_NaN_voxels(beta_contrast, voxels_without_NaN, 2);
logP_ols = fillin_NaN_voxels(logP_ols, voxels_without_NaN, 2); %#ok<NASGU>
contrast_variance = fillin_NaN_voxels(contrast_variance, voxels_without_NaN, 2); %#ok<NASGU>
residual = fillin_NaN_voxels(residual, voxels_without_NaN, 2);
beta_one_per_condition = fillin_NaN_voxels(beta_one_per_condition, voxels_without_NaN, 2); %#ok<NASGU>
beta_one_per_regressor = fillin_NaN_voxels(beta_one_per_regressor, voxels_without_NaN, 2); %#ok<NASGU>

% save
save(MAT_file, 'psc', 'null_response', ...
    'mean_signal', 'voxels_without_NaN', 'beta_contrast', ...
    'logP_ols', 'contrast_variance', 'beta_one_per_condition', ...
    'beta_one_per_regressor', 'df', 'P', 'residual', '-v7.3');

%% Permutation test

% number of permuted weights
% optionally perform permutation test
if I.n_perms > 0
    
    % estimate contrasts and residual from permutations
    beta_contrast_permtest = nan(...
        [I.n_perms, length(P.contrast_names), sum(voxels_without_NaN)]);
    residual_permtest = nan([I.n_perms, sum(voxels_without_NaN)]);
    for i = 1:I.n_perms
        [beta_contrast_permtest(i,:,:), ~, ~, ~, residual_permtest(i,:)] = ...
            regress_stats_ols(Y, W(randperm(n_trials),:), P.contrast_weights);
    end
    clear X;
    
    % convert to signed logP value
    logP_permtest = sig_via_null_gaussfit(...
        beta_contrast(:,voxels_without_NaN), beta_contrast_permtest);
    logP_residual_permtest = sig_via_null_gaussfit(...
        residual(:,voxels_without_NaN), residual_permtest);
    logP_residual_permtest = -logP_residual_permtest;
    
    % fill in NaN entries
    beta_contrast_permtest = fillin_NaN_voxels(beta_contrast_permtest, voxels_without_NaN, 3); %#ok<NASGU>
    residual_permtest = fillin_NaN_voxels(residual_permtest, voxels_without_NaN, 2); %#ok<NASGU>
    logP_permtest = fillin_NaN_voxels(logP_permtest, voxels_without_NaN, 2); %#ok<NASGU>
    logP_residual_permtest = fillin_NaN_voxels(logP_residual_permtest, voxels_without_NaN, 2); %#ok<NASGU>

    % save results to matfile
    save(MAT_file, 'beta_contrast_permtest', 'logP_permtest',...
        'residual_permtest', 'logP_residual_permtest', '-append');
    
end

% v2
