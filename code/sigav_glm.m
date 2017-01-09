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
% 
% 2016-09-21: Modified to deal with zero regressors and contrasts, Sam NH
% 
% 2016-12-21: Made it possible to whiten the data before performing regression
% analysis. Doing so causes the regression analysis to perform LCMV.

% optional arguments and defaults
I.onset_delay = 5;
I.offset_delay = 1;
I.n_perms = 0;
I.whiten = false;
I.remove_unspecified_trials = false; 
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
clear xi t nTR;

%% Measure null and convert to psc

% mean response to null
% -> 1 x voxels
xi = strcmp('NULL', T.conds);
assert(sum(xi) > 0);
null_response = mean(response(xi,:),1);
clear xi;

%% Optionally remove trials that did not have a corresponding conditionrather than setting them to zero

if I.remove_unspecified_trials
    xi = ismember(T.conds, P.condition_names);
    response  = response(xi,:);
    T.conds = T.conds(xi,:);
    T.condition_indices = T.condition_indices(xi,:);
    T.onsets = T.onsets(xi,:);
    T.durs = T.durs(xi,:);
    n_trials = sum(xi);
    clear xi;
end

%% convert to psc

% -> trials x voxels
psc = 100 * (response - repmat(null_response, n_trials, 1)) ...
    ./ repmat(null_response, n_trials, 1);  

%% Weights

% weights applied to each condition for each regressor
% -> trial x regressor
n_regressors = length(P.regressor_names);
W = zeros(n_trials,n_regressors);
for i = 1:n_trials
    xi = strcmp(T.conds{i}, P.condition_names);
    if any(xi)
        W(i,:) = P.regressor_weights(xi,:);
    end
end
nonzero_regressors = any(W~=0,1);
W = W(:,nonzero_regressors);

% weights with one regressor per condition
% -> trial x condition
n_conditions = length(P.condition_names);
W_one_per_condition = zeros(n_trials,n_conditions);
for i = 1:n_trials
    xi = strcmp(T.conds{i}, P.condition_names);
    W_one_per_condition(i,xi) = 1;
end
nonzero_conditions = any(W_one_per_condition~=0,1);
W_one_per_condition = W_one_per_condition(:,nonzero_conditions);

% contrasts with any nonzero weights
nonzero_contrasts = any(P.contrast_weights(nonzero_regressors,:)~=0,1);
C = P.contrast_weights(nonzero_regressors, nonzero_contrasts);

%% Check that zero-mean contrasts remain so when excluding zero regressors

n_contrasts = length(P.contrast_names);
for i = 1:n_contrasts
    if abs(sum(P.contrast_weights(:,i))) < 1e-10 ...
            && abs(sum(P.contrast_weights(nonzero_regressors,i))) > 1e-10
        error('Contrast "%s" is no longer zero mean when excluding zero regressors\n',...
            P.contrast_names{i});
    end
end

%% Optionally whiten data and apply whitening matrix to the regressors as well

if I.whiten
    Z = (psc * psc')^(-1/2);
    psc = Z * psc;
    W = Z * W;
    W_one_per_condition = Z * W_one_per_condition;
end

%% Regression analyses

% regress weighted event matrix
% -> contrasts x voxels
[beta_contrast, logP_ols, contrast_variance, df, residual] = ...
    regress_stats_ols(psc, W, C, 'demean', false); %#ok<ASGLU>

% separate beta weight for regressor
% redundant with previous analysis if contrast matrix is the identity
% -> regressors x voxels
beta_one_per_regressor = ...
    regress_stats_ols(psc, W, eye(size(W,2)), 'demean', false);

% separate beta weight for each condition
% -> condition x voxels
beta_one_per_condition = ...
    regress_stats_ols(psc, W_one_per_condition, ...
    eye(size(W_one_per_condition,2)), 'demean', false);

%% Fill in NaN entries and save

% set NaN for zero contrasts
beta_contrast = fillin_NaN(beta_contrast, nonzero_contrasts, 1);
logP_ols = fillin_NaN(logP_ols, nonzero_contrasts, 1);
contrast_variance = fillin_NaN(contrast_variance, nonzero_contrasts, 1);
beta_one_per_regressor = fillin_NaN(beta_one_per_regressor, nonzero_regressors, 1);
beta_one_per_condition = fillin_NaN(beta_one_per_condition, nonzero_conditions, 1);

% set NaN for voxels with NaN
beta_contrast = fillin_NaN(beta_contrast, voxels_without_NaN, 2);
logP_ols = fillin_NaN(logP_ols, voxels_without_NaN, 2); %#ok<NASGU>
contrast_variance = fillin_NaN(contrast_variance, voxels_without_NaN, 2); %#ok<NASGU>
beta_one_per_regressor = fillin_NaN(beta_one_per_regressor, voxels_without_NaN, 2); %#ok<NASGU>
beta_one_per_condition = fillin_NaN(beta_one_per_condition, voxels_without_NaN, 2); %#ok<NASGU>
residual = fillin_NaN(residual, voxels_without_NaN, 2);
psc = fillin_NaN(psc, voxels_without_NaN, 2);
null_response = fillin_NaN(null_response, voxels_without_NaN, 2); %#ok<NASGU>
mean_signal = fillin_NaN(mean_signal, voxels_without_NaN, 2); %#ok<NASGU>

% save
save(MAT_file, 'psc', 'null_response', ...
    'mean_signal', 'voxels_without_NaN', 'beta_contrast', ...
    'logP_ols', 'contrast_variance', 'beta_one_per_condition', ...
    'beta_one_per_regressor', 'df', 'P', 'residual', '-v7.3');

%% Permutation test

% number of permuted weights
% optionally perform permutation test
if I.n_perms > 0
    
    [parent_directory,perm_MAT_file] = fileparts(MAT_file);
    perm_MAT_file = [parent_directory '/' perm_MAT_file '_' num2str(I.n_perms) 'perms.mat'];
    clear parent_directory;
    
    % estimate contrasts and residual from permutations
    beta_contrast_permtest = nan([I.n_perms, sum(nonzero_contrasts), sum(voxels_without_NaN)]);
    residual_permtest = nan([I.n_perms, sum(voxels_without_NaN)]);
    for i = 1:I.n_perms
        [beta_contrast_permtest(i,:,:), ~, ~, ~, residual_permtest(i,:)] = ...
            regress_stats_ols(psc(:,voxels_without_NaN), ...
            W(randperm(n_trials),:), C);
    end
    
    % convert to signed logP value
    logP_permtest = sig_via_null_gaussfit(...
        beta_contrast(nonzero_contrasts,voxels_without_NaN), ...
        beta_contrast_permtest);
    logP_residual_permtest = sig_via_null_gaussfit(...
        residual(:,voxels_without_NaN), residual_permtest, 'tail', 'left');
    
    % set NaN for zero contrasts
    logP_permtest = fillin_NaN(logP_permtest, nonzero_contrasts, 1);
    beta_contrast_permtest = fillin_NaN(beta_contrast_permtest, nonzero_contrasts, 2);
    
    % set NaN for voxels with NaN
    logP_permtest = fillin_NaN(logP_permtest, voxels_without_NaN, 2); %#ok<NASGU>
    beta_contrast_permtest = fillin_NaN(beta_contrast_permtest, voxels_without_NaN, 3); %#ok<NASGU>
    logP_residual_permtest = fillin_NaN(logP_residual_permtest, voxels_without_NaN, 2); %#ok<NASGU>
    residual_permtest = fillin_NaN(residual_permtest, voxels_without_NaN, 2); %#ok<NASGU>

    % save results to matfile
    save(perm_MAT_file, 'beta_contrast_permtest', 'logP_permtest',...
        'residual_permtest', 'logP_residual_permtest');
    
end
