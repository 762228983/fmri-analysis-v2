function glm_second_level(data_matrix_files, para_files, ...
    MAT_file_second_level, MAT_files_first_level, varargin)

% 2016-08-26: Added residual statistics, Sam NH
%
% 2016-08-27: Modified how optional arguments are handled, current
% implementation requires the MAT files the results are saved to, to be
% specified
%
% 2016-09-09: Results of permutation test saved as a separate MAT file
% 
% 2016-09-09: No longer calls first level script, glm_event_regression.m

% handle optional inputs and defaults
I.n_perms = 0;
I.nuissance_regressor_files = cell(size(data_matrix_files));
I = parse_optInputs_keyvalue(varargin, I);

% number of total runs
n_runs = length(data_matrix_files);
assert(n_runs == length(para_files));
assert(n_runs == length(MAT_files_first_level));

%% Second level, ols stats

% load all runs
% beta_contrast_allruns: runs x contrast x voxel
for i = 1:n_runs
    
    % load first level analysis
    X = load(MAT_files_first_level{i},'beta_contrast',...
        'contrast_variance','P','df','residual');
    
    if i == 1
        beta_contrast_allruns = nan([n_runs, size(X.beta_contrast)]);
        contrast_variance_allruns = nan([n_runs, size(X.contrast_variance)]);
        residual_allruns = nan([n_runs, size(X.residual)]);
        dfs = nan(n_runs,1);
        P = X.P;  %#ok<NASGU>
    end
    
    beta_contrast_allruns(i,:,:) = X.beta_contrast;
    contrast_variance_allruns(i,:,:) = X.contrast_variance;
    residual_allruns(i,:,:) = X.residual;
    dfs(i) = X.df;
    
end
clear X;

% set beta contrast and residual to the average across runs
beta_contrast = squeeze_dims( mean(beta_contrast_allruns, 1), 1);
residual = squeeze_dims( mean(residual_allruns, 1), 1);

% ols stats across runs
[logP_fixed, contrast_variance_fixed] = ...
    fixed_effects(beta_contrast_allruns, contrast_variance_allruns, dfs); %#ok<ASGLU>
[logP_random, contrast_variance_random] = ...
    random_effects(beta_contrast_allruns); %#ok<ASGLU>

% save results
save(MAT_file_second_level, 'beta_contrast', 'logP_fixed', 'logP_random', ...
    'contrast_variance_fixed', 'contrast_variance_random', 'P', 'residual');

% clear variables no longer needed
clear beta_contrast_allruns contrast_variance_allruns ...
    logP_fixed contrast_variance_fixed logP_random contrast_variance_random;

%% Second level, permutation test

if I.n_perms == 0
    return;
end

% second level MAT files with permuted stats
[~,perm_MAT_file_second_level] = fileparts(MAT_file_second_level);
perm_MAT_file_second_level = [perm_MAT_file_second_level ...
    '_' num2str(I.n_perms) 'perms.mat'];

% first level MAT files with permuted stats
perm_MAT_files_first_level = cell(size(MAT_files_first_level));
for i = 1:n_runs
    [~,perm_MAT_files_first_level{i}] = fileparts(MAT_files_first_level{i});
    perm_MAT_files_first_level{i} = [MAT_files_first_level{i} ...
        '_' num2str(I.n_perms) 'perms.mat'];
end

% average across runs
for i = 1:n_runs
    X = load(perm_MAT_files_first_level{i}, ...
        'beta_contrast_permtest', 'residual_permtest');
    if i == 1
        beta_contrast_permtest = X.beta_contrast_permtest / n_runs;
        residual_permtest = X.residual_permtest / n_runs;
    else
        beta_contrast_permtest = beta_contrast_permtest ...
            + X.beta_contrast_permtest / n_runs;
        residual_permtest = residual_permtest + X.residual_permtest / n_runs;
    end
end

% estimate P value
load(MAT_file_second_level, 'beta_contrast', 'residual');
logP_permtest = sig_via_null_gaussfit(beta_contrast, beta_contrast_permtest); %#ok<NASGU>
logP_residual_permtest = -sig_via_null_gaussfit(residual, residual_permtest); %#ok<NASGU>

% save results
save(perm_MAT_file_second_level, 'logP_permtest', 'beta_contrast_permtest', ...
    'residual_permtest', 'logP_residual_permtest', '-append');


