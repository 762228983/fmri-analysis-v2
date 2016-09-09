function sigav_glm_second_level(...
    data_matrix_files, para_files, parameter_file, ...
    MAT_file_second_level, MAT_files_first_level, varargin)

% Combines results from applying sigav_glm.m to data from several different
% runs.
% 
% 2016-09-09 - Last Modified, Sam NH

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
        sigav_glm(data_matrix_files{i}, para_files{i}, ...
            parameter_file, MAT_files_first_level{i}, ...
            'onset_delay', I.onset_delay, 'offset_delay', I.offset_delay);
    end
    
end

if n_runs == 1
    return;
end

% load results form individual runs
if exist(MAT_file_second_level, 'file') && ~I.overwrite
    return;
end

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
matfile_vars = whos('-file', MAT_file_second_level);
varnames = {matfile_vars(:).name};

if ~any(strcmp('logP_permtest',varnames)) ...
        || ~any(strcmp('beta_contrast_permtest',varnames)) ...
        || ~any(strcmp('residual_permtest',varnames)) ...
        || ~any(strcmp('logP_residual_permtest',varnames))
    
    % average across runs
    for i = 1:n_runs
        X = load(MAT_files_first_level{i}, ...
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
    save(MAT_file_second_level, 'logP_permtest', 'beta_contrast_permtest', ...
        'residual_permtest', 'logP_residual_permtest', '-append');
    
end
