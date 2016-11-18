function glm_second_level(data_matrix_files, para_files, ...
    MAT_file_second_level, MAT_files_first_level, varargin)

% 2016-08-26: Added residual statistics, Sam NH
%
% 2016-08-27: Modified how optional arguments are handled, current
% implementation requires the MAT files the results are saved to, to be
% specified, Sam NH
%
% 2016-09-09: Results of permutation test saved as a separate MAT file, Sam NH
% 
% 2016-09-09: No longer calls first level script, glm_event_regression.m, Sam NH
% 
% 2016-09-10: Made residual permutation test one-tailed, Sam NH
% 
% 2016-09-21: Modified to deal with NaN values, Sam NH

% handle optional inputs and defaults
I.n_perms = 0;
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
        'contrast_variance','P','df','residual','beta_one_per_regressor');
    
    if i == 1
        beta_contrast_allruns = nan([n_runs, size(X.beta_contrast)]);
        beta_one_per_regressor_allruns = nan([n_runs, size(X.beta_one_per_regressor)]);
        contrast_variance_allruns = nan([n_runs, size(X.contrast_variance)]);
        residual_allruns = nan([n_runs, size(X.residual)]);
        dfs = nan(n_runs,1);
        P = X.P;  %#ok<NASGU>
    end
    try
    beta_contrast_allruns(i,:,:) = X.beta_contrast;
    beta_one_per_regressor_allruns(i,:,:) = X.beta_one_per_regressor;
    contrast_variance_allruns(i,:,:) = X.contrast_variance;
    residual_allruns(i,:,:) = X.residual;
    dfs(i) = X.df;
    catch
        keyboard
    end
end
clear X;

% set beta contrast and residual to the average across runs
beta_contrast = squeeze_dims( nanmean(beta_contrast_allruns, 1), 1);
beta_one_per_regressor = squeeze_dims( nanmean(beta_one_per_regressor_allruns, 1), 1); %#ok<NASGU>
residual = squeeze_dims( mean(residual_allruns, 1), 1);

% ols stats across runs
[logP_fixed, contrast_variance_fixed] = ...
    fixed_effects(beta_contrast_allruns, contrast_variance_allruns, dfs); %#ok<ASGLU>
[logP_random, contrast_variance_random] = ...
    random_effects(beta_contrast_allruns); %#ok<ASGLU>

% save results
save(MAT_file_second_level, 'beta_contrast', 'logP_fixed', 'logP_random', ...
    'contrast_variance_fixed', 'contrast_variance_random', 'P', 'residual', ...
    'beta_one_per_regressor');

% clear variables no longer needed
clear beta_contrast_allruns beta_one_per_regressor_allruns beta_one_per_regressor ...
    contrast_variance_allruns logP_fixed contrast_variance_fixed logP_random ...
    contrast_variance_random;

%% Second level, permutation test

if I.n_perms == 0
    return;
end

% second level MAT files with permuted stats
[parent_directory,perm_MAT_file_second_level] = fileparts(MAT_file_second_level);
perm_MAT_file_second_level = [parent_directory '/' perm_MAT_file_second_level ...
    '_' num2str(I.n_perms) 'perms.mat'];
clear parent_directory;

% first level MAT files with permuted stats
perm_MAT_files_first_level = cell(size(MAT_files_first_level));
for i = 1:n_runs
    [parent_directory,perm_MAT_files_first_level{i}] = fileparts(MAT_files_first_level{i});
    perm_MAT_files_first_level{i} = [parent_directory '/' perm_MAT_files_first_level{i} ...
        '_' num2str(I.n_perms) 'perms.mat'];
    clear parent_directory;
end

% average across runs
for i = 1:n_runs
    X = load(perm_MAT_files_first_level{i}, ...
        'beta_contrast_permtest', 'residual_permtest');
    if i == 1
        beta_contrast_permtest = X.beta_contrast_permtest;
        counts_beta_contrast_permtest = double(~isnan(beta_contrast_permtest));
        residual_permtest = X.residual_permtest;
        counts_residual_permtest = double(~isnan(residual_permtest));
    else
        beta_contrast_permtest = ...
            nanplus( beta_contrast_permtest, X.beta_contrast_permtest );
        xi = ~isnan(X.beta_contrast_permtest);
        counts_beta_contrast_permtest(xi) = counts_beta_contrast_permtest(xi)+1;
        clear xi;
        
        residual_permtest = ...
            nanplus(residual_permtest, X.residual_permtest);
        xi = ~isnan(X.residual_permtest);
        counts_residual_permtest(xi) = counts_residual_permtest(xi)+1;
        clear xi;
    end
end
beta_contrast_permtest = beta_contrast_permtest ./ counts_beta_contrast_permtest;
residual_permtest = residual_permtest ./ counts_residual_permtest;

% estimate P value
load(MAT_file_second_level, 'beta_contrast', 'residual');
logP_permtest = sig_via_null_gaussfit(beta_contrast, beta_contrast_permtest); %#ok<NASGU>
logP_residual_permtest = sig_via_null_gaussfit(residual, residual_permtest, ...
    'tail', 'left'); %#ok<NASGU>

% save results
save(perm_MAT_file_second_level, 'logP_permtest', 'beta_contrast_permtest', ...
    'residual_permtest', 'logP_residual_permtest');

function Z = nanplus(X,Y)

assert(all(size(X)==size(Y)))
Z = nan(size(X));

xi = ~isnan(X) & ~isnan(Y);
Z(xi) = X(xi) + Y(xi);

xi = ~isnan(X) & isnan(Y);
Z(xi) = X(xi);

xi = isnan(X) & ~isnan(Y);
Z(xi) = Y(xi);

