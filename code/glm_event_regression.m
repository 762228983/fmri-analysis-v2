function glm_event_regression(...
    data_matrix_file, para_file, parameter_file, MAT_file, varargin)

% Regression analysis with discrete events/blocks. Regressors are weighted
% events, and the beta weights from the regression analysis are multiplied by a
% contrast vector/matrix. Statistics are computed using ordinary least squares
% and a permutation test.
% 
% -- Variables saved to MAT_file --
% 
% beta_one_per_regressor: [NUM_REGRESSORS x NUM_VOXELS] matrix. Each row
% contains the beta weights for a single regressor across all voxels. NaN values
% indicate voxels that did not have values in the input data matrix.
% 
% beta_contrast: [NUM_CONTRASTS x NUM_VOXELS] matrix. Each row contains the
% values for a single contrast between the betas (i.e. a single column of the
% contrast_weights matrix). Equivalent to beta_one_per_regressor if the
% contrast_weight matrix was the identity matrix.
% 
% logP_ols: [NUM_CONTRASTS x NUM_VOXELS] matrix. Each row contains a
% significance measure (p-value) for a single contrast based on ordinary
% least-squares assumptions (independent, Gaussian-distributed errors). P-values
% are logarithmically transformed and signed based on whether the voxel
% responded positively or negatively to the contrast: sign(beta_contrast) .*
% -log10(p).
% 
% -- Variables saved to permutation file --
% 
% logP_permtest: [NUM_CONTRASTS x NUM_VOXELS] matrix. Each row contains a
% significance measure (p-value) for a single contrast based on a null
% distribution, computed by permuting the order of conditions in a run. P-values
% are logarithmically transformed and signed based on whether the beta contrast
% exceeded (positive) or fellow below (negative) the mean of the null
% distribution.
% 
% logP_residual_permtest: [1 x NUM_VOXELS] vector. A single row with a
% significance measure that evaluates the overall significance of the regression
% model. P-values are logarithmically transformed (-log10(p)) and are stricly
% positive.
% 
% 
% 2016-07-08: Generalized to have contrasts, Sam NH
% 
% 2016-08-26: Add residual statistics, Sam NH
% 
% 2016-08-27: Modified how optional arguments are handled, current
% implementation requires the MAT files the results are saved to, to be
% specified
% 
% 2016-09-09: Results are of permutation tests are saved as a separate MAT file
% 
% 2016-09-10: Made residual permutation test one-tailed
% 
% 2016-09-21: Modified to deal with zero regressors and contrasts, Sam NH


%% Outside directories

% handle optional inputs and defaults
I.n_perms = 0;
I.nuissance_regressor_file = [];
I = parse_optInputs_keyvalue(varargin, I);

% default glm parameters
P = load(parameter_file);
P = glm_default_parameters(P);

%% Format data matrix

% load data matrix
% -> time x voxel
load(data_matrix_file, 'data_matrix', 'TR');
Y = data_matrix;
clear data_matrix;

% remove voxels with NaN values
voxels_without_NaN = all(~isnan(Y));
Y = Y(:,voxels_without_NaN);

% re-scale voxels to have mean 100
% controls for errant scaling and ensures
% beta weights are by default in units of
% percent signal change
% -> time x voxel
Y = 100 * Y ./ ...
    repmat(mean(Y), [size(Y,1),1]);

% demean data matrix
% -> time x voxel
Y = Y - repmat(mean(Y), [size(Y,1),1]);

%% Condition matrix

% collection of "boxcars" with ones during times when a given condition was "on"
% -> time x event
boxcar_sr = 1;
[B, event_names] = boxcar_from_para(para_file, boxcar_sr);

% convolve boxcars with hrf
% -> time x event
B = convolve_with_hrf(B, 'fsfast-gamma-BOLD', boxcar_sr);

% interpolate condition matrix to the time-points at which
% the data were collected
% -> time x event
B = interp1( (0:length(B)-1)/boxcar_sr, B, (0:size(Y,1)-1)*TR );

%% Regressor weights for condition matrix

% weights applied to each condition for each regressor
% -> condition x regressor
n_events = size(B,2);
n_regressors = length(P.regressor_names);
W = zeros(n_events,n_regressors);
for i = 1:n_events
    xi = strcmp(event_names{i}, P.condition_names);
    if any(xi)
        W(i,:) = P.regressor_weights(xi,:);
    end
end
nonzero_regressors = any(W~=0,1);
W = W(:,nonzero_regressors);

% weights with one regressor per condition
n_conditions = length(P.condition_names);
W_one_per_condition = zeros(n_events,n_conditions);
for i = 1:n_events
    xi = strcmp(event_names{i}, P.condition_names);
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

%% Nuissance regressors

% load nuissance regressors from input file
if ~isempty(I.nuissance_regressor_file)
    load(I.nuissance_regressor_file, 'X_nuissance');
else
    X_nuissance = [];
end

% add linear trend regressor as a nuissance regressor
if P.linear_trend
    X_nuissance = [X_nuissance, zscore((1:size(Y,1))')];
end
n_nuissance = size(X_nuissance,2);

%% Regression

% regress weighted event matrix
if n_contrasts > 0
    [beta_contrast, logP_ols, contrast_variance, df, residual] = ...
        regress_stats_ols( Y, [B * W, X_nuissance], ...
        [C; zeros(n_nuissance, size(C,2))]); %#ok<ASGLU>
end

% separate beta weight for regressor
% redundant with previous analysis if contrast matrix is the identity
beta_one_per_regressor = ...
    regress_stats_ols(Y, [B * W, X_nuissance], ...
    [eye(size(W,2)); zeros(n_nuissance, size(W,2))]);

% separate beta weight for each condition
beta_one_per_condition = ...
    regress_stats_ols(Y, [B * W_one_per_condition, X_nuissance], ...
    [eye(size(W_one_per_condition,2)); zeros(n_nuissance, size(W_one_per_condition,2))]);

%% Save

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

save(MAT_file, 'voxels_without_NaN', 'beta_contrast', ...
    'logP_ols', 'contrast_variance', 'beta_one_per_condition', ...
    'beta_one_per_regressor', 'df', 'P', 'residual', '-v7.3');

%% Permutation test

% number of permuted weights
% optionally perform permutation test
if I.n_perms > 0
    
    [parent_directory,perm_MAT_file] = fileparts(MAT_file);
    perm_MAT_file = [parent_directory '/' perm_MAT_file '_' num2str(I.n_perms) 'perms.mat'];
    clear parent_directory;
    
    % can exclude null so that stimulus conditions only permuted
    % with respect to other stimulus conditions
    %     if optInputs(varargin, 'exclude-null')
    %         xi = ~strcmp('NULL', event_names);
    %         B = B(:,xi);
    %         W = W(xi,:);
    %         n_events = sum(xi);
    %     end
    
    % estimate contrasts and residusl from permutations
    beta_contrast_permtest = nan([I.n_perms, sum(nonzero_contrasts), sum(voxels_without_NaN)]);
    residual_permtest = nan([I.n_perms, sum(voxels_without_NaN)]);
    for i = 1:I.n_perms
        X = B * W(randperm(n_events), :);
        [beta_contrast_permtest(i,:,:), ~, ~, ~, residual_permtest(i,:)] ...
            = regress_stats_ols(Y, [X, X_nuissance], ...
            [C; zeros(n_nuissance, size(C,2))]);
    end
    clear X;
    
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


