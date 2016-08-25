function matfile = glm_event_regression(data_matrix_file, para_file, ...
    parameter_file, n_perms, matfile, nuissance_regressor_file, varargin)

% Regression analysis with discrete events/blocks. Regressors are weighted
% events, and the beta weights from the regression analysis are multiplied by a
% contrast vector/matrix. Statistics are computed using ordinary least squares
% and a permutation test.
% 
% 2016-07-08: Generalized to have contrasts.

%% Outside directories

% global root_directory
% if ~exist([root_directory '/general-analysis-code'], 'dir')
%     error('general-analysis-code not found in root_directory');
% end
% addpath([root_directory '/general-analysis-code']);

% default no permutation tests
if nargin < 4
    n_perms = 0;
end

% file to save results to
if nargin < 5
    matfile = [pwd '/' DataHash([data_matrix_file, para_file, ...
        parameter_file, n_perms])];
end

if nargin < 6
    nuissance_regressor_file = [];
end

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
% -> time x condition
boxcar_sr = 1;
[B, event_names] = boxcar_from_para(para_file, boxcar_sr);

% convolve boxcars with hrf
% -> time x condition
B = convolve_with_hrf(B, 'fsfast-gamma-BOLD', boxcar_sr);

% interpolate condition matrix to the time-points at which
% the data were collected
% -> time x condition
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

% weights with one regressor per condition
n_conditions = length(P.condition_names);
W_one_per_condition = zeros(n_events,n_conditions);
for i = 1:n_events
    xi = strcmp(event_names{i}, P.condition_names);
    W_one_per_condition(i,xi) = 1;
end

%% Nuissance regressors

% load nuissance regressors from input file
if ~isempty(nuissance_regressor_file)
    load(nuissance_regressor_file, 'X_nuissance');
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
[beta_contrast, logP_ols, contrast_variance, df] = ...
    regress_stats_ols( Y, [B * W, X_nuissance], ...
    [P.contrast_weights; zeros(n_nuissance, size(P.contrast_weights,2))]); %#ok<ASGLU>

% separate beta weight for regressor
% redundant with previous analysis if contrast matrix is the identity
beta_one_per_regressor = ...
    regress_stats_ols(Y, [B * W, X_nuissance], ...
    [eye(n_regressors); zeros(n_nuissance, n_regressors)]);

% separate beta weight for each condition
beta_one_per_condition = ...
    regress_stats_ols(Y, [B * W_one_per_condition, X_nuissance], ...
    [eye(n_conditions); zeros(n_nuissance, n_conditions)]);

% fill in NaN entries
beta_contrast = fillin_NaN_voxels(beta_contrast, voxels_without_NaN, 2);
logP_ols = fillin_NaN_voxels(logP_ols, voxels_without_NaN, 2); %#ok<NASGU>
contrast_variance = fillin_NaN_voxels(contrast_variance, voxels_without_NaN, 2); %#ok<NASGU>
beta_one_per_condition = fillin_NaN_voxels(beta_one_per_condition, voxels_without_NaN, 2); %#ok<NASGU>
beta_one_per_regressor = fillin_NaN_voxels(beta_one_per_regressor, voxels_without_NaN, 2); %#ok<NASGU>

% save
save(matfile, 'voxels_without_NaN', 'beta_contrast', ...
    'logP_ols', 'contrast_variance', 'beta_one_per_condition', ...
    'beta_one_per_regressor', 'df', 'P', '-v7.3');

%% Permutation test

% number of permuted weights
% optionally perform permutation test
if n_perms > 0
    
    % can exclude null so that stimulus conditions only permuted
    % with respect to other stimulus conditions
    if optInputs(varargin, 'exclude-null')
        xi = ~strcmp('NULL', event_names);
        B = B(:,xi);
        W = W(xi,:);
        n_events = sum(xi);
    end
    
    % estimate contrast from permutations
    beta_contrast_permtest = nan([n_perms, length(P.contrast_names), sum(voxels_without_NaN)]);
    for i = 1:n_perms
        X = B * W(randperm(n_events),:);
        beta_contrast_permtest(i,:,:) = regress_stats_ols(Y, X, P.contrast_weights);
    end
    clear X;
    
    % convert to signed logP value
    try
        logP_permtest = sig_via_null_gaussfit(...
            beta_contrast(:,voxels_without_NaN), beta_contrast_permtest);
    catch
        keyboard;
    end
    
    % fill in NaN entries
    beta_contrast_permtest = fillin_NaN_voxels(...
        beta_contrast_permtest, voxels_without_NaN, 3); %#ok<NASGU>
    logP_permtest = fillin_NaN_voxels(logP_permtest, voxels_without_NaN, 2); %#ok<NASGU>
        
    % save results to matfile
    save(matfile, 'beta_contrast_permtest', 'logP_permtest', '-append');
    
end

% v2



