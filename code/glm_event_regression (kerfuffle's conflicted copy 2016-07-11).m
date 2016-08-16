function matfile = glm_event_regression(data_matrix_file, para_file, ...
    parameter_file, n_perms, matfile, varargin)

% Regression analysis with discrete events/blocks. Regressors are weighted
% events, and the beta weights from the regression analysis are multiplied by a
% contrast vector/matrix. Statistics are computed using ordinary least squares
% and a permutation test.
% 
% 2016-07-08: Generalized to have contrasts.

keyboard;
%% Outside directories

global root_directory
if ~exist([root_directory '/general-analysis-code'], 'dir')
    error('general-analysis-code not found in root_directory');
end
addpath([root_directory '/general-analysis-code'])

% default no permutation tests
if nargin < 4
    n_perms = 0;
end

% file to save results to
if nargin < 5
    matfile = [pwd '/' DataHash([data_matrix_file, para_file, ...
        parameter_file, n_perms])];
end

    
%% Format data matrix

% load data matrix
% -> time x voxel
load(data_matrix_file, 'data_matrix', 'TR');
Y = data_matrix;
clear data_matrix;

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
[C, condition_names] = boxcar_from_para(para_file, boxcar_sr);

% convolve boxcars with hrf
% -> time x condition
C = convolve_with_hrf(C, 'fsfast-gamma-BOLD');

% interpolate condition matrix to the time-points at which
% the data were collected
% -> time x condition
C = interp1( (0:length(C)-1)/boxcar_sr, C, (0:size(Y,1)-1)*TR );

%% Regressor weights for condition matrix

% weights applied to each condition for each regressor
% -> condition x regressor
P = load(parameter_file);
n_conditions = size(C,2);
n_regressors = length(P.regressor_names);
W = zeros(n_conditions,n_regressors);
for i = 1:n_conditions
    xi = strcmp(condition_names{i}, P.condition_names);
    if ~isempty(xi)
        W(i,:) = P.regressor_weights(xi,:);
    end
end

%% Regression

[beta_contrast, logP_ols, contrast_variance] = ...
    regress_stats_ols(Y, C * W, P.contrast_weights);
save(matfile, 'beta_contrast', 'logP_ols', 'contrast_variance', '-v7.3');

%% Permutation test

% number of permuted weights
% optionally perform permutation test
if n_perms > 0
    
    % can exclude null so that stimulus conditions only permuted
    % with respect to other stimulus conditions
    if optInputs(varargin, 'exclude-null')
        xi = ~strcmp('NULL', condition_names);
        C = C(:,xi);
        W = W(xi,:);
        n_conditions = sum(xi);
    end
    
    % estimate contrast from permutations
    beta_contrast_permtest = nan([n_perms, size(beta_contrast)]);
    for i = 1:n_perms
        X = C * W(randperm(n_conditions),:);
        beta_contrast_permtest(i,:,:) = regress_stats_ols(Y, X, P.contrast_weights);
    end
    
    % convert to signed logP value
    logP_permtest = ...
        sig_via_null_gaussfit(beta_contrast, beta_contrast_permtest);
    
    % save results to matfile
    save(matfile, 'beta_contrast_permtest', 'logP_permtest', '-v7.3', '-append');
    
end

% v1



