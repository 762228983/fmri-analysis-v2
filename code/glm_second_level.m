function [matfile_second_level, matfile_first_level] = ...
    glm_second_level(data_matrix_files, para_files, ...
    parameter_file, n_perms, output_directory, varargin)

% global root_directory;
% if ~exist([root_directory '/general-analysis-code'], 'dir')
%     error('general-analysis-code not found in root_directory');
% end
% addpath([root_directory '/general-analysis-code']);

% optional arguments
if nargin < 4
    n_perms = 0;
end
if nargin < 5
    output_directory = pwd;
    fprintf('No output directory specified\nSaving to working directory\m')
end

% create output directory if not present
if ~exist(output_directory, 'dir')
    mkdir(output_directory);
end

%% First level / run-based analysis

% analyzing individual runs
n_runs = length(data_matrix_files);
matfile_first_level = cell(1,n_runs);
for i = 1:n_runs
    
    fprintf('Analyzing run %d\n',i);
    
    % first level analysis
    matfile_first_level{i} = ...
        [output_directory '/r' num2str(i) '_' num2str(n_perms) 'perms.mat'];
    
    % check if the output file already exists, if not perform analysis
    if ~exist(matfile_first_level{i}, 'file') ...
            || optInputs(varargin, 'overwrite') 
        
        % first level regression
        glm_event_regression(data_matrix_files{i}, para_files{i}, ...
            parameter_file, n_perms, matfile_first_level{i} );
    end
end

if n_runs == 1
    matfile_second_level = [];
    return;
end

%% Second level, ols stats

matfile_second_level = [output_directory '/allruns_' num2str(n_perms) 'perms.mat'];

% check if the output file already exists
if ~exist(matfile_second_level, 'file') || optInputs(varargin, 'overwrite') 
        
    % load all runs
    % beta_contrast_allruns: runs x contrast x voxel
    for i = 1:n_runs
        
        % load first level analysis
        X = load(matfile_first_level{i},'beta_contrast','contrast_variance','P','df');
        
        if i == 1
            beta_contrast_allruns = nan([n_runs, size(X.beta_contrast)]);
            contrast_variance_allruns = nan([n_runs, size(X.contrast_variance)]);
            dfs = nan(n_runs,1);
            P = X.P;
        end
        
        beta_contrast_allruns(i,:,:) = X.beta_contrast;
        contrast_variance_allruns(i,:,:) = X.contrast_variance;
        dfs(i) = X.df;
        
    end
    clear X;
    
    % set beta contrast to the average across runs
    beta_contrast = squeeze_dims(mean(beta_contrast_allruns),1);
    
    % ols stats across runs
    [logP_fixed, contrast_variance_fixed] = ...
        fixed_effects(beta_contrast_allruns, contrast_variance_allruns, dfs); %#ok<ASGLU>
    [logP_random, contrast_variance_random] = ...
        random_effects(beta_contrast_allruns); %#ok<ASGLU>
    
    % save results
    save(matfile_second_level, 'beta_contrast', 'logP_fixed', 'logP_random', ...
        'contrast_variance_fixed', 'contrast_variance_random', 'P');
    
    % clear variables no longer needed
    clear beta_contrast_allruns contrast_variance_allruns ...
        logP_fixed contrast_variance_fixed logP_random contrast_variance_random;

end

%% Second level, permutation test

if n_perms == 0
    return;
end
matfile_vars = whos('-file', matfile_second_level);
varnames = {matfile_vars(:).name};

if ~any(strcmp('logP_permtest',varnames)) ...
        || ~any(strcmp('beta_contrast_permtest',varnames))

    % average across runs
    for i = 1:n_runs
        X = load(matfile_first_level{i}, 'beta_contrast_permtest');
        if i == 1
            beta_contrast_permtest = X.beta_contrast_permtest / n_runs;
        else
            beta_contrast_permtest = beta_contrast_permtest ...
                + X.beta_contrast_permtest / n_runs;
        end
    end
        
    % estimate P value
    load(matfile_second_level, 'beta_contrast');
    logP_permtest = sig_via_null_gaussfit(beta_contrast, beta_contrast_permtest); %#ok<NASGU>
    
    % save results
    save(matfile_second_level, 'logP_permtest', 'beta_contrast_permtest', '-append');
    
end

