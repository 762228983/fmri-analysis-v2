function [matfile_second_level, matfile_first_level] = ...
    glm_second_level(data_matrix_files, para_files, parameter_file, varargin)

% 2016-08-26: Added residual statistics, Sam NH
% 
% 2016-08-27: Modified how optional arguments are handled

I.n_perms = 0;
I.output_directory = pwd;
I.nuissance_regressor_files = cell(size(data_matrix_files));
I.overwrite = false;
I = parse_optInputs_keyvalue(varargin, I);

% create output directory if not present
if ~exist(I.output_directory, 'dir')
    mkdir(I.output_directory);
end

%% First level / run-based analysis

% analyzing individual runs
n_runs = length(data_matrix_files);
matfile_first_level = cell(1,n_runs);
for i = 1:n_runs
    
    fprintf('Analyzing run %d\n',i);
    
    %     % hash value specific to the inputs of glm_event_regression
    %     args = {...
    %         data_matrix_files{i}, para_files{i}, ...
    %         parameter_file, I.n_perms, matfile_first_level{i}};
    %     hash_value = DataHash(args);
    
    % first level analysis
    matfile_first_level{i} = ...
        [I.output_directory '/r' num2str(i) ...
        '_' num2str(I.n_perms) 'perms.mat'];
    
    % check if the output file already exists, if not perform analysis
    if ~exist(matfile_first_level{i}, 'file') ...
            || I.overwrite
       
        % first level regression
        glm_event_regression(data_matrix_files{i}, para_files{i}, ...
            parameter_file, I.n_perms, matfile_first_level{i}, ...
            I.nuissance_regressor_files{i});
        
    end
end

if n_runs == 1
    matfile_second_level = [];
    return;
end

%% Second level, ols stats

matfile_second_level = [I.output_directory '/allruns_' num2str(I.n_perms) 'perms.mat'];

% check if the output file already exists
if ~exist(matfile_second_level, 'file') || I.overwrite
    
    % load all runs
    % beta_contrast_allruns: runs x contrast x voxel
    for i = 1:n_runs
        
        % load first level analysis
        X = load(matfile_first_level{i},'beta_contrast',...
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
    save(matfile_second_level, 'beta_contrast', 'logP_fixed', 'logP_random', ...
        'contrast_variance_fixed', 'contrast_variance_random', 'P', 'residual');
    
    % clear variables no longer needed
    clear beta_contrast_allruns contrast_variance_allruns ...
        logP_fixed contrast_variance_fixed logP_random contrast_variance_random;
    
end

%% Second level, permutation test

if I.n_perms == 0
    return;
end
matfile_vars = whos('-file', matfile_second_level);
varnames = {matfile_vars(:).name};

if ~any(strcmp('logP_permtest',varnames)) ...
        || ~any(strcmp('beta_contrast_permtest',varnames)) ...
        || ~any(strcmp('residual_permtest',varnames)) ...
        || ~any(strcmp('logP_residual_permtest',varnames))
    
    % average across runs
    for i = 1:n_runs
        X = load(matfile_first_level{i}, ...
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
    load(matfile_second_level, 'beta_contrast', 'residual');
    logP_permtest = sig_via_null_gaussfit(beta_contrast, beta_contrast_permtest); %#ok<NASGU>
    logP_residual_permtest = -sig_via_null_gaussfit(residual, residual_permtest); %#ok<NASGU>
    
    % save results
    save(matfile_second_level, 'logP_permtest', 'beta_contrast_permtest', ...
        'residual_permtest', 'logP_residual_permtest', '-append');
    
end

