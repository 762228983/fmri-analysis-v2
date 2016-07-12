function [weight_file_runaverage, weight_file_individual_runs] = ...
    component_voxel_weights(data_files, para_files, condition_weights_file, output_directory, n_perms, varargin)

if ~exist(output_directory, 'dir')
    mkdir(output_directory);
end

% analyzing individual runs
n_runs = length(data_files);
weight_file_individual_runs = cell(1,n_runs);
for i = 1:n_runs
    
    fprintf('Analyzing run %d\n',i);
    
    % name of the output file
    weight_file_individual_runs{i} = ...
        [output_directory '/r' num2str(i) '_nperms' num2str(n_perms) '.mat'];
    
    if ~exist(weight_file_individual_runs{i}, 'file') ...
            || optInputs(varargin, 'overwrite') % check if the output file already exists
        % perform the analysis if the output file doesn't already exist
        data_file = data_files{i};
        para_file = para_files{i};
        [component_weights, component_names, permuted_weights, logP] = ...
            glm_voxel_timecourses_weighted_boxcar(...
            data_file, para_file, condition_weights_file, ...
            'permutation-test', n_perms);
        
        save(weight_file_individual_runs{i}, 'component_weights', 'component_names', ...
            'permuted_weights', 'logP', 'data_file', 'para_file', '-v7.3');
    end
end

% clear temporary variables
clear component_weights component_names permuted_weights logP;

if n_runs == 1
    return;
end

% stats across runs
weight_file_runaverage = [output_directory '/all-runs_nperms' num2str(n_perms)  '.mat'];
if ~exist(weight_file_runaverage, 'file') || optInputs(varargin, 'overwrite') % check if the output file already exists
    
    % average across runs
    for i = 1:n_runs
        load(weight_file_individual_runs{i}, 'component_weights', 'component_names', 'permuted_weights');
        if i == 1
            component_weights_runaverage = zeros(size(component_weights));
            permuted_weights_runaverage = zeros(size(permuted_weights));
        end
        component_weights_runaverage =  component_weights_runaverage + component_weights / n_runs;
        permuted_weights_runaverage = permuted_weights_runaverage + permuted_weights / n_runs;
    end
    
    % compute average and standard deviation across permuted samples
    null_mean = squeeze(mean(permuted_weights_runaverage,1));
    null_std = squeeze(std(permuted_weights_runaverage,[],1));
    
    % use null mean and standard deviation to convert to a z-statistics
    component_weights_z = (component_weights_runaverage - null_mean) ./ null_std;
    
    % conver z to -log10[p]
    logP = -sign(component_weights_z) .* log10(2*normcdf(-abs(component_weights_z), 0, 1));
    
    % save
    component_weights = component_weights_runaverage;
    permuted_weights = permuted_weights_runaverage;
    save(weight_file_runaverage, 'component_weights', 'component_names', 'permuted_weights', 'logP', 'data_files', 'para_files');
    
end

