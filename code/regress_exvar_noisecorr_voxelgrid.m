function [Gcorr, corr_MAT_file] = regress_exvar_noisecorr_voxelgrid(...
    G, F, output_directory, varargin)

% Takes a feature matrix F (stim x feature) and voxel data structure G in the form
% of a grid, and computes predictions for each voxel's response
% 
% Alternatively, F and G can be mat files which contain the relevant
% matrix/structure. The variables in each file are assumed to named 'F' and 'D',
% such that load(F, 'F') and load(G, 'G') read in the appropriate variables. 
% 
% MAT_file is an optional MAT file that the results can be saved to.  

if ischar(F)
    load(F, 'F');
end

if ischar(G)
    load(G, 'G');
end

I.test_folds = 2;
I.train_folds = min(max(round((size(F,1) / I.test_folds) / 10),2),10);
I.regression_method = 'ridge';
I.K = [];
I.use_sbatch = false;
I.overwrite_predictions = false;
I.overwrite_correlations = false;
I = parse_optInputs_keyvalue(varargin, I);

% predictions
prediction_MAT_file = [output_directory '/predictions.mat'];
if ~exist(prediction_MAT_file, 'file') || I.overwrite_predictions
    
    % stimuli x rep x voxel
    D = grid2matrix(G);
    [n_stim, n_rep, n_vox] = size(D); %#ok<ASGLU>
    
    % remove zero or NaN voxels
    nonempty_voxels = ~(all(all(D < 1e-100,1),2) | all(all(isnan(D),1),2));
    D = D(:,:,nonempty_voxels);
    n_nonempty_vox = size(D,3);
    
    % stim x (rep * voxel)
    D = reshape(D, n_stim, n_rep  * n_nonempty_vox);
    
    if I.use_sbatch
        
        % directory to save results for each voxel
        directory_to_save_results_for_each_voxel = ...
            [output_directory '/predictions_individual_voxels'];
        if ~exist(directory_to_save_results_for_each_voxel, 'dir')
            mkdir(directory_to_save_results_for_each_voxel);
        end
        [Yh, folds] = regress_predictions_parallelize_with_slurm(...
            F, D, I.test_folds, I.regression_method, I.K, I.train_folds, ...
            directory_to_save_results_for_each_voxel);
        
    else
        
        Yh = nan(size(D));
        for i = 1:size(D,2)
            fprintf('Voxel %d\n',i); drawnow;
            [Yh(:,i), ~, ~, folds] = regress_predictions_from_3way_crossval(F, D(:,i), ...
                I.test_folds, I.regression_method, I.K, I.train_folds);
        end
        
    end
    
    save(prediction_MAT_file, 'Yh', 'D', 'folds', 'nonempty_voxels',...
        'n_stim', 'n_rep', 'n_vox', 'n_nonempty_vox');
    
end

% noise-corrected correlations
corr_MAT_file = [output_directory '/noise_corrected_correlations.mat'];
if ~exist(corr_MAT_file, 'file') || I.overwrite_correlations
    
    % load variables from previous step
    load(prediction_MAT_file, 'Yh', 'D', 'folds', 'nonempty_voxels', ...
        'n_stim', 'n_rep', 'n_vox', 'n_nonempty_vox');
    
    % separate out reps / voxels into different dimensions
    D = reshape(D, [n_stim, n_rep, n_nonempty_vox]);
    Yh = reshape(Yh, [n_stim, n_rep, n_nonempty_vox]);
    
    % compute correlation
    r = nan(n_nonempty_vox,1);
    for i = 1:size(D,3);
        r(i) = normalized_correlation_within_folds(D(:,:,i), Yh(:,:,i), folds); 
    end
    
    % add back in empty voxels as NaNs
    r = fillin_NaN(r, nonempty_voxels, 1);
    
    save(corr_MAT_file, 'r');
    
else
    
    load(corr_MAT_file, 'r');
    
end

% conver to grid
Gcorr = G;
Gcorr.grid_data{1} = nan(size(Gcorr.grid_data{1},1), size(Gcorr.grid_data{1},2));
Gcorr.grid_data{2} = nan(size(Gcorr.grid_data{2},1), size(Gcorr.grid_data{2},2));
Gcorr = matrix2grid(r, Gcorr);
save(corr_MAT_file, 'Gcorr');

