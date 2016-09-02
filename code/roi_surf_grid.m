function [psc, conditions, n_voxels_per_run_and_threshold] = ...
    roi_surf_grid(  us, grid_roi, grid_spacing_mm, ...
    localizer_info, test_info, fwhm, varargin )

% Primary top-level script for performing ROI analyses
%
% fprintf('v2\n'); drawnow;

% % -- Example Arguments --
% % subject
% us = 158;
%
% % roi constraint region and spacing of the grid
% grid_roi = 'hand-audctx';
% grid_spacing_mm = 2.8571/2;
%
% % amount the data was smoothed
% fwhm = 2.8571;
%
% % information about the test data
% clear test_info;
% test_info.exp = 'pitch_localizer_monkey';
% test_info.runtype = 'pitchloc2_combined_split';
%
% % information about the localizer data/contrasts
% clear localizer_info;
% localizer_info = [test_info, test_info];
%
% localizer_info(1).contrast = 'all_pitchloc2';
% localizer_info(1).contrast_sign = 1;
% localizer_info(1).threshold_type = 'absolute';
% localizer_info(1).thresholds = 3;
% localizer_info(1).max_runs_to_use = [];
%
% localizer_info(2).contrast = 'f0100-200_vs_f0800-1600';
% localizer_info(2).contrast_sign = 1;
% localizer_info(2).threshold_type = 'relative';
% localizer_info(2).thresholds = 0.05:0.05:0.3;
% localizer_info(2).max_runs_to_use = [];
%
% 2016-08-26: Created, Sam NH

I.verbose = true;
I.anatomical_mask = '';
I = parse_optInputs_keyvalue(varargin, I);

% default parameters
test_info = ...
    default_test_parameters(test_info, us);
localizer_info = ...
    default_localizer_parameters(localizer_info, us, test_info);

% number of thresholds for each localizer
n_localizers = length(localizer_info);
n_thresholds_per_localizer = nan(1, n_localizers);
for j = 1:n_localizers
    n_thresholds_per_localizer(j) = length(localizer_info(j).thresholds);
end

% initialize PSC matrix
% runs x conditions x thresholds
psc = nan([ length(test_info.runs), length(test_info.conditions), ...
    n_thresholds_per_localizer ]);

% number of voxels selected for different relative thresholds
n_voxels_per_run_and_threshold = ...
    nan([length(test_info.runs), n_thresholds_per_localizer]);

for k = 1:length(test_info.runs) % loop through runs
    
    if I.verbose
        % print information about this localizer
        fprintf('Test: %s, run %d\n', test_info.runtype, test_info.runs(k));
        drawnow;
    end
    
    % matrix of psc values for each voxels
    % condition x voxel matrix
    voxel_psc = test_psc(test_info, us, test_info.runs(k), ...
        fwhm, grid_spacing_mm, grid_roi);
        
    % create the mask
    n_voxels = size(voxel_psc,2);
    if isempty(I.anatomical_mask)
        mask = true(1,n_voxels);
    else
        mask = label2grid(I.anatomical_mask, grid_roi, grid_spacing_mm);
    end
    mask = mask > 0.99;
    
    % check there are no exactly zero values
    assert(~any(voxel_psc(:)==0));
    
    % read in the localizer contrast matrix
    localizer_contrast_stat_matrix = nan(n_localizers, n_voxels);
    for j = 1:n_localizers
        
        % ensure non-independence
        if strcmp(test_info.exp, localizer_info(j).exp) && strcmp(test_info.runtype, localizer_info(j).runtype)
            localizer_runs_to_use = setdiff(localizer_info(j).runs, test_info.runs(k));
        else
            localizer_runs_to_use = localizer_info(j).runs;
        end
        
        % subsample the localizer runs
        % (e.g. in order to match response reliability)
        if ~isinf(localizer_info(j).max_runs_to_use)
            
            % distance of localizer runs to test run
            dist = (localizer_runs_to_use - test_info.runs(k)).^2;
            
            % which runs to use based on distance
            switch localizer_info(j).use_nearest_or_farthest
                case 'nearest'
                    [~,xi] = sort(dist, 'ascend');
                case 'farthest'
                    fprintf('Using the farthest runs\n'); drawnow;
                    [~,xi] = sort(dist, 'descend');
                otherwise
                    error('nearest_or_farthest_runs cannot be %s', ...
                        localizer_info(j).use_nearest_or_farthest);
            end
            localizer_runs_to_use = ...
                localizer_runs_to_use(xi(1:localizer_info(j).max_runs_to_use));
            clear dist;
            
        end
        
        % check there is at least one usable localizer run
        if isempty(localizer_runs_to_use)
            error('There needs to be at least one usable localizer run.\n');
        end
        
        % print the runs being used
        if I.verbose
            fprintf('Localizer: %s, runs %s\n', ...
                localizer_info(j).contrast, sprintf('%d',localizer_runs_to_use));
            drawnow;
        end
        
        fprintf('Finding pstat\n\n'); drawnow;
        
        % stat to localizer voxels with
        localizer_contrast_stat_matrix(j,:) = localizer_stat(...
            localizer_info(j), us, localizer_runs_to_use, ...
            fwhm, grid_spacing_mm, grid_roi);
        
    end
    
    % check there are no exactly zero values
    assert(~any(localizer_contrast_stat_matrix(:)==0));
    % localizer_contrast_stat_matrix(localizer_contrast_stat_matrix==0) = NaN;
        
    % loop through all combinations of thresholds, selecting voxels, and
    % measuring mean PSC values
    for i = 1:prod(n_thresholds_per_localizer)
        
        % the index into the threshold vector for each localizer
        threshold_indices = cell(1,n_localizers);
        [threshold_indices{:}] = ind2sub(n_thresholds_per_localizer, i);
        
        % remove NaN voxels
        % logical all is applied to conditions for each voxel
        voxels_in_roi = find(...
            all(~isnan(voxel_psc),1) ...
            & all(~isnan(localizer_contrast_stat_matrix),1) ...
            & mask);
        
        % loop through the localizers
        for j = 1:n_localizers
            
            % relavent info for single localizer
            thresh = localizer_info(j).thresholds(threshold_indices{j});
            stat = localizer_contrast_stat_matrix(j,:) * sign(localizer_info(j).contrast_sign);
            
            % check all statistics are not NaN
            assert(all(~isnan(stat(voxels_in_roi(:)))));
            
            % select voxels from those left
            switch localizer_info(j).threshold_type
                case 'absolute'
                    xi = stat(voxels_in_roi) > thresh;
                    voxels_in_roi = voxels_in_roi(xi);
                case 'relative'
                    n_voxels_to_select = round( thresh * length(voxels_in_roi) );
                    [~,xi] = sort(stat(voxels_in_roi), 'descend');
                    voxels_in_roi = voxels_in_roi(xi(1:n_voxels_to_select));
                otherwise
                    error('Localizer "selection type" should be "absolute-threshold" or "relative-threshold" not %s', localizer_info(j).threshold_type);
            end
        end
        
        n_voxels_after_selection = length(voxels_in_roi);
        n_voxels_per_run_and_threshold(k, threshold_indices{:}) = n_voxels_after_selection;
        
        if isempty(voxels_in_roi)
            warning('No Voxels in ROI');
            continue;
        end
        
        voxel_psc_in_roi = voxel_psc(:, voxels_in_roi);
        
        % check voxels do not have NaN values
        assert(all(~isnan(voxel_psc_in_roi(:))));
        
        % average PSC values within selected voxels
        psc(k,:,threshold_indices{:}) = mean(voxel_psc_in_roi, 2);
        
    end
end

% remove single dimensions for the thresholds
n_thresh_dim = setdiff(n_thresholds_per_localizer,1);
if isempty(n_thresh_dim);
    n_thresh_dim = 1;
end
psc = reshape(psc, ...
    [length(test_info.runs), length(test_info.conditions), n_thresh_dim]);
n_voxels_per_run_and_threshold = reshape(...
    n_voxels_per_run_and_threshold, [length(test_info.runs), n_thresh_dim]);

% return conditions used
conditions = test_info.conditions;

% find/compute desired psc file
function voxel_psc = test_psc(test_info, us, test_run, ...
    fwhm, grid_spacing_mm, grid_roi)

switch test_info.psc_method
    case 'sigav'
        
        [~,MAT_file_first_level] = sigav_surf_grid(...
            test_info.exp, us, test_info.runtype, ...
            fwhm, grid_spacing_mm, grid_roi, ...
            test_info.condition_names_file, ...
            'runs', test_run, 'overwrite', test_info.overwrite, ...
            'para_prefix', test_info.para_prefix);
        assert(length(MAT_file_first_level)==1);
        load(MAT_file_first_level{1}, 'psc');
        voxel_psc = psc;
        
    case 'glm'
        
        [~,MAT_file_first_level] = glm_surf_grid(...
            test_info.exp, us, test_info.runtype, ...
            fwhm, test_info.analysis_name, ...
            grid_spacing_mm, grid_roi, test_info.n_perms, ...
            'runs', test_run, 'overwrite', test_info.overwrite, ...
            'plot_surf', false, 'plot_reliability', false, ...
            'para_prefix', test_info.para_prefix);
        assert(length(MAT_file_first_level)==1);
        load(MAT_file_first_level{1}, 'beta_one_per_regressor');
        voxel_psc = beta_one_per_regressor;
        
    otherwise
        
        error('No matching case');
        
end

fprintf('First level file\n%s\n', MAT_file_first_level{1}); drawnow;


% helper function that find the appropriate file with p-values
function loc_stat = localizer_stat(...
    localizer_info, us, localizer_runs_to_use, ...
    fwhm, grid_spacing_mm, grid_roi)

MAT_file_second_level = glm_surf_grid(...
    localizer_info.exp, us, localizer_info.runtype, ...
    fwhm, localizer_info.analysis_name, ...
    grid_spacing_mm, grid_roi, localizer_info.n_perms, ...
    'runs', localizer_runs_to_use, 'plot_surf', false,...
    'plot_reliability', false, 'overwrite', localizer_info.overwrite, ...
    'para_prefix', localizer_info.para_prefix);

fprintf('Second level file\n%s\n', MAT_file_second_level); drawnow;

if localizer_info.n_perms >= 100
    load(MAT_file_second_level, 'logP_permtest');
    loc_stat = logP_permtest;
else
    load(MAT_file_second_level, 'logP_fixed')
    loc_stat = logP_fixed;
end

% select the row of loc_stat with the desired contrast
load(MAT_file_second_level, 'P');
xi = strcmp(localizer_info.contrast, P.contrast_names);
assert(sum(xi)==1);
loc_stat = loc_stat(xi,:);

function localizer_info = default_localizer_parameters(...
    localizer_info, us, test_info)

n_localizers = length(localizer_info);

for j = 1:n_localizers
    if ~isfield(localizer_info(j), 'exp') || isempty(localizer_info(j).exp)
        localizer_info(j).exp = test_info.exp;
    end
    
    if ~isfield(localizer_info(j), 'runtype') || ...
            isempty(localizer_info(j).runtype)
        localizer_info(j).runtype = test_info.runtype;
    end
    
    if ~isfield(localizer_info(j), 'runs') || isempty(localizer_info(j).runs)
        localizer_info(j).runs = read_runs(...
            localizer_info(j).exp, us, localizer_info(j).runtype);
    end
    
    if ~isfield(localizer_info(j), 'max_runs_to_use') || ...
            isempty(localizer_info(j).max_runs_to_use)
        localizer_info(j).max_runs_to_use = inf;
    end
    
    if ~isfield(localizer_info(j), 'contrast_sign') || ...
            isempty(localizer_info(j).contrast_sign)
        localizer_info(j).contrast_sign = 1;
    end
    
    if ~isfield(localizer_info(j), 'use_nearest_or_farthest') || ...
            isempty(localizer_info(j).use_nearest_or_farthest)
        localizer_info(j).use_nearest_or_farthest = 'nearest';
    end
    
    if ~isfield(localizer_info(j), 'overwrite')
        localizer_info(j).overwrite = false;
    end
    
    if ~isfield(localizer_info(j), 'para_prefix')
        localizer_info(j).para_prefix = localizer_info(j).runtype;
    end

end

function test_info = default_test_parameters(test_info, us)

% runs to use
if ~isfield(test_info, 'runs')
    test_info.runs = read_runs(...
        test_info.exp, us, test_info.runtype);
end

% default use signal averaging to measure each voxel's response
if ~isfield(test_info, 'psc_method')
    test_info.psc_method = 'sigav';
end

% test info
if ~isfield(test_info, 'overwrite')
    test_info.overwrite = false;
end

% prefix to the para files
if ~isfield(test_info, 'para_prefix')
    test_info.para_prefix = test_info.runtype;
end

load(test_info.condition_names_file, 'condition_names');
test_info.conditions = condition_names;
