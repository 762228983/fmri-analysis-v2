function [comp_psc, condition_names, component_names] = ...
    component_localizer_surf_grid(  us, grid_roi, grid_spacing_mm, ...
    component_info, test_info, fwhm, varargin )

% 2016-09-09: Last modified, Sam NH

% optional arguments
I.verbose = true;
I.anatomical_mask = '';
I.thresh_logP_residual_permtest = -inf;
I = parse_optInputs_keyvalue(varargin, I);

% default parameters
test_info = ...
    default_test_parameters(test_info, us);
component_info = ...
    default_component_parameters(component_info, us);

% add a permutation test in order to evaluate significance
if I.thresh_logP_residual_permtest > -inf
    component_info.n_perms = 100;
end

for k = 1:length(test_info.runs) % loop through runs
    
    if I.verbose
        % print information about this localizer
        fprintf('Test: %s, run %d\n', test_info.runtype, test_info.runs(k));
        drawnow;
    end
    
    % matrix of psc values for each voxel
    % condition x voxel matrix
    [voxel_psc, condition_names] = psc_single_run(test_info, us, ...
        test_info.runs(k), fwhm, grid_spacing_mm, grid_roi);
    
    % create the mask
    n_voxels = size(voxel_psc,2);
    if isempty(I.anatomical_mask)
        mask = true(1,n_voxels);
    else
        mask = label2grid(I.anatomical_mask, grid_roi, grid_spacing_mm);
    end
    mask = mask > 0.99;
    
    % ensure non-independence
    if strcmp(test_info.exp, component_info.exp) ...
            && strcmp(test_info.runtype, component_info.runtype)
        localizer_runs_to_use = setdiff(...
            component_info.runs, test_info.runs(k));
    else
        localizer_runs_to_use = component_info.runs;
    end
    
    % component weights and significance values from a permutation test
    [comp_weights, logP_residual_permtest, component_names] = ...
        localizer_weights(component_info, us, localizer_runs_to_use, ...
        fwhm, grid_spacing_mm, grid_roi);
    
    % further select voxels based on significance values of permutation test
    if I.thresh_logP_residual_permtest > -inf;
        mask = mask & logP_residual_permtest > I.thresh_logP_residual_permtest;
    end
    
    if k == 1
        % component response matrix, initialization
        comp_psc = nan(length(test_info.runs), length(condition_names),...
            length(component_names));
    end
    
    % measure psc
    mask = mask & all(~isnan(voxel_psc),1) & all(~isnan(comp_weights),1);
    comp_psc(k,:,:) = voxel_psc(:,mask) * pinv(comp_weights(:,mask));
    
end

% helper function that find the appropriate file with p-values
function [beta_one_per_regressor, logP_residual_permtest, component_names] = ...
    localizer_weights(component_info, us, localizer_runs_to_use, ...
    fwhm, grid_spacing_mm, grid_roi) %#ok<STOUT>

[MAT_file_second_level, MAT_files_first_level, ...
    perm_MAT_file_second_level, perm_MAT_files_first_level, P] ...
    = glm_surf_grid(...
    component_info.exp, us, component_info.runtype, ...
    fwhm, component_info.analysis_name, ...
    grid_spacing_mm, grid_roi, ...
    'n_perms', component_info.n_perms, ...
    'analysis_type', component_info.analysis_type, ...
    'runs', localizer_runs_to_use, 'plot_surf', false,...
    'plot_reliability', false, 'overwrite', component_info.overwrite);

if length(localizer_runs_to_use)>1
    fprintf('Second level file\n%s\n', MAT_file_second_level); drawnow;
    load(MAT_file_second_level, 'beta_one_per_regressor');
    if component_info.n_perms > 0
        load(perm_MAT_file_second_level, 'logP_residual_permtest');
    else
        logP_residual_permtest = [];
    end
else
    fprintf('First level file\n%s\n', MAT_files_first_level{1}); drawnow;
    assert(length(MAT_files_first_level)==1);
    load(MAT_files_first_level{1}, 'beta_one_per_regressor');
    if component_info.n_perms > 0
        load(perm_MAT_files_first_level{1}, 'logP_residual_permtest');
    else
        logP_residual_permtest = [];
    end
end

component_names = P.regressor_names;

function component_info = default_component_parameters(...
    component_info, us)

% runs to use
if ~isfield(component_info, 'runs')
    component_info.runs = read_runs(...
        component_info.exp, us, component_info.runtype);
end

% test info
if ~isfield(component_info, 'overwrite')
    component_info.overwrite = false;
end


