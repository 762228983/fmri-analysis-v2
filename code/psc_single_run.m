function [voxel_psc, condition_names] = psc_single_run(...
    test_info, us, test_run, fwhm, grid_spacing_mm, grid_roi)

% Helper function fro roi_surf_grid.m and component_localizer_surf_grid.m
%
% 2016-09-02: Separated out from roi_surf_grid.m, Sam NH
% 
% 2016-09-09: Including sigav-glm and incorporating with standard GLM analysis

switch test_info.analysis_type
    case 'sigav'
        
        [~,MAT_file_first_level] = sigav_surf_grid(...
            test_info.exp, us, test_info.runtype, ...
            fwhm, grid_spacing_mm, grid_roi, ...
            test_info.condition_names_file, ...
            'runs', test_run, 'overwrite', test_info.overwrite);
        assert(length(MAT_file_first_level)==1);
        load(MAT_file_first_level{1}, 'psc');
        voxel_psc = psc;
        load(test_info.condition_names_file,  'condition_names');
        
    case {'sigav-glm','glm'}
        
        [~,MAT_file_first_level] = glm_surf_grid(...
            test_info.exp, us, test_info.runtype,...
            fwhm, test_info.analysis_name, ...
            grid_spacing_mm, grid_roi, ...
            'analysis_type', test_info.analysis_type, ...
            'runs', test_run, 'overwrite', test_info.overwrite, ...
            'plot_surf', false, 'plot_reliability', false);
        
        assert(length(MAT_file_first_level)==1);
        load(MAT_file_first_level{1}, 'beta_one_per_regressor','P');
        voxel_psc = beta_one_per_regressor;
        condition_names = P.regressor_names;
        
    otherwise
        
        error('No matching case');
        
end

fprintf('%s\n', MAT_file_first_level{1}); drawnow;