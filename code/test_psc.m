function voxel_psc = test_psc(test_info, us, test_run, ...
    fwhm, grid_spacing_mm, grid_roi)

% Helper function fro roi_surf_grid.m and component_localizer_surf_grid.m
% 
% 2016-09-2: Separated out from roi_surf_grid.m, Sam NH

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