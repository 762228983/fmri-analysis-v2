function test_info = default_test_parameters(test_info, us)

% Helper function fro roi_surf_grid.m and component_localizer_surf_grid.m
% 
% 2016-09-02: Separated out from roi_surf_grid.m, Sam NH
% 
% 2016-09-04: No longer defaults to older signal averaging analysis, condition
% names are not included in the defaults

% runs to use
if ~isfield(test_info, 'runs')
    test_info.runs = read_runs(...
        test_info.exp, us, test_info.runtype);
end

% test info
if ~isfield(test_info, 'overwrite')
    test_info.overwrite = false;
end

% prefix to the para files
if ~isfield(test_info, 'para_prefix')
    test_info.para_prefix = test_info.runtype;
end