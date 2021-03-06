function combine_runs_wrapper(exp, us, original_runtype, ...
    original_runs, fwhm, grid_roi, grid_spacing_mm, varargin)

% Wrapper for combine_runs_surf_grid, combine_runs_para_file, and
% combine_runs_scan_parameters. Executes all of the scripts necessary to create
% a new runtype from combinations of runs.
%
% 2016-01-09: Created, Sam NH

% optional input arguments
I.overwrite = false;
I.detrend = 0;
I.combined_runtype = [original_runtype '_combined'];
I.combined_runs = 1:length(original_runs);
I = parse_optInputs_keyvalue(varargin, I);

% combine grid files
combine_runs_surf_grid(exp, us, original_runtype, ...
    original_runs, fwhm, grid_roi, grid_spacing_mm, ...
    'combined_runtype', I.combined_runtype, ...
    'overwrite', I.overwrite, ...
    'combined_runs', I.combined_runs, ...
    'detrend', I.detrend);

% para files
combine_runs_para_file(exp, us, original_runtype, ...
    original_runtype, original_runs, ...
    'combined_para_prefix', I.combined_runtype, ...
    'overwrite', I.overwrite, ...
    'combined_runs', I.combined_runs);

% scan parameters
combine_runs_scan_parameters(exp, us, original_runtype, original_runs, ...
    'combined_runtype', I.combined_runtype, ...
    'overwrite', I.overwrite, ...
    'combined_runs', I.combined_runs)

% run orders
combine_runs_runorder_file(exp, us, original_runtype, original_runs, ...
    'combined_runtype', I.combined_runtype, ...
    'overwrite', I.overwrite, ...
    'combined_runs', I.combined_runs)

