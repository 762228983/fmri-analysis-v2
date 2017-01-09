function combine_runs_wrapper(exp, us, original_runtype, ...
    original_runs, fwhm, grid_roi, grid_spacing_mm, varargin)

% Wrapper for combine_runs_surf_grid, combine_runs_para_file, and
% combine_runs_scan_parameters. Executes all of the scripts necessary to create
% a new runtype from combinations of runs.
% 
% 2016-01-09: Created, Sam NH

% combine grid files
combine_runs_surf_grid(exp, us, original_runtype, ...
    original_runs, fwhm, grid_roi, grid_spacing_mm, varargin{:});

% para files
combine_runs_para_file(exp, us, original_runtype, original_runs, varargin{:});

% scan parameters
combine_runs_scan_parameters(exp, us, original_runtype, original_runs, varargin{:})

