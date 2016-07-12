function grid_mat_file = downsample_surface_timecourses(exp, us, runtype, r, fwhm, grid_spacing_mm, grid_roi, varargin)
% function downsample_surface_timecourses(exp, us, runtype, r, input_fname, grid_spacing_mm, plot_figures, varargin)
% 
% Downsamples preprocessed timecourse files to a grid on the flattened
% cortical surface. Useful because Freesurfer uses an excessively fine
% tessellation for most purposes. 
% 
% The function relies on the subfunction "interp_from_surface_to_grid.m"
% to perform the interpolation. 

global root_directory;

subjid = ['us' num2str(us)];

preprocessing_directory = [root_directory '/' exp '/analysis/preprocess/usub' num2str(us) '/' runtype '_r' num2str(r)];

% map to fsaverage template
if optInputs(varargin, 'monkey')
    surface_directory = [preprocessing_directory '/surf'];
else
    surface_directory = [preprocessing_directory '/myfsaverage'];
end

% flattened patches and constraint ROIs for monkeys and humans
if optInputs(varargin, 'monkey')
    freesurfer_directory = [root_directory '/freesurfer/' subjid];
else
    freesurfer_directory = [root_directory '/freesurfer/myfsaverage'];
end

patch_rh = [freesurfer_directory '/surf/rh.cortex.patch.flat'];
patch_lh = [freesurfer_directory '/surf/lh.cortex.patch.flat'];
roi_rh_label = [freesurfer_directory '/label/rh.' grid_roi '.label'];
roi_lh_label = [freesurfer_directory '/label/lh.' grid_roi '.label'];

if optInputs(varargin, 'plot')
    plot_figures = true;
else
    plot_figures = false;
end

% scan parameters
[TR, TA, nTR, n_disdaqs] = read_functional_scan_parameters(exp,us,runtype,r,varargin); %#ok<ASGLU>

% mat file with gridded values
grid_mat_file = [surface_directory '/' 'smooth-' num2str(fwhm) 'mm_grid-' num2str(grid_spacing_mm) 'mm_' grid_roi '.mat'];
if ~exist(grid_mat_file, 'file') || optInputs(varargin, 'overwrite')
    
    % surface files using Freesurfer's fine tessellation
    surface_rh = [surface_directory '/' 'rh.' 'smooth-' num2str(fwhm) 'mm.mgz'];
    surface_lh = [surface_directory '/' 'lh.' 'smooth-' num2str(fwhm) 'mm.mgz'];
    
    % interpolate to the grid
    G = interp_from_surface_to_grid(surface_rh, surface_lh, patch_rh, patch_lh, roi_rh_label, roi_lh_label, grid_spacing_mm, plot_figures); %#ok<NASGU>
    
    % save to the mat file
    save(grid_mat_file, 'G', 'TR');
    
end