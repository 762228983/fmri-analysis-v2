function [unwrapped_grid, G, output_file] = ...
    label2grid(label_name, grid_roi, grid_spacing_mm, varargin)

global root_directory;

% optional arguments and defaults
I.overwrite = false;
I = parse_optInputs_keyvalue(varargin, I);

% top-level freesurfer directory with fsaverage data
freesurfer_directory = [root_directory '/freesurfer/myfsaverage'];

% file to save grid to
output_directory = [freesurfer_directory '/grid_labels/' ...
    grid_roi '_' num2str(grid_spacing_mm) 'mm'];
output_file = [output_directory '/' label_name '.mat'];
output_directory
if ~exist(output_directory, 'dir')
    mkdir(output_directory);
end

% load previously computed results
if exist(output_file, 'file') && ~I.overwrite
    load(output_file, 'G');
    unwrapped_grid = grid2matrix(G); %#ok<NODEF>
    return;
end

% source label files
label_rh = [freesurfer_directory '/label/rh.' label_name '.label'];
label_lh = [freesurfer_directory '/label/lh.' label_name '.label'];

label_rh

% convert to surface file, right hemisphere
n_surfpts = 163842;
lb = read_label_SNH(label_rh);
surf_verts = zeros(n_surfpts,1);
surf_verts(lb.vnums+1) = 1;
surf_rh = [freesurfer_directory '/rh.' num2str(randi(2^52)) '.mgz'];
MRIwrite_surface(surf_verts, surf_rh, 'rh');
clear lb surf_verts;

% convert to surface file, left hemisphere
n_surfpts = 163842;
lb = read_label_SNH(label_lh);
surf_verts = zeros(n_surfpts,1);
surf_verts(lb.vnums+1) = 1;
surf_lh = [freesurfer_directory '/lh.' num2str(randi(2^52)) '.mgz'];
MRIwrite_surface(surf_verts, surf_lh, 'lh');
clear lb surf_verts;

% interpolate to the grid
plot_figures = false;
patch_rh = [freesurfer_directory '/surf/rh.cortex.patch.flat'];
patch_lh = [freesurfer_directory '/surf/lh.cortex.patch.flat'];
grid_rh_label = [freesurfer_directory '/label/rh.' grid_roi '.label'];
grid_lh_label = [freesurfer_directory '/label/lh.' grid_roi '.label'];
G = interp_from_surface_to_grid(surf_rh, surf_lh, patch_rh, patch_lh, ...
    grid_rh_label, grid_lh_label, grid_spacing_mm, plot_figures);

% delete temporary surface files
delete(surf_rh);
delete(surf_lh);

% save grid
save(output_file, 'G');

% unwrap grid to vector
unwrapped_grid = grid2matrix(G);
