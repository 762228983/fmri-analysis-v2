function [MAT_file_second_level, MAT_files_first_level,...
    analysis_directory, figure_directory] = ...
    sigav_surf_grid(exp, us, runtype, fwhm, ...
    grid_spacing_mm, grid_roi, condition_names_file, varargin)

% Signal averages responses measured on a surface grid
%
% 2016-08-24 - Created, Sam NH
%
% 2016-08-25 - Modified to work with sigav_second_level.m, Sam NH
%
% 2016-08-27: Modified how optional arguments are handled

global root_directory;

% optional arguments and defaults
I.overwrite = false;
I.onset_delay = 5;
I.offset_delay = 1;
I.remove_run_offsets = true;
I.runs = read_runs(exp, us, runtype);
I = parse_optInputs_keyvalue(varargin, I);

% analysis directory
analysis_directory = [...
    root_directory  '/' exp '/analysis/sigav'  ...
    '/fsaverage_smooth-' num2str(fwhm) 'mm' ...
    '_' 'grid-' num2str(grid_spacing_mm) 'mm' ...
    '_' grid_roi '_onsdelay' num2str(I.onset_delay) ...
    '_offdelay' num2str(I.offset_delay) '/usub' num2str(us) '/'];

% analysis directory
figure_directory = strrep(analysis_directory, 'analysis', 'figures');

% create this directories if not present
if ~exist(analysis_directory, 'dir')
    mkdir(analysis_directory);
end
if ~exist(figure_directory, 'dir')
    mkdir(figure_directory);
end

% create cell struction with para files and data files
fprintf('Converting surface files to data matrix...\n');
n_runs = length(I.runs);
para_files = cell(1, n_runs);
data_matrix_files = cell(1, n_runs);
MAT_files_first_level = cell(1, n_runs);
for i = 1:length(I.runs)
    
    r = I.runs(i);
    
    % read weighting file
    para_files{i} = [root_directory '/' exp '/data/para/usub' num2str(us) ...
        '/' runtype  '_r' num2str(r) '.par'];
    
    % TR
    TR = read_functional_scan_parameters(exp,us,runtype,r); %#ok<NASGU>
    
    % preprocessing directory with files in fsaverage space
    preproc_fsaverage_directory = [root_directory '/' exp '/analysis/preprocess' ...
        '/usub' num2str(us) '/' runtype '_r' num2str(r) '/myfsaverage'];
    
    % input surface grid
    grid_file = [preproc_fsaverage_directory '/' ...
        'smooth-' num2str(fwhm) 'mm' '_' ...
        'grid-' num2str(grid_spacing_mm) 'mm' '_' grid_roi '.mat'];
    
    % reformated data matrix file
    data_matrix_files{i} = ...
        strrep(grid_file, '.mat', '_unwrapped_data_matrix.mat');
    
    if ~exist(data_matrix_files{i}, 'file') || I.overwrite
        
        % reformat to voxel x datapoint/timepoint
        load(grid_file, 'G');
        data_matrix = grid2matrix(G); %#ok<NASGU>
        
        % save
        save(data_matrix_files{i}, 'data_matrix', 'TR', 'G', '-v7.3');
        
    end
    
    % file to save results of individual run analysis
    MAT_files_first_level{i} = ...
        [analysis_directory '/r' num2str(r) ...
        '_' num2str(n_perms) 'perms.mat'];
    
end

% file to save results of second level analysis pooling across runs
MAT_file_second_level = [analysis_directory '/r' sprintf('%d', I.runs) ...
        '_' num2str(n_perms) 'perms.mat'];

sigav_second_level(...
    data_matrix_files, para_files, condition_names_file, ...
    MAT_file_second_level, MAT_files_first_level, ...
    'onset_delay', I.onset_delay, 'offset_delay', I.offset_delay, ...
    'remove_run_offsets', I.remove_run_offsets, 'overwrite', I.overwrite);