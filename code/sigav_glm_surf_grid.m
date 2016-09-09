function [MAT_file_second_level, MAT_files_first_level,...
    analysis_directory, figure_directory] = ...
    sigav_glm_surf_grid(exp, us, runtype, fwhm, analysis_name, ...
    grid_spacing_mm, grid_roi, varargin)

% Pipeline wrapper for applying sigav_glm.m to data from a surface grid
% representation in the pipeline.
%
% 2016-09-09: Last modified, Sam NH

global root_directory;

% default parameters
I.overwrite = false;
I.onset_delay = 5;
I.offset_delay = 1;
I.para_prefix = runtype;

% default runs
I.runs = read_runs(exp, us, runtype);

% set optional inputs
I = parse_optInputs_keyvalue(varargin, I);

% analysis directory
analysis_directory = [...
    root_directory  '/' exp '/analysis/sigav-glm/' analysis_name  ...
    '/fsaverage_smooth-' num2str(fwhm) 'mm' ...
    '_' 'grid-' num2str(grid_spacing_mm) 'mm' ...
    '_' grid_roi '_onsdelay' num2str(I.onset_delay) ...
    '_offdelay' num2str(I.offset_delay) '/usub' num2str(us)];

% condition weight file
parameter_file = [root_directory '/' exp '/analysis/sigav-glm' ...
    '/' analysis_name '/parameters.mat'];

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
    
    % TR
    TR = read_functional_scan_parameters(exp,us,runtype,r); %#ok<NASGU>
    
    % para file
    para_files{i} = [root_directory '/' exp '/data/para/usub' num2str(us) ...
        '/' I.para_prefix  '_r' num2str(r) '.par'];
    
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
        [analysis_directory '/r' num2str(r) '.mat'];
    
end

% file to save results of second level analysis pooling across runs
MAT_file_second_level = [analysis_directory '/r' sprintf('%d', I.runs) '.mat'];

sigav_second_level(...
    data_matrix_files, para_files, parameter_file, ...
    MAT_file_second_level, MAT_files_first_level, ...
    'onset_delay', I.onset_delay, 'offset_delay', I.offset_delay, ...
    'overwrite', I.overwrite);