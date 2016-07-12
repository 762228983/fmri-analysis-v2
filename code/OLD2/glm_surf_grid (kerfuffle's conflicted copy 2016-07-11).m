function [matfile_second_level, matfile_first_level] = ...
    glm_surf_grid(exp, us, runtype, fwhm, analysis_name, ...
    grid_spacing_mm, grid_roi, n_perms, varargin)

global root_directory;

if nargin < 8
    n_perms = 10;
end

% analysis directory
analysis_directory = [root_directory  '/' exp '/analysis/glm/' analysis_name ...
    '/fsaverage_smooth-' num2str(fwhm) 'mm' ...
    '_' 'grid-' num2str(grid_spacing_mm) 'mm' ...
    '_' grid_roi '/usub' num2str(us) '/'];

% condition weight file
parameter_file = [root_directory '/' exp '/analysis/glm' ...
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
runs = read_runs(exp, us, runtype);
n_runs = length(runs);
para_files = cell(1, n_runs);
data_matrix_files = cell(1, n_runs);
for i = 1:length(runs)
    
    r = runs(i);
    
    % read weighting file
    para_files{i} = [root_directory '/' exp '/data/para/usub' num2str(us) ...
        '/' runtype  '_r' num2str(r) '.par'];
    
    % TR
    TR = read_functional_scan_parameters(exp,us,runtype,r,varargin); %#ok<NASGU>
    
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
    
    if ~exist(data_matrix_files{i}, 'file') || optInputs(varargin, 'overwrite')
                
        % reformat to voxel x datapoint/timepoint
        load(grid_file, 'G');
        data_matrix = grid2matrix(G); %#ok<NASGU,NODEF>
        
        % save
        save(data_matrix_files{i}, 'data_matrix', ...
            'TR', 'G', 'dims_rh', 'dims_lh', '-v7.3');
        
    end
end
    
% perform the second level analysis
[matfile_second_level, matfile_first_level] = ...
    glm_second_level(data_matrix_files, para_files, parameter_file, ...
    n_perms, analysis_directory, varargin{:});

% plot reliability across runs
glm_contrast_map_reliability(matfile_first_level,...
    analysis_directory, figure_directory, varargin{:});

if optInputs(varargin, 'noplot')
    return;
end

% select first or second level analysis
if length(data_matrix_files) > 1
    matfile = matfile_second_level;
else
    matfile = matfile_first_level{1};
end

keyboard;

% grid and interpolate a stat of interest to the surface
stat_to_plot = 'logP_permtest';
load(grid_file, 'G');
X = load(matfile, stat_to_plot, 'P'); %#ok<NASGU>
G = matrix2grid(X.(stat_to_plot)', G);
surf = grid2surface(G);

% plot surface
color_range = [-5 5];
n_contrasts = size(surf,3);
for i = 1:n_contrasts
    hemis = {'rh','lh'};
    for q = 1:2
        
        figure_file = [figure_directory '/' 'pmap_' P.contrast_names{i} '_' hemis{q} '_central95.png'];
        
        if ~exist(figure_file, 'file') || optInputs(varargin, 'overwrite')
            
            if ~exist('figh', 'var') % Plot figures with outlines of standard ROIs
                close all;
                figh = figure;
                pos = get(figh,'Position');
                set(figh, 'Position', [pos(1:2), 800 800]);
            end
           
            % plot surface map, relative threshold, central 95% of the distribution
            plot_fsaverage_1D_overlay(surf(:,q,i),hemis{q},'parula',color_range,figh);
            export_fig(figure_file,'-png','-r100','-nocrop');
        end
    end
end


% % binary map with relative cutoff
% i = 5;
% [Nx,x] = hist(logP(:,i),100);
% Cx = cumsum(Nx/sum(Nx));
% [~,xi] = unique(Cx);
% x = x(xi);
% Cx = Cx(xi);
% color_range = interp1(Cx,x,[0.025 0.975]);
% 
% % plot
% figure;
% subplot(1,2,1);
% imagesc(flipud(rot90(G.grid_data{1}(:,:,i))));
% colormap('parula');
% colorbar;
% title('Right Hemi');
% subplot(1,2,2);
% imagesc(fliplr(flipud(rot90(G.grid_data{2}(:,:,i))))); %#ok<FLUDLR>
% colormap('parula');
% title('Left Hemi');
% colorbar;


% % grid stats
% stats_to_grid_and_resample = {'logP_permtest', 'logP_fixed'};
% matfile_gridded = strrep(matfile,  '.mat', '_gridded.mat');
% for i = 1:length(stats_to_grid_and_resample)
%     X = load(matfile, stats_to_grid_and_resample{i}); %#ok<NASGU>
%     eval([stats_to_grid_and_resample{i} ...
%         ' = matrix2grid(X.' stats_to_grid_and_resample{i} ''', G);']);   
%     if i == 1
%         save(matfile_gridded, stats_to_grid_and_resample{i});
%     else
%         save(matfile_gridded, '-append', stats_to_grid_and_resample{i});
%     end    
% end
% clear X;