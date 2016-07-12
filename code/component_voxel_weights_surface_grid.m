function [weight_file_runaverage, weight_file_individual_runs] = ...
    component_voxel_weights_surface_grid(exp, us, runtype, fwhm, ...
    grid_spacing_mm, grid_roi, component_analysis_name, n_perms, varargin)

global root_directory;

% condition weight file
condition_weights_file = [root_directory '/' exp '/analysis/' component_analysis_name '/condition-weights.mat'];

% analysis directory
analysis_directory = [root_directory  '/' exp '/analysis/' component_analysis_name '/voxel-weights/' 'usub' num2str(us) ...
    '/' 'fsaverage_smooth-' num2str(fwhm) 'mm_grid-' num2str(grid_spacing_mm) 'mm_' grid_roi];

% analysis directory
figure_directory = [root_directory  '/' exp '/figures/' component_analysis_name '/voxel-weights/' 'usub' num2str(us) ...
    '/' 'fsaverage_smooth-' num2str(fwhm) 'mm_grid-' num2str(grid_spacing_mm) 'mm_' grid_roi];    
    
% create this directory if it doesn't exist
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
data_files = cell(1, n_runs);
for i = 1:length(runs)
    
    r = runs(i);
    
    % read weighting file
    para_files{i} = [root_directory '/' exp '/data/para/usub' num2str(us) '/' runtype  '_r' num2str(r) '.par'];
    
    % TR
    TR = read_functional_scan_parameters(exp,us,runtype,r,varargin);
    
    % input data files
    data_files{i} = [root_directory '/' exp '/analysis/preprocess/usub' num2str(us) '/' runtype '_r' num2str(r) '/myfsaverage' ...
        '/' 'smooth-' num2str(fwhm) 'mm_grid-' num2str(grid_spacing_mm) 'mm_' grid_roi '_unwrapped_data_matrix.mat'];
    
    if ~exist(data_files{i}, 'file') || optInputs(varargin, 'overwrite')
        
        % read input surface grid and form data matrix
        grid_file = [root_directory '/' exp '/analysis/preprocess/usub' num2str(us) '/' runtype '_r' num2str(r) '/myfsaverage' ...
            '/' 'smooth-' num2str(fwhm) 'mm_grid-' num2str(grid_spacing_mm) 'mm_' grid_roi '.mat'];
        load(grid_file, 'G');
        dims_rh = size(G.grid_data{1});
        dims_lh = size(G.grid_data{2});
        data_matrix = [reshape(G.grid_data{1}, [prod(dims_rh(1:2)),dims_rh(3)])', reshape(G.grid_data{2}, [prod(dims_lh(1:2)),dims_lh(3)])'];
        
        save(data_files{i}, 'data_matrix', 'TR', 'G', 'dims_rh', 'dims_lh');
        
    end
end
    
% perform the analysis
[weight_file_runaverage, weight_file_individual_runs] = ...
    component_voxel_weights(data_files, para_files, condition_weights_file, ...
    analysis_directory, n_perms, varargin{:});

% plot reliability
component_voxel_weights_reliability(...
    weight_file_individual_runs, analysis_directory, ...
    figure_directory, varargin{:});

% load and reformat the computed significance maps
load([analysis_directory '/all-runs_nperms' num2str(n_perms)  '.mat'], 'logP', 'component_names');
n_components = size(logP,2);
load(data_files{1}, 'G', 'dims_rh', 'dims_lh');
G_logP = G;
G_logP.grid_data{1} = reshape(logP(1:prod(dims_rh(1:2)),:), [dims_rh(1:2), n_components]);
G_logP.grid_data{2} = reshape(logP(prod(dims_rh(1:2))+1:end,:), [dims_lh(1:2), n_components]);

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
% imagesc(flipud(rot90(G_logP.grid_data{1}(:,:,i))));
% colormap('parula');
% colorbar;
% title('Right Hemi');
% subplot(1,2,2);
% imagesc(fliplr(flipud(rot90(G_logP.grid_data{2}(:,:,i))))); %#ok<FLUDLR>
% colormap('parula');
% title('Left Hemi');
% colorbar;

% create surface maps
nsurfpts = 163842;
interpolated_surface_map = nan(1, nsurfpts);
for i = 1:n_components
    hemis = {'rh','lh'};
    for q = 1:2
        surface_file = [analysis_directory '/' hemis{q} '.' component_names{i} '.mgz'];
        if ~exist(surface_file, 'file') || optInputs(varargin, 'overwrite')
            interpolated_surface_map(G_logP.vi{q}+1) = interp2(G_logP.grid_x{q},G_logP.grid_y{q} ,G_logP.grid_data{q}(:,:,i), G_logP.vras{q}(:,1), G_logP.vras{q}(:,2), 'linear');
            switch hemis{q}
                case 'rh'
                    MRIwrite_surface(interpolated_surface_map, surface_file, 'rh');
                case 'lh'
                    MRIwrite_surface(interpolated_surface_map, surface_file, 'lh');
            end
        end
    end
end


color_range = [-5 5];

for i = 1:n_components
    for q = 1:2
        
        figure_file = [figure_directory '/' 'pmap_' component_names{i} '_' hemis{q} '_central95.png'];
        
        if ~exist(figure_file, 'file') || optInputs(varargin, 'overwrite')
            
            if ~exist('figh', 'var') % Plot figures with outlines of standard ROIs
                close all;
                figh = figure;
                pos = get(figh,'Position');
                set(figh, 'Position', [pos(1:2), 800 800]);
            end
            
            surface_file = [analysis_directory '/' hemis{q} '.' component_names{i} '.mgz'];
            surf = MRIread(surface_file);
            
            %         % by default plots the central 95% of the distribution
            %         [Nx,x] = hist(logP(:,i),100);
            %         Cx = cumsum(Nx/sum(Nx));
            %         [~,xi] = unique(Cx);
            %         x = x(xi);
            %         Cx = Cx(xi);
            %         color_range = interp1(Cx,x,[0.025 0.975]);
            
            % plot surface map, relative threshold, central 95% of the distribution
            plot_fsaverage_1D_overlay(surf.vol,hemis{q},'parula',color_range,figh);
            export_fig(figure_file,'-png','-r100','-nocrop');
        end
    end
end
