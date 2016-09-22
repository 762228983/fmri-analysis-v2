function [corr_test_retest, data_set_sizes, quantiles] = glm_voxel_reliability(...
    MAT_files_first_level, analysis_directory, figure_directory, varargin)

% Reliability of regressors across voxels
% 
% 2016-09-13: Created, Sam NH

% optional arguments and defaults
I.voxel_stat = 'beta_one_per_regressor';
I.overwrite = false;
I.group_runs = {};
I = parse_optInputs_keyvalue(varargin, I);

% quantiles of the distribution to measure
quantiles = 0.1:0.1:0.9;
n_quantiles = length(quantiles);

% measure test-retest correlation for different data set sizes
n_smps = 1e3;
mat_file = [analysis_directory '/voxel_reliability_' I.voxel_stat '_' num2str(n_smps) 'smps.mat'];
if ~exist(mat_file, 'file') || I.overwrite
        
    n_runs = length(MAT_files_first_level);
    
    % load weights from each run
    for i = 1:n_runs
        
        % load weights for single run
        X = load(MAT_files_first_level{i}, I.voxel_stat, 'P');
        
        % initialize weights
        if i == 1
            [n_betas, n_voxels] = size(X.(I.voxel_stat));
            stat_all_runs = nan([n_betas, n_voxels, n_runs]);
        end
        
        % assign weights for this run
        stat_all_runs(:,:,i) = X.(I.voxel_stat);
        
    end
    clear X;
    
    if ~isempty(I.group_runs)
        n_run_groups = length(I.group_runs);
        stat_group_average = nan([n_betas,n_voxels,n_run_groups]);
        for i = 1:n_run_groups
            stat_group_average(:,:,i) = nanmean(stat_all_runs(:,:,I.group_runs{i}),3);
        end
        stat_all_runs = stat_group_average;
        n_runs = n_run_groups;
    end
    
    % voxels without any NaNs
    voxels_without_NaNs = all(all(~isnan(stat_all_runs),1),3);
    
    % differt numbers of runs to test when measuring split-half reliability
    data_set_sizes = [2.^(0:1:log2(floor(n_runs/2))), floor(n_runs/2)];
    data_set_sizes = unique(data_set_sizes);
    n_data_set_sizes = length(data_set_sizes);
    
    fprintf('Sampling voxel reliability...\n');
    corr_test_retest = nan(n_smps, n_data_set_sizes, n_quantiles);
    for i = 1:n_data_set_sizes
        for j = 1:n_smps
            % random permutation of runs
            xi = randperm(n_runs);
            xi_split1 = xi(1:data_set_sizes(i));
            xi_split2 = xi(data_set_sizes(i)+1:2*data_set_sizes(i));
            r = fastcorr( ...
                mean(stat_all_runs(:,voxels_without_NaNs,xi_split1),3),...
                mean(stat_all_runs(:,voxels_without_NaNs,xi_split2),3));
            corr_test_retest(j,i,:) = quantile(r(:), quantiles);
        end
    end
    save(mat_file, 'corr_test_retest', 'n_runs', 'data_set_sizes', 'n_data_set_sizes');
else
    load(mat_file, 'corr_test_retest', 'n_runs', 'data_set_sizes', 'n_data_set_sizes');
end
clear xi xi_split1 xi_split2 r;

% infer test-retest reliability for all runs using spearman brown
corr_allruns = correlation_power_analysis(...
    corr_test_retest(:,n_data_set_sizes,:), n_runs/floor(n_runs/2), 0);

% add to test-retest matrix
corr_test_retest = cat(2, corr_test_retest, corr_allruns);
data_set_sizes = [data_set_sizes, n_runs];
n_data_set_sizes = n_data_set_sizes+1; %#ok<NASGU>

% plot
figure_fname = [figure_directory '/voxel_reliability_' I.voxel_stat '_' num2str(n_smps) 'smps.pdf'];
if ~exist(figure_fname, 'file') || I.overwrite
    figure;
    cols = parula(n_quantiles);
    h = errorbar_plot_from_samples(corr_test_retest, log2(data_set_sizes));
    for i = 1:length(h)
        set(h(i), 'Color', cols(i,:));
    end
    set(gca, 'XTick', log2(data_set_sizes), 'XTickLabel', data_set_sizes);
    xL = xlim;
    hold on; plot(xL, [0 0], 'k--');
    xlabel('Number of Runs'); ylabel('Test-Retest Correlation');
    h = legend(cellstr(num2str(quantiles')),'Location','EastOutside');
    set(get(h, 'title'), 'string', 'Quantiles');
    box off;
    title('Contrast Map Reliability');
    export_fig(figure_fname,'-pdf','-transparent');
end