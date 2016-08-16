function corr_test_retest = component_voxel_weights_reliability(...
    weight_file_individual_runs, analysis_directory, figure_directory, varargin)

global root_directory;

% add analysis code directory, throw error if it doesn't exist
if ~exist([root_directory '/general-analysis-code'], 'dir')
    error('general-analysis-code directory is not in the root_directory');
end
addpath([root_directory '/general-analysis-code']);

% load weights from each run
n_runs = length(weight_file_individual_runs);
for i = 1:n_runs
    
    % load weights for single run
    load(weight_file_individual_runs{i}, 'component_weights', 'component_names');
    
    % initialize weights
    if i == 1
        [n_voxels, n_components] = size(component_weights);
        weights_allruns = nan([n_voxels, n_runs, n_components]);
    end
    
    % assign weights for this run
    weights_allruns(:,i,:) = component_weights;
    
end

% voxels without any NaNs
vox_without_NaNs = all(all(~isnan(weights_allruns),2),3);

% differt numbers of runs to test when measuring split-half reliability
data_set_sizes = [2.^(0:1:log2(floor(n_runs/2))), floor(n_runs/2)];
data_set_sizes = unique(data_set_sizes);
n_data_set_sizes = length(data_set_sizes);

% measure test-retest correlation for different data set sizes
n_smps = 1000;
mat_file = [analysis_directory '/test-retest-reliability-' ...
    num2str(n_smps) 'smps-' DataHash(weight_file_individual_runs)  '.mat'];
if ~exist(mat_file, 'file') || optInputs(varargin, 'overwrite')
    fprintf('Sampling test-retest reliability...\n');
    corr_test_retest = nan(n_smps, n_data_set_sizes, n_components);
    for i = 1:n_data_set_sizes
        for k = 1:n_components
            weights_single_component = weights_allruns(vox_without_NaNs,:,k);                  
            for j = 1:n_smps
                % random permutation of runs
                xi = randperm(n_runs);
                xi_split1 = xi(1:data_set_sizes(i));
                xi_split2 = xi(data_set_sizes(i)+1:2*data_set_sizes(i));
                split1 = weights_single_component(:,xi_split1);
                split2 = weights_single_component(:,xi_split2);
                corr_test_retest(j,i,k) = corr(mean(split1,2), mean(split2,2));
            end
        end
    end
    save(mat_file, 'corr_test_retest');
else
    load(mat_file, 'corr_test_retest');
end
clear xi xi_split1 xi_split2 split1 split2;

% infer test-retest reliability for all runs using spearman brown
corr_allruns = correlation_power_analysis(...
    corr_test_retest(:,n_data_set_sizes,:), n_runs/floor(n_runs/2), 0);

% add to test-retest matrix
corr_test_retest = cat(2, corr_test_retest, corr_allruns);
data_set_sizes = [data_set_sizes, n_runs];
n_data_set_sizes = n_data_set_sizes+1; %#ok<NASGU>

% plot
figure_fname = [figure_directory '/test-retest-reliability-' num2str(n_smps) 'smps.pdf'];
if ~exist(figure_fname, 'file') || optInputs(varargin, 'overwrite')
    figure;
    errorbar_plot_from_samples(corr_test_retest, log2(data_set_sizes));
    set(gca, 'XTick', log2(data_set_sizes), 'XTickLabel', data_set_sizes);
    xL = xlim;
    hold on; plot(xL, [0 0], 'k--');
    xlabel('Number of Runs'); ylabel('Test-Retest Correlation');
    legend(component_names,'Location','EastOutside');
    box off;
    export_fig(figure_fname,'-pdf','-transparent');
end





