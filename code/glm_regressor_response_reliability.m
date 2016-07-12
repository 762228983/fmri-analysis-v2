function corr_test_retest = glm_regressor_response_reliability(...
    matfile_first_level, analysis_directory, figure_directory, varargin)

global root_directory;

% add analysis code directory, throw error if it doesn't exist
if ~exist([root_directory '/general-analysis-code'], 'dir')
    error('general-analysis-code directory is not in the root_directory');
end
addpath([root_directory '/general-analysis-code']);

% load weights from each run
n_runs = length(matfile_first_level);
for i = 1:n_runs
        
    % load weights for single run
    load(matfile_first_level{i}, ...
        'beta_one_per_regressor', 'beta_one_per_condition', 'P');
    
    % initialize weights
    if i == 1
        [n_regressors, n_voxels] = size(beta_one_per_regressor);
        n_conditions = size(beta_one_per_condition,1);
        R = nan(n_conditions, n_voxels, n_runs);
        W = nan(n_regressors, n_voxels, n_runs);
    end
    
    % assign weights for this run
    W(:,:,i) = beta_one_per_regressor;
    R(:,:,i) = beta_one_per_condition;
    
end

% voxels without any NaNs
voxels_without_NaNs = all(all(~isnan(R),1),3);

% differt numbers of runs to test when measuring split-half reliability
data_set_sizes = [2.^(0:1:log2(floor(n_runs/2))), floor(n_runs/2)];
data_set_sizes = unique(data_set_sizes);
n_data_set_sizes = length(data_set_sizes);

% measure test-retest correlation for different data set sizes
n_smps = 1e3;
mat_file = [analysis_directory ...
    '/regressor-response-reliability_' num2str(n_smps) 'smps.mat'];
if ~exist(mat_file, 'file') || optInputs(varargin, 'overwrite')
    fprintf('Sampling regressor response reliability...\n');
    corr_test_retest = nan(n_smps, n_data_set_sizes, n_regressors);
    for i = 1:n_data_set_sizes
        for j = 1:n_smps
            % random permutation of runs
            xi = randperm(n_runs);
            xi_split1 = xi(1:data_set_sizes(i));
            xi_split2 = xi(data_set_sizes(i)+1:2*data_set_sizes(i));
            
            R_split1 = mean(R(:,voxels_without_NaNs,xi_split1),3);
            R_split2 = mean(R(:,voxels_without_NaNs,xi_split2),3);
            W_split1 = mean(W(:,voxels_without_NaNs,xi_split1),3);
            W_split2 = mean(W(:,voxels_without_NaNs,xi_split2),3);
            
            R1 = R_split1 * pinv(W_split2);
            R2 = R_split2 * pinv(W_split1);
            
            corr_test_retest(j,i,:) = fastcorr(R1, R2);
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
figure_fname = [figure_directory '/' ...
    'regressor-response-reliability_' num2str(n_smps) 'smps.pdf'];
if ~exist(figure_fname, 'file') || optInputs(varargin, 'overwrite')
    figure;
    errorbar_plot_from_samples(corr_test_retest, log2(data_set_sizes));
    set(gca, 'XTick', log2(data_set_sizes), 'XTickLabel', data_set_sizes);
    xL = xlim;
    hold on; plot(xL, [0 0], 'k--');
    xlabel('Number of Runs'); ylabel('Test-Retest Correlation');
    legend(P.regressor_names,'Location','EastOutside');
    box off;
    title('Regressor Response Reliability');
    export_fig(figure_fname,'-pdf','-transparent');
end





