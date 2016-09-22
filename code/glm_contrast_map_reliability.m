function corr_test_retest = glm_contrast_map_reliability(...
    matfile_first_level, analysis_directory, figure_directory, varargin)

% global root_directory;
% 
% 2016-08-27: Modified how optional arguments are handled

% optional arguments and defaults
I.overwrite = false;
I = parse_optInputs_keyvalue(varargin, I);

% differt numbers of runs to test when measuring split-half reliability
n_runs = length(matfile_first_level);
data_set_sizes = [2.^(0:1:log2(floor(n_runs/2))), floor(n_runs/2)];
data_set_sizes = unique(data_set_sizes);
n_data_set_sizes = length(data_set_sizes);
    
% measure test-retest correlation for different data set sizes
n_smps = 1e3;
mat_file = [analysis_directory '/map-reliability_' num2str(n_smps) 'smps.mat'];
if ~exist(mat_file, 'file') || I.overwrite
    
    
    % load weights from each run
    for i = 1:n_runs
        
        % load weights for single run
        load(matfile_first_level{i}, 'beta_contrast', 'P');
        
        % initialize weights
        if i == 1
            [n_contrasts, n_voxels] = size(beta_contrast);
            contrast_allruns = nan([n_voxels, n_runs, n_contrasts]);
        end
        
        % assign weights for this run
        contrast_allruns(:,i,:) = beta_contrast';
        
    end
    
    % voxels without any NaNs
    voxels_without_NaNs = all(all(~isnan(contrast_allruns),2),3);
    

    fprintf('Sampling contrast map reliability...\n');
    corr_test_retest = nan(n_smps, n_data_set_sizes, n_contrasts);
    for i = 1:n_data_set_sizes
        for k = 1:n_contrasts
            single_contrast = contrast_allruns(voxels_without_NaNs,:,k);                  
            for j = 1:n_smps
                % random permutation of runs
                xi = randperm(n_runs);
                xi_split1 = xi(1:data_set_sizes(i));
                xi_split2 = xi(data_set_sizes(i)+1:2*data_set_sizes(i));
                split1 = single_contrast(:,xi_split1);
                split2 = single_contrast(:,xi_split2);
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
figure_fname = [figure_directory '/map-reliability_' num2str(n_smps) 'smps.pdf'];
if ~exist(figure_fname, 'file') || I.overwrite
    figure;
    errorbar_plot_from_samples(corr_test_retest, log2(data_set_sizes));
    set(gca, 'XTick', log2(data_set_sizes), 'XTickLabel', data_set_sizes);
    xL = xlim;
    hold on; plot(xL, [0 0], 'k--');
    xlabel('Number of Runs'); ylabel('Test-Retest Correlation');
    legend(P.contrast_names,'Location','EastOutside');
    box off;
    title('Contrast Map Reliability');
    export_fig(figure_fname,'-pdf','-transparent');
end





