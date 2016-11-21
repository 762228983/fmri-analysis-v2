function P = regress_predictions_voxelgrid(G, F, output_directory, varargin)

% Takes a feature matrix F (stim x feature) and voxel data structure G in the form
% of a grid, and computes predictions for each voxel's response
% 
% Alternatively, F and G can be mat files which contain the relevant
% matrix/structure. The variables in each file are assumed to named 'F' and 'D',
% such that load(F, 'F') and load(G, 'G') read in the appropriate variables. 
% 
% MAT_file is an optional MAT file that the results can be saved to.  

if ischar(F)
    load(F, 'F');
end

if ischar(G)
    load(G, 'G');
end

assert(all(size(D) == [1,2]));

addpath([root_directory '/general-analysis-code']);

I.test_folds = 2;
I.train_folds = min(max(round((size(F,1) / I.test_folds) / 10),2),10);
I.regression_method = 'ridge';
I.K = [];
I.use_sbatch = false;
I = parse_optInputs_keyvalue(varargin, I);

if I.use_sbatch
    
    Yh = regress_predictions_parallelize_with_slurm(...
        F, grid2matrix(G), test_folds, method, K, train_folds, output_directory);
    
else
    
    D = grid2matrix(G);
    Yh = size(D);
    for i = size(D,2)
        Yh(:,i) = regress_predictions_from_3way_crossval(F, D(:,i), ...
            I.test_folds, I.method, I.K, I.train_folds, output_directory);
    end
end

P = matrix2grid(Yh);

fprintf('Regression against feature vectors\n');
hemis = {'rh', 'lh'};

for q = 1:length(hemis)
    
    fprintf('%s\n', hemis{q}); drawnow;
    
    xdim = size(D{q},2);
    ydim = size(D{q},3);
    weights{q} = nan(n_features, xdim, ydim);
    bestK{q} = nan(xdim, ydim);
    mse_bestK{q} = nan(xdim, ydim);

    fprintf('xdim: ');
    for i = 1:xdim
        fprintf('%d ', i); drawnow;
        if mod(i,10)==0
            fprintf('\n');
        end
        for j = 1:ydim
            if all(~isnan(D{q}(:,i,j)))
                [b, bestK{q}(i,j), ~, ~, mse_bestK{q}(i,j)] = ...
                    regress_weights_from_2way_crossval(...
                    F, D{q}(:,i,j), folds, regression_type, K);
                weights{q}(:,i,j) = b(2:end);
            end
        end
    end
end

if nargin > 2
    save(MAT_file, 'weights', 'mse_bestK', 'bestK');
end