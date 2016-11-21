function regress_voxelgrid(D, F)

% Takes a feature matrix F (stim x feature) and voxel data structure in the form
% of a grid (X x Y x stim), computes weights and mean-squared errors using
% cross-validated leave-one-out regression.
% 
% Alternatively, F and D can be mat files which contain the relevant
% matrix/structure. The variables in each file are assumed to named 'F' and 'D',
% such that load(F, 'F') and load(D, 'D') read in the appropriate variables. 
% 
% MAT_file is an optional MAT file that the results can be saved to.  

if ischar(F)
    load(F, 'F');
end

if ischar(D)
    load(D, 'D');
end

assert(all(size(D) == [1,2]));

addpath([root_directory '/general-analysis-code']);

regression_type = 'ridge';
K = 2.^(-30:50);

% change number of folds
folds = 10;
if optInputs(varargin, 'folds')
    folds = varargin{optInputs(varargin, 'folds')+1};
end

fprintf('Regression against feature vectors\n');

hemis = {'rh', 'lh'};
weights = cell(1,2);
bestK = cell(1,2);
mse_bestK = cell(1,2);
n_features = size(F,2);

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