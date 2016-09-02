function Y = fillin_NaN_voxels(X, voxels_without_NaNs, DIM)

% Helper function for glm_event_regression.m and sigav.m

assert(size(X,DIM)==sum(voxels_without_NaNs));
dims = size(X);
dims(DIM) = length(voxels_without_NaNs);

xi = cell(1,length(dims));
for i = 1:length(dims)
    if i == DIM
        xi{i} = find(voxels_without_NaNs);
    else
        xi{i} = 1:dims(i);
    end
end

Y = nan(dims);
Y(xi{:}) = X;