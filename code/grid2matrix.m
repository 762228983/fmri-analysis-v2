function X = grid2matrix(G)

% unwrap a grid file to a matrix;

dims_rh = size(G.grid_data{1});
dims_lh = size(G.grid_data{2});

if length(dims_rh) < 3
    dims_rh = [dims_rh, 1];
end

if length(dims_lh) < 3
    dims_lh = [dims_lh, 1];
end
    
X = [...
    reshape(G.grid_data{1}, [prod(dims_rh(1:2)),dims_rh(3:end)])', ...
    reshape(G.grid_data{2}, [prod(dims_lh(1:2)),dims_lh(3:end)])'];