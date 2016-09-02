function P = glm_default_parameters(P)

if ~isfield(P, 'linear_trend')
    P.linear_trend = false;
end

if ~isfield(P, 'n_whitematter_PCs')
    P.n_whitematter_PCs = 0;
end

% check to make sure there is not an errant parameter specified
% e.g. due to a mispelling of an intended parameter
possible_parameters = {'condition_names', 'regressor_names', ...
    'regressor_weights', 'contrast_names', 'contrast_weights', ...
    'linear_trend', 'n_whitematter_PCs'};
all_parameters = fieldnames(P);
for i = 1:length(all_parameters)
    if ~any(strcmp(all_parameters{i}, possible_parameters))
        error('%s not a recognized parameter.', all_parameters{i});
    end
end