function P = glm_default_parameters(P)

if ~isfield(P, 'detrend')
    P.detrend = 1;
end

if ~isfield(P, 'n_whitematter_PCs')
    P.n_whitematter_PCs = 0;
end

% check to make sure there is not an errant parameter specified
% e.g. due to a mispelling of an intended parameter
possible_parameters = {'condition_names', 'regressor_names', ...
    'regressor_weights', 'contrast_names', 'contrast_weights', ...
    'detrend', 'n_whitematter_PCs','para_prefix'};
all_parameters = fieldnames(P);
for i = 1:length(all_parameters)
    if ~any(strcmp(all_parameters{i}, possible_parameters))
        error('%s not a recognized parameter.', all_parameters{i});
    end
end