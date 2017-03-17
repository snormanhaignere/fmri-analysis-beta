function [W, W_one_per_condition, C, ...
    nonzero_regressors, nonzero_contrasts, nonzero_conditions] = ...
    sigav_weights_and_contrasts(P, T)

% weights applied to each condition for each regressor
% -> trial x regressor
n_trials = length(T.onsets);
n_regressors = length(P.regressor_names);
W = zeros(n_trials,n_regressors);
for i = 1:n_trials
    xi = strcmp(T.conds{i}, P.condition_names);
    if any(xi)
        W(i,:) = P.regressor_weights(xi,:);
    end
end
nonzero_regressors = any(W~=0,1);
W = W(:,nonzero_regressors);

% weights with one regressor per condition
% -> trial x condition
n_conditions = length(P.condition_names);
W_one_per_condition = zeros(n_trials,n_conditions);
for i = 1:n_trials
    xi = strcmp(T.conds{i}, P.condition_names);
    W_one_per_condition(i,xi) = 1;
end
nonzero_conditions = any(W_one_per_condition~=0,1);
W_one_per_condition = W_one_per_condition(:,nonzero_conditions);

% contrasts with any nonzero weights
nonzero_contrasts = any(P.contrast_weights(nonzero_regressors,:)~=0,1);
C = P.contrast_weights(nonzero_regressors, nonzero_contrasts);

%% Check that zero-mean contrasts remain so when excluding zero regressors

n_contrasts = length(P.contrast_names);
for i = 1:n_contrasts
    if abs(sum(P.contrast_weights(:,i))) < 1e-10 ...
            && abs(sum(P.contrast_weights(nonzero_regressors,i))) > 1e-10
        sprintf(['Contrast "%s" is no longer zero mean '...
            'when excluding zero regressors\n'], P.contrast_names{i});
        keyboard;
    end
end