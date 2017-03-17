function [psc, mean_signal, T, voxels_without_NaN, null_response] = ...
    sigav_data(data_matrix_file, para_file, parameter_file, varargin)

% Calculates signal averaged responses from a data matrix
% 
% Calculates percent signal change for a set of events/stimuli by signal
% averaging a fixed number of time-points after event onset. A set of predictors
% is then regressed against the signal averaged values, and the beta weights
% from this analysis are contrasted. Stats are computed using vanilla OLS
% equations, assuming independent, Gaussian errors, and a permutation test is
% used to compute stats by shuffling the order of conditions.
% 
% 2017-01-12: Created, Sam NH
% 
% 2017-02-04: Unspecified conditions are always removed, except for NULL

% optional arguments and defaults
I.onset_delay = 5;
I.offset_delay = 1;
I.n_perms = 0;
I = parse_optInputs_keyvalue(varargin, I);

% analysis parameters
P = load(parameter_file);
P = glm_default_parameters(P);

%% Format data matrix

% load data matrix
% -> time x voxel
load(data_matrix_file, 'data_matrix', 'TR');
Y = data_matrix;
clear data_matrix;

% remove voxels with NaN values
voxels_without_NaN = all(~isnan(Y));
Y = Y(:,voxels_without_NaN);
n_voxels_without_NaN = sum(voxels_without_NaN);

% average signal for each voxel
% 1 x voxels
mean_signal = mean(Y,1); 

%% Optionally detrend

if P.detrend > 0
    Y = detrend_poly(Y, P.detrend, 'restore_mean', true);
end

%% Average response to each event

% timing information about each event
T = read_para(para_file);

n_trials = length(T.onsets);
n_TR = size(Y,1);
t = (0:n_TR-1)*TR;
response = nan(n_trials, n_voxels_without_NaN);
for i = 1:n_trials
    xi = t >= T.onsets(i) + I.onset_delay ...
        & t <= T.onsets(i) + T.durs(i) + I.offset_delay;
    response(i,:) = mean(Y(xi,:),1);
end
clear xi t nTR;

%% Measure null and convert to psc

% mean response to null
% -> 1 x voxels
xi = strcmp('NULL', T.conds);
assert(sum(xi) > 0);
null_response = mean(response(xi,:),1);
clear xi;

%% Remove trials without corresponding condition, except for NULL trials

xi = ismember(T.conds, [P.condition_names(:); 'NULL']);
response  = response(xi,:);
T.conds = T.conds(xi,:);
T.condition_indices = T.condition_indices(xi,:);
T.onsets = T.onsets(xi,:);
T.durs = T.durs(xi,:);
n_trials = sum(xi);
clear xi;

%% convert to psc

% -> trials x voxels
psc = 100 * (response - repmat(null_response, n_trials, 1)) ...
    ./ repmat(null_response, n_trials, 1);

