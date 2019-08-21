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
% 
% 2018-05-28: Exclude voxels with mean timecourse near zero.
% 
% 2018-05-30: Add an optional tsnr threshold
% 
% 2019-08-20: Added ability to measure multiple timepoints in between onset
% and offset so as to measure a signal averaged timecourse

% optional arguments and defaults
I.onset_delay = 5;
I.offset_delay = 1;
I.n_perms = 0;
I.tsnr_threshold = 30;
I.n_tps = 1;
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

% temporal SNR
tsnr = mean(Y,1) ./ std(Y,[],1);

% remove voxels with NaN values
assert(max(mean(Y,1))>10);
voxels_without_NaN = all(~isnan(Y)) & (mean(Y,1) > 1e-10) & (tsnr > I.tsnr_threshold);
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
response = nan(n_trials, n_voxels_without_NaN, I.n_tps);
for i = 1:n_trials
    if I.n_tps == 1
        xi = t >= T.onsets(i) + I.onset_delay ...
            & t <= T.onsets(i) + T.durs(i) + I.offset_delay;
        response(i,:) = mean(Y(xi,:),1);
    else
        win = T.onsets(i) + linspace(I.onset_delay, T.durs(i) + I.offset_delay, I.n_tps)';
        assert(win(1)>=t(1) && win(2)<=t(end));
        response(i,:,:) = interp1(t, Y, win, 'linear')';
    end
end
clear xi t nTR;

%% Measure null and convert to psc

% mean response to null
% -> 1 x voxels x time
xi = strcmp('NULL', T.conds);
assert(sum(xi) > 0);
null_response = mean(response(xi,:,:),1);
clear xi;

%% Remove trials without corresponding condition, except for NULL trials

xi = ismember(T.conds, [P.condition_names(:); 'NULL']);
response  = response(xi,:, :);
T.conds = T.conds(xi,:);
T.condition_indices = T.condition_indices(xi,:);
T.onsets = T.onsets(xi,:);
T.durs = T.durs(xi,:);
n_trials = sum(xi);
clear xi;

%% convert to psc

% -> trials x voxels x timepoints
null_response_replicated = repmat(null_response, [n_trials, 1, 1]);
psc = 100 * (response - null_response_replicated) ./ null_response_replicated;

