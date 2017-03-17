function [TR, TA, nTR, n_disdaqs, fname] = ...
    read_functional_scan_parameters(exp,us,runtype,r,varargin)

% function [TR, TA, nTR, n_disdaqs] = read_functional_scan_parameters(exp,us,runtype,r)
% 
% Returns the TR, acquisition time / TA, number of TRs and the number of discarded acquisitions (disdaqs)
% for a functional run
% 
% 2016-04-6: Last modified by Sam NH
% 
% 2017-01-27: Made it possible to have multiple runs to account for combinations
% of runs.
% 
% 2017-02-23: Now returns the file it read from

global root_directory;

% file with parameters for a specific run and subject
fname_run_specific = [root_directory '/' exp '/data/brain/scan-parameters/' ...
    runtype '_us' num2str(us) '_r' sprintf('%d',r) '.txt'];

% general-purpose file if subject specific file is not found
fname_general = [root_directory '/' exp '/data/brain/scan-parameters/' runtype '.txt'];

% read runtypes from file, throw error if file not found
if exist(fname_run_specific, 'file')
    fname = fname_run_specific;
elseif exist(fname_general, 'file')
    fname = fname_general;
else
    error('No runtype file found. One of following two files should exist:\n%s\n%s', ...
        fname_run_specific, fname_general);
end

fid = fopen(fname);

% read and assign variable names
x = textscan(fid, '%s%f');
variable_name = x{1};
variable_value = x{2};

xi = strcmp('TR', variable_name);
TR = variable_value(xi);

xi = strcmp('TA', variable_name);
TA = variable_value(xi);

xi = strcmp('nTR', variable_name);
nTR = variable_value(xi);

xi = strcmp('nDISDAQs', variable_name);
n_disdaqs = variable_value(xi);

% close matlab file identifier
fclose(fid);

