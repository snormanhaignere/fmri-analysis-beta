function recon_all(us,varargin)

% function recon_all(s)
% 
% runs recon-all on a subject, should only be called once per subject
% the command is run with nohup, so it is fine to exit the terminal
% after starting the reconstruction.
% 
% 2016-08-27: Modified how optional arguments are handled

global root_directory;

% optional arguments and defaults
I.overwrite = false;
I.freesurfer_version = '5.3.0';
I.localGI = false;
I.runanyways = false;
I = parse_optInputs_keyvalue(varargin, I);

freesurfer_directory = [root_directory '/freesurfer'];
if ~exist(freesurfer_directory,'dir')
    mkdir(freesurfer_directory); 
end

log_directory = [freesurfer_directory '/log'];
if ~exist(log_directory,'dir')
    mkdir(log_directory); 
end

% subject id used by freesurfer
subjid = ['us' num2str(us)];

input_file = [root_directory '/anatomicals/' subjid '/struct.nii.gz'];

freesurfer_subject_directory = [freesurfer_directory '/' subjid];
if exist(freesurfer_subject_directory,'dir') && I.overwrite
    unix(['rm -f ' freesurfer_subject_directory]);
end

if ~exist(freesurfer_subject_directory,'dir')
    unix_freesurfer_version(I.freesurfer_version, ['recon-all -i ' input_file ' -subjid ' subjid]);
    fprintf('Creating Freesurfer directory:\n\n%s\n\n', freesurfer_subject_directory);
end

flags = ' ';
if I.localGI
    flags = [flags '-localGI '];
end

if ~exist([freesurfer_subject_directory '/surf/rh.inflated'],'file') || I.overwrite || I.runanyways
    log_file = [log_directory '/recon-all_' subjid '.txt'];
    unix_freesurfer_version(I.freesurfer_version, ['nohup recon-all -all' flags '-subjid ' subjid ' > ' log_file ' &']);
    fprintf('Surface reconstruction begun\n\nCheck log-file for progress:\n\n%s\n', log_file);
end
