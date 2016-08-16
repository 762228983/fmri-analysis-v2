function recon_all(us,varargin)

% function recon_all(s)
% 
% runs recon-all on a subject, should only be called once per subject
% the command is run with nohup, so it is fine to exit the terminal
% after starting the reconstruction.

global root_directory;

freesurfer_directory = [root_directory '/freesurfer'];
if ~exist(freesurfer_directory,'dir'); 
    mkdir(freesurfer_directory); 
end

log_directory = [freesurfer_directory '/log'];
if ~exist(log_directory,'dir'); 
    mkdir(log_directory); 
end

% subject id used by freesurfer
subjid = ['us' num2str(us)];

input_file = [root_directory '/anatomicals/' subjid '/struct.nii.gz'];

freesurfer_subject_directory = [freesurfer_directory '/' subjid];
if exist(freesurfer_subject_directory,'dir') && optInputs(varargin, 'overwrite')
    unix(['rm -f ' freesurfer_subject_directory]);
end

if ~exist(freesurfer_subject_directory,'dir')
    unix_freesurfer_version('5.3.0', ['recon-all -i ' input_file ' -subjid ' subjid]);
    fprintf('Creating Freesurfer directory:\n\n%s\n\n', freesurfer_subject_directory);
end

if ~exist([freesurfer_subject_directory '/surf/rh.inflated'],'file') || optInputs(varargin, 'overwrite');
    log_file = [log_directory '/recon-all_' subjid '.txt'];
    unix_freesurfer_version('5.3.0', ['nohup recon-all -all -subjid ' subjid ' > ' log_file ' &']);
    fprintf('Surface reconstruction begun\n\nCheck log-file for progress:\n\n%s\n', log_file);
end
