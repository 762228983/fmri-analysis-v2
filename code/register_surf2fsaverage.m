function register_surf2fsaverage(exp,us,runtype,r,hemi,input_fname,varargin)

global root_directory;

% FSL and freesurfer versions
[~, freesurfer_version] = read_version_info(exp,varargin{:});

subjid = ['us' num2str(us)];

% volume file to resample
preprocessing_directory = ...
    [root_directory '/' exp '/analysis/preprocess/' ...
    'usub' num2str(us) '/' runtype '_r' num2str(r)];

% surface_file
fsaverage_directory = [preprocessing_directory '/surf'];
surface_file = [fsaverage_directory '/' hemi '.' input_fname '.mgz'];

% map to fsaverage
fsaverage_directory = [preprocessing_directory '/myfsaverage'];
fsaverage_file = [fsaverage_directory '/' hemi '.' input_fname '.mgz'];
if ~exist(fsaverage_directory,'dir');
    mkdir(fsaverage_directory);
end

if ~exist(fsaverage_file,'file') || optInputs(varargin,'overwrite')
  unix_freesurfer_version(freesurfer_version, ...
      ['mri_surf2surf --srcsubject ' subjid ' --sval ' surface_file  ...
      ' --trgsubject myfsaverage --tval ' fsaverage_file  ' --hemi ' hemi]);
end