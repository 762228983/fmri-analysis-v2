function register_func2highres_flirt(exp, us, runtype, r, varargin)

% function register_func2highres_flirt(exp, us, runtype, r, varargin)
% registers a functional image to a highres anatomical using FSL's
% linear registration algorithm, flirt

global root_directory;

% degrees of freedom, 6 = rigid, 12 = affine
dof = '6';

% range of search
search_range = '-180 180';
if optInputs(varargin, '90')
    search_range = '-90 90';
end

% addpath(genpath('/software/Freesurfer/5.3.0/matlab'));
% fsl version for this experiment
[fsl_version, freesurfer_version] = read_version_info(exp,varargin{:});

% directory to save registration files to
preprocessing_directory = [root_directory '/' exp '/analysis/preprocess/usub' num2str(us) '/' runtype '_r' num2str(r)];
registration_directory = [preprocessing_directory '/reg_flirt'];
if ~exist(registration_directory, 'dir')
  mkdir(registration_directory);
end

% directory with freesurfer reconstruction
subjid = ['us' num2str(us)];
freesurfer_directory = [root_directory '/freesurfer/' subjid];

% brain-masked functional image used for motion correction
exfunc = [preprocessing_directory '/example_func_brain_masked.nii.gz']; 

% brain-masked structural image in custom MGH format
highres_mgz = [freesurfer_directory '/mri/brainmask.mgz']; 

% brain-masked structural image converted to nifti format for use by FSL
highres_nifti = [freesurfer_directory '/mri/brainmask.nii.gz']; 
if ~exist(highres_nifti, 'file') || optInputs(varargin, 'overwrite')
    unix_freesurfer_version(freesurfer_version, ['mri_convert ' highres_mgz ' ' highres_nifti]);
end

% calculate functional to anatomical registration matrix
exfunc2highres_regmat_fsl = [registration_directory '/example_func2highres.mat']; % functional to anatomical registration matrix
if ~exist(exfunc2highres_regmat_fsl, 'file') || optInputs(varargin, 'overwrite')
    unix_fsl(fsl_version, ['flirt -in ' exfunc ' -ref ' highres_nifti ' -cost corratio -dof ' dof ' -searchrx ' search_range ' -searchry ' search_range ' -searchrz ' search_range ' -omat ' exfunc2highres_regmat_fsl]);
end

% invert registration
highres2exfunc_regmat_fsl = [registration_directory '/highres2example_func.mat']; % anatomical to functional registration matrix
if ~exist(highres2exfunc_regmat_fsl, 'file') || optInputs(varargin, 'overwrite')
    unix_fsl(fsl_version, ['convert_xfm -omat ' highres2exfunc_regmat_fsl ' -inverse ' exfunc2highres_regmat_fsl]);
end

% convert to freesurfer style matrix
exfunc2highres_regmat_freesurf = [registration_directory '/example_func2highres.dat'];
if ~exist(exfunc2highres_regmat_freesurf, 'file') || optInputs(varargin, 'overwrite')
    unix_freesurfer_version( freesurfer_version, ['tkregister2 --s ' subjid ' --mov ' exfunc ' --targ ' highres_mgz ' --reg ' exfunc2highres_regmat_freesurf ' --fsl ' exfunc2highres_regmat_fsl ' --noedit ']);
end

% interpolate functional to anatomical to evaluate success
exfunc2highres_nifti = [registration_directory '/example_func2highres.nii.gz']; % functional image interpolated to anatomical
if ~exist(exfunc2highres_nifti, 'file') || optInputs(varargin, 'overwrite')
    unix_fsl(fsl_version, ['flirt -interp trilinear -in ' exfunc ' -ref ' highres_nifti ' -applyxfm -init ' exfunc2highres_regmat_fsl ' -out ' exfunc2highres_nifti]);
end

% view results using freeview
if optInputs(varargin, 'freeview')
  unix_freesurfer_version(freesurfer_version,['freeview ' highres_mgz ':grayscale=10,150 ' exfunc2highres_nifti ':grayscale=0,2000 &'])
end
  
% view results using tkregister2
if optInputs(varargin, 'tkregister2')
  tkregister_title = ['us' num2str(us) '-r' num2str(r) '-flirt'];
  unix_freesurfer_version( freesurfer_version, ['tkregister2 --s ' subjid ' --mov ' exfunc ' --targ ' highres_mgz ' --reg ' exfunc2highres_regmat_freesurf ' --title ' tkregister_title ' --surf &']);
end