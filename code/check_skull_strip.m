function check_skull_strip(exp, us, runtype, runs, varargin)

% function check_skull_strip(exp, us, runtype, runs, varargin)
% useful function for check success of skull-stripping algorithm

global root_directory;

freesurfer_version = read_freesurfer_version(exp,varargin{:});
for i = 1:length(runs)
    preprocessing_direcotry = [root_directory '/' exp '/analysis/preprocess/usub' num2str(us) '/' runtype '_r' num2str(runs(i)) '/'];
    unix_freesurfer_version(freesurfer_version, ['freeview ' preprocessing_direcotry 'example_func.nii.gz:grayscale=0,2000 ' preprocessing_direcotry 'brain_mask.nii.gz:colormap=heat:opacity=0.5 &']);
end