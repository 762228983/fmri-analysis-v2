function bet_struct(us,favalue,varargin)

% function bet(usubs,varargin)
%
% bet skull-strips the reoriented data using
% fsl's skull-stripping utility
% fa values must be specified for each new subject

global root_directory

input_file = [root_directory '/anatomicals/us' num2str(us) '/struct.nii.gz'];
bet_file = [root_directory '/anatomicals/us' num2str(us) '/bet.nii.gz'];
bet_mask_file = [root_directory '/anatomicals/us' num2str(us) '/bet_mask.nii.gz'];

if ~exist(bet_file,'file') ||  ~exist(bet_mask_file,'file') ||optInputs(varargin, 'overwrite')
    unix_fsl('5.0', ['bet2 ' input_file ' ' strrep(bet_file,'.nii.gz','') ' -f ' num2str(favalue) ' -m']);
end

if optInputs(varargin, 'check')
    unix_freesurfer_version('5.3.0', ['freeview ' input_file ':grayscale=0,2000 ' bet_mask_file ':colormap=heat:opacity=0.5 &']);
end