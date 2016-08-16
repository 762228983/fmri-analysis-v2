function [fsl_version, freesurfer_version] = read_version_info(exp, varargin)

fsl_version = '5.0';
freesurfer_version = '5.3.0';

% global root_directory
% 
% fname = [root_directory '/' exp '/code/version-info.txt'];
% 
% % throw an error if the file doesn't exist
% if ~exist(fname, 'file');
%     error('%s file does not exist\n');
% end
% 
% % read the contents of the file
% fid = fopen(fname, 'r');
% x = textscan(fid, '%s%s');
% fclose(fid);
% 
% % extract fsl and freesurfer version from contents
% software_names = x{1};
% software_versions = x{2};
% fclose all;
% 
% xi = strcmp('fsl-version', software_names);
% fsl_version = software_versions{xi};
% 
% xi = strcmp('freesurfer-version', software_names);
% freesurfer_version = software_versions{xi};