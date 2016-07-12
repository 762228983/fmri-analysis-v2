function runtypes = read_runtypes(exp,us,varargin)

% function runtypes = read_runtypes(exp,us,varargin)
% 
% Returns the runtypes for a particular subject.
% First looks for a subject-specific file with runs the particular subject of interest.
% If it can find the subject-specific file, it then looks for a generic file with runtypes for all
% subjects.
% 
% 2016-04-6: Last modified by Sam NH

global root_directory;

% file with runtypes for a specific subject
fname_subject_specific = [root_directory '/' exp '/data/runtypes/runtypes_for_us' num2str(us) '.txt'];

% file to use if there is no subject-specific file 
fname_allsubjects = [root_directory '/' exp '/data/runtypes/runtypes_for_all_subjects.txt'];

% read runtypes from file, throw error if file not found
if exist(fname_subject_specific, 'file')
    fid = fopen(fname_subject_specific, 'r');
    x = textscan(fid,'%s'); fclose(fid);
    runtypes = x{1};
elseif exist(fname_allsubjects, 'file')
    fid = fopen(fname_allsubjects, 'r');
    x = textscan(fid,'%s'); fclose(fid);
    runtypes = x{1};
else
    error('No runtype file found. One of following two files should exist:\n%s\n%s', fname_subject_specific, fname_allsubjects)
end
