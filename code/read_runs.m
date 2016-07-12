function [runnum, sequenceid, scanid] = read_runs(exp, us, runtype, varargin)

% Returns the run numbers and corresponding sequence and scan ids for a particular experiment,
% subject and runtype. 
% 
% 2016-04-6: Last modified by Sam NH


global root_directory;

% file name with run order information
fname = [root_directory '/' exp '/data/brain/runorders/' runtype '_us' num2str(us) '.txt'];

% test that file exists
if ~exist(fname,'file');
    error('Run order file doesn''t exist:\n%s\n', fname);
end

% read file
fid = fopen(fname,'r');
x = textscan(fid,'%d%d%d'); 
fclose(fid);
runnum = x{1};
sequenceid = x{2};
scanid = x{3};